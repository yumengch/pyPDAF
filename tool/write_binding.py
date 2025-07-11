import re
from pathlib import Path
import textwrap


WRAP_WIDTH = 80
TYPE_MAP = {
  'integer':      ('INTEGER',    'c_int'),
  'real':         ('REAL',       'c_double'),
  'logical':      ('LOGICAL',    'c_bool'),
  'complex':      ('COMPLEX',    'c_double_complex'),
  'character':    ('CHARACTER',  'c_char')
}


def split_outside_parens(s: str):
    """
    Split s on commas that are not enclosed in parentheses.
    E.g. "a(:,:), b, c(:)" → ["a(:,:)", " b", " c(:)"]
    """
    parts = []
    buf = ""
    depth = 0
    for ch in s:
        if ch == "(":
            depth += 1
            buf += ch
        elif ch == ")":
            depth -= 1
            buf += ch
        elif ch == "," and depth == 0:
            parts.append(buf)
            buf = ""
        else:
            buf += ch
    if buf:
        parts.append(buf)
    return parts


def _process_buffer(buffer, buffer_comments, arg_names, decls, comments):
    """
    Helper: split `buffer` on '::', then split RHS on commas outside parens,
    match to arg_names, record decl_code and join buffer_comments into one.
    """
    try:
        left, right = buffer.split("::", 1)
    except ValueError:
        return

    decl_code = left.strip().lower()
    var_list  = right.strip()
    if len(buffer_comments) == 0:
        buffer_comments = ["" for _ in split_outside_parens(var_list)]

    for chunk, comment in zip(split_outside_parens(var_list), buffer_comments):
        name_with_shape = chunk.strip()
        # match base name + optional shape
        m = re.match(r'(\w+)\s*(\([^\)]*\))?', name_with_shape)
        if not m:
            continue
        base = m.group(1).lower()
        shape = m.group(2) or ""
        if base in arg_names:
            decls.append((base, decl_code, shape))
            if comment:
                comments[base] = comment


def wrap_comma_list(prefix, items, subsequent_indent="  ", suffix=""):
    """
    Wrap a comma-separated list under WRAP_WIDTH, inserting '&'
    at end of every line except the last, and appending `suffix`
    (e.g. ' bind(c)') to that last line.
    """
    lines = []
    current = prefix
    for i, item in enumerate(items):
        sep = ", " if i < len(items) - 1 else ")"
        addition = item + sep
        # if adding this would exceed limit (allowing space for ' &' on wrapped lines)
        if len(current) + len(addition) + (3 if sep==", " else len(suffix)) > WRAP_WIDTH:
            lines.append(current + " &")
            current = subsequent_indent + item + sep
        else:
            current += addition
    # if items were empty, we might have just the prefix
    if current == prefix:
        current += ")"

    # append suffix to the final line
    if suffix:
        current = current + suffix
    lines.append(current)
    return lines


def extract_subroutine_blocks(lines):
    """
    Return a list of [start_line_idx, end_line_idx, block_lines] for each
    SUBROUTINE … END SUBROUTINE found **outside** any INTERFACE / ABSTRACT INTERFACE.
    """
    blocks = []
    interface_depth = 0
    i = 0
    n = len(lines)

    while i < n:
        line = lines[i]

        # Entering an INTERFACE or ABSTRACT INTERFACE
        if re.match(r'\s*(abstract\s+)?interface\b', line, re.I):
            interface_depth += 1
            i += 1
            continue

        # Leaving an INTERFACE block
        if re.match(r'\s*end\s+interface\b', line, re.I):
            if interface_depth > 0:
                interface_depth -= 1
            i += 1
            continue

        # Only consider SUBROUTINE if we're *not* inside an interface
        if interface_depth == 0 and re.match(r'\s*subroutine\s+', line, re.I):
            start = i
            # find the matching END SUBROUTINE
            i += 1
            while i < n and not re.match(r'\s*end\s+subroutine\b', lines[i], re.I):
                i += 1
            # include the END SUBROUTINE line
            if i < n:
                end = i
                blocks.append(lines[start:end+1])
            i += 1
            continue

        i += 1

    return blocks


def extract_multiline_signature(block):
    sig_lines = []
    inside = False
    start_idx = None
    for idx, line in enumerate(block):
        if not inside and re.match(r'\s*subroutine\s+', line, re.I):
            inside = True
            start_idx = idx
        if inside:
            sig_lines.append(line.strip())
            if ")" in line:
                break
    if not sig_lines:
        raise ValueError("No subroutine signature found")

    # Join and clean
    full_sig = " ".join(sig_lines).replace("&", "")
    match = re.match(r'\s*subroutine\s+(\w+)\s*\(([^)]*)\)', full_sig, re.I)
    if not match:
        raise ValueError(f"Could not parse subroutine signature from:\n{full_sig}")

    name = match.group(1)
    args = [x.strip().lower() for x in match.group(2).split(',') if x.strip()]
    return name, args, start_idx


def extract_declarations_by_args(block, start_idx, arg_names):
    """
    From block[start_idx:] find all “::” declarations, including
    multi-line EXTERNAL/PROCEDURE and type declarations.
    Returns:
      decls   = [(var_name, decl_code), …]
      comments = { var_name: comment_str, … }
    """
    decls = []
    comments = {}

    buffer = ""           # holds the code part across continuations
    buffer_comments = []   # list of comment fragments

    collecting = False
    for line in block[start_idx:]:
        # only start after IMPLICIT NONE
        if not collecting:
            if re.match(r'\s*implicit\s+none', line, re.I):
                collecting = True
            continue

        # stop at end of subroutine
        if re.match(r'\s*end\s+subroutine', line, re.I):
            break

        raw = line.rstrip("\n")
        # 1) split off any comment
        if "!" in raw:
            code_part, comment_part = raw.split("!", 1)
            comment_part = comment_part.replace("<", "", 1).strip()
        else:
            code_part, comment_part = raw, None

        # --- Start a new declaration if we see '::' and we're not in one ---
        if not buffer and "::" in code_part:
            buffer = code_part.strip()
            if comment_part:
                buffer_comments = [comment_part]
            else:
                buffer_comments = []
            # if this line ends with '&', stay in continuation
            if buffer.strip().endswith("&"):
                continuing_decl = True
                buffer = buffer.rstrip("&").rstrip()
            else:
                continuing_decl = False
                # single-line decl: process immediately
                _process_buffer(buffer, buffer_comments, arg_names, decls, comments)
                buffer = ""
                buffer_comments = []
            continue

        # --- If we're in a continuation, absorb blank/comment lines too ---
        if buffer and continuing_decl:
            # If there's no code at all, just collect comment and keep going
            if code_part.strip() == "":
                if comment_part:
                    buffer_comments[-1] = buffer_comments[-1] + comment_part
                continue

            # Otherwise, append code & comment
            buffer += " " + code_part.strip()
            if comment_part:
                buffer_comments.append(comment_part)

            # Still continuing?
            if buffer.strip().endswith("&"):
                continuing_decl = True
                buffer = buffer.rstrip("&").rstrip()
                continue
            else:
                continuing_decl = False
                # now we've got the full multi-line decl
                _process_buffer(buffer, buffer_comments, arg_names, decls, comments)
                buffer = ""
                buffer_comments = []
            continue

    return decls, comments


def compare_signature_and_declarations(arg_list, decls):
    """
    Given:
      arg_list : ['u_collect_state', 'u_obs_op', 'u_init_obs', ...]
      decls    : [('u_collect_state', 'external'), ('u_obs_op', 'external'), ...]
    Returns:
      missing : [ names in arg_list but not in decls ]
      extra   : [ names in decls    but not in arg_list ]
    """
    decl_vars = {var for var, _, _ in decls}
    missing = [arg for arg in arg_list if arg not in decl_vars]
    extra   = [var for var, _, _ in decls   if var not in arg_list]
    return missing, extra


def generate_bindc_wrapper(name, arg_list, decls, comments):
    """
    name      : original subroutine name (string)
    arg_list  : list of argument names in order, already lowercased
    decls     : list of (var, decl_code) from extract_declarations_by_args
    comments  : dict var->comment_str
    """
    wrapper_name = f"c__{name}"

    # detect obs_f / obs_l arguments
    obs_args = [
        var for var, code, _ in decls
        if re.match(r"type\s*\(\s*obs_f\s*\)", code, re.I)
        or re.match(r"type\s*\(\s*obs_l\s*\)", code, re.I)
    ]

    # All other args stay in the C-interface
    non_obs = [
            v for v in arg_list
            if v not in obs_args
        ]

    # Build the C-binding argument list:
    c_args = []
    if obs_args:
        c_args.append("i_obs")
    c_args = c_args + non_obs

    # Wrap signature with continuations and bind(c)
    sig_lines = wrap_comma_list(
        prefix=f"SUBROUTINE c__{name}(",
        items=c_args,
        subsequent_indent="   ",
        suffix=" bind(c)"
    )

    lines = []
    lines.extend(sig_lines)

    if len(c_args) > 0:
        lines.append("   use iso_c_binding")
        lines.append("")

        # categorize args
        typed_args    = []
        external_args = []
        decl_map = { var: (code, shape) for var, code, shape in decls }

        # declare i_obs if needed
        if obs_args:
            lines.append("   ! index into observation arrays")
            lines.append("   INTEGER(c_int), INTENT(in) :: i_obs")
            lines.append("")

        for var in non_obs:
            code, _ = decl_map.get(var, "")
            if code.startswith("external") or "procedure" in code:
                external_args.append(var)
            else:
                typed_args.append(var)

        # emit typed args
        for var in typed_args:
            code, shape = decl_map[var]
            # figure out fortran base type → (F_TYPE, C_KIND)
            base_type = re.match(r'(integer|real|logical|complex|character)', code).group(1).lower()
            f_type, c_kind = TYPE_MAP[base_type]
            # pointer if it has 'pointer' in the code
            pointer = ", POINTER" if "pointer" in code.lower() else ""
            # preserve shape if present
            shape_attr = f", DIMENSION{shape}" if shape else ""
            # keep intent
            intent = re.search(r'intent\s*\(\s*(in|out|inout)\s*\)', code, re.I)
            intent_str = f", INTENT({intent.group(1).lower()})" if intent else ""
            # comment
            if comments.get(var):
                lines.append(f"   ! {comments[var]}")
            else:
                lines.append(f"   !")
            lines.append(f"   {f_type}({c_kind}){pointer}{shape_attr}{intent_str} :: {var}")
        lines.append("")

        # emit procedure args
        for var in external_args:
            # comment
            if comments.get(var):
                lines.append(f"   ! {comments[var]}")
            else:
                lines.append(f"   !")
            proc_type = f"c__{var}_pdaf"
            lines.append(f"   procedure({proc_type}) :: {var}")
        lines.append("")

    # build the CALL to the original:
    call_args = []
    for v in arg_list:
        if v in obs_args:
            call_args.append(f"{v}(i_obs)")
        else:
            call_args.append(v)

    call_lines = wrap_comma_list(
        f"   call {name}(",
        call_args,
        subsequent_indent="      "
    )
    lines.extend(call_lines)
    # 8) end
    lines.append(f"END SUBROUTINE {wrapper_name}")
    lines.append("")  # blank after call

    return "\n".join(lines)


def process_file(filepath, output_dir):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    blocks = extract_subroutine_blocks(lines)

    if not blocks:
        print (f"No subroutines found in {filepath}")
        return  # No subroutines, skip

    for block in blocks:
        try:
            name, arg_list, start_idx = extract_multiline_signature(block)
            decls, comments = extract_declarations_by_args(block, start_idx, arg_list)
            missing, extra = compare_signature_and_declarations(arg_list, decls)
            if missing:
                print(f"!!! Missing declarations for: {missing} in {name} in {filepath}")
            if extra:
                print(f"!!! Declared but not in signature: {extra} in {name} in {filepath}")
            # print (decls)
            wrapper = generate_bindc_wrapper(name, arg_list, decls, comments)

            if 'pdafomi_assim' in name.lower():
                output_path = Path(output_dir) / 'pdafomi_assim_c_binding.f90'
            elif 'pdafomi_put' in name.lower():
                output_path = Path(output_dir) / 'pdafomi_put_c_binding.f90'
            elif '_cb' in name.lower():
                output_path = Path(output_dir) / 'pdaf_callback_c_binding.f90'
            elif 'pdafomi_set' in name.lower():
                output_path = Path(output_dir) / 'pdafomi_set_c_binding.f90'
            elif 'pdafomi' in name.lower():
                output_path = Path(output_dir) / 'pdafomi_c_binding.f90'
            elif 'pdaflocalomi_assim' in name.lower():
                output_path = Path(output_dir) / 'pdaflocalomi_assim_c_binding.f90'
            elif 'pdaflocalomi_put' in name.lower():
                output_path = Path(output_dir) / 'pdaflocalomi_put_c_binding.f90'
            elif 'pdaflocalomi' in name.lower():
                output_path = Path(output_dir) / 'pdaflocalomi_c_binding.f90'
            elif 'pdaflocal_assim' in name.lower():
                output_path = Path(output_dir) / 'pdaflocal_assim_c_binding.f90'
            elif 'pdaflocal_put' in name.lower():
                output_path = Path(output_dir) / 'pdaflocal_put_c_binding.f90'
            elif 'pdaflocal' in name.lower():
                output_path = Path(output_dir) / 'pdaflocal_c_binding.f90'
            elif 'pdaf_set' in name.lower():
                output_path = Path(output_dir) / 'pdaf_set_c_binding.f90'
            elif 'pdaf_get' in name.lower():
                output_path = Path(output_dir) / 'pdaf_get_c_binding.f90'
            elif 'pdaf_assim' in name.lower():
                output_path = Path(output_dir) / 'pdaf_assim_c_binding.f90'
            elif 'pdaf_put' in name.lower():
                output_path = Path(output_dir) / 'pdaf_put_c_binding.f90'
            elif 'pdaf_diag' in name.lower():
                output_path = Path(output_dir) / 'pdaf_diag_c_binding.f90'
            elif 'pdaf_iau' in name.lower():
                output_path = Path(output_dir) / 'pdaf_iau_c_binding.f90'
            elif 'pdaf3_assim' in name.lower():
                output_path = Path(output_dir) / 'pdaf3_assim_c_binding.f90'
            elif 'pdaf3_put' in name.lower():
                output_path = Path(output_dir) / 'pdaf3_put_c_binding.f90'
            else:
                output_path = Path(output_dir) / f"pdaf_c_binding.f90"
            output_path.parent.mkdir(parents=True, exist_ok=True)
            with open(output_path, 'a') as f:
                f.write(wrapper + "\n")
        except Exception as e:
            print(f"Skipping block in {filepath}: {e}")


def process_directory(src_dir, dst_dir):
    for path in Path(src_dir).rglob("*"):
        if path.suffix.lower() in (".f90", ".f", ".f95", ".f03", ".f08", ".for", ".ftn"):
            if path.name.lower() in ('pdaf_cb_procedures.f90', 'pdaf_assim_interfaces.f90',
                                     'pdaf_mod_core.f90', 'pdafomi.f90', 'pdaf.f90'):
                continue
            process_file(path, dst_dir)


if __name__ == "__main__":
    process_directory('PDAF-PDAF_V3.0beta/src', 'bindc_files')
