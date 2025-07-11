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

    for chunk in split_outside_parens(var_list):
        name_with_shape = chunk.strip()
        # match base name + optional shape
        m = re.match(r'(\w+)\s*(\([^\)]*\))?', name_with_shape)
        if not m:
            continue
        base = m.group(1).lower()
        m = re.search(r'dimension\s*\(([^)]+)\)', decl_code, re.IGNORECASE)
        shape = m.group(1) if m else ""
        if base in arg_names:
            decls.append((base, decl_code, shape))
            comments[base] = buffer_comments


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
        else:
            code_part, comment_part = raw, None

        # --- Start a new declaration if we see '::' and we're not in one ---
        if comment_part and buffer_comments == []:
            buffer_comments = [comment_part.strip()]
            comment_part = ""

        if not buffer and "::" in code_part:
            buffer = code_part.strip()
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

        if not buffer and code_part.strip() == "":
            if comment_part:
                buffer_comments.append( comment_part.strip() )
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


def get_pyx_arg(arg_list, decl_map):
    pyx_in_arg_list = []
    pyx_out_arg_list = []
    for arg in arg_list:
        code, _ = decl_map.get(arg, ("", ""))
        match = re.search(r'intent\(in\)|intent\(inout\)', code, re.IGNORECASE)
        if match:
            pyx_in_arg_list.append(arg)
        match = re.search(r'intent\(out\)|intent\(inout\)', code, re.IGNORECASE)
        if match:
            pyx_out_arg_list.append(arg)
    return pyx_in_arg_list, pyx_out_arg_list


def get_c_def(src_dir):
    filepath = Path(src_dir) / Path('pdaf_c_cb_interface.f90')
    with open(filepath, 'r') as f:
        lines = f.readlines()

    blocks = extract_subroutine_blocks(lines)

    if not blocks:
        print (f"No subroutines found in {filepath}")
        return  # No subroutines, skip

    cb_interface = {}
    for block in blocks:
        # try:
        name, arg_list, start_idx = extract_multiline_signature(block)
        decls, comments = extract_declarations_by_args(block, start_idx, arg_list)
        missing, extra = compare_signature_and_declarations(arg_list, decls)
        if missing:
            print(f"!!! Missing declarations for: {missing} in {name} in {filepath}")
        if extra:
            print(f"!!! Declared but not in signature: {extra} in {name} in {filepath}")

        name = name.lower()
        cb_interface[name] = {
            'args': arg_list,
            'decls': decls,
            'comments': comments
        }
        # except Exception as e:
        #     print(f"Skipping block in {filepath}: {e}")

    return cb_interface


def get_pyx_def(src_dir):
    cb_interface = get_c_def(src_dir)
    for proc in cb_interface:
        decls_map = {var: (code, shape) for var, code, shape in cb_interface[proc]['decls']}
        pyx_in_args, pyx_out_args = get_pyx_arg(cb_interface[proc]['args'], decls_map)
        cb_interface[proc]['pyx_in_args'] = pyx_in_args
        cb_interface[proc]['pyx_out_args'] = pyx_out_args
    return cb_interface

