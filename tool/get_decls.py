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

def extract_dimension_shape(decl_code: str) -> str:
    # Find the DIMENSION keyword
    match = re.search(r'dimension\s*\(', decl_code, re.IGNORECASE)
    if not match:
        return ""

    start = match.end()  # position after '('
    depth = 1
    i = start
    while i < len(decl_code):
        if decl_code[i] == '(':
            depth += 1
        elif decl_code[i] == ')':
            depth -= 1
            if depth == 0:
                return decl_code[start:i].strip()
        i += 1
    return ""  # if unmatched

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
        # m = re.search(r'dimension\s*\(([^)]+)\)', decl_code, re.IGNORECASE)
        # shape = m.group(1) if m else ""
        shape = extract_dimension_shape(decl_code) 
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

    for line in block[start_idx:]:
        # stop at end of subroutine
        if re.match(r'\s*end\s+subroutine', line, re.I):
            break

        if 'implicit none' in line:
            continue


        raw = line.rstrip("\n")
        # 1) split off any comment
        if "!" in raw:
            code_part, comment_part = raw.split("!", 1)
        else:
            code_part, comment_part = raw, None

        if 'use ' in code_part.lower():
            continue

        # --- Start a new declaration if we see '::' and we're not in one ---
        if comment_part and buffer_comments == []:
            buffer_comments = [comment_part.strip()]
            comment_part = ""

        if not buffer and "::" in code_part:
            buffer = code_part.strip()
            # single-line decl: process immediately
            _process_buffer(buffer, buffer_comments, arg_names, decls, comments)
            buffer = ""
            buffer_comments = []
            continue

        if not buffer and code_part.strip() == "":
            if comment_part:
                buffer_comments.append(comment_part.strip())
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
