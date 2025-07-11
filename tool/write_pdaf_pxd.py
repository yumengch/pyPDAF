import re
from pathlib import Path
import get_decls
import get_decls_cb

WRAP_WIDTH = 80

TYPE_MAP = {
  'integer':      'int* ',
  'real':         'double* ',
  'logical':      'bint* ',
  'character':    'char* '
}



def wrap_comma_list(prefix, arg_list, decl_map, cb_interface, subsequent_indent="    ", suffix=""):
    """
    Wrap a comma-separated list under WRAP_WIDTH, inserting '&'
    at end of every line except the last, and appending `suffix`
    (e.g. ' bind(c)') to that last line.
    """
    lines = []
    current = prefix
    for i, arg in enumerate(arg_list):
        sep = ", " if i < len(arg_list) - 1 else ")"
        code, _ = decl_map.get(arg, ("", ""))
        # check user-supplied procedures
        match = re.search(r'procedure\((.*?)\)', code)
        if match:
            if current != subsequent_indent:
                lines.append(current)
            current = subsequent_indent
            addition = ""

            proc_type = match.group(1)
            assert proc_type in cb_interface, f'{proc_type} is not in cb_interface'
            # if this is a C binding, we need to use the C interface name
            decl_map_cb = { var: (code, shape) for var, code, shape in cb_interface[proc_type]['decls'] }
            cb_prefix = subsequent_indent + f'void (*{proc_type.lower()})('
            addition_list = wrap_comma_list(
                prefix=cb_prefix,
                arg_list=cb_interface[proc_type]['args'],
                decl_map=decl_map_cb,
                cb_interface=cb_interface,
                subsequent_indent=" "*len(cb_prefix),
                suffix=sep
            )
            lines.extend(addition_list)
        elif 'pointer' in code or 'dimension(:' in code.lower():
            addition = f'CFI_cdesc_t* ' + arg.lower() + sep
        else:
            base_type = re.match(r'(integer|real|logical|character)', code).group(1).lower()
            c_type = TYPE_MAP[base_type]
            addition = c_type + arg.lower() + sep

        # if adding this would exceed limit (allowing space for ' &' on wrapped lines)
        if len(current) + len(addition) + (4 if sep==", " else len(suffix)) > WRAP_WIDTH:
            lines.append(current)
            current = subsequent_indent + addition
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


def generate_extern_pxd(name, arg_list, decls, cb_interface):
    """
    name      : original subroutine name (string)
    arg_list  : list of argument names in order, already lowercased
    decls     : list of (var, decl_code) from extract_declarations_by_args
    comments  : dict var->comment_str
    """
    decl_map = { var: (code, shape) for var, code, shape in decls }

    lines = wrap_comma_list(
                prefix=f"cdef extern void {name.lower()}(",
                arg_list=arg_list,
                decl_map=decl_map,
                cb_interface=cb_interface,
                suffix=" noexcept nogil;"
            )
    lines.append("")
    return "\n".join(lines)

def process_file(filepath, dst_dir):

    cb_interface = get_decls_cb.get_c_def(filepath.parent)

    with open(filepath, 'r') as f:
        lines = f.readlines()

    blocks = get_decls.extract_subroutine_blocks(lines)

    if not blocks:
        print (f"No subroutines found in {filepath}")
        return  # No subroutines, skip

    for block in blocks:
        name, arg_list, start_idx = get_decls.extract_multiline_signature(block)
        decls, _ = get_decls.extract_declarations_by_args(block, start_idx, arg_list)
        missing, extra = get_decls.compare_signature_and_declarations(arg_list, decls)
        if missing:
            print(f"!!! Missing declarations for: {missing} in {name} in {filepath}")
        if extra:
            print(f"!!! Declared but not in signature: {extra} in {name} in {filepath}")

        extern_decls = generate_extern_pxd(name, arg_list, decls, cb_interface)

        filename = filepath.stem + '.pxd'
        output_path = Path(dst_dir) / filename
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'a') as f:
            f.write(extern_decls + "\n")

def process_directory(src_dir, dst_dir):
    for path in Path(src_dir).rglob("*"):
        if path.suffix.lower() in (".f90"):
            if path.stem == "pdaf_c_cb_interface":
                continue
            process_file(path, dst_dir)

if __name__ == "__main__":
    # process_directory('/home/users/ia923171/pyPDAF_dev/pyPDAF/src/fortran', '')
    process_directory('/home/users/ia923171/pyPDAF_dev/pyPDAF/src/fortran', '/home/users/ia923171/pyPDAF_dev/pyPDAF/src/pyPDAF')