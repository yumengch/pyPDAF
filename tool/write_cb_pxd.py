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

def wrap_comma_list(prefix, arg_list, decl_map, subsequent_indent="    ", suffix=""):
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
        if 'pointer' in code or 'dimension(:' in code.lower():
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


def generate_extern_pxd(name, arg_list, decls):
    """
    name      : original subroutine name (string)
    arg_list  : list of argument names in order, already lowercased
    decls     : list of (var, decl_code) from extract_declarations_by_args
    comments  : dict var->comment_str
    """
    decl_map = { var: (code, shape) for var, code, shape in decls }
    lines = wrap_comma_list(
                prefix=f"cdef void {name.lower()}(",
                arg_list=arg_list,
                decl_map=decl_map,
                suffix=" noexcept nogil;"
            )
    lines.append(f"cdef void* {name[3:].lower()} = NULL;")
    lines.append("")
    return "\n".join(lines)

def process_file(src_dir, dst_dir):

    cb_interface = get_decls_cb.get_c_def(src_dir)

    for name in cb_interface:
        cb_decls = generate_extern_pxd(name, arg_list=cb_interface[name]['args'],
                                       decls=cb_interface[name]['decls'])
        filename = 'pdaf_c_cb_interface.pxd'
        output_path = Path(dst_dir) / filename
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'a') as f:
            f.write(cb_decls + "\n")

        # except Exception as e:
        #     print(f"Skipping block in {filepath}: {e}")

if __name__ == "__main__":
    # process_directory('/home/users/ia923171/pyPDAF_dev/pyPDAF/src/fortran', '')
    process_file('/home/users/ia923171/pyPDAF_dev/pyPDAF/src/fortran', '/home/users/ia923171/pyPDAF_dev/pyPDAF/src/pyPDAF')
