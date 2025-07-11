import math
import re
from pathlib import Path

import numpy as np
import get_decls
import get_decls_cb
from docstring import docstrings


WRAP_WIDTH = 80

TYPE_MAP = {
  'integer':      'int',
  'real':         'double',
  'logical':      'bool',
  'character':    'str'
}

TYPE_MAP_ARR = {
  'integer':      'np.intc',
  'real':         'np.float64',
  'logical':      'np.bool',
  'character':    'np.str_'
}

TYPE_MAP_CNP = {
  'integer':      'cnp.int32_t',
  'real':         'cnp.float64_t',
}



DUMMY_DOC = "Checking the corresponding PDAF documentation in https://pdaf.awi.de\n" \
            "    For internal subroutines checking corresponding PDAF comments."


def write_header():
    lines = []
    lines.append("import sys")
    lines.append("import numpy as np")

    return lines


def get_pyx_arg(arg_list, decl_map):
    pyx_in_arg_list = []
    pyx_out_arg_list = []
    for arg in arg_list:
        code, _ = decl_map.get(arg, ("", ""))
        match = re.search(r'procedure\((.*?)\)', code)
        if match:
            proc_type = match.group(1).lower()
            pyx_in_arg_list.append(proc_type.replace("c__", "py__"))
        else:
            match = re.search(r'intent\(in\)|intent\(inout\)', code, re.IGNORECASE)
            if match:
                pyx_in_arg_list.append(arg)
            match = re.search(r'intent\(out\)|intent\(inout\)', code, re.IGNORECASE)
            if match:
                pyx_out_arg_list.append(arg)
    return pyx_in_arg_list, pyx_out_arg_list


def wrap_comma_list(prefix, items, subsequent_indent="    ", suffix=""):
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
        if len(current) + len(addition) + (4 if sep==", " else len(suffix)) > WRAP_WIDTH:
            lines.append(current)
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


def write_c_signature(pyx_name, arg_list, decl_map):
    pyx_def_arg_list = []
    for arg in arg_list:
        code, _ = decl_map.get(arg, ("", ""))
        if 'dimension(:' in code.lower():
            pyx_def_arg_list.append(f'CFI_cdesc_t* {arg}')
        else:
            base_type = re.match(r'(integer|real|logical|character)', code).group(1).lower()
            pyx_def_arg_list.append(f'{TYPE_MAP[base_type]}* {arg}')

    return wrap_comma_list(
                prefix=f"cdef void {pyx_name}(",
                items=pyx_def_arg_list,
                suffix=" noexcept with gil:"
            )


def write_memory_view(arg_list, decl_map, indent="    "):
    """Convert 2D arrays to fortran contiguous arrays
    """
    # convert np arrays to memoryview
    lines = []
    for arg in arg_list:
        code, shape = decl_map.get(arg, ("", ""))
        if 'dimension' not in code.lower():
            # not a Fortran array, skip
            continue

        if 'dimension(:' in code.lower():
            n_dim = len(shape.split(','))
            s = f'cdef CFI_index_t {arg}_subscripts[{n_dim}]'
            s = f'cdef size_t {arg}_dim[{n_dim}]'
            lines.append(indent + s)
            for i in range(n_dim):
                s = f'{arg}_subscripts[{i}] = 0'
                lines.append(indent + s)
                s = f'{arg}_dim[{i}] = {arg}.dim[{i}].extent'
                lines.append(indent + s)

            s = f'cdef double *{arg}_ptr = <{TYPE_MAP[base_type]} *>CFI_address({arg}, {arg}_subscripts)'
            lines.append(indent + s)
            s = f'cdef {TYPE_MAP[base_type]}[::1,:] {arg}_np = '
            n_colons = len(shape.split(','))
            s_colons = [f':{arg}_dim[{i}]' for i in range(1, n_colons)]
            s_colons = f':{arg}_dim[0]:1,' if n_colons == 1 else \
                ",".join(s_colons)  # remove last comma
            s += f'np.asarray(<{TYPE_MAP[base_type]}[{s_colons}]> {arg}_ptr, order="F")'
            lines.append(indent + s)

        else:
            colons = shape.split(',')
            n_colons = len(colons)
            s_colons = ':,' * (n_colons - 1)
            s_colons = "::1" if n_colons == 1 else \
                "::1," + s_colons[:-1]  # remove last comma
            base_type = re.match(r'(integer|real|logical|character)', code).group(1).lower()
            s = f'cdef {TYPE_MAP[base_type]}[{s_colons}] '
            s_colons = [f':{colons[i]}[0]' for i in range(1, n_colons)]
            s_colons = f':{colons[0]}[0]:1,' if n_colons == 1 else \
                ",".join(s_colons)  # remove last comma
            s += f'{arg}_np = np.asarray(<{TYPE_MAP[base_type]}[{s_colons}]> {arg}, order="F")'
            lines.append(indent + s)
    return lines


def write_func_call(name, arg_list, decl_map, indent="    "):
    """write the C function call in Cython"""
    # get argument list of the subroutine
    in_args = []
    for arg in arg_list:
        code, shape = decl_map.get(arg, ("", ""))
        if 'dimension' in code.lower():
            in_args.append(f'{arg}_np.base')
        else:
            in_args.append(f'{arg}[0]')

    out_args = []
    for arg in arg_list:
        code, _ = decl_map.get(arg, ("", ""))
        match = re.search(r'intent\(out\)|intent\(inout\)', code, re.IGNORECASE)
        if match:
            s = f'{arg}_np' if 'dimension' in code.lower() else f'{arg}'
            out_args.append(s)
    prefix = indent + ','.join(out_args) + f' = (<object>{name[3:].lower()})('
    lines = wrap_comma_list(prefix=prefix,
                            items=in_args,
                            subsequent_indent=len(prefix)*' ')

    return lines


def write_assert(arg_list, decl_map, indent="    "):
    lines = []
    for arg in arg_list:
        code, shape = decl_map.get(arg, ("", ""))
        if 'dimension' not in code.lower():
            continue

        s = f'cdef double[::1,:] {arg}_new'
        lines.append(indent + s)
        n_colons = len(shape.split(','))
        s_colons = '0,' * n_colons
        s_colons = s_colons[:-1]  # remove last comma
        s = f'if {arg} != &{arg}_np[{s_colons}]:'
        lines.append(indent + s)
        colons = shape.split(',')
        n_colons = len(colons)
        s_colons = [f':{colons[i]}[0]' for i in range(1, n_colons)]
        s_colons = f':{colons[0]}[0]:1,' if n_colons == 1 else \
            ",".join(s_colons)  # remove last comma
        if 'dimension(:' in code.lower():
            s_colons = [f':{arg}_dim[{i}]' for i in range(1, n_colons)]
            s_colons = f':{arg}_dim[0]:1,' if n_colons == 1 else \
                ",".join(s_colons)  # remove last comma
        s = f'{arg}_new = np.asarray(<double[{s_colons}]> {arg}, order="F")'
        lines.append(2*indent + s)
        s = f'{arg}_new[...] = {arg}_np'
        lines.append(2*indent + s)
        s = f'warnings.warn("The memory address of {arg} is changed in c__add_obs_err_pdaf."'
        lines.append(2*indent + s)
        s = '"The values are copied to the original Fortran array, and can slow-down the system.", RuntimeWarning)'
        lines.append(3*indent + s)

    return lines


def generate_pyx(name, arg_list, decls, comments):
    """
    name      : original subroutine name (string)
    arg_list  : list of argument names in order, already lowercased
    decls     : list of (var, decl_code) from extract_declarations_by_args
    comments  : dict var->comment_str
    """

    lines = []
    # replace procedure argument names with their procedure type names
    decl_map = { var: (code, shape) for var, code, shape in decls }
    pyx_in_arg_list, pyx_out_arg_list = get_pyx_arg(arg_list, decl_map)
    for arg in arg_list:
        code, _ = decl_map.get(arg, ("", ""))
        match = re.search(r'procedure\((.*?)\)', code)
        if match:
            proc_type = match.group(1).lower()
            decl_map[proc_type] = decl_map.pop(arg)
            comments[proc_type] = comments.pop(arg)

    # write the function signature
    lines_sig = write_c_signature(name, arg_list, decl_map)
    lines.extend(lines_sig)

    # write the function body
    # memory views
    lines_mv = write_memory_view(arg_list, decl_map)
    lines.extend(lines_mv)

    lines.append("")

    lines_call = write_func_call(name, arg_list, decl_map)
    lines.extend(lines_call)

    lines.append("")

    lines_assert = write_assert(pyx_out_arg_list, decl_map)
    lines.extend(lines_assert)

    return "\n".join(lines)


def process_file(filepath, dst_dir):
    lines = write_header()

    filename = filepath.stem + '.pyx'
    output_path = Path(dst_dir) / filename
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'a') as f:
        f.write( "\n".join(lines))
        f.write( "\n\n")

    cb_interface = get_decls_cb.get_c_def(filepath.parent)
    for name in cb_interface:
        arg_list = cb_interface[name]['args']
        decls = cb_interface[name]['decls']
        comments = cb_interface[name]['comments']

        extern_decls = generate_pyx(name, arg_list, decls, comments)

        with open(output_path, 'a') as f:
            f.write(extern_decls + "\n\n\n")

        # except Exception as e:
        #     print(f"Skipping block in {filepath}: {e}")


def process_directory(src_dir, dst_dir):
    for path in Path(src_dir).rglob("*"):
        if path.suffix.lower() in (".f90"):
            if path.stem == "pdaf_c_cb_interface":
                process_file(path, dst_dir)


if __name__ == "__main__":
    # process_directory('/home/users/ia923171/pyPDAF_dev/pyPDAF/src/fortran', '')
    process_directory('/home/users/ia923171/pyPDAF_dev/pyPDAF/src/fortran', '/home/users/ia923171/pyPDAF_dev/pyPDAF/src/pyPDAF')