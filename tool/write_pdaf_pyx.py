import re
from pathlib import Path
import get_decls
import get_decls_cb
from docstring import docstrings


WRAP_WIDTH = 80

TYPE_MAP = {
  'integer':      'int ',
  'real':         'double ',
  'logical':      'bint ',
  'character':    'char '
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
    lines.append("cimport numpy as cnp")
    lines.append("from . cimport pdaf_c_cb_interface as pdaf_cb")
    lines.append("from .cfi_binding cimport CFI_cdesc_t, CFI_address, CFI_index_t, CFI_establish")
    lines.append("from .cfi_binding cimport CFI_attribute_other, CFI_type_double, CFI_type_int")
    lines.append("from .cfi_binding cimport CFI_cdesc_rank1, CFI_cdesc_rank2, CFI_cdesc_rank3")
    lines.append("")
    lines.append("try:")
    lines.append("    import mpi4py")
    lines.append("    mpi4py.rc.initialize = False")
    lines.append("except ImportError:")
    lines.append("    pass")
    lines.append("")
    lines.append("# Global error handler")
    lines.append("def global_except_hook(exctype, value, traceback):")
    lines.append("    from traceback import print_exception")
    lines.append("    try:")
    lines.append("        import mpi4py.MPI")
    lines.append("")
    lines.append("        if mpi4py.MPI.Is_initialized():")
    lines.append("            try:")
    lines.append("                sys.stderr.write('Uncaught exception was '"
                 "'detected on rank {}.\\n'.format(")
    lines.append("                    mpi4py.MPI.COMM_WORLD.Get_rank()))")

    lines.append("                print_exception(exctype, value, traceback)")
    lines.append("                sys.stderr.write(\"\\n\")")
    lines.append("                sys.stderr.flush()")
    lines.append("            finally:")
    lines.append("                try:")
    lines.append("                    mpi4py.MPI.COMM_WORLD.Abort(1)")
    lines.append("                except Exception as e:")
    lines.append("                    sys.stderr.write('MPI Abort failed, this process will hang.\\n')")
    lines.append("                    sys.stderr.flush()")
    lines.append("                    raise e")
    lines.append("        else:")
    lines.append("            sys.__excepthook__(exctype, value, traceback)")
    lines.append("    except ImportError:")
    lines.append("        sys.__excepthook__(exctype, value, traceback)")
    lines.append("")
    lines.append("sys.excepthook = global_except_hook")

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
                # currently only pdafomi_diag has intent(inout) pointers
                # we do not take input pointers
                if 'pointer' not in code.lower():
                    pyx_in_arg_list.append(arg)
            match = re.search(r'intent\(out\)', code, re.IGNORECASE)
            if match and 'dimension(:' in code.lower() and 'pointer' not in code.lower():
                # if it's an assumed shape output array, we need to handle it differently
                # it will need an input array to establish a CFI
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


def write_pyx_signature(pyx_name, arg_list, decl_map):
    pyx_def_arg_list = []
    for arg in arg_list:
        code, shape = decl_map.get(arg, ("", ""))
        if 'dimension' in code.lower():
            n_colons = len(shape.split(',')) - 1
            s_colons_t = ':,' * n_colons
            s_colons = "::1" if n_colons == 0 else \
                "::1," + s_colons_t[:-1]  # remove last comma
            base_type = re.match(r'(integer|real|logical|character)', code).group(1).lower()
            pyx_def_arg_list.append(f'{TYPE_MAP[base_type]}[{s_colons}] {arg}')
        elif 'procedure' in code.lower():
            pyx_def_arg_list.append(arg)
        else:
            base_type = re.match(r'(integer|real|logical|character)', code).group(1).lower()
            pyx_def_arg_list.append(f'{TYPE_MAP[base_type]} {arg}')

    return wrap_comma_list(
                prefix=f"def {pyx_name}(",
                items=pyx_def_arg_list,
                suffix=":"
            )


def get_docstring_type(code, shape):
    if 'procedure' in code.lower():
        return 'Callable'
    elif 'dimension' in code.lower():
        base_type = re.match(r'(integer|real|logical|character)', code).group(1).lower()
        return f'ndarray[{TYPE_MAP_ARR[base_type]}, ndim={len(shape.split(","))}]'
    else:
        base_type = re.match(r'(integer|real|logical|character)', code).group(1).lower()
        return TYPE_MAP[base_type]


def write_docstring_args(arg_in_list, arg_out_list, cb_interface, decl_map, comments,
                         callback="", indent="    "):
    lines = []
    if len(arg_in_list) > 0:
        lines.append('')
        header = f'{callback} Parameters'.strip()
        lines.append(indent + header)
        lines.append(indent+'-'*len(header))
        for arg in arg_in_list:
            code, shape = decl_map.get(arg, ("", ""))
            arg_type = get_docstring_type(code, shape)
            lines.append(indent + f'{arg} : {arg_type}')
            lines.append(2*indent + f'\n{2*indent}'.join(comments[arg]))
            if 'dimension' in code.lower():
                lines.append(2*indent + f'Array shape: ({shape})')
            if arg_type == 'Callable':
                decl_cb_map = {arg_cb: (code, shape)
                               for arg_cb, code, shape in cb_interface[arg]['decls']}
                com_lines = write_docstring_args(cb_interface[arg]['pyx_in_args'],
                                                 cb_interface[arg]['pyx_out_args'],
                                                 cb_interface,
                                                 decl_cb_map,
                                                 cb_interface[arg]['comments'],
                                                 indent=2*indent,
                                                 callback="Callback")
                lines.extend(com_lines)
                lines.append('')

    if len(arg_in_list) > 0:
        lines.append('')
        header = f'{callback} Returns'.strip()
        lines.append(indent + header)
        lines.append(indent+'-'*len(header))
        for arg in arg_out_list:
            code, shape = decl_map.get(arg, ("", ""))
            arg_type = get_docstring_type(code, shape)
            lines.append(indent + f'{arg} : {arg_type}')
            lines.append(2*indent + f'\n{2*indent}'.join(comments[arg]))
            if 'dimension' in code.lower():
                lines.append(2*indent + f'Array shape: ({shape})')

    return lines


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
            s = f'cdef CFI_cdesc_rank{n_dim} {arg}_cfi'
            lines.append(indent + s)
            s = f'cdef CFI_cdesc_t *{arg}_ptr = <CFI_cdesc_t *> &{arg}_cfi'
            lines.append(indent + s)
            if 'pointer' in code.lower():
                continue

            s = f'cdef size_t {arg}_nbytes = {arg}.nbytes'
            lines.append(indent + s)
            s = f'cdef CFI_index_t {arg}_extent[{n_dim}]'
            lines.append(indent + s)
            for i in range(n_dim):
                s = f'{arg}_extent[{i}] = {arg}.shape[{i}]'
                lines.append(indent + s)

        n_colons = len(shape.split(',')) - 1
        s_colons_t = ':,' * n_colons
        s_colons = "::1" if n_colons == 0 else \
            "::1," + s_colons_t[:-1]  # remove last comma
        base_type = re.match(r'(integer|real|logical|character)', code).group(1).lower()
        s = f'cdef cnp.ndarray[{TYPE_MAP_CNP[base_type]}, ndim={n_colons+1}, mode="fortran", negative_indices=False, cast=False] '
        if 'intent(out)' in code.lower() and 'dimension(:' not in code.lower():
            s += f'{arg}_np = np.zeros(({shape}), dtype={TYPE_MAP_ARR[base_type]}, order="F")'
            lines.append(indent + s)
            s = f'cdef {TYPE_MAP[base_type]}[{s_colons}] {arg} = {arg}_np'
            lines.append(indent + s)
        else:
            s += f'{arg}_np = np.asarray({arg}, dtype={TYPE_MAP_ARR[base_type]}, order="F")'
            lines.append(indent + s)
    return lines


def write_input_cfi(arg_list, decl_map, indent="    "):
    """Convert 2D arrays to fortran contiguous arrays
    """
    # convert np arrays to memoryview
    lines = []
    for arg in arg_list:
        code, shape = decl_map.get(arg, ("", ""))
        if 'dimension' not in code.lower():
            # not a Fortran array, skip
            continue

        if 'intent(inout)' in code.lower():
            continue
        if 'intent(out)' in code.lower():
            continue

        if 'dimension(:' in code.lower():
            n_dim = len(shape.split(','))
            s = f'cdef CFI_cdesc_rank{n_dim} {arg}_cfi'
            lines.append(indent + s)
            s = f'cdef CFI_cdesc_t *{arg}_ptr = <CFI_cdesc_t *> &{arg}_cfi'
            lines.append(indent + s)
            if 'pointer' in code.lower():
                continue
            s = f'cdef size_t {arg}_nbytes = {arg}.nbytes'
            lines.append(indent + s)
            s = f'cdef CFI_index_t {arg}_extent[{n_dim}]'
            lines.append(indent + s)
            for i in range(n_dim):
                s = f'{arg}_extent[{i}] = {arg}.shape[{i}]'
                lines.append(indent + s)

    return lines


def write_user_cython(arg_list, decl_map, indent='    '):
    """convert python function to Cython function"""
    lines = []
    for arg in arg_list:
        code, _ = decl_map.get(arg, ("", ""))
        if 'procedure' in code.lower():
            lines.append(indent+ f'pdaf_cb.{arg.replace("py__", "")} = <void*>{arg}')
    return lines


def write_return_def(arg_list, decl_map, indent='    '):
    lines = []
    for arg in arg_list:
        code, _ = decl_map.get(arg, ("", ""))
        if 'intent(inout)' in code.lower():
            # only define pure return arguments
            continue
        if 'dimension' in code.lower():
            # array and pointers are defined in memory view
            continue
        base_type = re.match(r'(integer|real|logical|character)', code).group(1).lower()
        s = f'cdef {TYPE_MAP[base_type]} {arg}'
        lines.append(indent + s)
    return lines


def write_func_call(name, arg_list, decl_map, indent="    "):
    """write the C function call in Cython"""
    # special treatment for init subroutine
    lines = [indent + 'with nogil:']
    # establish CFI for assumed shape arrays
    for arg in arg_list:
        code, shape = decl_map.get(arg, ("", ""))
        if 'dimension' not in code.lower():
            # not a Fortran array, skip
            continue
        if 'pointer' in code.lower():
            # pointer arrays do not need CFI_establish with CFI
            continue

        if 'dimension(:' in code.lower():
            n_dim = len(shape.split(','))
            s_zeros = '0,' * n_dim
            s_zeros = s_zeros[:-1]
            base_type = re.match(r'(integer|real|logical|character)', code).group(1).lower()
            s = f'CFI_establish({arg}_ptr, &{arg}[{s_zeros}], CFI_attribute_other,'
            lines.append(2*indent + s)
            s = f'CFI_type_{TYPE_MAP[base_type]}, {arg}_nbytes, {n_dim}, {arg}_extent)'
            lines.append(2*indent + len('CFI_establish(')*' ' + s)
            lines.append("")

    # get argument list of the subroutine
    c_args = []
    for arg in arg_list:
        code, shape = decl_map.get(arg, ("", ""))
        if 'dimension(:' in code.lower():
            c_args.append(f'{arg}_ptr')
        elif 'dimension' in code.lower():
            n_zeros = len(shape.split(','))
            s_zeros = '0,' * n_zeros
            s_zeros = s_zeros[:-1]
            c_args.append(f'&{arg}[{s_zeros}]')
        elif 'procedure' in code.lower():
            proc_type = re.search(r'procedure\((.*?)\)', code).group(1).lower()
            c_args.append(f'pdaf_cb.{proc_type.replace("py__", "c__")}')
        else:
            c_args.append(f'&{arg}')
    prefix = 2*indent + f'{name.lower()}('
    lines_call = wrap_comma_list(prefix=prefix,
                                 items=c_args,
                                 subsequent_indent=len(prefix)*' ')

    lines.extend(lines_call)
    lines.append("")

    return lines


def write_return(arg_list, decl_map, indent="    "):
    if len(arg_list) == 0:
        return []
    lines = []

    for arg in arg_list:
        code, shape = decl_map.get(arg, ("", ""))
        if 'pointer' in code.lower():
            n_dim = len(shape.split(','))
            s = f'cdef CFI_index_t {arg}_subscripts[{n_dim}]'
            lines.append(indent + s)
            for i in range(n_dim):
                s = f'{arg}_subscripts[{i}] = 0'
                lines.append(indent + s)
            base_type = re.match(r'(integer|real|logical|character)', code).group(1).lower()
            s = f'cdef {TYPE_MAP[base_type]}* {arg}_ptr_np'
            lines.append(indent + s)
            s = f'{arg}_ptr_np = <{TYPE_MAP[base_type]}*>CFI_address({arg}_ptr, {arg}_subscripts)'
            lines.append(indent + s)
            s = f'cdef cnp.ndarray[{TYPE_MAP_CNP[base_type]}, ndim={n_dim}, mode="fortran", negative_indices=False, cast=False] '
            s_colons = f':{arg}_ptr.dim[0].extent:1'
            s_colons += ',' if n_dim > 1 else ''
            s_colons += ','.join([f':{arg}_ptr.dim[{i}].extent' for i in range(1, n_dim)])
            s += f'{arg}_np = np.asarray(<{TYPE_MAP[base_type]}[{s_colons}]> {arg}_ptr_np, order="F")'
            lines.append(indent + s)

    s = 'return '
    for arg in arg_list:
        code, shape = decl_map.get(arg, ("", ""))
        s += f'{arg}_np, ' if 'dimension' in code.lower() else f'{arg}, '

    lines.append(indent + s[:-2])
    return lines


def generate_pyx(name, arg_list, decls, comments, cb_interface, is_internal):
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
    c_arg_list = []
    for arg in arg_list:
        code, _ = decl_map.get(arg, ("", ""))
        match = re.search(r'procedure\((.*?)\)', code)
        if match:
            proc_type = match.group(1).lower()
            decl_map[proc_type.replace("c__", "py__")] = decl_map.pop(arg)
            comments[proc_type.replace("c__", "py__")] = comments.pop(arg)
            cb_interface[proc_type.replace("c__", "py__")] = cb_interface.pop(proc_type)
            c_arg_list.append(proc_type.replace('c__','py__'))
        else:
            c_arg_list.append(arg)

    # write the function signature
    pyx_name = re.sub(r'c__pdaf_', '_' if is_internal else '', name.lower())
    pyx_name = re.sub(r'c__pdafomi_', '_' if is_internal else '', pyx_name.lower())
    pyx_name = re.sub(r'c__pdaflocal_', '_' if is_internal else '', pyx_name.lower())
    pyx_name = re.sub(r'c__pdaflocalomi_', '_' if is_internal else '', pyx_name.lower())
    pyx_name = re.sub(r'c__pdaf3_', '_' if is_internal else '', pyx_name.lower())
    pyx_name = re.sub(r'c__pdaf', '_' if is_internal else '', pyx_name.lower())
    lines_sig = write_pyx_signature(pyx_name, pyx_in_arg_list, decl_map)
    lines.extend(lines_sig)

    # write the docstrings
    lines.append('    """' + docstrings.docstrings.get(pyx_name, DUMMY_DOC))
    lines_args = write_docstring_args(pyx_in_arg_list, pyx_out_arg_list,
                                      cb_interface, decl_map, comments)
    lines.extend(lines_args)
    lines.append('    """')

    # write the function body
    # memory views
    lines_mv = write_memory_view(pyx_out_arg_list, decl_map)
    lines.extend(lines_mv)

    lines_cfi = write_input_cfi(pyx_in_arg_list, decl_map)
    lines.extend(lines_cfi)

    lines_cb = write_user_cython(pyx_in_arg_list, decl_map)
    lines.extend(lines_cb)

    lines_return_def = write_return_def(pyx_out_arg_list, decl_map)
    lines.extend(lines_return_def)

    lines_call = write_func_call(name, c_arg_list, decl_map)
    lines.extend(lines_call)

    lines_return = write_return(pyx_out_arg_list, decl_map)
    lines.extend(lines_return)

    return "\n".join(lines)


def process_file(filepath, dst_dir):

    with open(filepath, 'r') as f:
        lines = f.readlines()

    blocks = get_decls.extract_subroutine_blocks(lines)

    if not blocks:
        print (f"No subroutines found in {filepath}")
        return  # No subroutines, skip
    lines = write_header()

    filename = filepath.stem + '.pyx'
    output_path = Path(dst_dir) / filename
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'a') as f:
        f.write( "\n".join(lines))
        f.write( "\n\n")

    for block in blocks:
        # try:
        cb_interface = get_decls_cb.get_pyx_def(filepath.parent)
        name, arg_list, start_idx = get_decls.extract_multiline_signature(block)
        decls, comments = get_decls.extract_declarations_by_args(block, start_idx, arg_list)
        missing, extra = get_decls.compare_signature_and_declarations(arg_list, decls)
        if missing:
            print(f"!!! Missing declarations for: {missing} in {name} in {filepath}")
        if extra:
            print(f"!!! Declared but not in signature: {extra} in {name} in {filepath}")

        extern_decls = generate_pyx(name, arg_list, decls, comments, cb_interface, is_internal='internal' in filepath.stem)

        with open(output_path, 'a') as f:
            f.write(extern_decls + "\n\n\n")

        # except Exception as e:
        #     print(f"Skipping block in {filepath}: {e}")


def process_directory(src_dir, dst_dir):
    for path in Path(src_dir).rglob("*"):
        if path.suffix.lower() in (".f90"):
            if path.stem == "pdaf_c_cb_interface":
                continue
            process_file(path, dst_dir)

if __name__ == "__main__":
    # process_directory('/home/users/ia923171/pyPDAF_dev/pyPDAF/src/fortran', '')
    process_directory('/home/users/ia923171/pyPDAF_dev/pyPDAF/src/fortran', '/home/users/ia923171/pyPDAF_dev/pyPDAF/src/pyPDAF')
