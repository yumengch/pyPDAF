"""This is a simple script for converting C binding
fortran routines in pyPDAF/fortran to Cython definition and
implementation files.
"""
import typing
import re

# todo: type for pointer should be matching its corresponding fortran pointer type
# here double* is used because PDAF only uses double
conv = {'integer' : 'int', 'logical': 'bint', 'real': 'double', 'procedure':'void', 'character':'CFI_cdesc_t', 'type':'double*'}
pyconv = {'integer' : 'int', 'logical': 'bool', 'real': 'float', 'procedure':'void', 'character':'CFI_cdesc_t', 'type':'double*'}
npyconv = {'integer' : 'np.intc', 'real': 'np.float64'}

def extract_dimension_name(s:str) -> str | None:
    """extract dimension names from dimension
    """
    # Match one or more word characters (letters, digits, and underscores), but exclude digits
    match = re.search(r'[a-zA-Z_]+', s)
    if match:
        return match.group()
    return None


def get_pyx_arg_list(subroutine_name:str, arg_info: dict[str, dict[str, str|bool|None|list[str]]]) -> tuple[list[str], list[str]]:
    """Retrive dimension of input array arguments of a subroutine.
    If the dimension is an input argument of the subroutine, we do
    not write it explicitly as they can be obtained from
    numpy array shape. However, in the callback functions,
    these functions will be called by Fortran, all arguments must
    stay in the argument list.

    Parameters
    ----------
    subroutine_name : str
        The name of the subroutine
    arg_info : dict
        The dictionary of argument information

    Returns
    -------
    arg_list : list
        The list of arguments to be written in the Cython file
    """
    ArrayDims: list[str] = []
    for _, info in arg_info.items():
        # we need to keep the dimension of output arrays
        # when the dimension is only used by output arrays
        if info['intent'] == 'out': continue
        if info['array']:
            assert type(info['dimension']) == list, f"dimension in info {info['dimension']} is not a list"
            for dim in info['dimension']:
                dimsize = extract_dimension_name(dim)
                if dimsize in arg_info:
                    if dimsize not in ArrayDims:
                        ArrayDims.append(dimsize)

    arg_list : list[str] = []
    for arg_name, info in arg_info.items():
        # we collect all arguments that is not input array dimensions
        # however, we collect all arguments whatsoever for callback functions
        if arg_name not in ArrayDims or ('_cb' in subroutine_name):
            arg_list.append(arg_name)

    return arg_list, ArrayDims


def write_func_def(f:typing.TextIO, subroutine_name:str,
                   arg_info:dict[str, dict[str, str|bool|None|list[str]]],
                   arg_list: list[str]) -> None:
    """write def function_name(arg1, arg2, ...) in Cython

    Parameters
    ----------
    f : file
        The file to write the function definition
    subroutine_name : str
        The name of the subroutine
    arg_info : dict
        The dictionary of argument information
    arg_list : list
        The list of arguments to be written in the Cython file
    """
    # define function
    if subroutine_name[3:11] == 'pdafomi_':
        funcname = subroutine_name[7:]
    elif subroutine_name[3:8] == 'pdaf_':
        funcname = subroutine_name[8:]
    else:
        funcname = subroutine_name[3:]

    s = f'def {funcname} ('
    indent = ' '*len(s)
    count = 0
    # only arguments with input values
    for argname in arg_list:
        info = arg_info[argname]
        assert  type(info['type']) is str, f"type in info {info['type']} is not a str"

        # we do not write purely return arguments
        # callback functions are exceptions
        if type(info['intent']) is str:
            if ('in' not in info['intent']) and ('_cb' not in funcname):
                continue

        # write procedure
        if info['type'] == 'procedure':
            assert  type(info['kind']) is str, f"kind in info {info['kind']} is not a str"
            s += f'py{info["kind"][1:]}'

        elif info['array']:
            assert type(info['dimension']) is list, f"dimension in info {info['dimension']} is not a list"
            # array are passed as memory views
            s += conv[info['type']]+'[:'
            # all callback functiosn should be fortran contiguous
            if '_cb' in subroutine_name or len(info['dimension']) == 1:
                s += ':1'
            for _ in range(len(info['dimension']) - 1):
                s += ',:'
            s += f'] {argname}'
        else:
            s += f'{conv[info["type"]]} {argname}'

        s += f',\n{indent}'
        count += 1

    if count > 0:
        n = len(f',\n{indent}')
        s = s[:-n] + f'\n{indent[:-1]}):\n'
    else:
        s += '):\n'
    f.write(s)


def write_docstring(f:typing.TextIO, subroutine_name:str,
                    arg_info: dict[str, dict[str, str|bool|None|list[str]]], arg_list:list[str]) -> None:
    """write docstring based on Fortran arguments comments
    """
    s = '    \"\"\"'
    # todo: this needs to be done better
    s += f'See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/{subroutine_name[3:]} or PDAF source files \n\n'
    f.write(s)

    indent = ' '*4
    s = indent + 'Parameters\n'
    s += indent + '----------\n'
    count = 0
    for arg in arg_list:
        info = arg_info[arg]
        assert type(info['kind']) is str, f"kind in info {info['kind']} is not a str"
        assert type(info['type']) is str, f"type in info {info['type']} is not a str"
        assert type(info['comment']) is str, f"comment in info {info['comment']} is not a str"
        if type(info['intent']) is str:
            if 'in' not in info['intent'] \
                and ('_cb' not in subroutine_name):
                continue
        if info['type'] == 'procedure':
            s += indent + f'py{info["kind"][1:]} : func\n'
        elif info['array']:
            s += indent+arg +f' : ndarray[{pyconv[info["type"]]}]\n'
        else:
            s += indent + f'{arg} : {pyconv[info["type"]]}\n'
        s += 2*indent + info['comment'] +'\n'
        count += 1
    if count > 0:
        f.write(s)

    s = '\n'
    s += indent + 'Returns\n'
    s += indent + '-------\n'
    count = 0
    for arg, info in arg_info.items():
        assert type(info['type']) is str, f"type in info {info['type']} is not a str"
        assert type(info['comment']) is str, f"comment in info {info['comment']} is not a str"
        if info['intent'] is None:
            continue
        if type(info['intent']) is str:
            if 'out' not in info['intent']:
                continue

        if info['array']:
            s += indent+arg +f' : ndarray[{pyconv[info["type"]]}]\n'
        elif info['type'] == 'type':
            s += indent+arg + f' : ndarray[float]\n'
        else:
            s += indent + f'{arg} : {pyconv[info["type"]]}\n'
        s += 2*indent + info['comment'] +'\n'
        count += 1
    if count > 0:
        f.write(s)

    s = '    \"\"\"\n'
    f.write(s)


def write_memory_view(f:typing.TextIO, subroutine_name:str, arg_info: dict[str, dict[str, str|bool|None|list[str]]]):
    """Convert 2D arrays to fortran contiguous arrays
    """
    # convert np arrays to memoryview
    for argname, info in arg_info.items():
        if not info['array']:
            continue

        if type(info['intent']) is str:
            if 'in' not in info['intent'] and '_cb' not in subroutine_name:
                continue
        assert type(info['dimension']) is list, f"dimension in info {info['dimension']} is not a list"
        assert  type(info['type']) is str, f"type in info {info['type']} is not a str"
        if len(info['dimension']) <= 1: continue
        s = ' '*4 + f'cdef {conv[info["type"]]}[::1'
        s += f'] {argname}_f = np.asfortranarray({argname}).ravel(order="F")\n'
        f.write(s)


def write_dims(f, arg_info, ArrayDims):
    """get dimensions of input array arguments.
    """

    # define dimension of arrays
    s = '    cdef int '
    count = 0
    for dim in ArrayDims:
        s += f'{dim}, '
        count += 1
    s = s[:-2]+'\n' if count > 0 else ''
    f.write(s)

    # sort the dimension of arrays
    d = sorted(arg_info.items(), key=lambda item: -1 if not item[1]['array'] else len(item[1]['dimension']), reverse=True)
    ndim = {k: v['dimension'] for k, v in d if v['array'] and 'in' in v['intent']}

    for argname, dims in ndim.items():
        dimnames = []
        adjusts = []
        for i, dim in enumerate(dims):
            dimname, adjust = dim, ''
            for sign in ['+', '-']:
                if sign in dim:
                    dimname, adjust = dim.split(sign)
                    adjust = ' '.join([sign, adjust])

            dimname = dimname if dimname in ArrayDims else '_'
            if dimname in ArrayDims:
                ArrayDims.remove(dimname)

            dimnames.append(dimname)
            adjusts.append(adjust)

        if all([dimname == '_' for dimname in dimnames]):
            continue

        s = ''
        sadjust = '    '
        for i, (dimname, adjust) in enumerate(zip(dimnames, adjusts)):
            s += ' '*4 + f'{dimname} = {argname}.shape[{i}]\n'
            for sign, nsign in zip(['+', '-'], ['-', '+']):
                if sign in adjust:
                    sadjust += f'{dimname} = {dimname} {adjust.replace(sign, nsign)}\n'
        f.write(s)
        if not sadjust.isspace():
            f.write(sadjust)
    f.write('\n' if count > 0 else '')


def write_user_Cython(f:typing.TextIO, arg_info: dict[str, dict[str, str|bool|None|list[str]]]) -> None:
    # convert python function to Cython function
    count = 0
    for argname, info in arg_info.items():
        assert type(info['kind']) is str, f"kind in info {info['kind']} is not a str"
        if info['type'] == 'procedure':
            s = ' '*4 + f'c__PDAFcython.{info["kind"][3:]} = <void*>py{info["kind"][1:]}\n'
            f.write(s)
            if argname == 'u_init_ens':
                s = ' '*4 + 'c__PDAFcython.init_ens_pdaf_single_member = <void*>py__init_ens_pdaf\n'
                f.write(s)
        count += 1
    f.write('\n' if count > 0 else '')


def write_return_def(f, arg_info):
    indent = ' '*4
    count = 0
    # define return variables
    for argname, info in arg_info.items():
        if info['intent'] != 'out':
            continue
        if info['array']:
            s = indent + f'cdef {conv[info["type"]]} [::1] {argname}'
            s += f' = np.zeros(({', '.join(info["dimension"])}), dtype={npyconv[info["type"]]}'
            s += ').ravel()\n'
        else:
            s = indent + f'cdef {conv[info["type"]]} {argname}\n'
        f.write(s)
        count += 1
    f.write('\n' if count > 0 else '')


def write_func_call(f, subroutine_name, arg_info):
    # special treatment for init subroutine
    # because init_ens_pdaf uses a different interface
    if subroutine_name == 'c__pdaf_init':
        s = ' '*4 + f'if (filtertype == 0) or (filtertype == 200 and subtype == 0):\n'
        s += ' '*8 + 'c__pdaf_init (&filtertype, &subtype, &stepnull,\n'
        s += ' '*8 + '              &param_int[0], &dim_pint,\n'
        s += ' '*8 + '              &param_real[0], &dim_preal,\n'
        s += ' '*8 + '              &comm_model, &comm_filter, &comm_couple,\n'
        s += ' '*8 + '              &task_id, &n_modeltasks, &in_filterpe,\n'
        s += ' '*8 + '              c__PDAFcython.c__init_ens_pdaf_single_member,\n'
        s += ' '*8 + '              &in_screen, &flag)\n'
        s += ' '*4 + 'else:\n'
        f.write(s)
    indent = ' '*4 if subroutine_name != 'c__pdaf_init' else ' '*8
    # call the actual subroutine
    s = indent + f'{subroutine_name} ('
    indent = ' '*len(s)
    for argname, info in arg_info.items():
        if info['array']:
            if len(info["dimension"]) > 1 and ('in' in info['intent']
              or '_cb' in subroutine_name):
                s += f'&{argname}_f[0]'
            else:
                s += f'&{argname}[0]'
        elif info['type'] == 'procedure':
            s += f'c__PDAFcython.{info["kind"]}'
        else:
            s += f'&{argname}'
        s += ',\n' + indent

    if len(arg_info) > 0:
        n = len(f',\n{indent}')
        s = s[:-n] + f'\n{indent[:-1]})\n'
    else:
        s += ')\n'
    f.write(s + '\n')


def write_returns(f, arg_info):
    s = '    return '
    ss = ''
    count = 0
    for argname, info in arg_info.items():
        if info['intent'] is None:
            continue
        if 'out' not in info['intent']:
            continue
        # todo: hard coded pointer to np array conversion
        if argname == 'dims':
            ss = '    dims = np.asarray(dims)\n'
            continue

        if info['array']:
            s += f'np.asarray({argname}).reshape(({', '.join(info["dimension"])}), order=''\'F\'''), '
        elif info['type'] == 'type':
            s += f'np.asarray(<double[:np.prod(dims)]> {argname}).reshape(dims, order=''\'F\'''), \\\n           '
        else:
            s += f'{argname}, '
        count += 1
    f.write(ss)
    if count > 0:
        s = s[:-2]
        f.write(s + '\n\n')


def writeCallBackDef(f, name, routine):
    # define function
    if name[3:11] == 'pdafomi_':
        funcname = name[11:]
    elif name[3:8] == 'pdaf_':
        funcname = name[8:]
    else:
        funcname = name[3:]

    s = f'def {funcname} ():\n'
    indent = ' '*4
    s += indent + f'return {name}\n'
    f.write(s+'\n\n')


def write_PDAF_calls(filename:str, func_info: dict[str, dict[str, dict[str, str|bool|None|list[str]]]]) -> None:
    """write the PDAF interface calls"""
    with open(filename, 'w') as f:
        # MPI exception handling
        s:str  = 'import sys\n'
        s += '\n'
        s += 'import pyPDAF.UserFunc as PDAFcython\n'
        s += 'cimport pyPDAF.UserFunc as c__PDAFcython\n'
        s += '\n'
        s += 'import numpy as np\n'
        s += 'cimport numpy as cnp\n'
        s += '\n'
        s += 'try:\n'
        s += '    import mpi4py\n'
        s += '    mpi4py.rc.initialize = False\n'
        s += 'except ImportError:\n'
        s += '    pass\n'
        s += '\n'
        s += '# Global error handler\n'
        s += 'def global_except_hook(exctype, value, traceback):\n'
        s += '    from traceback import print_exception\n'
        s += '    try:\n'
        s += '        import mpi4py.MPI\n'
        s += '\n'
        s += '        if mpi4py.MPI.Is_initialized():\n'
        s += '            try:\n'
        s += '                sys.stderr.write("Uncaught exception was detected on rank {}. \\n".format(\n'
        s += '                    mpi4py.MPI.COMM_WORLD.Get_rank()))\n'
        s += '\n'
        s += '                print_exception(exctype, value, traceback)\n'
        s += '                sys.stderr.write("\\n")\n'
        s += '                sys.stderr.flush()\n'
        s += '            finally:\n'
        s += '                try:\n'
        s += '                    mpi4py.MPI.COMM_WORLD.Abort(1)\n'
        s += '                except Exception as e:\n'
        s += '                    sys.stderr.write("MPI Abort failed, this process will hang.\\n")\n'
        s += '                    sys.stderr.flush()\n'
        s += '                    raise e\n'
        s += '        else:\n'
        s += '            sys.__excepthook__(exctype, value, traceback)\n'
        s += '    except ImportError:\n'
        s += '        sys.__excepthook__(exctype, value, traceback)\n'
        s += '\n'
        s += 'sys.excepthook = global_except_hook\n'
        s += '\n\n'
        f.write(s)

        for subroutine_name in func_info:
            arg_list, ArrayDims = get_pyx_arg_list(subroutine_name.lower(), func_info[subroutine_name])
            write_func_def(f, subroutine_name.lower(), func_info[subroutine_name], arg_list)
            write_docstring(f, subroutine_name, func_info[subroutine_name], arg_list)
            write_memory_view(f, subroutine_name.lower(), func_info[subroutine_name])
            if '_cb' not in subroutine_name:
                write_dims(f, func_info[subroutine_name], ArrayDims)
                write_user_Cython(f, func_info[subroutine_name])
                write_return_def(f, func_info[subroutine_name])
            write_func_call(f, subroutine_name.lower(), func_info[subroutine_name])
            write_returns(f, func_info[subroutine_name])


if __name__ == '__main__':
    import get_interface_info
    import write_pxd
    user_func_info = get_interface_info.get_func_info(['../pyPDAF/fortran/U_PDAF_interface_c_binding.F90'])
    PDAF_func_info = get_interface_info.get_func_info(['../pyPDAF/fortran/PDAF_c_binding.F90', '../pyPDAF/fortran/PDAFomi_obs_c_binding.F90'])
    write_pxd.write_Pxd_file('PDAF.pxd', PDAF_func_info, user_func_info)
    write_PDAF_calls('PDAF.pyx', PDAF_func_info)






