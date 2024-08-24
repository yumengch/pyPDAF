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
mypy_conv = {'integer' : 'int', 'logical': 'bool', 'real': 'float', 'character':'str',}

def get_current_column(file:typing.TextIO) -> int:
    current_pos = file.tell()
    
    # Start reading backwards from the current position
    line_start_pos = current_pos
    
    while line_start_pos > 0:
        line_start_pos -= 1
        file.seek(line_start_pos)
        char = file.read(1)
        if char == '\n':
            line_start_pos += 1
            break
    
    # The column number is the difference between the current position and the last newline position
    column = current_pos - line_start_pos
    return column

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


def write_func_def(f:typing.TextIO, subroutine_name:str, user_func_info:dict[str, dict[str, dict[str, str|bool|None|list[str]]]], arg_info:dict[str, dict[str, str|bool|None|list[str]]], arg_list: list[str]) -> None:
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
    if subroutine_name[3:11].lower() == 'pdafomi_':
        funcname = subroutine_name[7:]
    elif subroutine_name[3:8].lower() == 'pdaf_':
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
            s += f'py{info["kind"][1:]} : '
            # print (info["kind"][1:])
            user_arg_info = user_func_info[f'c{info["kind"][1:]}']
            s += 'Callable['
            indent_user = ' '*len(s.split('\n')[-1])
            for i, (_, u_arg_info) in enumerate(user_arg_info.items()):
                if u_arg_info['array']:
                    s_dims:str = 'int, '*len(u_arg_info['dimension'])
                    s_dims = s_dims[:-2]
                    s_dtype: str = npyconv[u_arg_info['type'].lower()]
                    s += f'np.ndarray[tuple[{s_dims}], np.dtype[{s_dtype}]], '
                else:
                    s += mypy_conv[u_arg_info['type']] + ', '
                if len(s.split('\n')[-1]) > 70 and i < len(user_arg_info) - 1:
                    s += '\n'
                    s += indent_user
            s = s[:-2]
            s += ']'
        elif info['array']:
            assert type(info['dimension']) is list, f"dimension in info {info['dimension']} is not a list"
            s_dims:str = 'int, '*len(info['dimension'])
            s_dims = s_dims[:-2]
            s_dtype: str = npyconv[info['type'].lower()]
            s += f'{argname}: np.ndarray[tuple[{s_dims}], np.dtype[{s_dtype}]]'
        else:
            s += f'{argname}: {conv[info["type"].lower()]}'

        s += f','
        if len(s.split('\n')[-1]) > 70 or info['type'] == 'procedure':
            s += '\n'
            s += indent
        count += 1


    n = len(f',\n{indent}')
    if s[-n:] == f',\n{indent}':
        s = s[:-n] + f'\n{indent[:-1]})'
    else:
        s = s[:-1] if count > 0 else s
        s += ')'
    f.write(s)


def write_returns(f:typing.TextIO, arg_info:dict[str, dict[str, str|bool|None|list[str]]]) -> None:

    count = 0
    indent = ' '*get_current_column(f)
    s = ' -> tuple['
    for argname, info in arg_info.items():
        if info['intent'] is None:
            continue
        if 'out' not in info['intent']:
            continue
        # todo: hard coded pointer to np array conversion
        if argname == 'dims':
            continue

        if info['array']:
            s_dims:str = 'int, '*len(info['dimension'])
            s_dims = s_dims[:-2]
            s_dtype: str = npyconv[info['type'].lower()]
            s += f'np.ndarray[tuple[{s_dims}], np.dtype[{s_dtype}]], '
        elif info['type'] == 'type':
            s_dims:str = 'int, '*len(arg_info['dims']['dimension'])
            s += f'np.ndarray[tuple[{s_dims}], np.dtype[float]], '
        else:
            s += f'{mypy_conv[info["type"]]}, '
        if len(s.split('\n')[-1]) > 70:
            s += f'\n{indent}'
        count += 1

    if count > 1:
        n = len(f',\n{indent}')
        if s[-n:] == f',\n{indent}':
            s = s[:-n-3] + ']: ...'
        else:
            s = s[:-2] + ']: ...'
        f.write(s + '\n\n')
    elif count == 1:
        s = ' -> ' + s[10:-2] + ': ...'
        f.write(s + '\n\n')
    else:
        f.write(' -> None: ... \n\n')

def write_PDAF_calls(filename:str, user_func_info:dict[str, dict[str, dict[str, str|bool|None|list[str]]]], func_info: dict[str, dict[str, dict[str, str|bool|None|list[str]]]]) -> None:
    """write the PDAF interface calls"""
    with open(filename, 'w+') as f:
        s = 'import numpy as np\n'
        s += 'from typing import Callable'
        s += '\n\n'
        f.write(s)

        for subroutine_name in func_info:
            arg_list, ArrayDims = get_pyx_arg_list(subroutine_name.lower(), func_info[subroutine_name])
            write_func_def(f, subroutine_name, user_func_info, func_info[subroutine_name], arg_list)
            write_returns(f, func_info[subroutine_name])

if __name__ == '__main__':
    import get_interface_info
    user_func_info = get_interface_info.get_func_info(['../src/pyPDAF/fortran/U_PDAF_interface_c_binding.F90'])
    PDAF_func_info = get_interface_info.get_func_info(['../src/pyPDAF/fortran/PDAF_c_binding.F90', '../src/pyPDAF/fortran/PDAFomi_obs_c_binding.F90'])
    write_PDAF_calls('PDAF.pyi', user_func_info, PDAF_func_info)
