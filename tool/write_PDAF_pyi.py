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
mypy_conv = {'integer' : 'int', 'logical': 'bool', 'real': 'float', 'character':'str', 'logical': 'bool',}

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


def get_pyx_arg_list(subroutine_name:str, arg_info: dict[str, dict[str, str|list[str]]]) -> tuple[list[str], list[str]]:
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
        if len(info['dimension']) > 0:
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


def write_func_def(f:typing.TextIO, subroutine_name:str, user_func_info:dict[str, dict[str, dict[str, str|list[str]]]], arg_info:dict[str, dict[str, str|list[str]]], arg_list: list[str]) -> None:
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
    elif subroutine_name[3:13].lower() == 'pdaflocal_':
        funcname = subroutine_name[7:]
    elif subroutine_name[3:16].lower() == 'pdaflocalomi_':
        funcname = subroutine_name[7:]
    else:
        funcname = subroutine_name[3:]

    s:str = f'def {funcname} ('
    indent:str = ' '*len(s)
    count:int = 0
    count_callback_returns:int = 0
    len_callbakc_returns:int = 0
    s_dims: str
    s_dtype: str
    # only arguments with input values
    for argname in arg_list:
        info = arg_info[argname]
        assert  type(info['type']) is str, f"type in info {info['type']} is not a str"

        # we do not write purely return arguments
        # callback functions are exceptions
        if (info['intent'] != '') and ('in' not in info['intent']) and ('_cb' not in funcname):
            continue

        # write procedure
        if info['type'] == 'procedure':
            assert  type(info['kind']) is str, f"kind in info {info['kind']} is not a str"
            s += f'py{info["kind"][1:]} : '
            user_arg_info = user_func_info[f'c{info["kind"][1:]}']
            s += 'Callable[['
            indent_user = ' '*len(s.split('\n')[-1])
            # write the arguments of the callback fuction
            for i, (_, u_arg_info) in enumerate(user_arg_info.items()):
                assert type(u_arg_info['type']) is str, f"type in info {u_arg_info['type']} is not a str"
                if len(u_arg_info['dimension']) > 0:
                    s_dims = 'int, '*len(u_arg_info['dimension'])
                    s_dims = s_dims[:-2]
                    s_dtype = npyconv[u_arg_info['type'].lower()]
                    s += f'np.ndarray[tuple[{s_dims}], np.dtype[{s_dtype}]], '
                else:
                    s += mypy_conv[u_arg_info['type']] + ', '
                if len(s.split('\n')[-1]) > 70 and i < len(user_arg_info) - 1:
                    s += '\n'
                    s += indent_user
            n = len(f', \n{indent_user}')
            if s[-n:] == f', \n{indent_user}':
                s = s[:-n]
            else:
                s = s[:-2]
            s += '], '
            count_callback_returns = 0
            len_callbakc_returns = len(s)
            s += 'tuple['
            # write the returns of the callback function
            for i, (_, u_arg_info) in enumerate(user_arg_info.items()):
                assert type(u_arg_info['type']) is str, f"type in info {u_arg_info['type']} is not a str"
                if 'out' in u_arg_info['intent']:
                    count_callback_returns += 1
                    if len(u_arg_info['dimension']) > 0:
                        s_dims = 'int, '*len(u_arg_info['dimension'])
                        s_dims = s_dims[:-2]
                        s_dtype = npyconv[u_arg_info['type'].lower()]
                        s += f'np.ndarray[tuple[{s_dims}], np.dtype[{s_dtype}]], '
                    else:
                        s += mypy_conv[u_arg_info['type']] + ', '
                    if len(s.split('\n')[-1]) > 70 and i < len(user_arg_info) - 1:
                        s += '\n'
                        s += indent_user
            if count_callback_returns > 0:
                n = len(f', \n{indent_user}')
                if s[-n:] == f', \n{indent_user}':
                    s = s[:-n]
                else:
                    s = s[:-2]
                if count_callback_returns == 1:
                    s = s[:len_callbakc_returns] + s[len_callbakc_returns+len('tuple['):]
                    s += ']'
                else:
                    s += ']]'


            else:
                s = s[:-6]
                s += 'None]'
        elif len(info['dimension']) > 0:
            s_dims = 'int, '*len(info['dimension'])
            s_dims = s_dims[:-2]
            s_dtype = npyconv[info['type'].lower()]
            s += f'{argname}: np.ndarray[tuple[{s_dims}], np.dtype[{s_dtype}]]'
        else:
            s += f'{argname}: {mypy_conv[info["type"].lower()]}'

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


def write_returns(f:typing.TextIO, arg_info:dict[str, dict[str, str|list[str]]]) -> None:
    s_dims : str
    s_dtype : str
    count = 0
    indent = ' '*get_current_column(f)
    s = ' -> tuple['
    for argname, info in arg_info.items():
        if 'out' not in info['intent']:
            continue
        # todo: hard coded pointer to np array conversion
        if argname == 'dims':
            continue

        assert type(info['type']) is str, f"type in info {info['type']} is not a str"
        if len(info['dimension']) > 0:
            s_dims = 'int, '*len(info['dimension'])
            s_dims = s_dims[:-2]
            s_dtype = npyconv[info['type'].lower()]
            s += f'np.ndarray[tuple[{s_dims}], np.dtype[{s_dtype}]], '
        elif info['type'] == 'type':
            s_dims = 'int, '*len(arg_info['dims']['dimension'])
            s += f'np.ndarray[tuple[{s_dims}], np.dtype[float]], '
        else:
            s += f'{mypy_conv[info["type"]]}, '
        if len(s.split('\n')[-1]) > 70:
            s += f'\n{indent}'
        count += 1

    if count > 1:
        n = len(f',\n{indent}')
        if s[-n:] == f',\n{indent}':
            s = s[:-n-3] + ']:'
        else:
            s = s[:-2] + ']:'
        f.write(s + '\n')
    elif count == 1:
        s = ' -> ' + s[10:-2] + ':'
        f.write(s + '\n')
    else:
        f.write(' -> None:\n')


def write_docstring(f:typing.TextIO, subroutine_name:str, user_func_info:dict[str, dict[str, dict[str, str|list[str]]]],
                    arg_info: dict[str, dict[str, str|list[str]]], arg_list:list[str]) -> None:
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
    s_dims:str
    s_dtype: str
    for arg in arg_list:
        info = arg_info[arg]
        assert type(info['kind']) is str, f"kind in info {info['kind']} is not a str"
        assert type(info['type']) is str, f"type in info {info['type']} is not a str"
        assert type(info['comment']) is str, f"type in info {info['type']} is not a str"

        if (info['intent'] != '') and ('in' not in info['intent']) and ('_cb' not in subroutine_name):
            continue

        if info['type'] == 'procedure':
            user_arg_info = user_func_info[f'c{info["kind"][1:]}']
            s += indent + f'py{info["kind"][1:]} : '
            s += 'Callable['
            for u_arg, u_arg_info in user_arg_info.items():
                assert type(u_arg_info['type']) is str, f"type in u_arg_info {u_arg_info['type']} is not a str"
                if len(u_arg_info['dimension']) > 0:
                    s_dims = ', '.join(u_arg_info['dimension'])
                    s_dtype = npyconv[u_arg_info['type']]
                    s += f'{u_arg}:ndarray[tuple[{s_dims}], {s_dtype}], '
                else:
                    s += u_arg + ':' + mypy_conv[u_arg_info['type']] + ', '
            s = s[:-2]
            s += ']\n'
            s += 2*indent + info['comment'] +'\n\n'
            s += 2*indent + '**Callback Parameters**\n\n'
            for u_arg, u_arg_info in user_arg_info.items():
                assert type(u_arg_info['type']) is str, f"type in u_arg_info {u_arg_info['type']} is not a str"
                assert type(u_arg_info['comment']) is str, f"type in u_arg_info {u_arg_info['type']} is not a str"
                if len(u_arg_info['dimension']) > 0:
                    s_dims = ', '.join(u_arg_info['dimension'])
                    s_dtype = npyconv[u_arg_info['type']]
                    s += 2*indent + f'* {u_arg} : ndarray[tuple[{s_dims}], {s_dtype}]\n\n'
                else:
                    s += 2*indent + f'* {u_arg}' + ' : ' + mypy_conv[u_arg_info['type']] + '\n\n'
                s += 3*indent + '* ' + u_arg_info['comment'] +'\n\n'
 
            s += 2*indent + '**Callback Returns**\n\n'
            for u_arg, u_arg_info in user_arg_info.items():
                if 'out' not in u_arg_info['intent']: continue
                assert type(u_arg_info['type']) is str, f"type in u_arg_info {u_arg_info['type']} is not a str"
                assert type(u_arg_info['comment']) is str, f"type in u_arg_info {u_arg_info['type']} is not a str"
                if len(u_arg_info['dimension']) > 0:
                    s_dims = ', '.join(u_arg_info['dimension'])
                    s_dtype = npyconv[u_arg_info['type']]
                    s += 2*indent + f'* {u_arg} : ndarray[tuple[{s_dims}], {s_dtype}]\n\n'
                else:
                    s += 2*indent + f'* {u_arg}' + ':' + mypy_conv[u_arg_info['type']] + '\n\n'
                s += 3*indent+ '* ' + u_arg_info['comment'] +'\\n\n'
            s += '\n'
        elif len(info['dimension']) > 0:
            s += indent+arg +f' : '
            s_dims = ', '.join(info['dimension'])
            s_dtype = npyconv[info['type'].lower()]
            s += f'ndarray[tuple[{s_dims}], {s_dtype}]\n'
            s += 2*indent + info['comment'] +'\n'
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

        if 'out' not in info['intent']:
            continue

        if len(info['dimension']) > 0:
            s += indent+arg +f' : '
            s_dims = ', '.join(info['dimension'])
            s_dtype = npyconv[info['type'].lower()]
            s += f'ndarray[tuple[{s_dims}], {s_dtype}]\n '
        elif info['type'] == 'type':
            s += indent+arg + f' : ndarray[float]\n'
        else:
            s += indent + f'{arg} : {pyconv[info["type"]]}\n'

        s += 2*indent + info['comment'] +'\n'
        count += 1
    if count > 0:
        f.write(s)

    s = '    \"\"\"\n'
    s += '    ...\n\n'
    f.write(s)


def write_PDAF_calls(filename:str, user_func_info:dict[str, dict[str, dict[str, str|list[str]]]], func_info: dict[str, dict[str, dict[str, str|list[str]]]]) -> None:
    """write the PDAF interface calls"""
    with open(filename, 'w+') as f:
        s = 'import numpy as np\n'
        s += 'from typing import Callable'
        s += '\n\n'
        f.write(s)

        for subroutine_name in func_info:
            arg_list, _ = get_pyx_arg_list(subroutine_name.lower(), func_info[subroutine_name])
            write_func_def(f, subroutine_name, user_func_info, func_info[subroutine_name], arg_list)
            write_returns(f, func_info[subroutine_name])
            write_docstring(f, subroutine_name, user_func_info, func_info[subroutine_name], arg_list)

if __name__ == '__main__':
    import get_interface_info
    user_func_info = get_interface_info.get_func_info(['../src/fortran/U_PDAF_interface_c_binding.F90'])
    PDAF_func_info = get_interface_info.get_func_info(['../src/fortran/PDAF_c_binding.F90',
                                                       '../src/fortran/PDAFomi_obs_c_binding.F90',
                                                       '../src/fortran/PDAFlocal_c_binding.F90', ])
    write_PDAF_calls('PDAF.pyi', user_func_info, PDAF_func_info)
