"""This is a simple script for converting C binding
fortran routines in pyPDAF/fortran to Cython definition and
implementation files.
"""
import os
import typing
import re
import docstrings

# todo: type for pointer should be matching its corresponding fortran pointer type
# here double* is used because PDAF only uses double
conv = {'integer': 'int', 'logical': 'bint', 'real': 'double',
        'procedure': 'void', 'character': 'CFI_cdesc_t', 'type': 'double*'}
pyconv = {'integer': 'int', 'logical': 'bool', 'real': 'float',
          'procedure': 'void', 'character': 'CFI_cdesc_t', 'type': 'double*'}
npyconv = {'integer': 'np.intc', 'real': 'np.float64'}
mypy_conv = {'integer': 'int', 'logical': 'bool',
             'real': 'float', 'character': 'str', }


def make_ordinal(n) -> str:
    '''
    Convert an integer into its ordinal representation::

    This code is from `https://stackoverflow.com/a/50992575`_.
        make_ordinal(0)   => '0th'
        make_ordinal(3)   => '3rd'
        make_ordinal(122) => '122nd'
        make_ordinal(213) => '213th'
    '''
    n = int(n)
    if 11 <= (n % 100) <= 13:
        suffix = 'th'
    else:
        suffix = ['th', 'st', 'nd', 'rd', 'th'][min(n % 10, 4)]
    return str(n) + suffix


def extract_dimension_name(s: str) -> str | None:
    """extract dimension names from dimension
    """
    # Match one or more word characters (letters, digits, and underscores), but exclude digits
    match = re.search(r'[a-zA-Z_]+', s)
    if match:
        return match.group()
    return None


def get_pyx_arg_list(subroutine_name: str,
                     arg_info: dict[str, dict[str, str | list[str]]]
                     ) -> tuple[list[str], list[str]]:
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
    array_dims: list[str] = []
    for _, info in arg_info.items():
        # we need to keep the dimension of output arrays
        # when the dimension is only used by output arrays
        if info['intent'] == 'out':
            continue
        if len(info['dimension']) == 0:
            continue

        for dim in info['dimension']:
            dimsize = extract_dimension_name(dim)
            if dimsize in arg_info:
                if dimsize not in array_dims:
                    array_dims.append(dimsize)

    arg_list: list[str] = []
    for arg_name, info in arg_info.items():
        # we collect all arguments that is not input array dimensions
        # however, we collect all arguments whatsoever for callback functions
        if arg_name not in array_dims or ('_cb' in subroutine_name):
            arg_list.append(arg_name)

    return arg_list, array_dims


def write_func_def(f: typing.TextIO, subroutine_name: str,
                   arg_info: dict[str, dict[str, str | list[str]]],
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

    s = f'def {funcname} ('
    indent = ' '*len(s)
    count = 0
    # only arguments with input values
    for argname in arg_list:
        info = arg_info[argname]
        assert isinstance(info['type'], str), \
            f"type in info {info['type']} is not a str"
        # we do not write purely return arguments
        # callback functions are exceptions
        if (info['intent'] != '') and ('in' not in info['intent']) and ('_cb' not in funcname):
            continue

        # write procedure
        if info['type'] == 'procedure':
            assert isinstance(info['kind'], str), \
                f"kind in info {info['kind']} is not a str"
            s += f'py{info["kind"][1:]}'
        elif len(info['dimension']) > 0:
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


def write_docstring(
        f: typing.TextIO, subroutine_name: str,
        user_func_info: dict[str, dict[str, dict[str, str | list[str]]]],
        arg_info: dict[str, dict[str, str | list[str]]],
        arg_list: list[str]) -> None:
    """write docstring based on Fortran arguments comments
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
    s = '    r\"\"\"'
    # todo: this needs to be done better
    s += docstrings.docstrings[funcname]+'\n\n'
    f.write(s)

    indent = ' '*4
    s = indent + 'Parameters\n'
    s += indent + '----------\n'
    count = 0
    s_dims: str
    s_dtype: str
    for arg in arg_list:
        info = arg_info[arg]
        assert isinstance(info['kind'], str), \
            f"kind in info {info['kind']} is not a str"
        assert isinstance(info['type'], str), \
            f"type in info {info['type']} is not a str"
        assert isinstance(info['comment'], str), \
            f"type in info {info['type']} is not a str"

        if (info['intent'] != ''
            ) and ('in' not in info['intent']
                   ) and ('_cb' not in subroutine_name):
            continue

        if info['type'] == 'procedure':
            user_arg_info = user_func_info[f'c{info["kind"][1:]}']
            s += indent + f'py{info["kind"][1:]} : '
            s += 'Callable['
            for u_arg, u_arg_info in user_arg_info.items():
                assert isinstance(u_arg_info['type'], str), \
                    f"type in u_arg_info {u_arg_info['type']} is not a str"
                if len(u_arg_info['dimension']) > 0:
                    s_dims = ', '.join(u_arg_info['dimension'])
                    s_dtype = npyconv[u_arg_info['type']]
                    s += f'{u_arg} : ndarray[tuple[{s_dims}], {s_dtype}], '
                else:
                    s += u_arg + ':' + mypy_conv[u_arg_info['type']] + ', '
            s = s[:-2]
            s += ']\n'
            s += 2*indent + info['comment'].replace(
                '\\n', f'\n{(2*indent)[:-1]}') + '\n\n'
            s += 2*indent + '**Callback Parameters**\n'
            for u_arg, u_arg_info in user_arg_info.items():
                assert isinstance(u_arg_info['type'], str), \
                    f"type in u_arg_info {u_arg_info['type']} is not a str"
                assert isinstance(u_arg_info['comment'], str), \
                    f"type in u_arg_info {u_arg_info['type']} is not a str"
                if len(u_arg_info['dimension']) > 0:
                    s_dims = ', '.join(u_arg_info['dimension'])
                    s_dtype = npyconv[u_arg_info['type']]
                    s += 3*indent + f'* **{u_arg}** :' \
                        f' ndarray[tuple[{s_dims}], {s_dtype}]\n'
                else:
                    s += 3*indent + f'* **{u_arg}** : ' + \
                        mypy_conv[u_arg_info['type']] + '\n'
                s += 4 * indent + '* ' + u_arg_info['comment'].replace(
                    '\n', f'\n{2 * indent}   ') + '\n'

            s += '\n'
            s += 2*indent + '**Callback Returns**\n'
            for u_arg, u_arg_info in user_arg_info.items():
                if 'out' not in u_arg_info['intent']:
                    continue
                assert isinstance(u_arg_info['type'], str), \
                    f"type in u_arg_info {u_arg_info['type']} is not a str"
                assert isinstance(u_arg_info['comment'], str), \
                    f"type in u_arg_info {u_arg_info['type']} is not a str"
                if len(u_arg_info['dimension']) > 0:
                    s_dims = ', '.join(u_arg_info['dimension'])
                    s_dtype = npyconv[u_arg_info['type']]
                    s += 3*indent + f'* **{u_arg}** :' \
                        f' ndarray[tuple[{s_dims}], {s_dtype}]\n'
                else:
                    s += 3*indent + f'* **{u_arg}** : ' + \
                        mypy_conv[u_arg_info['type']] + '\n'
                s += 4 * indent + '* ' + u_arg_info['comment'].replace(
                    '\n', f'\n{2 * indent}   ') + '\n'
        elif len(info['dimension']) > 0:
            s += indent+arg + ' : '
            s_dims = ', '.join(info['dimension'])
            s_dtype = npyconv[info['type'].lower()]
            s += f'ndarray[tuple[{s_dims}], {s_dtype}]\n'
            # comments of the array
            s += 2*indent + info['comment'] + '\n\n'
            # adding documentation for array dimensions
            if len(info['dimension']) > 1:
                i = 0
                for dim in info['dimension']:
                    dim0 = dim.split('+')[0].split('-')[0].replace(' ', '')
                    if dim in arg_info:
                        i += 1
                        s += 2*indent
                        s += 'The' if i == 1 else 'the'
                        s += f'{make_ordinal(i)}-th dimension {dim0}' \
                            f' is {arg_info[dim0]["comment"]}'
                        s += ';\n' if i == 1 else '.\n'
            else:
                dim = info['dimension'][0]
                dim0 = dim.split('+')[0].split('-')[0].replace(' ', '')
                if dim in arg_info:
                    s += 2*indent + f'The array dimension `{dim0}`' \
                        f' is {arg_info[dim0]["comment"]}.\n'

        else:
            s += indent + f'{arg} : {pyconv[info["type"]]}\n'
            s += 2*indent + info['comment'] + '\n'

        count += 1

    if count > 0:
        f.write(s)

    s = '\n'
    s += indent + 'Returns\n'
    s += indent + '-------\n'
    count = 0
    for arg, info in arg_info.items():
        assert isinstance(info['type'], str), \
            f"type in info {info['type']} is not a str"
        assert isinstance(info['comment'], str), \
            f"comment in info {info['comment']} is not a str"

        if 'out' not in info['intent']:
            continue

        if len(info['dimension']) > 0:
            s += indent+arg + ' : '
            s_dims = ', '.join(info['dimension'])
            s_dtype = npyconv[info['type'].lower()]
            s += f'ndarray[tuple[{s_dims}], {s_dtype}]\n'
        elif info['type'] == 'type':
            s += indent+arg + ' : ndarray[float]\n'
        else:
            s += indent + f'{arg} : {pyconv[info["type"]]}\n'

        s += 2*indent + info['comment'] + '\n'

        # add comments for array dimensions
        if len(info['dimension']) > 0:
            s += '\n'
            if len(info['dimension']) > 1:
                i = 0
                for dim in info['dimension']:
                    dim0 = dim.split('+')[0].split('-')[0].replace(' ', '')
                    if dim in arg_info:
                        i += 1
                        s += 2*indent + f'The {make_ordinal(i)}-th ' \
                            f'dimension {dim0} is {arg_info[dim0]["comment"]}\n'
            else:
                dim = info['dimension'][0]
                dim0 = dim.split('+')[0].split('-')[0].replace(' ', '')
                if dim in arg_info:
                    s += 2*indent + f'The array dimension `{dim0}`' \
                        f' is {arg_info[dim0]["comment"]}\n'
        count += 1
    if count > 0:
        f.write(s)

    s = '    \"\"\"\n\n'
    f.write(s)


def write_memory_view(
        f: typing.TextIO, subroutine_name: str,
        arg_info: dict[str, dict[str, str | list[str]]]):
    """Convert 2D arrays to fortran contiguous arrays
    """
    # convert np arrays to memoryview
    for argname, info in arg_info.items():
        if len(info['dimension']) <= 1:
            continue

        if 'in' not in info['intent'] and '_cb' not in subroutine_name:
            continue

        assert isinstance(info['type'], str), \
            f"type in info {info['type']} is not a str"
        s = ' '*4 + f'cdef {conv[info["type"]]}[::1'
        s += f'] {argname}_f = np.asfortranarray({argname}).ravel(order="F")\n'
        f.write(s)


def write_dims(
        f: typing.TextIO, arg_info: dict[str, dict[str, str | list[str]]],
        array_dims: list[str]) -> None:
    """get dimensions of input array arguments.
    """
    # define dimension of arrays
    s = '    cdef int '
    count = 0
    for dim in array_dims:
        s += f'{dim}, '
        count += 1
    s = s[:-2]+'\n' if count > 0 else ''
    f.write(s)

    # sort the dimension of arrays
    d = sorted(arg_info.items(),
               key=lambda item: -1
               if len(item[1]['dimension']) == 0 else
               len(item[1]['dimension']), reverse=True)
    ndim = {k: v['dimension'] for k, v in d if len(
        v['dimension']) > 0 and 'in' in v['intent']}

    for argname, dims in ndim.items():
        dimnames = []
        adjusts = []
        for i, dim in enumerate(dims):
            dimname, adjust = dim, ''
            for sign in ['+', '-']:
                if sign in dim:
                    dimname, adjust = dim.split(sign)
                    adjust = ' '.join([sign, adjust])

            dimname = dimname if dimname in array_dims else '_'
            if dimname in array_dims:
                array_dims.remove(dimname)

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
                    sadjust += f'{dimname} = {dimname} {
                        adjust.replace(sign, nsign)}\n'
        f.write(s)
        if not sadjust.isspace():
            f.write(sadjust)
    f.write('\n' if count > 0 else '')


def write_user_cython(f: typing.TextIO,
                      arg_info: dict[str, dict[str, str | list[str]]]) -> None:
    """convert python function to Cython function"""
    count = 0
    for _, info in arg_info.items():
        assert isinstance(info['kind'], str), f"kind in info {
            info['kind']} is not a str"
        if info['type'] == 'procedure':
            s = ' '*4 + f'c__PDAFcython.{info["kind"][3:]}' \
                f' = <void*>py{info["kind"][1:]}\n'
            f.write(s)
            count += 1
    f.write('\n' if count > 0 else '')


def write_return_def(f: typing.TextIO,
                     arg_info: dict[str, dict[str, str | list[str]]]) -> None:
    indent = ' '*4
    count = 0
    # define return variables
    for argname, info in arg_info.items():
        if info['intent'] != 'out':
            continue

        assert isinstance(info['type'], str), \
            f"type in info {info['type']} is not a str"
        if len(info['dimension']) > 0:
            s = indent + f'cdef {conv[info["type"]]} [::1] {argname}'
            s += f' = np.zeros(({', '.join(info["dimension"])}),' \
                f' dtype={npyconv[info["type"]]}'
            s += ').ravel()\n'
        else:
            s = indent + f'cdef {conv[info["type"]]} {argname}\n'
        f.write(s)
        count += 1
    f.write('\n' if count > 0 else '')


def write_func_call(f: typing.TextIO, subroutine_name: str,
                    arg_info: dict[str, dict[str, str | list[str]]]) -> None:
    """write the C function call in Cython"""
    # special treatment for init subroutine
    s = '    with nogil:\n'
    f.write(s)
    indent = ' '*8
    # call the actual subroutine
    s = indent + f'{subroutine_name.lower()} ('
    indent = ' '*len(s)
    for argname, info in arg_info.items():
        if len(info['dimension']) > 0:
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
        if 'out' not in info['intent']:
            continue
        # todo: hard coded pointer to np array conversion
        if argname == 'dims':
            ss = '    dims = np.asarray(dims)\n'
            continue

        if len(info['dimension']) > 0:
            s += f'np.asarray({argname}).reshape(('\
                f'{', '.join(info["dimension"])}), order=''\'F\'''), '
        elif info['type'] == 'type':
            s += f'np.asarray(<double[:np.prod(dims)]> {argname}' \
                ').reshape(dims, order=''\'F\'''), \\\n           '
        else:
            s += f'{argname}, '
        count += 1
    f.write(ss)
    if count > 0:
        s = s[:-2]
        f.write(s + '\n\n')


def write_pdaf_calls(
        filename: str, user_func_info,
        func_info: dict[str, dict[str, dict[str, str | list[str]]]]) -> None:
    """write the PDAF interface calls"""
    with open(filename, 'w', encoding="utf-8") as f:
        # MPI exception handling
        s: str = 'import sys\n'
        s += 'import numpy as np\n'
        s += 'from . cimport UserFunc as c__PDAFcython\n'
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
        s += '                sys.stderr.write('
        s += '"Uncaught exception was detected on rank {}. \\n".format(\n'
        s += '                    mpi4py.MPI.COMM_WORLD.Get_rank()))\n'
        s += '\n'
        s += '                print_exception(exctype, value, traceback)\n'
        s += '                sys.stderr.write("\\n")\n'
        s += '                sys.stderr.flush()\n'
        s += '            finally:\n'
        s += '                try:\n'
        s += '                    mpi4py.MPI.COMM_WORLD.Abort(1)\n'
        s += '                except Exception as e:\n'
        s += '                    sys.stderr.write('
        s += '"MPI Abort failed, this process will hang.\\n")\n'
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
            arg_list, array_dims = get_pyx_arg_list(
                subroutine_name, func_info[subroutine_name])
            write_func_def(
                f, subroutine_name, func_info[subroutine_name],
                arg_list)
            write_docstring(f, subroutine_name, user_func_info,
                            func_info[subroutine_name], arg_list)
            write_memory_view(f, subroutine_name, func_info[subroutine_name])
            if '_cb' not in subroutine_name:
                write_dims(f, func_info[subroutine_name], array_dims)
                write_user_cython(f, func_info[subroutine_name])
                write_return_def(f, func_info[subroutine_name])
            write_func_call(f, subroutine_name, func_info[subroutine_name])
            write_returns(f, func_info[subroutine_name])


if __name__ == '__main__':
    import get_interface_info
    import write_pxd
    user_func_info = get_interface_info.get_func_info(
        [os.path.join('..', 'src', 'fortran', 'U_PDAF_interface_c_binding.F90')])
    PDAF_func_info = get_interface_info.get_func_info(
        [os.path.join('..', 'src', 'fortran', 'PDAF_c_binding.F90'),
         os.path.join('..', 'src', 'fortran', 'PDAFomi_obs_c_binding.F90'),
         os.path.join('..', 'src', 'fortran', 'PDAFlocal_c_binding.F90'),])
    write_pxd.write_Pxd_file('PDAF.pxd', PDAF_func_info, user_func_info)
    write_pdaf_calls('PDAF.pyx', user_func_info, PDAF_func_info)
