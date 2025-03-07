import os
import typing
import re

pyconv: dict[str, str] = {'integer': 'int',
                          'logical': 'bool', 'real': 'float', 'character': 'str'}
conv = {'integer': 'int', 'logical': 'bint',
        'real': 'double', 'character': 'CFI_cdesc_t'}
special_functions = ['c__init_ens_pdaf', 'c__prepoststep_pdaf']


def extract_dimension_name(s: str) -> str | None:
    """extract dimension names from dimension
    """
    # Match one or more word characters (letters, digits, and underscores), but exclude digits
    match = re.search(r'[a-zA-Z_]+', s)
    if match:
        return match.group()
    return None


def write_user_def(f: typing.TextIO, subroutine_name: str,
                   arg_info: dict[str, dict[str, str | list[str]]]) -> None:
    """write the interface of Python user-defined function to the file

    Parameters
    ----------
    f : TextIO
        the file to write the interface
    subroutine_name : str
        the name of the subroutine
    arg_info : dict
        a dictionary containing information about the function arguments
    """
    s = f'def py{subroutine_name[1:]}('
    for argname, one_arg_info in arg_info.items():
        assert isinstance(one_arg_info["type"], str), \
            f'Unknown arg_info variable type: {one_arg_info["type"]}'
        s += f'{argname}:{pyconv[one_arg_info["type"]]}, '
    s = s[:-2] + '):\n'
    f.write(s)


def write_user_docstring(f: typing.TextIO,
                         arg_info: dict[str, dict[str, str | list[str]]]
                         ) -> None:
    """write the docstring of the Python user-defined function to the file

    Parameters
    ----------
    f : TextIO
        the file to write the docstring
    arg_info : dict
        a dictionary containing information about the function arguments
    """
    s = '    \"\"\"'
    # We will need to write better docstring for each function
    s += 'See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ \n\n'
    f.write(s)

    indent = ' '*4
    s = indent + 'Parameters\n'
    s += indent + '----------\n'
    count = 0
    for arg, info in arg_info.items():
        assert isinstance(info["type"], str), \
            f'Unknown arg_info variable type: {info["type"]}'
        assert isinstance(info["comment"], str), \
            f'Unknown comment variable type: {info["comment"]}'
        if len(info['dimension']) > 0:
            s += indent+arg + f' : ndarray[{pyconv[info["type"]]}]\n'
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
        if 'out' not in info['intent']:
            continue

        assert isinstance(info["type"], str), \
            f'Unknown arg_info variable type: {
            info["type"]}'
        assert isinstance(info["comment"], str), \
            f'Unknown comment variable type: {
            info["comment"]}'

        if len(info['dimension']) > 0:
            s += indent+arg + f' : ndarray[{pyconv[info["type"]]}]\n'
        elif info['type'] == 'type':
            s += indent+arg + ' : ndarray[float]\n'
        else:
            s += indent + f'{arg} : {pyconv[info["type"]]}\n'
        s += 2*indent + info['comment'] + '\n'
        count += 1

    if count > 0:
        f.write(s+'\n')

    s = indent + 'Raises\n'
    s += indent + '------\n'
    s += indent + 'RuntimeError\n'
    s += 2*indent + 'No user-supplied function\n'
    f.write(s)

    s = '    \"\"\"\n'
    f.write(s)


def write_c_user_def(
        f: typing.TextIO, subroutine_name: str,
        arg_info: dict[str, dict[str, str | list[str]]]) -> None:
    """write the definition of the C user-defined function to the file

    Parameters
    ----------
    f : TextIO
        the file to write the definition
    subroutine_name : str
        the name of the subroutine
    arg_info : dict
        a dictionary containing information about the function arguments
    """
    s = f'cdef void {subroutine_name} ('
    f.write(s)
    for arg, info in arg_info.items():
        assert isinstance(info["type"], str), \
            f'Unknown arg_info variable type: {info["type"]}'
        s = f'{conv[info["type"]]}* {arg}'
        s += '' if arg == list(arg_info)[-1] else ', '
        f.write(s)
    s = ') noexcept with gil:'
    f.write(s+'\n\n')


def write_array_conversion(f: typing.TextIO,
                           arg_info: dict[str, dict[str, str | list[str]]]
                           ) -> dict[str, list[str]]:
    """Convert input C pointers to numpy arrays and get the array size for each arguments

    Parameters
    ----------
    f : TextIO
        the file to write the conversion
    arg_info : dict
        a dictionary containing information about the function arguments
    """
    count: int = 0
    dim_ens_assert: bool = False
    indent: str = ' '*4
    arg_dims: dict[str, list[str]] = {}
    for argname, info in arg_info.items():
        if len(info['dimension']) == 0:
            continue

        # replace the dimension name by C input arguments
        dims: list[str] = []
        for d in info['dimension']:
            # some dimension delcaration contains calculations
            dimname = extract_dimension_name(d)
            assert dimname is not None, f'Unknown dimension name: {d}'
            dim: str
            if dimname != d:
                dim = f'({d})'
                # dim_ens[0] - 1 should be handled elsewhere.
                # we nevertheless add a check for dim_ens[0]
                if 'dim_ens-1' in d.replace(' ', '') and not dim_ens_assert:
                    s = '\n' + indent + \
                        'assert dim_ens[0] > 1, '\
                        '"ensemble size must be > 1 for ensemble filters."\n\n'
                    f.write(s)
                    dim_ens_assert = True
            else:
                dim = d
            dims.append(dim.replace(dimname, f'{dimname}[0]'))

        # we only need to define the array size once
        if dims in arg_dims.values():
            for argname_exist, dims_exist in arg_dims.items():
                if dims_exist == dims:
                    break
            arg_dims[argname] = [f'__USE_ARG_{argname_exist}', ]
        else:
            arg_dims[argname] = dims

        # remove brackets from dims
        dims = [d.replace('(', '') for d in dims]
        dims = [d.replace(')', '') for d in dims]
        # convert the C pointer to numpy array
        assert isinstance(info["type"], str), \
            f'Unknown arg_info variable type: {info["type"]}'
        if len(info['dimension']) == 1:
            s = indent + f'cdef {conv[info["type"]]}[::1]' \
                f' {argname}_np = np.asarray(' \
                f'<{conv[info["type"]]}[:{dims[0]}]> {argname})\n'
            f.write(s)
        else:
            s = indent + f'cdef {conv[info["type"]]}[::1'
            for _ in range(len(info['dimension']) - 1):
                s += ',:'
            s += f'] {argname}_np = np.asarray(' \
                f'<{conv[info["type"]]}[:{dims[0]}:1,'
            for dim in dims[1:]:
                s += f':{dim},'
            s = s[:-1]
            s += f']> {argname}, order=''\'F\''')\n'

            f.write(s)
        count += 1
    f.write('\n' if count > 0 else '')

    return arg_dims


def write_c_u_calls(f: typing.TextIO, subroutine_name: str,
                    arg_info: dict[str, dict[str, str | list[str]]]) -> None:
    """write the Python user-defined function calls to the file

    Parameters
    ----------
    f : TextIO
        the file to write the calls
    subroutine_name : str
        the name of the subroutine
    arg_info : dict
        a dictionary containing information about the function arguments
    """
    indent: str = ' '*4
    count: int = 0
    s: str = indent
    # write the function returns
    for arg, info in arg_info.items():
        if 'out' not in info['intent']:
            continue

        if len(info['dimension']) > 0:
            s += arg + '_np, '
        else:
            s += f'{arg}[0], '

        count += 1

    if count > 0:
        s = s[:-2]
        s += ' = '
    # call Python function using C pointers
    s += f'(<object>{subroutine_name[3:]})('
    count = 0
    # write the input function arguments
    for arg, info in arg_info.items():
        if len(info['dimension']) > 0:
            s += arg + '_np.base, '
        else:
            s += f'{arg}[0], '
        count += 1
    s = s[:-2] if count > 0 else s
    s += ')\n'
    f.write(s)


def check_output_array_memory(
        f: typing.TextIO, subroutine_name: str, arg_dims: dict
        [str, list[str]],
        arg_info: dict[str, dict[str, str | list[str]]]) -> None:
    """Check the user-supplied output arrays does not change the address of the input C pointers

    Parameters
    ----------
    f : TextIO
        the file to write the checks
    subroutine_name : str
        the name of the subroutine
    arg_dims : dict
        a dictionary containing the array size of the function arguments
    arg_info : dict
        a dictionary containing information about the function arguments
    """
    indent: str = ' '*4
    for argname, info in arg_info.items():
        if 'out' not in info['intent']:
            continue

        if len(info['dimension']) == 0:
            continue

        assert isinstance(info["type"], str), \
            f'Unknown arg_info variable type: {info["type"]}'
        s = '\n'
        s += indent + f'cdef {conv[info["type"]]}[::1'
        for i in range(len(info['dimension']) - 1):
            s += ',:'
        s += f'] {argname}_new\n'

        # check if the memory address of the numpy array is different from the input C pointer
        dim_0: str = '0,'*len(info['dimension'])
        dim_0 = dim_0[:-1]
        s += indent + f'if {argname} != &{argname}_np[{dim_0}]:\n'

        # if it is not, we assign values of the output array to
        # the input C pointer and raise a warning
        if len(info['dimension']) == 1:
            s += indent * 2 + f'{argname}_new = '\
                f'np.asarray(<{conv[info["type"]]} '\
                f'[: {info["dimension"][0]} [0]]> {argname})\n'
        else:
            s += indent*2 + f'{argname}_new = np.asarray(<{
                conv[info["type"]]}['
            argname_exist = argname if '_USE_ARG_' not in \
                arg_dims[argname][0] else arg_dims[argname][0].replace(
                    '__USE_ARG_', '')
            s += f':{arg_dims[argname_exist][0]}:1,'
            for dim in arg_dims[argname_exist][1:]:
                s += f':{dim},'
            s = s[:-1]
            s += f']> {argname}, order=''\'F\''')\n'

        s += indent*2 + f'{argname}_new[...] = {argname}_np\n'
        s += indent*2 + f'warnings.warn("The memory address of {
            argname} is changed in {subroutine_name}." \n '
        s += indent*2 + '"The values are copied to the original ' \
            'Fortran array, and can slow-down the system.", RuntimeWarning)\n'
        f.write(s)


def write_special_functions(f: typing.TextIO, subroutine_name: str) -> None:
    f_code: typing.TextIO = open(os.path.join('special_functions',
                                              f'{subroutine_name}.pyx'),
                                 'r'
                                 )
    for line in f_code:
        f.write(line)
    f_code.close()


def writeUserCalls(
        filename: str,
        func_info: dict[str, dict[str, dict[str, str | list[str]]]]) -> None:
    with open(filename, 'w') as f:
        s = 'import numpy as np\n'
        s += 'import warnings\n'
        s += '\n\n'
        f.write(s)
        # for subroutine_name in func_info:
        #     write_user_def(f, subroutine_name.lower(), func_info[subroutine_name])
        #     # write_user_docstring(f, func_info[subroutine_name])
        #     s = '    raise RuntimeError(\'...Wrong '
        #     s += f'py{subroutine_name[1:]} is called!!!...\')'
        #     f.write(s+'\n\n\n')

        for subroutine_name in func_info:
            if subroutine_name in special_functions:
                write_special_functions(f, subroutine_name)
                f.write('\n\n')
                continue

            write_c_user_def(f, subroutine_name,
                             func_info[subroutine_name])
            arg_dims = write_array_conversion(
                f, func_info[subroutine_name])
            write_c_u_calls(f, subroutine_name,
                            func_info[subroutine_name])
            f.write('\n')
            check_output_array_memory(f, subroutine_name, arg_dims,
                                      func_info[subroutine_name])
            f.write('\n\n')


if __name__ == '__main__':
    import get_interface_info
    import write_pxd
    user_func_info = get_interface_info.get_func_info(
        ['../src/fortran/U_PDAF_interface_c_binding.F90'])
    write_pxd.write_Pxd_file('UserFunc.pxd', user_func_info, user_func_info)
    writeUserCalls('UserFunc.pyx', user_func_info)
