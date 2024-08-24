conv:dict[str, str] = {'integer' : 'int', 'logical': 'bint', 'real': 'double', 'procedure':'void', 'character':'CFI_cdesc_t', 'type':'double*'}

def write_Pxd_file(filename:str,
                   func_info:dict[str, dict[str, dict[str, str|bool|None|list[str]]]],
                   user_func_info:dict[str, dict[str, dict[str, str|bool|None|list[str]]]]) -> None:
    """
    A function that writes the content of the func_info dictionary to a .pxd file with specific formatting.

    Parameters
    ----------
    filename: str
        the name of the .pxd file to be written
    func_info : dict
        a dictionary containing information about the function
    user_func_info : dict
        a dictionary containing additional user-defined function information
    """
    s : str
    indent : str
    routine_name : str
    argname : str
    with open(filename, 'w') as f:
        for routine_name in func_info:
            # Declare pointer to Python user-supplied callback functions
            if func_info == user_func_info:
                s = f'cdef void*  {routine_name[3:]} = NULL;\n'
                f.write(s)

            if func_info == user_func_info:
                # Decalre Cython user-supplied callback functions
                s = f'cdef void {routine_name} ('
            else:
                # Decalre Fortran subroutine interfaces
                s = f'cdef extern void {routine_name.lower()} ('
            indent = ' '*len(s)
            f.write(s)

            # write the function arguments
            for argname, arginfo in func_info[routine_name].items():
                if arginfo["type"] == 'procedure':
                    # write the user-supplied function interface
                    user_arg_info = [user_func_info[user_func_name]
                            for user_func_name in user_func_info if user_func_name == arginfo['kind']][0]
                    s = f'void (*{arginfo['kind']})('
                    indent2 = indent + ' '*len(s)
                    for _, u_info in user_arg_info.items():
                        assert type(u_info["type"]) is str, f'Unknown u_info variable type: {u_info["type"]}'
                        s += f'{conv[u_info["type"]]}*'
                        s += f',\n{indent2}'

                    n = len(f',\n{indent2}')
                    s = s[:-n] + f'\n{indent2[:-1]})'
                else:
                    # simply write the argument type and name
                    assert type(arginfo["type"]) is str, f'Unknown arginfo variable type: {arginfo["type"]}'
                    s = f'{conv[arginfo["type"]]}* {argname}'

                s += f'\n{indent[:-1]}' if argname == list(func_info[routine_name])[-1] else f',\n{indent}'
                f.write(s)

            s = ') noexcept;'
            f.write(s+'\n')

