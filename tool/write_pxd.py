conv = {'integer' : 'int', 'logical': 'bint', 'real': 'double', 'procedure':'void', 'character':'CFI_cdesc_t', 'type':'double*'}

def write_Pxd_file(filename, func_info, user_func_info):
    """
    A function that writes the content of the func_info dictionary to a .pxd file with specific formatting.
    Parameters:
        - filename: the name of the .pxd file to be written
        - func_info: a dictionary containing information about the function
        - user_func_info: a dictionary containing additional user-defined function information
    """
    with open(filename, 'w') as f:
        for routine_name in func_info:
            if func_info == user_func_info:
                s = f'cdef void {routine_name.lower()} ('
            else:
                s = f'cdef extern void {routine_name.lower()} ('
            indent = ' '*len(s)
            f.write(s)
            for argname, arginfo in func_info[routine_name].items():
                if arginfo["type"] == 'procedure':
                    user_arg_info = [user_func_info[user_func_name] 
                            for user_func_name in user_func_info if user_func_name.lower() == arginfo['kind']][0]
                    s = f'void (*{arginfo['kind']})('
                    indent2 = indent + ' '*len(s)
                    for _, u_info in user_arg_info.items():
                        s += f'{conv[u_info["type"]]}*'
                        s += f',\n{indent2}'

                    n = len(f',\n{indent2}')
                    s = s[:-n] + f'\n{indent2[:-1]})'
                else:
                    s = f'{conv[arginfo["type"]]}* {argname}'

                s += f'\n{indent[:-1]}' if argname == list(func_info[routine_name])[-1] else f',\n{indent}'
                f.write(s)

            s = ') noexcept;'
            f.write(s+'\n')
