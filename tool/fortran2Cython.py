"""This is a simple script for converting C binding
fortran routines in pyPDAF/fortran to Cython definition and
implementation files.
"""
# todo: type for pointer should be matching its corresponding fortran pointer type
# here double* is used because PDAF only uses double
conv = {'integer' : 'int', 'logical': 'bint', 'real': 'double', 'procedure':'void', 'character':'CFI_cdesc_t', 'type':'double*'}
pyconv = {'integer' : 'int', 'logical': 'bool', 'real': 'float', 'procedure':'void', 'character':'CFI_cdesc_t', 'type':'double*'}
cnpyconv = {'integer' : 'cnp.int32_t', 'real': 'cnp.float64_t'}
npyconv = {'integer' : 'np.intc', 'real': 'np.float64'}

def write_PDAF_calls(filename, func_info):
    with open(filename, 'w') as f:
        # MPI exception handling
        s = 'import pyPDAF.UserFunc as PDAFcython\n'
        s += 'cimport pyPDAF.UserFunc as c__PDAFcython\n\n'
        s += 'import numpy as np\n'
        s += 'cimport numpy as cnp\n'
        s += 'import sys\n'
        s += 'from traceback import print_exception\n'
        s += 'import mpi4py.MPI as MPI\n'
        s += '# Global error handler\n'
        s += 'def global_except_hook(exctype, value, traceback):\n'
        s += '    \n'
        s += '    try:\n'
        s += '        \n'
        s += '        sys.stderr.write("\\n*****************************************************\\n")\n'
        s += '        sys.stderr.write("Uncaught exception was detected on rank {}. \\n".format(\n'
        s += '            MPI.COMM_WORLD.Get_rank()))\n'
        s += '        \n'
        s += '        print_exception(exctype, value, traceback)\n'
        s += '        sys.stderr.write("*****************************************************\\n\\n\\n")\n'
        s += '        sys.stderr.write("\\n")\n'
        s += '        sys.stderr.write("Calling MPI_Abort() to shut down MPI processes...\\n")\n'
        s += '        sys.stderr.flush()\n'
        s += '    finally:\n'
        s += '        try:\n'
        s += '            MPI.COMM_WORLD.Abort(1)\n'
        s += '        except Exception as e:\n'
        s += '            sys.stderr.write("*****************************************************\\n")\n'
        s += '            sys.stderr.write("Sorry, we failed to stop MPI, this process will hang.\\n")\n'
        s += '            sys.stderr.write("*****************************************************\\n")\n'
        s += '            sys.stderr.flush()\n'
        s += '            raise e\n'
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

def get_pyx_arg_list(subroutine_name, arg_info):
    """Retrive dimension of input array arguments of a subroutine.
    If the dimension is an input argument of the subroutine, we do
    not write it explicitly as they can be obtained from
    numpy array shape. However, in the callback functions,
    these functions will be called by Fortran, all arguments must
    stay in the argument list.
    """
    ArrayDims = []
    for _, info in arg_info.items():
        # we need to keep the dimension of output arrays
        # when the dimension is only used by output arrays
        if info['intent'] == 'out': continue
        if info['array']:
            dims = info['dimension']
            for dim in dims:
                dimsize = dim.replace('-', '+').split('+')[0]
                if dimsize in arg_info:
                    ArrayDims.append(dimsize)
    ArrayDims = list(set(ArrayDims))

    arg_list = []
    for arg_name, info in arg_info.items():
        if arg_name not in ArrayDims or ('_cb' in subroutine_name):
            arg_list.append(arg_name)

    return arg_list, ArrayDims

def write_func_def(f, subroutine_name, arg_info, arg_list):
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
        if info['intent'] is not None:
            if ('in' not in info['intent']) \
                and ('_cb' not in funcname):
                continue

        if info['type'] == 'procedure':
            s += f'py{info["kind"][1:]}'
        elif info['array']:
            s += f'cnp.ndarray[{cnpyconv[info["type"]]}, ndim={len(info["dimension"])}] {argname}'
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

def write_docstring(f, subroutine_name, arg_info, arg_list):
    s = '    \"\"\"'
    s += f'See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/{subroutine_name[3:]} or PDAF source files \n\n'
    f.write(s)

    indent = ' '*4
    s = indent + 'Parameters\n'
    s += indent + '----------\n'
    count = 0
    for arg in arg_list:
        info = arg_info[arg]
        if info['intent'] is not None:
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
        if info['intent'] is None:
            continue
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

def write_memory_view(f, subroutine_name, arg_info):
    # convert np arrays to memoryview
    for argname, info in arg_info.items():
        if not info['array']:
            continue
        if ('in' not in info['intent']) \
            and ('_cb' not in subroutine_name):
            continue
        s = ' '*4 + f'cdef {conv[info["type"]]}[::1] '
        s += f'{argname}_view = np.array({argname}'
        if info['type'] == 'integer':
            s += ', dtype=np.intc'
        s += ').ravel(order=''\'F\''')\n'
        f.write(s)

def write_dims(f, arg_info, ArrayDims):
    """Write dimensions for array arguments.
    """
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

        # s = ' '*4
        s = ''
        sadjust = '    '
        for i, (dimname, adjust) in enumerate(zip(dimnames, adjusts)):
            s += ' '*4 + f'{dimname} = {argname}.shape[{i}]\n'
            for sign, nsign in zip(['+', '-'], ['-', '+']):
                if sign in adjust:
                    sadjust += f'{dimname} = {dimname} {adjust.replace(sign, nsign)}\n'

        # s += f' = {argname}.shape\n'
        f.write(s)
        if not sadjust.isspace():
            f.write(sadjust)
    f.write('\n' if count > 0 else '')

def write_user_Cython(f, arg_info):
    # convert python function to Cython function
    count = 0
    for _, info in arg_info.items():
        if info['type'] == 'procedure':
            s = ' '*4 + f'PDAFcython.py{info["kind"][1:]} = py{info["kind"][1:]}\n'
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
            s = indent + f'cdef {conv[info["type"]]} [::1] '
            s += f'{argname}_view'
            s += f' = np.zeros(({', '.join(info["dimension"])}), dtype={npyconv[info["type"]]}'
            s += ').ravel()\n'
        else:
            s = indent + f'cdef {conv[info["type"]]} {argname}\n'
        f.write(s)
        count += 1
    f.write('\n' if count > 0 else '')

def write_func_call(f, subroutine_name, arg_info):
    # call the actual subroutine
    s = ' '*4 + f'{subroutine_name} ('
    indent = ' '*len(s)
    for argname, info in arg_info.items():
        if info['array']:
            s += f'&{argname}_view[0]'
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
            ss = '    dims = np.asarray(dims_view)\n'
            continue

        if info['array']:
            s += f'np.asarray({argname}_view).reshape(({', '.join(info["dimension"])}), order=''\'F\'''), '
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


def write_user_def(f, subroutine_name, arg_info):
    s = f'def py{subroutine_name[1:]}('
    for argname, _ in arg_info.items():
        s += argname + ', '
    s = s[:-2] + '):\n'
    f.write(s)


def write_user_docstring(f, arg_info):
    s = '    \"\"\"'
    s += 'See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ \n\n'
    f.write(s)

    indent = ' '*4
    s = indent + 'Parameters\n'
    s += indent + '----------\n'
    count = 0
    for arg, info in arg_info.items():
        if info['array']:
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
        if info['intent'] is None:
            continue
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
        f.write(s+'\n')

    s = indent + 'Raises\n'
    s += indent + '------\n'
    s += indent + 'RuntimeError\n'
    s += 2*indent + 'No user-supplied function\n'
    f.write(s)

    s = '    \"\"\"\n'
    f.write(s)


def write_C_user_def(f, subroutine_name, arg_info):
    s = f'cdef void {subroutine_name} ('
    indent = ' '*len(s)
    f.write(s)
    for arg, info in arg_info.items():
        s = f'{conv[info["type"]]}* {arg}'
        s += '' if arg == list(arg_info)[-1] else ', '
        f.write(s)
    s = ') noexcept:'
    f.write(s+'\n')


def write_array_conversion(f, arg_info):
    count = 0
    for argname, info in arg_info.items():
        if not info['array']:
            continue
        dims = ''
        for d in info['dimension']:
            have_sign = False
            for sign in ['+', '-']:
                if sign in d:
                    have_sign = True
                    dimname, adjust = d.split(sign)
                    dims += dimname+'[0]' + sign + ''.join(adjust)
            if not have_sign:
                dims += d+'[0]'
            dims += ', '
        dims = dims[:-2]
        indent = ' '*4
        if len(set(dims.replace(" ", "").split(','))) == 1 and dims.replace(" ", "").split(',')[0] == 'dim_ens[0]-1':
            s = indent + 'if dim_ens[0] > 1:\n'
            s += indent + ' '*4 + f'{argname}_np = np.asarray(<{conv[info["type"]]}[:np.prod(({dims}))]> \n{arglen+indent}{argname})'
            s += f'.reshape(({dims}), order=''\'F\''')\n'
            s += indent + 'else:\n'
            s += indent + ' '*4 + f'{argname}_np = None\n'
        else:
            arglen = ' '*len(f'{argname}_np = np.asarray(')
            s = indent + f'{argname}_np = np.asarray(<{conv[info["type"]]}[:np.prod(({dims}))]> \n{arglen+indent}{argname})'
            s += f'.reshape(({dims}), order=''\'F\''')\n'
        f.write(s)
        count += 1
    f.write('\n' if count > 0 else '')


def write_C_U_calls(f, subroutine_name, arg_info):
    s = ' '*4
    count = 0
    for arg, info in arg_info.items():
        if info['intent'] is None:
            continue
        if 'out' not in info['intent']:
            continue

        if info['array']:
            s += arg +'_np_tmp, '
        else:
            s += f'{arg}[0], '
        count += 1
    if count > 0:
        s = s[:-2]
        s += ' = '

    s += f'py{subroutine_name[1:]}('
    count = 0
    for arg, info in arg_info.items():
        if info['array']:
            s += arg +'_np, '
        else:
            s += f'{arg}[0], '
        count += 1
    s = s[:-2] if count > 0 else s
    s += ')\n'
    f.write(s)


def convertToArrays(f, subroutine_name, arg_info):
    indent = ' '*4
    s = ''
    for argname, info in arg_info.items():
        if info['intent'] is None:
            continue
        if 'out' not in info['intent']:
            continue

        if info['array']:
            if len(set(info['dimension'])) == 1 and info['dimension'][0] == 'dim_ens-1':
                s += indent + f'cdef double[::1] {argname}_view\n'
                s += indent + 'if dim_ens[0] > 1:\n'
                s += indent + ' '*4 + f'{argname}_np[:] = {argname}_np_tmp[:]\n'
                s += indent + ' '*4 + f'{argname}_view = {argname}_np.ravel(order=''\'F\''')\n'
                s += indent + ' '*4 + f'assert {argname} == &{argname}_view[0], '\
                               f'\'reference (memory address) of {argname} has changed in {subroutine_name}.\'\n'
            else:
                s += indent + f'{argname}_np[:] = {argname}_np_tmp[:]\n'
                s += indent + f'cdef {conv[info["type"]]}[::1] '
                s += f'{argname}_view = {argname}_np.ravel(order=''\'F\''')\n'
                s += indent + f'assert {argname} == &{argname}_view[0], '\
                               f'\'reference (memory address) of {argname} has changed in {subroutine_name}.\'\n'
    f.write(s)


def writeUserCalls(filename, func_info):
    with open(filename, 'w') as f:
        s = 'import numpy as np\n'
        s += 'import mpi4py.MPI as MPI\n'
        s += 'import sys\n\n'
        s += 'from traceback import print_exception\n'
        s += '# Global error handler\n'
        s += 'def global_except_hook(exctype, value, traceback):\n'
        s += '    \n'
        s += '    try:\n'
        s += '        \n'
        s += '        sys.stderr.write("\\n*****************************************************\\n")\n'
        s += '        sys.stderr.write("Uncaught exception was detected on rank {}. \\n".format(\n'
        s += '            MPI.COMM_WORLD.Get_rank()))\n'
        s += '        \n'
        s += '        print_exception(exctype, value, traceback)\n'
        s += '        sys.stderr.write("*****************************************************\\n\\n\\n")\n'
        s += '        sys.stderr.write("\\n")\n'
        s += '        sys.stderr.write("Calling MPI_Abort() to shut down MPI processes...\\n")\n'
        s += '        sys.stderr.flush()\n'
        s += '    finally:\n'
        s += '        try:\n'
        s += '            MPI.COMM_WORLD.Abort(1)\n'
        s += '        except Exception as e:\n'
        s += '            sys.stderr.write("*****************************************************\\n")\n'
        s += '            sys.stderr.write("Sorry, we failed to stop MPI, this process will hang.\\n")\n'
        s += '            sys.stderr.write("*****************************************************\\n")\n'
        s += '            sys.stderr.flush()\n'
        s += '            raise e\n'
        s += '\n'
        s += 'sys.excepthook = global_except_hook\n'
        s += '\n\n'
        f.write(s)
        for subroutine_name in func_info:
            write_user_def(f, subroutine_name.lower(), func_info[subroutine_name])
            write_user_docstring(f, func_info[subroutine_name])
            s = '    raise RuntimeError(\'...Wrong '
            s += f'py{subroutine_name[1:]} is called!!!...\')'
            f.write(s+'\n\n\n')

        for subroutine_name in func_info:
            write_C_user_def(f, subroutine_name.lower(), func_info[subroutine_name])
            write_array_conversion(f, func_info[subroutine_name])
            write_C_U_calls(f, subroutine_name.lower(), func_info[subroutine_name])
            f.write('\n')
            convertToArrays(f, subroutine_name.lower(), func_info[subroutine_name])
            f.write('\n\n')

if __name__ == '__main__':
    import get_interface_info
    import write_pxd
    user_func_info = get_interface_info.get_func_info(['../pyPDAF/fortran/U_PDAF_interface_c_binding.F90'])
    PDAF_func_info = get_interface_info.get_func_info(['../pyPDAF/fortran/PDAF_c_binding.F90', '../pyPDAF/fortran/PDAFomi_obs_c_binding.F90'])
    write_pxd.write_Pxd_file('PDAF.pxd', PDAF_func_info, user_func_info)
    write_pxd.write_Pxd_file('UserFunc.pxd', user_func_info, user_func_info)
    write_PDAF_calls('PDAF.pyx', PDAF_func_info)
    writeUserCalls('UserFunc.pyx', user_func_info)






