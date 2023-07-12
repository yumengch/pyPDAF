"""This is a simple script for converting C binding
fortran routines in pyPDAF/fortran to Cython definition and
implementation files.
"""

# todo: type for pointer should be matching its corresponding fortran pointer type
# here double* is used because PDAF only uses double
conv = {'integer' : 'int', 'logical': 'bint', 'real': 'double', 'procedure':'void', 'character':'CFI_cdesc_t', 'type':'double*'}
pyconv = {'integer' : 'int', 'logical': 'bool', 'real': 'float', 'procedure':'void', 'character':'CFI_cdesc_t', 'type':'double*'}

# todo: array dimension can only handle one + or -

def mergeLine(f, line):
    # remove all white lines
    skip_condition = line.isspace()
    line = line.strip().replace('\n', '')
    # ignore lines starting with use as they are module import
    skip_condition = skip_condition or line[:3] == 'use'
    # ignore implicit none lines
    skip_condition = skip_condition or line.replace(' ', '') == 'implicitnone'
    if skip_condition:
        return ' '

    # merge all concatenated lines into one string
    while '&' in line.split('!')[0].lower():
        line = line.replace('&', '')
        line = line + f.readline().lower().strip().replace('\n', '')
    return line


def getArgs(string):
    args = []
    left = []
    s = ''
    for c in string:
        if c == '(':
            left.append(c)
        if c == ',' and left == []:
            args.append(s.strip())
            s = ''
            continue
        if c == ')':
            left.pop()
        s += c
    args.append(s.strip())
    return args


def getArgInfo(ArgInfo, line, comment):
    declaration = line.replace(' ', '').split('::')
    args = getArgs(declaration[1])
    attrs = declaration[0].replace(')', '(').split('(')

    for arg in args:
        varname = arg.split('(')[0]
        if varname not in ArgInfo:
            continue

        ArgInfo[varname]['type'] =  attrs[0]
        ArgInfo[varname]['kind'] = attrs[1].split('kind=')[-1]
        ArgInfo[varname]['intent'] = attrs[attrs.index(',intent') + 1] if ',intent' in attrs else None
        ArgInfo[varname]['array'] = ',dimension' in attrs or '(' in arg
        ArgInfo[varname]['size'] = attrs[attrs.index(',dimension')+1] if ',dimension' in attrs else None
        ArgInfo[varname]['size'] = arg[arg.index('(')+1:-1] if '(' in arg else ArgInfo[arg]['size']
        ArgInfo[varname]['comment'] = comment.strip()

    return ArgInfo


def getArguments(lines):
    # get the arg list and subroutine name
    line = lines.pop(0)
    arg_list = line.replace(' ', '').replace(')', '(').split('(')[1].split(',')
    ArgInfo = {'name':line.split('(')[0].split()[-1]}
    for arg in arg_list:
        if arg == '':
            continue
        ArgInfo[arg] = {}

    # get information in the arg list for C declaration
    comment = ''
    for line in lines:
        # todo: assume all variable declaration is separated by ::
        # this is not the case for older fortran standards
        if '!' in line:
            comment = ''.join([comment, line.strip().replace('!', '\n        ')])

        if '::' not in line:
            continue

        assert '!' not in line, 'if ! is in  the same line as ::, '\
                                'it doesn\'t fit the assumption that '\
                                'comment and arguments declaration are separated.'
        ArgInfo = getArgInfo(ArgInfo, line, comment)

        comment = ''
    return ArgInfo


def getFuncInfo(filename):
    with open(filename, 'r') as f:
        line = mergeLine(f, f.readline().lower())
        lines = []
        # using start to exclude unwanted content
        start = False
        while line:
            line = mergeLine(f, f.readline().lower())

            if not line.isspace() and start:
                lines.append(line)

            if line == 'contains' or line == 'abstract interface':
                start = True

            if line == 'end interface':
                start = False

    FuncInfo = []
    subroutine = []
    for line in lines:
        subroutine.append(line)

        if 'endsubroutine' in line.replace(' ', ''):
            assert 'subroutine' in subroutine[0], f'{subroutine[0]}'
            FuncInfo.append(getArguments(subroutine))
            subroutine = []

    return FuncInfo


def writePxdFile(filename, UserFuncInfo, PDAFInfo):
    with open(filename, 'w') as f:
        for routine in PDAFInfo:
            name = routine.pop('name')
            s = f'cdef extern void {name} ('
            indent = ' '*len(s)
            f.write(s)
            for arg, info in routine.items():
                if info["type"] == 'procedure':
                    info = [UserFunc for UserFunc in UserFuncInfo if UserFunc['name'] == info['kind']][0]
                    u_name = info['name']
                    s = f'void (*{u_name})('
                    indent2 = indent + ' '*len(s)
                    for u_arg, u_info in info.items():
                        if u_arg == 'name':
                            continue
                        s += f'{conv[u_info["type"]]}*'
                        s += f',\n{indent2}'

                    n = len(f',\n{indent2}')
                    s = s[:-n] + f'\n{indent2[:-1]})'
                else:
                    s = f'{conv[info["type"]]}* {arg}'

                s += f'\n{indent[:-1]}' if arg == list(routine)[-1] else f',\n{indent}'
                f.write(s)

            s = ');'
            f.write(s+'\n')
            routine['name'] = name


def getPyxArgList(name, routine):
    ArrayDims = []
    for arg, info in routine.items():
        if arg == 'name':
            continue
        if info['intent'] is not None:
            if ('in' not in info['intent']) \
                and ('_cb' not in name):
                continue
        if info['array']:
            dims = info['size'].split(',')
            for dim in dims:
                dimsize = dim.replace('-', '+').split('+')[0]
                if dimsize in routine:
                    ArrayDims.append(dimsize)
    ArrayDims = list(set(ArrayDims))

    arg_list = []
    for arg, info in routine.items():
        if arg == 'name':
            continue
        if info['intent'] is not None:
            if ('in' not in info['intent']) \
                and ('_cb' not in name):
                continue
        if arg not in ArrayDims \
            or ('_cb' in name):
            arg_list.append(arg)

    return arg_list, ArrayDims


def writeFuncDef(f, name, routine, arg_list):
    # define function
    if name[3:11] == 'pdafomi_':
        funcname = name[11:]
    elif name[3:8] == 'pdaf_':
        funcname = name[8:]
    else:
        funcname = name[3:]

    s = f'def {funcname} ('
    indent = ' '*len(s)
    count = 0
    for arg in arg_list:
        info = routine[arg]
        if info['intent'] is not None:
            if ('in' not in info['intent']) \
                and ('_cb' not in funcname):
                continue

        if info['type'] == 'procedure':
            s += f'py{info["kind"][1:]}'
        elif info['array']:
            s += arg
        else:
            s += f'{conv[info["type"]]} {arg}'

        s += f',\n{indent}'
        count += 1

    if count > 0:
        n = len(f',\n{indent}')
        s = s[:-n] + f'\n{indent[:-1]}):\n'
    else:
        s += '):\n'
    f.write(s)


def writeDocString(f, name, routine, arg_list):
    s = '    \"\"\"'
    s += 'See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ \n\n'
    f.write(s)

    indent = ' '*4
    s = indent + 'Parameters\n'
    s += indent + '----------\n'
    count = 0
    for arg in arg_list:
        info = routine[arg]
        if info['intent'] is not None:
            if 'in' not in info['intent'] \
                and ('_cb' not in name):
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
    for arg, info in routine.items():
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


def writeMemoryView(f, name, routine):
    # convert np arrays to memoryview
    for arg, info in routine.items():
        if not info['array']:
            continue
        if ('in' not in info['intent']) \
            and ('_cb' not in name):
            continue
        s = ' '*4 + f'cdef {conv[info["type"]]}[::1] '
        s += f'{arg}_view = np.array({arg}'
        if info['type'] == 'integer':
            s += ', dtype=np.intc'
        s += ').ravel(order=''\'F\''')\n'
        f.write(s)


def writeDims(f, routine, ArrayDims):

    s = '    cdef int '
    count = 0
    for dim in ArrayDims:
        s += f'{dim}, '
        count += 1
    s = s[:-2]+'\n' if count > 0 else ''
    f.write(s)

    d = sorted(routine.items(), key=lambda item: -1 if not item[1]['array'] else len(item[1]['size'].split(',')), reverse=True)
    ndim = {k: v['size'].split(',') for k, v in d if v['array'] and 'in' in v['intent']}

    for arg, dims in ndim.items():
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

        s = '    '
        sadjust = '    '
        for dimname, adjust in zip(dimnames, adjusts):
            s += f'{dimname}, '
            for sign, nsign in zip(['+', '-'], ['-', '+']):
                if sign in adjust:
                    sadjust += f'{dimname} = {dimname} {adjust.replace(sign, nsign)}\n'

        s += f' = {arg}.shape\n'
        f.write(s)
        if not sadjust.isspace():
            f.write(sadjust)
    f.write('\n' if count > 0 else '')


def writeUserCython(f, routine):
    # convert python function to Cython function
    count = 0
    for arg, info in routine.items():
        if info['type'] == 'procedure':
            s = ' '*4 + f'PDAFcython.py{info["kind"][1:]} = py{info["kind"][1:]}\n'
            f.write(s)
        count += 1
    f.write('\n' if count > 0 else '')


def writeReturnDef(f, routine):
    indent = ' '*4
    count = 0
    # define return variables
    for arg, info in routine.items():
        if info['intent'] != 'out':
            continue
        if info['array']:
            s = indent + f'cdef {conv[info["type"]]} [::1] '
            s += f'{arg}_view'
            s += f' = np.zeros(({info["size"]})'
            if info['type'] == 'integer':
                s += ', dtype=np.intc'
            s += ').ravel()\n'
        else:
            s = indent + f'cdef {conv[info["type"]]} {arg}\n'
        f.write(s)
        count += 1
    f.write('\n' if count > 0 else '')


def writeFuncCall(f, name, routine):
    # call the actual subroutine
    s = ' '*4 + f'{name} ('
    indent = ' '*len(s)
    for arg, info in routine.items():
        if info['array']:
            s += f'&{arg}_view[0]'
        elif info['type'] == 'procedure':
            s += f'c__PDAFcython.{info["kind"]}'
        else:
            s += f'&{arg}'
        s += ',\n' + indent

    if len(routine) > 0:
        n = len(f',\n{indent}')
        s = s[:-n] + f'\n{indent[:-1]})\n'
    else:
        s += ')\n'
    f.write(s + '\n')


def writeReturns(f, routine):
    s = '    return '
    ss = ''
    count = 0
    for arg, info in routine.items():
        if info['intent'] is None:
            continue
        if 'out' not in info['intent']:
            continue
        # todo: hard coded pointer to np array conversion
        if arg == 'dims':
            ss = '    dims = np.asarray(dims_view)\n'
            continue

        if info['array']:
            s += f'np.asarray({arg}_view).reshape(({info["size"]}), order=''\'F\'''), '
        elif info['type'] == 'type':
            s += f'np.asarray(<double[:np.prod(dims)]> {arg}).reshape(dims, order=''\'F\'''), \\\n           '
        else:
            s += f'{arg}, '
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


def writePDAFcalls(filename, PDAFInfo):
    with open(filename, 'w') as f:
        s = 'import pyPDAF.UserFunc as PDAFcython\n'
        s += 'cimport pyPDAF.UserFunc as c__PDAFcython\n\n'
        s += 'import numpy as np\n'
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
        for routine in PDAFInfo:
            name = routine.pop('name')
            arg_list, ArrayDims = getPyxArgList(name, routine)
            writeFuncDef(f, name, routine, arg_list)
            writeDocString(f, name, routine, arg_list)
            writeMemoryView(f, name, routine)
            if '_cb' not in name:
                writeDims(f, routine, ArrayDims)
                writeUserCython(f, routine)
                writeReturnDef(f, routine)
            writeFuncCall(f, name, routine)
            writeReturns(f, routine)
            routine['name'] = name


def writeUserPxdFile(filename, UserFuncInfo):
    with open(filename, 'w') as f:
        for routine in UserFuncInfo:
            name = routine.pop('name')
            s = f'cdef void {name} ('
            indent = ' '*len(s)
            f.write(s)
            for arg, info in routine.items():
                s = f'{conv[info["type"]]}* {arg}'
                s += f'\n{indent[:-1]}' if arg == list(routine)[-1] else f',\n{indent}'
                f.write(s)

            s = ');'
            f.write(s+'\n')
            routine['name'] = name


def writeUserDef(f, name, routine):
    s = f'def py{name[1:]}('
    for arg, info in routine.items():
        s += arg + ', '
    s = s[:-2] + '):\n'
    f.write(s)


def writeUserDocString(f, routine):
    s = '    \"\"\"'
    s += 'See detailed explanation of the routine in https://pdaf.awi.de/trac/wiki/ \n\n'
    f.write(s)

    indent = ' '*4
    s = indent + 'Parameters\n'
    s += indent + '----------\n'
    count = 0
    for arg, info in routine.items():
        # if info['intent'] is not None:
        #     if 'in' not in info['intent']:
        #         continue
        if info['array']:
            s += indent+arg +f' : ndarray[{pyconv[info["type"]]}]\n'
        else:
            s += indent + f'{arg} : {pyconv[info["type"]]}\n'
        s += 2*indent + info['comment'] +'\n'
        if info['array']:
            s += 2*indent + f'shape is ({info["size"]})\n'
        count += 1
    if count > 0:
        f.write(s)

    s = '\n'
    s += indent + 'Returns\n'
    s += indent + '-------\n'
    count = 0
    for arg, info in routine.items():
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


def writeCUserDef(f, name, routine):
    s = f'cdef void {name} ('
    indent = ' '*len(s)
    f.write(s)
    for arg, info in routine.items():
        s = f'{conv[info["type"]]}* {arg}'
        s += '' if arg == list(routine)[-1] else ', '
        f.write(s)
    s = '):'
    f.write(s+'\n')


def writeArrayConversion(f, routine):
    count = 0
    for arg, info in routine.items():
        if not info['array']:
            continue

        dims = ''
        for d in info['size'].split(','):
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
            s += indent + ' '*4 + f'{arg}_np = np.asarray(<{conv[info["type"]]}[:np.prod(({dims}))]> \n{arglen+indent}{arg})'
            s += f'.reshape(({dims}), order=''\'F\''')\n'
            s += indent + 'else:\n'
            s += indent + ' '*4 + f'{arg}_np = None\n'
        else:
            arglen = ' '*len(f'{arg}_np = np.asarray(')
            s = indent + f'{arg}_np = np.asarray(<{conv[info["type"]]}[:np.prod(({dims}))]> \n{arglen+indent}{arg})'
            s += f'.reshape(({dims}), order=''\'F\''')\n'
        f.write(s)
        count += 1
    f.write('\n' if count > 0 else '')


def writeCUCalls(f, name, routine):
    s = ' '*4
    count = 0
    for arg, info in routine.items():
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

    s += f'py{name[1:]}('
    count = 0
    for arg, info in routine.items():
        if info['array']:
            s += arg +'_np, '
        else:
            s += f'{arg}[0], '
        count += 1
    s = s[:-2] if count > 0 else s
    s += ')\n'
    f.write(s)


def convertToArrays(f, name, routine):
    indent = ' '*4
    s = ''
    for arg, info in routine.items():
        if info['intent'] is None:
            continue
        if 'out' not in info['intent']:
            continue

        if info['array']:
            if len(set(info['size'].split(','))) == 1 and info['size'].split(',')[0] == 'dim_ens-1':
                s += indent + f'cdef double[::1] {arg}_view\n'
                s += indent + 'if dim_ens[0] > 1:\n'
                s += indent + ' '*4 + f'{arg}_np[:] = {arg}_np_tmp[:]\n'
                s += indent + ' '*4 + f'{arg}_view = {arg}_np.ravel(order=''\'F\''')\n'
                s += indent + ' '*4 + f'assert {arg} == &{arg}_view[0], '\
                               f'\'reference (memory address) of {arg} has changed in {name}.\'\n'
            else:
                s += indent + f'{arg}_np[:] = {arg}_np_tmp[:]\n'
                s += indent + f'cdef {conv[info["type"]]}[::1] '
                s += f'{arg}_view = {arg}_np.ravel(order=''\'F\''')\n'
                s += indent + f'assert {arg} == &{arg}_view[0], '\
                               f'\'reference (memory address) of {arg} has changed in {name}.\'\n'
    f.write(s)


def writeUserCalls(filename, UserFuncInfo):
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
        for routine in UserFuncInfo:
            name = routine.pop('name')
            writeUserDef(f, name, routine)
            writeUserDocString(f, routine)
            s = '    raise RuntimeError(\'...Wrong '
            s += f'py{name[1:]} is called!!!...\')'
            f.write(s+'\n\n\n')
            routine['name'] = name

        for routine in UserFuncInfo:
            name = routine.pop('name')
            writeCUserDef(f, name, routine)
            writeArrayConversion(f, routine)
            writeCUCalls(f, name, routine)
            f.write('\n')
            convertToArrays(f, name, routine)
            f.write('\n\n')
            routine['name'] = name

if __name__ == '__main__':
    UserFuncInfo = getFuncInfo('../pyPDAF/fortran/U_PDAF_interface_c_binding.F90')
    writeUserPxdFile('U_PDAF_interface_c_binding.pxd', UserFuncInfo)
    writeUserCalls('U_PDAF_interface_c_binding.pyx', UserFuncInfo)

    PDAFInfo = getFuncInfo('../pyPDAF/fortran/PDAF_c_binding.F90')
    writePxdFile('PDAF_c_binding.pxd', UserFuncInfo, PDAFInfo)
    writePDAFcalls('PDAF_c_binding.pyx', PDAFInfo)
    PDAFInfo = getFuncInfo('../pyPDAF/fortran/PDAFomi_obs_c_binding.F90')
    writePxdFile('PDAFomi_obs_c_binding.pxd', UserFuncInfo, PDAFInfo)
    writePDAFcalls('PDAFomi_obs_c_binding.pyx', PDAFInfo)




