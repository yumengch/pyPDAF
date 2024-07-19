"""This module reads .F90 files and returns a list of functions

This module can only be used with pyPDAF/fortran files or PDAFc interface files.

This module has no intention to be generalised for general parser of Fortran files.
"""
import typing
import warnings
import re

def extract_outermost_bracket_content(s:str) -> str|None:
    stack : list[int] = []
    outer_start = None
    outer_end = None

    for i, char in enumerate(s):
        if char == '(':
            if not stack:
                outer_start = i
            stack.append(i)
        elif char == ')':
            stack.pop()
            if not stack:
                outer_end = i
                break

    if outer_start is not None and outer_end is not None:
        return s[outer_start + 1 : outer_end]
    else:
        return None

def merge_brackets(string_list : list[str]) -> list[str]:
    string_list_new : list[str] = []
    # merge arguments if arrays definitions are spliited.
    merge:bool = False
    current_merged:str = ""
    n_left_bracket:int = 0
    for string in string_list:
        if "(" in string:
            merge = True
            current_merged += string+','
            n_left_bracket += string.count('(')
            if ")" in string:
                n_left_bracket -= string.count(')')
                if n_left_bracket == 0:
                    string_list_new.append(current_merged[:-1])
                    merge = False
                    current_merged = ""
        elif merge:
            current_merged += string+','
            if ")" in string:
                n_left_bracket -= string.count(')')
                if n_left_bracket == 0:
                    string_list_new.append(current_merged[:-1])
                    merge = False
                    current_merged = ""
        else:
            string_list_new.append(string)

    return string_list_new

def merge_line(f: typing.TextIO, line:str) -> str:
    """
    A function to merge lines concatenate by ampersand,
    removing white spaces, handling module imports, and concatenating lines.

    Parameters:
        f : file object
        line : str

    Returns:
        str
    """
    # remove all white lines
    skip_condition:bool = line.isspace()
    # remove line break
    line = line.strip().replace('\n', '')
    # ignore lines starting with use as they are module import
    skip_condition = skip_condition or line[:3].lower() == 'use'
    # ignore implicit none lines
    skip_condition = skip_condition or line.replace(' ', '').lower() == 'implicitnone'
    if skip_condition:
        return ' '

    # merge all concatenated lines into one string
    # Here we assume no ampersand in the start of the line
    while '&' in line.split('!')[0].lower():
        line = line.replace('&', '')
        line = line + f.readline().lower().strip().replace('\n', '')
    return line

def get_subroutine_name(line:str) -> str:
    """
    Extracts the subroutine name from the given line.

    Args:
        line (str): The line of code containing the subroutine declaration.

    Returns:
        str: The name of the subroutine.
    """
    name:str = re.split("subroutine", line.split('(')[0], flags=re.IGNORECASE)[-1].replace(' ', '')
    return name

def get_args_attr(args_list:str) -> list[str]:
    """get argument names from a list of arguments

    Args:
        args_list (str): list of arguments
    """
    args_split : list[str] = args_list.split(',')
    return merge_brackets(args_split)

def get_arg_info(arg_info:dict[str, dict[str, str|bool|None|list[str]]], line:str, comment:str):
    declaration = line.replace(' ', '').split('::')
    args = get_args_attr(declaration[1])
    attrs = get_args_attr(declaration[0].lower())

    for arg in args:
        argname = arg.split('(')[0]
        # skip variables that are not in the argument list
        if argname not in arg_info:
            warnings.warn(f'{argname} is not in f{list(arg_info.keys())}')
            continue

        arg_info[argname]['type'] =  attrs[0].split('(')[0]
        arg_info[argname]['kind'] = attrs[0].split('(')[1].split(')')[0]
        intent_string = [s for s in attrs if 'intent' in s]
        arg_info[argname]['intent'] = intent_string[0].replace(' ', '')[7:-1] if len(intent_string) != 0 else None
        dimension_string : list[str|None] = [s for s in attrs if 'dimension' in s]
        if len(dimension_string) == 0:
            if '(' in arg:
                dimension_string =  [extract_outermost_bracket_content(arg),]
        else:
            dimension_string =  [d[10:-1] for d in dimension_string]
        arg_info[argname]['array'] = len(dimension_string) > 0
        arg_info[argname]['dimension'] = merge_brackets(dimension_string[0].split(',')) if arg_info[argname]['array'] else None
        arg_info[argname]['comment'] = comment.strip()
        if arg_info[argname]['array']:
            arg_info[argname]['comment'] = arg_info[argname]['comment'] + '; Dimension: ' + dimension_string[0]
    return arg_info

def get_arguments(subroutine:list[str]):
    # get the arg list and subroutine name
    line:str = subroutine[0]
    arg_list:list[str] = line.replace(' ', '').replace(')', '(').split('(')[1].split(',')
    arg_info:dict[str, dict[str, str|bool|None|list[str]]] = {}
    for arg in arg_list:
        if arg == '':
            continue
        arg_info[arg] = {}

    # get information in the arg list for C declaration
    comment:str = ''
    for line in subroutine[1:]:
        if '!' in line:
            comment = ''.join([comment, line.split('!')[1].strip()])

        if '::' not in line:
            continue


        assert '!' not in line, 'if ! is in  the same line as ::, '\
                                'it doesn\'t fit the assumption that '\
                                'comment and arguments declaration are separated.'

        arg_info = get_arg_info(arg_info, line, comment)

        comment = ''

    return arg_info

def get_subroutine_list(filename:str) -> dict[str, list[str]]:
    # Read the source file
    with open(filename, 'r') as f:
        subroutines:dict[str, list[str]] = {}
        # The first line is either comment or program/module name
        line:str = merge_line(f, f.readline().lower())
        # using start to exclude unwanted content
        # when start is true, the following lines are the content we want
        do_append:bool = False
        while line:
            line_orig = merge_line(f, f.readline())
            line = line_orig.lower()
            # get the line without comments
            no_comment_line:str = line.split('!')[0].replace(' ', '')
            # when one subroutine starts
            start_condition:bool = 'subroutine' in no_comment_line
            start_condition = start_condition and '(' in no_comment_line
            start_condition = start_condition and ')' in no_comment_line
            # when one subroutine ends
            end_condition:bool = 'endsubroutine' in no_comment_line
            # get the subroutine name and start construct the string list
            if start_condition:
                do_append = True
                subroutine_name:str = get_subroutine_name(line_orig)
                subroutines[subroutine_name] = []
                subroutines[subroutine_name].append(line)
            # append the line to the string list of current subroutine
            if not line.isspace() and do_append:
                subroutines[subroutine_name].append(line)
            if end_condition:
                do_append = False
    return subroutines

def get_func_info(filenames:list[str]) -> dict[str, dict[str,list[dict[str,str|bool|None|list[str]]]]]:
    func_info : dict[str, dict[str,list[dict[str,str|bool|None|list[str]]]]] = {}
    for filename in filenames:
        subroutines = get_subroutine_list(filename)

        for subroutine_name, subroutine in subroutines.items():
            arg_info = get_arguments(subroutine)
            assert subroutine_name not in func_info, 'Error: same subroutine name in files'
            func_info[subroutine_name] = arg_info

    return func_info


if __name__ == '__main__':
    get_func_info(['../pyPDAF/fortran/PDAF_c_binding.F90', ])
    get_func_info(['../pyPDAF/fortran/PDAFomi_obs_c_binding.F90', ])
    get_func_info(['../pyPDAF/fortran/U_PDAF_interface_c_binding.F90', ])