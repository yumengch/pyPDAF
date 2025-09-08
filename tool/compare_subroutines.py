import os
import re

def preprocess_fortran_file(file_path):
    with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()

    merged_lines = []
    current_line = ""

    for line in lines:
        stripped = line.strip()
        if "!" in stripped:
            code_part, comment_part = stripped.split("!", 1)
            stripped = code_part.strip()
        if stripped.endswith("&"):
            current_line += stripped[:-1] + " "
        elif stripped.startswith("&"):
            current_line += stripped[1:] + " "
        else:
            current_line += stripped
            merged_lines.append(current_line)
            current_line = ""

    return merged_lines

def extract_subroutines_from_dir(path):
    subroutine_signatures = set()

    for root, _, files in os.walk(path):
        for file in files:
            if file.lower().endswith(('.f90', '.f', '.f95', '.f03', '.f08')):
                if file == 'pdaf_c_f_interface.f90':
                    continue
                full_path = os.path.join(root, file)
                lines = preprocess_fortran_file(full_path)

                interface_depth = 0

                for line in lines:
                    # Entering an INTERFACE or ABSTRACT INTERFACE
                    if re.match(r'\s*(abstract\s+)?interface\b', line, re.I):
                        interface_depth += 1
                        continue

                    # Leaving an INTERFACE block
                    if re.match(r'\s*end\s+interface\b', line, re.I):
                        if interface_depth > 0:
                            interface_depth -= 1
                        continue

                    if interface_depth != 0: continue
                    match = re.search(r'\bsubroutine\s+([a-zA-Z0-9_]+)\s*\(([^)]*)\)', line, re.IGNORECASE)
                    if match:
                        name = match.group(1).lower()
                        args = match.group(2).replace(" ", "").lower()
                        signature = f"{name}({args})"
                        # removing c__ in pyPDAF c binding files
                        if name[:3] == 'c__':
                            signature = signature[3:]
                        # replace thisobs_l and thisobs with i_obs
                        if 'thisobs_l,thisobs' in signature:
                            signature = signature.replace('thisobs_l,thisobs', 'i_obs')
                        if 'thisobs,thisobs_l' in signature:
                            signature = signature.replace('thisobs,thisobs_l', 'i_obs')
                        if 'thisobs_l' in signature:
                            signature = signature.replace('thisobs_l', 'i_obs')
                        if 'thisobs' in signature:
                            signature = signature.replace('thisobs', 'i_obs')
                        subroutine_signatures.add(signature)

    return subroutine_signatures

def compare_subroutines(dir_old, dir_new):
    old_subs = extract_subroutines_from_dir(dir_old)
    new_subs = extract_subroutines_from_dir(dir_new)
    # new_subs = set()

    # print (old_subs)

    # for sub in new_subs:
    #     if 'pdaf3_assim_offline(' in sub:
    #         print(f"New subroutine found: {sub}")

    # for sub in old_subs:
    #     if 'pdaf3_assim_offline(' in sub:
    #         print(f"Old subroutine found: {sub}")

    added = new_subs - old_subs
    removed = old_subs - new_subs
    unchanged = old_subs & new_subs

    return added, removed, unchanged

if __name__ == "__main__":
    added, removed, unchanged = compare_subroutines('src/fortran',
                                                    'PDAF/src')

    print("\n--- Added Subroutines ---")
    for sig in sorted(added):
        print(sig)

    print("\n--- Removed Subroutines ---")
    for sig in sorted(removed):
        print(sig)

    print(f"\n--- Summary ---")
    print(f"Added:   {len(added)}")
    print(f"Removed: {len(removed)}")
    print(f"Unchanged: {len(unchanged)}")
