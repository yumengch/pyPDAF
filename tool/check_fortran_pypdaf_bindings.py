"""Compare C-bind Fortran wrappers with pyPDAF Cython bindings.

The script checks both directions:

* each ``c__...`` subroutine in ``src/fortran`` has a matching Cython
  ``cdef extern`` declaration and is called from a ``.pyx`` file;
* each pyPDAF ``c__...`` extern declaration or call site refers to a C-bind
  wrapper that exists in ``src/fortran``.

The parser is intentionally simple, but it handles wrapped Fortran subroutine
headers by collecting lines until ``bind(c)`` appears.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Iterable


IGNORE_DEFAULT = {
    "c__pdaf_parse_int",
    "c__pdaf_parse_real",
    "c__pdaf_parse_string",
    "c__pdaf_parse_logical",
}


FORTRAN_SUB_RE = re.compile(r"\bsubroutine\s+(c__\w+)\b", re.IGNORECASE)
CYTHON_EXTERN_RE = re.compile(
    r"\bcdef\s+(?:extern\s+)?[A-Za-z_]\w*(?:\s*\*)?\s+(c__\w+)\b",
    re.IGNORECASE,
)
CYTHON_CALL_RE = re.compile(r"\b(c__\w+)\s*\(")
PY_DEF_RE = re.compile(r"^\s*def\s+([A-Za-z_]\w*)\s*\(", re.MULTILINE)

PYTHON_PREFIXES = (
    "c__pdaflocalomi_",
    "c__pdaflocalomi",
    "c__pdaflocal_",
    "c__pdaflocal",
    "c__pdafomi_",
    "c__pdafomi",
    "c__pdaf3_",
    "c__pdaf_",
    "c__pdaf",
)

CALLBACK_BRIDGE_FILES = {
    "pdaf_c_cb_interface.f90",
    "pdaf_c_cb_interface.pxd",
    "pdaf_c_cb_interface.pyx",
}


def norm(name: str) -> str:
    return name.lower()


def strip_comment(line: str) -> str:
    if "!" not in line:
        return line
    return line.split("!", 1)[0]


def logical_lines(path: Path) -> Iterable[str]:
    current = ""
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        stripped = strip_comment(line).strip()
        if not stripped:
            continue

        if stripped.endswith("&"):
            current += stripped[:-1].strip() + " "
            continue

        if stripped.startswith("&"):
            current += stripped[1:].strip()
        else:
            current += stripped

        yield " ".join(current.split())
        current = ""

    if current:
        yield " ".join(current.split())


def collect_fortran_symbols(root: Path) -> dict[str, Path]:
    symbols: dict[str, Path] = {}
    for path in sorted(root.rglob("*.f90")):
        for line in logical_lines(path):
            if re.match(r"^\s*end\s+subroutine\b", line, re.IGNORECASE):
                continue
            if "bind(c" not in line.lower():
                continue
            match = FORTRAN_SUB_RE.search(line)
            if not match:
                continue
            symbols[norm(match.group(1))] = path
    return symbols


def collect_cython_externs(root: Path) -> dict[str, Path]:
    symbols: dict[str, Path] = {}
    for path in sorted(root.rglob("*.pxd")):
        text = path.read_text(encoding="utf-8", errors="ignore")
        for match in CYTHON_EXTERN_RE.finditer(text):
            symbols[norm(match.group(1))] = path
    return symbols


def collect_cython_calls(root: Path) -> dict[str, Path]:
    symbols: dict[str, Path] = {}
    for path in sorted(root.rglob("*.pyx")):
        text = path.read_text(encoding="utf-8", errors="ignore")
        for match in CYTHON_CALL_RE.finditer(text):
            symbols.setdefault(norm(match.group(1)), path)
    return symbols


def collect_python_defs(root: Path) -> set[str]:
    defs: set[str] = set()
    for path in sorted(root.rglob("*.pyx")):
        text = path.read_text(encoding="utf-8", errors="ignore")
        defs.update(match.group(1).lower() for match in PY_DEF_RE.finditer(text))
    return defs


def expected_python_name(symbol: str) -> str:
    """Return the expected Python wrapper name for a C-bind symbol.

    Most pyPDAF wrappers drop the module prefix, e.g. ``c__pdafomi_init`` maps
    to ``init``. Some PDAF filter names are fused to the ``pdaf`` prefix in the
    Fortran wrapper, e.g. ``c__pdafenkf_update`` maps to ``enkf_update``.
    Names that would start with a digit are exposed in Python with a leading
    underscore, e.g. ``c__pdaf_3dvar_init`` and ``c__pdaf3dvar_update`` map to
    ``_3dvar_init`` and ``_3dvar_update``.
    """
    name = norm(symbol)
    for prefix in PYTHON_PREFIXES:
        if name.startswith(prefix):
            name = name[len(prefix) :]
            break
    else:
        name = name.removeprefix("c__")

    if name and name[0].isdigit():
        return f"_{name}"
    return name


def without_callback_bridges(symbols: set[str], locations: dict[str, Path]) -> set[str]:
    return {
        name
        for name in symbols
        if locations.get(name, Path()).name.lower() not in CALLBACK_BRIDGE_FILES
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fortran-root", default="src/fortran")
    parser.add_argument("--python-root", default="src/pyPDAF")
    parser.add_argument(
        "--ignore",
        action="append",
        default=[],
        help="Additional c__ symbol to ignore; can be passed more than once.",
    )
    parser.add_argument(
        "--include-callback-bridges",
        action="store_true",
        help="Include callback bridge symbols from pdaf_c_cb_interface.f90.",
    )
    args = parser.parse_args()

    fortran_root = Path(args.fortran_root)
    python_root = Path(args.python_root)
    ignored = IGNORE_DEFAULT | {norm(item) for item in args.ignore}

    fortran_symbols = collect_fortran_symbols(fortran_root)
    cython_externs = collect_cython_externs(python_root)
    cython_calls = collect_cython_calls(python_root)
    python_defs = collect_python_defs(python_root)

    checked = {name for name in fortran_symbols if name not in ignored}
    if not args.include_callback_bridges:
        checked = without_callback_bridges(checked, fortran_symbols)
    cython_extern_names = set(cython_externs) - ignored
    cython_call_names = set(cython_calls) - ignored
    if not args.include_callback_bridges:
        cython_extern_names = without_callback_bridges(cython_extern_names, cython_externs)
        cython_call_names = without_callback_bridges(cython_call_names, cython_calls)

    missing_extern = sorted(checked - cython_extern_names)
    extra_extern = sorted(cython_extern_names - checked)
    missing_call = sorted(checked - cython_call_names)
    extra_call = sorted(cython_call_names - checked)
    missing_python_def = sorted(
        name for name in checked if expected_python_name(name) not in python_defs
    )

    print(f"Fortran C-bind subroutines: {len(fortran_symbols)}")
    print(f"Ignored symbols: {len(ignored & set(fortran_symbols))}")
    print(f"Checked symbols: {len(checked)}")
    print(f"Cython extern declarations: {len(cython_externs)}")
    print(f"Cython call sites: {len(cython_calls)}")
    print()

    def show(title: str, items: list[str], locations: dict[str, Path] | None = None) -> None:
        print(f"{title}: {len(items)}")
        for item in items:
            suffix = f"  [{locations[item]}]" if locations and item in locations else ""
            print(f"  {item}{suffix}")
        print()

    show("Fortran symbols missing from pyPDAF .pxd extern declarations", missing_extern, fortran_symbols)
    show("pyPDAF .pxd extern declarations missing from src/fortran", extra_extern, cython_externs)
    show("Fortran symbols missing from pyPDAF .pyx calls", missing_call, fortran_symbols)
    show("pyPDAF .pyx calls missing from src/fortran", extra_call, cython_calls)
    show("Missing same-name Python def wrappers", missing_python_def, fortran_symbols)

    if missing_extern or extra_extern or missing_call or extra_call:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
