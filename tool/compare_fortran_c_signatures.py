"""Compare PDAF Fortran signatures with pyPDAF C-wrapper signatures.

The older ``compare_subroutines.py`` script answers "which routine names are
new or removed?". This script answers the more specific binding-maintenance
question: for routines that exist in both ``PDAF/src`` and ``src/fortran``, do
the normalized argument lists still match?

The comparison is intentionally name based. For C wrappers the leading
``c__`` is removed. For OMI wrappers the derived-type arguments ``thisobs`` and
``thisobs_l`` are normalized to ``i_obs`` because pyPDAF wrappers pass the
observation type index and look up the derived-type instances internally.
"""
from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


SKIP_FILES_DEFAULT = {
    "pdaf_c_f_interface.f90",
}

PROC_RE = re.compile(
    r"^\s*"
    r"(?:(?:pure|elemental|recursive|module)\s+)*"
    r"(?:(?:integer|real|logical|complex|double\s+precision|character"
    r"|type\s*\([^)]*\)|class\s*\([^)]*\))"
    r"(?:\s*\([^)]*\))?\s+)?"
    r"(?P<kind>subroutine|function)\s+"
    r"(?P<name>[a-z_]\w*)\s*"
    r"(?:\((?P<args>[^)]*)\))?"
    r"(?:\s+result\s*\(\s*(?P<result>[a-z_]\w*)\s*\))?",
    re.IGNORECASE,
)

INTERFACE_START_RE = re.compile(r"^\s*(abstract\s+)?interface\b", re.IGNORECASE)
INTERFACE_END_RE = re.compile(r"^\s*end\s+interface\b", re.IGNORECASE)


@dataclass(frozen=True)
class Procedure:
    name: str
    kind: str
    args: tuple[str, ...]
    result: str | None
    path: Path
    line: int
    raw_name: str
    raw_args: tuple[str, ...]

    @property
    def signature(self) -> str:
        return f"{self.name}({','.join(self.args)})"

    @property
    def location(self) -> str:
        return f"{self.path}:{self.line}"


def strip_comment(line: str) -> str:
    """Remove Fortran comments, keeping the simple generated-code case sane."""
    if "!" not in line:
        return line
    return line.split("!", 1)[0]


def logical_lines(path: Path) -> Iterable[tuple[int, str]]:
    """Yield continued Fortran logical lines with the start line number."""
    current = ""
    start_line = 0
    for lineno, line in enumerate(path.read_text(encoding="utf-8", errors="ignore").splitlines(), 1):
        stripped = strip_comment(line).strip()
        if not stripped:
            continue

        if not current:
            start_line = lineno

        if stripped.endswith("&"):
            current += stripped[:-1].strip() + " "
            continue

        if stripped.startswith("&"):
            current += stripped[1:].strip() + " "
        else:
            current += stripped

        yield start_line, " ".join(current.split())
        current = ""

    if current:
        yield start_line, " ".join(current.split())


def split_args(args: str | None) -> tuple[str, ...]:
    if not args:
        return ()
    return tuple(arg.strip().lower() for arg in args.split(",") if arg.strip())


def normalize_name(name: str) -> str:
    name = name.lower()
    if name.startswith("c__"):
        name = name[3:]
    return name


def normalize_args(args: Iterable[str]) -> tuple[str, ...]:
    normalized = tuple(arg.lower() for arg in args)
    joined = ",".join(normalized)

    # C wrappers generally expose an observation type index instead of the
    # PDAF-OMI derived-type objects stored in module arrays.
    replacements = (
        ("thisobs_l,thisobs", "i_obs"),
        ("thisobs,thisobs_l", "i_obs"),
        ("thisobs_l", "i_obs"),
        ("thisobs", "i_obs"),
    )
    for old, new in replacements:
        joined = joined.replace(old, new)

    return tuple(arg for arg in joined.split(",") if arg)


def collect_procedures(root: Path, skip_files: set[str]) -> dict[str, list[Procedure]]:
    found: dict[str, list[Procedure]] = {}
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.suffix.lower() not in {".f90", ".f", ".f95", ".f03", ".f08"}:
            continue
        if path.name.lower() in skip_files:
            continue

        interface_depth = 0
        for line_no, line in logical_lines(path):
            if INTERFACE_START_RE.match(line):
                interface_depth += 1
                continue
            if INTERFACE_END_RE.match(line):
                interface_depth = max(0, interface_depth - 1)
                continue
            if interface_depth:
                continue

            # Avoid END SUBROUTINE / END FUNCTION lines.
            if re.match(r"^\s*end\s+(subroutine|function)\b", line, re.IGNORECASE):
                continue

            match = PROC_RE.search(line)
            if not match:
                continue

            raw_name = match.group("name")
            raw_args = split_args(match.group("args"))
            proc = Procedure(
                name=normalize_name(raw_name),
                kind=match.group("kind").lower(),
                args=normalize_args(raw_args),
                result=match.group("result").lower() if match.group("result") else None,
                path=path,
                line=line_no,
                raw_name=raw_name,
                raw_args=raw_args,
            )
            found.setdefault(proc.name, []).append(proc)

    return found


def choose_signature(procs: list[Procedure]) -> Procedure:
    """Prefer a public upstream definition over duplicates with same name."""
    return sorted(procs, key=lambda proc: (str(proc.path), proc.line))[0]


def format_proc(proc: Procedure) -> str:
    raw = f"{proc.raw_name}({','.join(proc.raw_args)})"
    if proc.kind == "function" and proc.result:
        raw += f" result({proc.result})"
    return f"{proc.signature}  [{proc.location}; raw: {raw}]"


def signatures_match(upstream: Procedure, wrapper: Procedure) -> bool:
    """Return whether a wrapper preserves the callable argument contract.

    For upstream subroutines the normalized arguments should match exactly.
    For upstream functions, C bindings are written as subroutines with one
    trailing output argument containing the function return value.
    """
    if upstream.args == wrapper.args:
        return True
    if upstream.kind == "function" and wrapper.kind == "subroutine":
        return wrapper.args[: len(upstream.args)] == upstream.args and (
            len(wrapper.args) == len(upstream.args) + 1
        )
    return False


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare PDAF/src routine signatures with src/fortran C wrappers.",
    )
    parser.add_argument("--pdaf-src", default="PDAF/src", help="Path to upstream PDAF Fortran sources.")
    parser.add_argument("--c-wrapper-src", default="src/fortran", help="Path to C-wrapper Fortran sources.")
    parser.add_argument(
        "--skip-file",
        action="append",
        default=[],
        help="File name to skip, case-insensitive. Can be passed more than once.",
    )
    parser.add_argument(
        "--show-matches",
        action="store_true",
        help="Print routines whose normalized argument lists match.",
    )
    args = parser.parse_args()

    skip_files = SKIP_FILES_DEFAULT | {name.lower() for name in args.skip_file}
    pdaf = collect_procedures(Path(args.pdaf_src), skip_files)
    wrappers = collect_procedures(Path(args.c_wrapper_src), skip_files)

    pdaf_names = set(pdaf)
    wrapper_names = set(wrappers)
    common_names = sorted(pdaf_names & wrapper_names)

    missing_wrappers = sorted(pdaf_names - wrapper_names)
    extra_wrappers = sorted(wrapper_names - pdaf_names)
    mismatches: list[tuple[Procedure, Procedure]] = []
    matches: list[tuple[Procedure, Procedure]] = []

    for name in common_names:
        upstream = choose_signature(pdaf[name])
        wrapper = choose_signature(wrappers[name])
        if signatures_match(upstream, wrapper):
            matches.append((upstream, wrapper))
        else:
            mismatches.append((upstream, wrapper))

    print(f"PDAF/src procedures: {sum(len(items) for items in pdaf.values())} ({len(pdaf)} unique names)")
    print(f"src/fortran C-wrapper procedures: {sum(len(items) for items in wrappers.values())} ({len(wrappers)} unique names)")
    print(f"Common names: {len(common_names)}")
    print()

    print(f"Signature mismatches: {len(mismatches)}")
    for upstream, wrapper in mismatches:
        print(f"  {upstream.name}")
        print(f"    PDAF/src:    {format_proc(upstream)}")
        print(f"    src/fortran: {format_proc(wrapper)}")
    print()

    print(f"Missing wrappers for PDAF/src names: {len(missing_wrappers)}")
    for name in missing_wrappers:
        print(f"  {format_proc(choose_signature(pdaf[name]))}")
    print()

    print(f"Extra wrapper names not found in PDAF/src: {len(extra_wrappers)}")
    for name in extra_wrappers:
        print(f"  {format_proc(choose_signature(wrappers[name]))}")
    print()

    if args.show_matches:
        print(f"Matching signatures: {len(matches)}")
        for upstream, wrapper in matches:
            print(f"  {upstream.signature}")

    return 1 if mismatches else 0


if __name__ == "__main__":
    raise SystemExit(main())
