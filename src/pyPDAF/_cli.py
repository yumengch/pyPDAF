"""Command-line helpers for pyPDAF discovery functions."""
from __future__ import annotations

import argparse


def print_version_main() -> None:
    """Print the linked PDAF version."""
    from pyPDAF import print_version

    print_version()


def print_filter_types_main() -> None:
    """Print available PDAF filter types."""
    from pyPDAF import print_filter_types

    print_filter_types(verbose=3)


def options_filters_main() -> None:
    """Print option information for a PDAF filter type."""
    parser = argparse.ArgumentParser(
        description="Print PDAF option information for one filter type."
    )
    parser.add_argument(
        "filtertype",
        type=int,
        help="PDAF filter type identifier.",
    )
    args = parser.parse_args()

    from pyPDAF import options_filters

    options_filters(args.filtertype)
