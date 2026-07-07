"""Partial type stubs for :mod:`pyPDAF.PDAF.internal`."""
from typing import Any

def time_temp(timerid: int) -> float:
    """Return the last elapsed interval for a PDAF timer.

    The value is the duration, in seconds, of the most recent completed timing
    interval for ``timerid``. Timing intervals are controlled with
    :func:`timeit`, typically by starting a timer with ``operation="new"`` and
    stopping it with ``operation="old"``.

    Parameters
    ----------
    timerid : int
        PDAF timer identifier.

    Returns
    -------
    time_temp : float
        Elapsed time in seconds for the last completed interval of this timer.
    """

def time_tot(timerid: int) -> float:
    """Return the accumulated elapsed time for a PDAF timer.

    The value is the total time, in seconds, accumulated over all completed
    timing intervals for ``timerid`` since the timer arrays were initialized.
    Timing intervals are controlled with :func:`timeit`.

    Parameters
    ----------
    timerid : int
        PDAF timer identifier.

    Returns
    -------
    time_tot : float
        Accumulated elapsed time in seconds for this timer.
    """

def memcount_get(id: int, munit: str) -> float:
    """Read a process-local PDAF memory counter.

    The counter value is returned in the requested memory unit. PDAF supports
    four one-character unit codes: ``"B"`` for bytes, ``"K"`` for kibibytes,
    ``"M"`` for mebibytes, and ``"G"`` for gibibytes. Lowercase unit codes are
    also accepted by PDAF.

    Parameters
    ----------
    id : int
        Identifier of the memory counter to read.
    munit : str
        Unit code for the returned value. Supported values are ``"B"``,
        ``"K"``, ``"M"``, and ``"G"``.

    Returns
    -------
    memcount_value : float
        Process-local memory count converted to the requested unit. PDAF
        returns ``0.0`` for unsupported unit codes.
    """

def memcount_get_global(id: int, munit: str, comm: int) -> float:
    """Read a globally reduced PDAF memory counter.

    This function first reads the process-local memory counter and then uses
    ``comm`` to compute the global value with MPI. The counter value is
    returned in the requested memory unit. PDAF supports four one-character
    unit codes: ``"B"`` for bytes, ``"K"`` for kibibytes, ``"M"`` for
    mebibytes, and ``"G"`` for gibibytes. Lowercase unit codes are also
    accepted by PDAF.

    Parameters
    ----------
    id : int
        Identifier of the memory counter to read.
    munit : str
        Unit code for the returned value. Supported values are ``"B"``,
        ``"K"``, ``"M"``, and ``"G"``.
    comm : int
        Fortran MPI communicator used for the global reduction.

    Returns
    -------
    memcount_value : float
        Globally reduced memory count converted to the requested unit. PDAF
        returns ``0.0`` for unsupported unit codes.
    """

def __getattr__(name: str) -> Any: ...
