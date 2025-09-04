# pylint: disable=invalid-name, no-name-in-module, wrong-import-position
import os
import sys

try:
    import mpi4py
    mpi4py.rc.initialize = False
except ImportError:
    pass

def _append_to_sharedlib_load_path():
    """Ensure the shared libraries in this package can be loaded on Windows.

    Windows lacks a concept equivalent to RPATH: Python extension modules
    cannot find DLLs installed outside the DLL search path. This function
    ensures that the location of the shared libraries distributed inside this
    Python package is in the DLL search path of the process.

    The Windows DLL search path includes the path to the object attempting
    to load the DLL: it needs to be augmented only when the Python extension
    modules and the DLLs they require are installed in separate directories.
    Cygwin does not have the same default library search path: all locations
    where the shared libraries are installed need to be added to the search
    path.

    This function is very similar to the snippet inserted into the main
    ``__init__.py`` of a package by ``delvewheel`` when it vendors external
    shared libraries.

    .. note::

        `os.add_dll_directory` is only available for Python 3.8 and later, and
        in the Conda ``python`` packages it works as advertised only for
        version 3.10 and later. For older Python versions, pre-loading the DLLs
        with `ctypes.WinDLL` may be preferred.
    """
    basedir = os.path.dirname(__file__)
    if os.name == 'nt':
        os.add_dll_directory(basedir)
    elif sys.platform == 'cygwin':
        os.environ['PATH'] = os.pathsep.join((os.environ['PATH'], basedir))

# Global error handler
def global_except_hook(exctype, value, traceback):
    from traceback import print_exception
    try:
        import mpi4py.MPI

        if mpi4py.MPI.Is_initialized():
            try:
                sys.stderr.write('Uncaught exception was ''detected on rank {}.\n'.format(
                    mpi4py.MPI.COMM_WORLD.Get_rank()))
                print_exception(exctype, value, traceback)
                sys.stderr.write("\n")
                sys.stderr.flush()
            finally:
                try:
                    mpi4py.MPI.COMM_WORLD.Abort(1)
                except Exception as e:
                    sys.stderr.write('MPI Abort failed, this process will hang.\n')
                    sys.stderr.flush()
                    raise e
        else:
            sys.__excepthook__(exctype, value, traceback)
    except ImportError:
        sys.__excepthook__(exctype, value, traceback)

sys.excepthook = global_except_hook

_append_to_sharedlib_load_path()

from pyPDAF.PDAF3 import init, init_forecast, set_parallel
from pyPDAF.PDAF3 import assimilate, assim_offline
from pyPDAF.PDAF3 import assimilate_local_nondiagr, assimilate_global_nondiagr, \
                         assimilate_lnetf_nondiagr, assimilate_lknetf_nondiagr, \
                         assimilate_enkf_nondiagr, assimilate_nonlin_nondiagr
from pyPDAF.PDAF3 import  assim_offline_local_nondiagr, \
                        assim_offline_global_nondiagr, \
                        assim_offline_lnetf_nondiagr, \
                        assim_offline_lknetf_nondiagr, \
                        assim_offline_enkf_nondiagr, \
                        assim_offline_lenkf_nondiagr, \
                        assim_offline_nonlin_nondiagr
from pyPDAF.PDAF3 import assimilate_3dvar_all, assim_offline_3dvar_all
from pyPDAF.PDAF3 import assimilate_3dvar_nondiagr, assimilate_en3dvar_estkf_nondiagr, \
                         assimilate_en3dvar_lestkf_nondiagr, \
                         assimilate_hyb3dvar_estkf_nondiagr, \
                         assimilate_hyb3dvar_lestkf_nondiagr
from pyPDAF.PDAF3 import assim_offline_3dvar_nondiagr, assim_offline_en3dvar_estkf_nondiagr, \
                         assim_offline_en3dvar_lestkf_nondiagr, \
                         assim_offline_hyb3dvar_estkf_nondiagr, \
                         assim_offline_hyb3dvar_lestkf_nondiagr
from pyPDAF.PDAF3 import generate_obs, generate_obs_offline
from pyPDAF.PDAF import deallocate
