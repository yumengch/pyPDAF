import sys
import numpy as np
cimport numpy as cnp
from . cimport pdaf_c_cb_interface as pdaf_cb
from .cfi_binding cimport CFI_cdesc_t, CFI_address, CFI_index_t, CFI_establish
from .cfi_binding cimport CFI_attribute_other, CFI_type_double, CFI_type_int
from .cfi_binding cimport CFI_cdesc_rank1, CFI_cdesc_rank2, CFI_cdesc_rank3

try:
    import mpi4py
    mpi4py.rc.initialize = False
except ImportError:
    pass

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

def set_doassim(int  i_obs, int  doassim):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int 
        index into observation arrays
    doassim : int 
        Flag whether to assimilate this observation type

    Returns
    -------
    """
    with nogil:
        c__pdafomi_set_doassim(&i_obs, &doassim)



def set_disttype(int  i_obs, int  disttype):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int 
        index into observation arrays
    disttype : int 
        Index of distance type

    Returns
    -------
    """
    with nogil:
        c__pdafomi_set_disttype(&i_obs, &disttype)



def set_ncoord(int  i_obs, int  ncoord):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int 
        index into observation arrays
    ncoord : int 
        Number of coordinates

    Returns
    -------
    """
    with nogil:
        c__pdafomi_set_ncoord(&i_obs, &ncoord)



def set_obs_err_type(int  i_obs, int  obs_err_type):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int 
        index into observation arrays
    obs_err_type : int 
        Type of observation error

    Returns
    -------
    """
    with nogil:
        c__pdafomi_set_obs_err_type(&i_obs, &obs_err_type)



def set_use_global_obs(int  i_obs, int  use_global_obs):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int 
        index into observation arrays
    use_global_obs : int 
        Set whether to use global full observations

    Returns
    -------
    """
    with nogil:
        c__pdafomi_set_use_global_obs(&i_obs, &use_global_obs)



def set_inno_omit(int  i_obs, double  inno_omit):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int 
        index into observation arrays
    inno_omit : double 
        Set observation omission error variance level

    Returns
    -------
    """
    with nogil:
        c__pdafomi_set_inno_omit(&i_obs, &inno_omit)



def set_inno_omit_ivar(int  i_obs, double  inno_omit_ivar):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int 
        index into observation arrays
    inno_omit_ivar : double 
        Value of inverse variance to omit observation

    Returns
    -------
    """
    with nogil:
        c__pdafomi_set_inno_omit_ivar(&i_obs, &inno_omit_ivar)



def set_id_obs_p(int  i_obs, int  nrows, int  dim_obs_p, int [::1,:] id_obs_p):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int 
        index into observation arrays
    nrows : int 
        number of rows required in observation operator
    dim_obs_p : int 
        number of process local observations
    id_obs_p : ndarray[np.intc, ndim=2]
        Indices of process-local observed field in state vector
        Array shape: (nrows, dim_obs_p)

    Returns
    -------
    """
    with nogil:
        c__pdafomi_set_id_obs_p(&i_obs, &nrows, &dim_obs_p, &id_obs_p[0,0])



def set_icoeff_p(int  i_obs, int  nrows, int  dim_obs_p, 
    double [::1,:] icoeff_p):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int 
        index into observation arrays
    nrows : int 
        number of rows required in observation operator
    dim_obs_p : int 
        number of process local observations
    icoeff_p : ndarray[np.float64, ndim=2]
        Interpolation coeffs. for obs. operator
        Array shape: (nrows, dim_obs_p)

    Returns
    -------
    """
    with nogil:
        c__pdafomi_set_icoeff_p(&i_obs, &nrows, &dim_obs_p, &icoeff_p[0,0])



def set_domainsize(int  i_obs, int  ncoord, double [::1] domainsize):
    """Checking the corresponding PDAF documentation in https://pdaf.awi.de
    For internal subroutines checking corresponding PDAF comments.

    Parameters
    ----------
    i_obs : int 
        index into observation arrays
    ncoord : int 
        number of coordinates considered for localizations
    domainsize : ndarray[np.float64, ndim=1]
        Size of domain for periodicity (<=0 for no periodicity)
        Array shape: (ncoord)

    Returns
    -------
    """
    with nogil:
        c__pdafomi_set_domainsize(&i_obs, &ncoord, &domainsize[0])



