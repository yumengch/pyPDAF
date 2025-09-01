import sys
import numpy as np
cimport numpy as cnp
from pyPDAF cimport pdaf_c_cb_interface as pdaf_cb
from pyPDAF.cfi_binding cimport CFI_cdesc_t, CFI_address, CFI_index_t, CFI_establish
from pyPDAF.cfi_binding cimport CFI_attribute_other, CFI_type_double, CFI_type_int
from pyPDAF.cfi_binding cimport CFI_cdesc_rank1, CFI_cdesc_rank2, CFI_cdesc_rank3

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

def iau_init(int  type_iau_in, int  nsteps_iau_in):
    """Initialise parameters for incremental analysis updates, IAU.

    It further allocates the array in which the ensemble increments are stored.
    This array exists on all processes that are part of model tasks.

    This is usually called after :func:`pyPDAF.PDAF.init`.
    It is important that this is called by all model processes because all these
    processes need the information on the IAU configuration and
    need to allocate the increment array.

    Parameters
    ----------
    type_iau_in : int
        - (0) no IAU
        - (1) constant increment weight 1/nsteps_iau
        - (2) Linear IAU weight with maximum in middle of IAU period. This weight
              linearly increases and decreases. This type is usually used if
              the IAU is applied in re-running over the previous observation period.
        - (3) Zero weights for null mode (can be used to apply IAU on user side).
              This stores the increment information, but does not apply the increment.
              One can use :func:`pyPDAF.PDAF.iau_set_pointer` to access the increment
              array.
    nsteps_iau_in : int
        number of time steps in IAU

    Returns
    -------
    flag : int
        Status flag
    """
    cdef int  flag
    with nogil:
        c__pdaf_iau_init(&type_iau_in, &nsteps_iau_in, &flag)

    return flag


def iau_reset(int  type_iau_in, int  nsteps_iau_in):
    """Modify the IAU type and the number of IAU time steps during a run

    While :func:`pyPDAF.PDAF.iau_init` sets the IAU type and number of IAU steps
    initially, one can change these settings during a model run.

    This function has to be called by all processes that are model processes.
    A common place is to call the function after an analysis step or
    in `distribute_state_pdaf`.

    Parameters
    ----------
    type_iau_in : int
        - (0) no IAU
        - (1) constant increment weight 1/nsteps_iau
        - (2) Linear IAU weight with maximum in middle of IAU period. This weight
              linearly increases and decreases. This type is usually used if
              the IAU is applied in re-running over the previous observation period.
        - (3) Zero weights for null mode (can be used to apply IAU on user side).
              This stores the increment information, but does not apply the increment.
              One can use :func:`pyPDAF.PDAF.iau_set_pointer` to access the increment
              array.
    nsteps_iau_in : int
        number of time steps in IAU

    Returns
    -------
    flag : int
        Status flag
    """
    cdef int  flag
    with nogil:
        c__pdaf_iau_reset(&type_iau_in, &nsteps_iau_in, &flag)

    return flag


def iau_set_weights(int  iweights, double [::1] weights):
    """Provide a user-specified vector of increment weights.

    While :func:`pyPDAF.PDAF.iau_init` allows to choose among
    pre-defined weight functions, one might like to use a different function
    and the corresponding weights can be set here.

    All model processes must call the routine.
    A common place is to call the function after an analysis step or in `distribute_state_pdaf`.

    Parameters
    ----------
    iweights : int
        Length of weights input vector
        If iweights is different from the number of IAU steps set in
        :func:`pyPDAF.PDAF.iau_init` or :func:`pyPDAF.PDAF.iau_reset`,
        only the minimum of iweights and the set IAU steps is filled with
        the provided weights vector.
    weights : ndarray[np.float64, ndim=1]
        Input weight vector
        Array shape: (iweights)
    """
    with nogil:
        c__pdaf_iau_set_weights(&iweights, &weights[0])



def iau_set_pointer():
    """Set a pointer to the ensemble increments array.

    This gives direct access to the increment array,
    e.g. to analyze it or to write it into a file for restarting.

    If it is called by each single process, but it only provides a pointer to
    the process-local part of the increment array.

    For domain-decomposed models, this array only includes the state vector
    part for the process domain. In addition, it usually only contains a
    sub-ensemble unless one uses the flexible parallelization mode with a
    single model task. For the fully parallel mode, the process(es) of a
    single model task only hold a single ensemble state.

    Returns
    -------
    iau_ptr_np : np.ndarray
        The increment array (process-local part)
    flag : int
        Status flag
    """
    cdef CFI_cdesc_rank2 iau_ptr_cfi
    cdef CFI_cdesc_t *iau_ptr_ptr = <CFI_cdesc_t *> &iau_ptr_cfi
    cdef int  flag
    with nogil:
        c__pdaf_iau_set_pointer(iau_ptr_ptr, &flag)

    cdef CFI_index_t iau_ptr_subscripts[2]
    iau_ptr_subscripts[0] = 0
    iau_ptr_subscripts[1] = 0
    cdef double * iau_ptr_ptr_np
    iau_ptr_ptr_np = <double *>CFI_address(iau_ptr_ptr, iau_ptr_subscripts)
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] iau_ptr_np = np.asarray(<double [:iau_ptr_ptr.dim[0].extent:1,:iau_ptr_ptr.dim[1].extent]> iau_ptr_ptr_np, order="F")
    return iau_ptr_np, flag


def iau_init_inc(int  dim_p, int  dim_ens_l, double [::1,:] ens_inc):
    """Fill the process-local increment array.

    A common use is when the IAU should be applied from the initial time of a run,
    for example if the increment was computed in a previous assimilation run
    and then stored for restarting. Since PDAF can only compute an increment
    in an analysis step, the user needs to provide the increment.

    The function is called after :func:`pyPDAF.PDAF.iau_init`.
    It has to be called by all processes that are model processes and one needs
    to provide the task-local ensemble
    (i.e. with local ensemble size `dim_ens_l=1` for the fully parallel mode,
    and usually `dim_ens_l>1` for the flexible parallelization mode).
    The function cannot be called in `init_ens_pdaf` since this function is only
    executed by filter processes instead of all model processes.

    Parameters
    ----------
    dim_p : int
        PE-local dimension of model state
    dim_ens_l : int
        Number of ensemble members that are run by a loop for each model task
    ens_inc : ndarray[np.float64, ndim=2]
        PE-local increment ensemble
        Array shape: (dim_p, dim_ens_l)

    Returns
    -------
    flag : int
        Status flag
    """
    cdef int  flag
    with nogil:
        c__pdaf_iau_init_inc(&dim_p, &dim_ens_l, &ens_inc[0,0], &flag)

    return flag


def iau_add_inc(py__collect_state_pdaf, py__distribute_state_pdaf):
    """Apply IAU to model forecasts in flexible parallel mode.

    PDAF automatically apply IAU in fully parallel model. However, one has to
    apply IAU with this function in flexible parallel mode.

    This has to be used in each model time stepping.

    Parameters
    ----------
    py__collect_state_pdaf : Callable
        Collect a state vector for PDAF

    py__distribute_state_pdaf : Callable
        Distribute a state vector for PDAF
    """
    pdaf_cb.collect_state_pdaf = <void*>py__collect_state_pdaf
    pdaf_cb.distribute_state_pdaf = <void*>py__distribute_state_pdaf
    with nogil:
        c__pdaf_iau_add_inc(pdaf_cb.c__collect_state_pdaf,
                            pdaf_cb.c__distribute_state_pdaf)

def iau_set_ens_pointer():
    """Set a pointer to the ensemble increments array.

    This is the same as :func:`pyPDAF.PDAF.iau_set_pointer`.

    This gives direct access to the increment array,
    e.g. to analyze it or to write it into a file for restarting.

    If it is called by each single process, but it only provides a pointer to
    the process-local part of the increment array.

    For domain-decomposed models, this array only includes the state vector
    part for the process domain. In addition, it usually only contains a
    sub-ensemble unless one uses the flexible parallelization mode with a
    single model task. For the fully parallel mode, the process(es) of a
    single model task only hold a single ensemble state.

    Returns
    -------
    iau_ptr_np : np.ndarray
        The increment array (process-local part)
    flag : int
        Status flag
    """
    cdef CFI_cdesc_rank2 iau_ptr_cfi
    cdef CFI_cdesc_t *iau_ptr_ptr = <CFI_cdesc_t *> &iau_ptr_cfi
    cdef int  flag
    with nogil:
        c__pdaf_iau_set_ens_pointer(iau_ptr_ptr, &flag)

    cdef CFI_index_t iau_ptr_subscripts[2]
    iau_ptr_subscripts[0] = 0
    iau_ptr_subscripts[1] = 0
    cdef double * iau_ptr_ptr_np
    iau_ptr_ptr_np = <double *>CFI_address(iau_ptr_ptr, iau_ptr_subscripts)
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="fortran", negative_indices=False, cast=False] iau_ptr_np = np.asarray(<double [:iau_ptr_ptr.dim[0].extent:1,:iau_ptr_ptr.dim[1].extent]> iau_ptr_ptr_np, order="F")
    return iau_ptr_np, flag

def iau_set_state_pointer():
    """Set a pointer to the state increments array.

    This gives direct access to the increment array used in the ensemble
    optimal interpolation mode. This can be used
    e.g. to analyze it or to write it into a file for restarting.

    If it is called by each single process, but it only provides a pointer to
    the process-local part of the increment array.

    For domain-decomposed models, this array only includes the state vector
    part for the process domain. In addition, it usually only contains a
    sub-ensemble unless one uses the flexible parallelization mode with a
    single model task. For the fully parallel mode, the process(es) of a
    single model task only hold a single ensemble state.

    Returns
    -------
    iau_ptr_np : np.ndarray
        The increment array (process-local part)
    flag : int
        Status flag
    """
    cdef CFI_cdesc_rank1 iau_x_ptr_cfi
    cdef CFI_cdesc_t *iau_x_ptr_ptr = <CFI_cdesc_t *> &iau_x_ptr_cfi
    cdef int  flag
    with nogil:
        c__pdaf_iau_set_state_pointer(iau_x_ptr_ptr, &flag)

    cdef CFI_index_t iau_x_ptr_subscripts[1]
    iau_x_ptr_subscripts[0] = 0
    cdef double * iau_x_ptr_ptr_np
    iau_x_ptr_ptr_np = <double *>CFI_address(iau_x_ptr_ptr, iau_x_ptr_subscripts)
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="fortran", negative_indices=False, cast=False] iau_x_ptr_np = np.asarray(<double [:iau_x_ptr_ptr.dim[0].extent:1]> iau_x_ptr_ptr_np, order="F")
    return iau_x_ptr_np, flag