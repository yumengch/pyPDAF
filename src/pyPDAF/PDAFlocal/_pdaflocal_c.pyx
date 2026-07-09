import sys
import numpy as np
cimport numpy as cnp
from pyPDAF cimport pdaf_c_cb_interface as pdaf_cb
from pyPDAF.cfi_binding cimport CFI_cdesc_t, CFI_address, CFI_index_t, CFI_establish
from pyPDAF.cfi_binding cimport CFI_attribute_other, CFI_type_double, CFI_type_int
from pyPDAF.cfi_binding cimport CFI_cdesc_rank1, CFI_cdesc_rank2, CFI_cdesc_rank3

def set_indices(int  dim_l, int [::1] map):
    r"""set_indices(dim_l: int, map: np.ndarray) -> None

    Set index vector to map local state vector to global state vectors.

    This is called in the user-supplied function `py__init_dim_l_pdaf`.
    This function only sets the mapping for given domain index `domain_p`
    between local state vector and global state vector.
    Each element of map is an index of the global state vector. The index starts from 1.

    E.g., map[0] = 2 means that the first element of local state vector is the
    2rd element of the global state vector.

    Parameters
    ----------
    dim_l : int
        Dimension of local state vector
    map : ndarray[np.intc, dim=1]
        Index array for mapping between local and global state vector
        shape: (dim_l,)
    """
    with nogil:
        c__pdaflocal_set_indices(&dim_l, &map[0])



def set_increment_weights(int  dim_l, double [::1] weights):
    r"""set_increment_weights(dim_l: int, weights: np.ndarray) -> None

    Initialises a PDAF_internal local array of increment weights.

    This is called in the user-supplied function `py__init_dim_l_pdaf`.

    The weights are applied in in :func:`pyPDAF.PDAFlocal.l2g_cb` where the local state vector
    is weighted by given weights. These can e.g. be used to apply a vertical localisation.

    In vertical localization, the local state vector is a full vertical column
    of the model grid. In this case, one can make the increment weight depending
    on the height (or depth) of a grid point.

    Another application is to implement weakly-coupled assimilation in which
    the local state vector contains all variables, but only a subset of them is updated.

    This is achieved by givening those element that should not be updated the weight 0.

    Parameters
    ----------
    dim_l : int
        Dimension of local state vector
    weights : ndarray[np.float64, dim=1]
        Weights array. Shape: (dim_l,)
    """
    with nogil:
        c__pdaflocal_set_increment_weights(&dim_l, &weights[0])



def clear_increment_weights():
    r"""clear_increment_weights() -> None

    Deallocates the local increment weight vector in :func:`pyPDAF.PDAFlocal.set_increment_weights`.
    """
    with nogil:
        c__pdaflocal_clear_increment_weights()



