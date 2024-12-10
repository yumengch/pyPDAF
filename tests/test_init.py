import numpy as np
import mpi4py.MPI as MPI
import importlib
import pytest


def test_dim_ens_1(filter_type, subtype):
    import pyPDAF.PDAF as PDAF
    """Test the initialisation when #
    dimension of `uinv` is (dim_ens - 1, dim_ens - 1).
    """
    def init_ens_pdaf(filtertype:int, dim_p:int, dim_ens:int,
                      state_p:np.ndarray, uinv:np.ndarray,
                      ens_p:np.ndarray, status_pdaf:int
                      ) -> tuple[np.ndarray, np.ndarray,
                                 np.ndarray, int]:
        assert uinv.shape == (dim_ens - 1, dim_ens - 1), \
            f"Expected uinv.shape to be {(dim_ens - 1, dim_ens - 1)}"\
                f" but got {uinv.shape}"
        return state_p, uinv, ens_p, status_pdaf

    dim_p = 20
    dim_ens = 4
    forget = 1.0

    param_int = np.array([dim_p, dim_ens, ], dtype=np.intc)
    param_float = np.array([forget, ])
    if filter_type == 200:
        param_int = np.array([dim_p, dim_ens, 1, dim_p, dim_p], dtype=np.intc)

    _, _, status = PDAF.init(filter_type,
                            subtype,
                            0,
                            param_int, param_float,
                            MPI.COMM_WORLD.py2f(),
                            MPI.COMM_WORLD.py2f(),
                            MPI.COMM_WORLD.py2f(),
                            1, 1, True,
                            init_ens_pdaf, 0)


def test_dim_ens(filter_type, subtype):
    import pyPDAF.PDAF as PDAF
    """Test the initialisation when #
    dimension of `uinv` is (dim_ens, dim_ens).
    """
    def init_ens_pdaf(filtertype:int, dim_p:int, dim_ens:int,
                      state_p:np.ndarray, uinv:np.ndarray,
                      ens_p:np.ndarray, status_pdaf:int
                      ) -> tuple[np.ndarray, np.ndarray,
                                 np.ndarray, int]:
        assert uinv.shape == (dim_ens, dim_ens), \
            f"Expected uinv.shape to be {(dim_ens, dim_ens)}"\
                f" but got {uinv.shape}"
        return state_p, uinv, ens_p, status_pdaf

    dim_p = 20
    dim_ens = 4
    forget = 1.0

    param_int = np.array([dim_p, dim_ens, ], dtype=np.intc)
    param_float = np.array([forget, ])
    if filter_type == 0:
        param_float = np.array([forget, 1e-5])

    _, _, status = PDAF.init(filter_type,
                            subtype,
                            0,
                            param_int,
                            param_float,
                            MPI.COMM_WORLD.py2f(),
                            MPI.COMM_WORLD.py2f(),
                            MPI.COMM_WORLD.py2f(),
                            1, 1, True,
                            init_ens_pdaf, 0)


def test_dim_1(filter_type, subtype):
    import pyPDAF.PDAF as PDAF
    """Test the initialisation when #
    dimension of `uinv` is (1, 1).
    """
    def init_ens_pdaf(filtertype:int, dim_p:int, dim_ens:int,
                      state_p:np.ndarray, uinv:np.ndarray,
                      ens_p:np.ndarray, status_pdaf:int
                      ) -> tuple[np.ndarray, np.ndarray,
                                 np.ndarray, int]:
        assert uinv.shape == (1, 1), \
            f"Expected uinv.shape to be {(1, 1)}"\
                f" but got {uinv.shape}"
        return state_p, uinv, ens_p, status_pdaf

    dim_p = 20
    dim_ens = 4
    if filter_type == 200 and subtype == 0:
        dim_ens = 1
    forget = 1.0

    param_int = np.array([dim_p, dim_ens, ], dtype=np.intc)
    param_float = np.array([forget, ])
    if filter_type == 0:
        param_float = np.array([forget, 1e-5])

    _, _, status = PDAF.init(filter_type,
                            subtype,
                            0,
                            param_int,
                            param_float,
                            MPI.COMM_WORLD.py2f(),
                            MPI.COMM_WORLD.py2f(),
                            MPI.COMM_WORLD.py2f(),
                            1, 1, True,
                            init_ens_pdaf, 0)

import sys
filter_type = int(sys.argv[1])
subtype = int(sys.argv[2])
if filter_type in [1, 3, 6, 7, 200, 200, 200, 200]:
    if filter_type == 200 and subtype != 0:
        test_dim_ens_1(filter_type, subtype)
if filter_type in [0, 4, 5, 9, 10, 11]:
    test_dim_ens(filter_type, subtype)
if filter_type in [2, 8, 12, 100, 200]:
    if filter_type == 200 and subtype == 0:
        test_dim_1(filter_type, subtype)

