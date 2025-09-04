# pylint: disable=no-name-in-module
"""Namespace for pyPDAF.PDAF package
"""
from pyPDAF.PDAF._pdaf_c import correlation_function, deallocate, eofcovar, \
                                force_analysis, gather_dim_obs_f, gather_obs_f, \
                                gather_obs_f2, gather_obs_f_flex, \
                                gather_obs_f2_flex, init, init_forecast, \
                                local_weight, local_weights, print_filter_types, \
                                print_da_types, print_info, reset_forget, sample_ens, \
                                get_fcst_info
from pyPDAF.PDAF.setter import set_comm_pdaf, set_debug_flag, set_ens_pointer, \
                               set_iparam, set_memberid, set_offline_mode, \
                               set_rparam, set_seedset, set_smoother_ens
from pyPDAF.PDAF.get import get_assim_flag, get_localfilter, \
                            get_local_type, get_memberid, get_obsmemberid, \
                            get_smoother_ens
from pyPDAF.PDAF.iau import iau_init, iau_reset, iau_set_pointer
from pyPDAF.PDAF.diag import diag_ensmean, \
                             diag_stddev_nompi, \
                             diag_stddev, \
                             diag_variance_nompi, \
                             diag_variance, \
                             diag_rmsd_nompi, \
                             diag_rmsd, \
                             diag_crps_mpi, \
                             diag_crps_nompi, \
                             diag_effsample, \
                             diag_ensstats, \
                             diag_compute_moments, \
                             diag_histogram, \
                             diag_reliability_budget
