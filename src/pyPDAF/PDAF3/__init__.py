"""Namespace for PDAF3 module."""
from ._pdaf3_c import init, init_forecast, set_parallel
from .assim import assimilate, assim_offline, \
                   assimilate_3dvar_all, assim_offline_3dvar_all, \
                   assimilate_local_nondiagr, assimilate_global_nondiagr, \
                   assimilate_lnetf_nondiagr, assimilate_lknetf_nondiagr, \
                   assimilate_enkf_nondiagr, assimilate_nonlin_nondiagr, \
                   assimilate_3dvar_nondiagr, assimilate_en3dvar_estkf_nondiagr, \
                   assimilate_en3dvar_lestkf_nondiagr, \
                   assimilate_hyb3dvar_estkf_nondiagr, \
                   assimilate_hyb3dvar_lestkf_nondiagr, \
                   assim_offline_local_nondiagr, \
                   assim_offline_global_nondiagr, \
                   assim_offline_lnetf_nondiagr, \
                   assim_offline_lknetf_nondiagr, \
                   assim_offline_enkf_nondiagr, \
                   assim_offline_lenkf_nondiagr, \
                   assim_offline_nonlin_nondiagr, \
                   assim_offline_3dvar_nondiagr, \
                   assim_offline_en3dvar_estkf_nondiagr, \
                   assim_offline_en3dvar_lestkf_nondiagr, \
                   assim_offline_hyb3dvar_estkf_nondiagr, \
                   assim_offline_hyb3dvar_lestkf_nondiagr, \
                   generate_obs, generate_obs_offline

