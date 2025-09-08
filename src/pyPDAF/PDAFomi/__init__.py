"""PDAFomi module"""
from ._pdafomi_c import init, init_local, check_error, gather_obs, \
                        gather_obsstate, get_interp_coeff_tri, \
                        get_interp_coeff_lin1d, get_interp_coeff_lin, \
                        init_dim_obs_l_iso, init_dim_obs_l_noniso, \
                        init_dim_obs_l_noniso_locweights, obs_op_gridpoint, \
                        obs_op_gridavg, obs_op_extern, obs_op_interp_lin, \
                        obs_op_adj_gridpoint, obs_op_adj_gridavg, \
                        obs_op_adj_interp_lin, observation_localization_weights, \
                        set_debug_flag, set_dim_obs_l, set_localization, \
                        set_localization_noniso, set_localize_covar_iso, \
                        set_localize_covar_noniso, \
                        set_localize_covar_noniso_locweights, \
                        set_obs_diag, set_domain_limits, \
                        get_domain_limits_unstr, store_obs_l_index, \
                        store_obs_l_index_vdist
from .diag import diag_dimobs, diag_get_hx, diag_get_hxmean, \
                  diag_get_ivar, diag_get_obs, diag_nobstypes, \
                  diag_obs_rmsd, diag_stats
from .setter import set_doassim, set_disttype, set_ncoord, \
                    set_obs_err_type, set_use_global_obs, \
                    set_inno_omit, set_inno_omit_ivar, \
                    set_id_obs_p, set_icoeff_p, set_domainsize, set_name