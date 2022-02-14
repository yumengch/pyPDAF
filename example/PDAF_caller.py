import numpy as np
import pyPDAF.PDAF as PDAF
import functools
import U_PDAF
import U_PDAFomi


class init_pdaf:
    def __init__(self, assim_dim, infl, filter_options, 
                    localization, model, pe, obs, screen):

        if (filter_options.filtertype == 2):
            # EnKF with Monte Carlo init
            self.setEnKFOptions(6, 2, assim_dim, infl, filter_options)
        else:
            # All other filters
            self.setETKFOptions(7, 2, assim_dim, infl, filter_options)

        U_init_ens_pdaf = functools.partial(
                            U_PDAF.init_ens_pdaf, 
                                model, pe, assim_dim)

        status = PDAF.init(filter_options.filtertype, 
                                 filter_options.subtype, 
                                 self.filter_param_i,  
                                 self.filter_param_r, 
                                 pe.COMM_model.py2f(), 
                                 pe.COMM_filter.py2f(), 
                                 pe.COMM_couple.py2f(), pe.task_id, 
                                 pe.n_modeltasks, pe.filterpe, 
                                 U_init_ens_pdaf, screen)
        try:
            assert status == 0, \
                f'ERROR {status} \
                in initialization of PDAF - stopping! \
                (PE f{pe.mype_world})'
        except AssertionError:
            pe.abort_parallel()

        U_next_observation_pdaf = functools.partial(
                                        U_PDAF.next_observation_pdaf,
                                            model, pe, obs[0].delt_obs)
        U_distribute_state_pdaf = functools.partial(
                                        U_PDAF.distribute_state_pdaf,
                                            model)
        U_prepoststep_ens_pdaf = functools.partial(
                                        U_PDAF.prepoststep_ens_pdaf,
                                            assim_dim, model, pe, obs)

        PDAF.get_state(U_next_observation_pdaf,
                        U_distribute_state_pdaf,
                        U_prepoststep_ens_pdaf)

        localization.set_lim_coords(model.nx_p, pe)

    def setEnKFOptions(self, dim_pint, dim_preal, 
                             assim_dim, infl, filter_options):
        self.filter_param_i = np.zeros(dim_pint, dtype=int)
        self.filter_param_r = np.zeros(dim_preal)

        self.filter_param_i[0] = assim_dim.dim_state_p
        self.filter_param_i[1] = assim_dim.dim_ens
        self.filter_param_i[2] = filter_options.rank_analysis_enkf
        self.filter_param_i[3] = filter_options.incremental
        self.filter_param_i[4] = 0

        self.filter_param_r[0] = infl.forget

    def setETKFOptions(self, dim_pint, dim_preal, 
                             assim_dim, infl, filter_options):
        self.filter_param_i = np.zeros(dim_pint, dtype=int)
        self.filter_param_r = np.zeros(dim_preal)

        self.filter_param_i[0] = assim_dim.dim_state_p
        self.filter_param_i[1] = assim_dim.dim_ens
        self.filter_param_i[2] = 0
        self.filter_param_i[3] = filter_options.incremental
        self.filter_param_i[4] = infl.type_forget
        self.filter_param_i[5] = filter_options.type_trans
        self.filter_param_i[6] = filter_options.type_sqrt

        self.filter_param_r[0] = infl.forget


class assimilate_pdaf:
    def __init__(self, model, obs, pe, assim_dim, localization, filtertype):
        localfilter = PDAF.get_localfilter()

        U_collect_state_pdaf = functools.partial(
                                        U_PDAF.collect_state_pdaf,
                                        model, assim_dim) 
        U_next_observation_pdaf = functools.partial(
                                        U_PDAF.next_observation_pdaf,
                                            model, pe, obs[0].delt_obs)
        U_distribute_state_pdaf = functools.partial(
                                        U_PDAF.distribute_state_pdaf,
                                            model)
        U_prepoststep_ens_pdaf = functools.partial(
                                        U_PDAF.prepoststep_ens_pdaf,
                                            assim_dim, model, pe, obs)
        U_init_dim_obs_PDAFomi = functools.partial(
                                        U_PDAFomi.init_dim_obs_pdafomi, 
                                        obs, 
                                        localization.local_range, 
                                        pe.mype_filter, 
                                        model.nx, 
                                        model.nx_p)

        U_obs_op_PDAFomi = functools.partial(
                                        U_PDAFomi.obs_op_pdafomi, 
                                        obs)

        if (localfilter==1):
            U_init_n_domains_pdaf = functools.partial(
                                        localization.init_n_domains_pdaf,
                                                        assim_dim)
            U_init_dim_l_pdaf = functools.partial(
                                        localization.init_dim_l_pdaf, 
                                        model.nx_p, pe.mype_filter)
            U_init_dim_obs_l_pdafomi = functools.partial(
                                        U_PDAFomi.init_dim_obs_l_pdafomi, 
                                                    obs, localization)
            status = PDAF.PDAFomi_assimilate_local(U_collect_state_pdaf,
                             U_distribute_state_pdaf,
                             U_init_dim_obs_PDAFomi,
                             U_obs_op_PDAFomi,
                             U_prepoststep_ens_pdaf,
                             U_init_n_domains_pdaf,
                             U_init_dim_l_pdaf,
                             U_init_dim_obs_l_pdafomi,
                             localization.g2l_state_pdaf,
                             localization.l2g_state_pdaf,
                             U_next_observation_pdaf)
        else:
            if filtertype!=8:
                status = PDAF.PDAFomi_assimilate_global(U_collect_state_pdaf,
                                              U_distribute_state_pdaf,
                                              U_init_dim_obs_PDAFomi,
                                              U_obs_op_PDAFomi,
                                              U_prepoststep_ens_pdaf,
                                              U_next_observation_pdaf)
            else:
                U_localize_covar_pdafomi = functools.partial(
                                                U_PDAFomi.localize_covar_pdafomi,
                                                obs, localization, 
                                                mype_filter, nx_p)
                status = PDAF.PDAFomi_assimilate_lenkf(U_collect_state_pdaf,
                                              U_distribute_state_pdaf,
                                              U_init_dim_obs_PDAFomi,
                                              U_obs_op_PDAFomi,
                                              U_prepoststep_ens_pdaf,
                                              U_localize_covar_pdafomi,
                                              U_next_observation_pdaf)

        if status != 0:
            print(('ERROR ', status, 
                    ' in PDAF_put_state - stopping! (PE ', 
                        pe.mype_world,')'))
            pe.abort_parallel()