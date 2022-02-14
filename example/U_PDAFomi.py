
def init_dim_obs_pdafomi(list_of_obs, local_range, 
                          mype_filter, nx, nx_p, step, dim_obs):
    dim_obs = 0
    for obs in list_of_obs:
        if(obs.doassim):
            obs.init_dim_obs(step, dim_obs, local_range, 
                                mype_filter, nx, nx_p)
            dim_obs += obs.dim_obs
    return dim_obs


def obs_op_pdafomi(list_of_obs, step, state_p, ostate):
    for obs in list_of_obs:
        obs.obs_op(step, state_p, ostate)

def init_dim_obs_l_pdafomi(list_of_obs, localization, 
                                domain_p, step, dim_obs, dim_obs_l):
    for obs in list_of_obs:
        obs.init_dim_obs_l(localization, 
                                domain_p, step, dim_obs, dim_obs_l)


def localize_covar_pdafomi(list_of_obs, localization, 
                            mype_filter, nx_p, HP_p, HPH):
    dim_p = HPH.shape[0]
    coords_p = np.zeros((2, dim_p))
    offset = mype_filter*nx_p

    coords_p[0] = np.where(np.ones(nx_p))[1] + offset
    coords_p[1] = np.where(np.ones(nx_p))[0]

    for i_obs_f in list_of_obs_f:
        i_obs_f.localize_covar(localization, HP_p, HPH, coords_p)


def deallocate_obs_pdafomi(list_of_obs, step):
    for obs in list_of_obs:
        obs.deallocate_obs(step)