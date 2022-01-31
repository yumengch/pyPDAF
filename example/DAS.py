import numpy as np
import pyPDAF.PDAF.PDAFomi as PDAFomi
from Localization import Localization
import PDAF_caller


class FilterOptions:
    def __init__(self, filtertype, subtype):
        # Select filter algorithm
        self.filtertype  = filtertype
        # Subtype of filter algorithm
        self.subtype = subtype

    def setTransformTypes(self, type_trans, type_sqrt, incremental, 
                            covartype, rank_analysis_enkf):
        # Type of ensemble transformation 
        self.type_trans = type_trans
        # Type of transform matrix square-root
        # (0) symmetric square root, (1) Cholesky decomposition
        self.type_sqrt = type_sqrt
        # (1) to perform incremental updating (only in SEIK/LSEIK!)
        self.incremental = incremental
        # Definition of factor in covar. matrix used in SEIK
        # (0) for dim_ens^-1 (old SEIK)
        # (1) for (dim_ens-1)^-1 (real ensemble covariance matrix)
        # This parameter has also to be set internally in PDAF_init.
        self.covartype = covartype
        # rank to be considered for inversion of HPH in analysis of EnKF;
        # (0) for analysis w/o eigendecomposition
        # if set to >=ensemble size, it is reset to ensemble size - 1
        self.rank_analysis_enkf = rank_analysis_enkf


class Inflation:
    def __init__(self, type_forget, forget):
        # Type of forgetting factor
        # (0) fixed
        # (1) global adaptive
        # (2) local adaptive for LSEIK/LETKF/LESTKF
        self.type_forget = type_forget
        # forgeting factor
        self.forget  = forget


class AssimilationDimensions:
    def __init__(self, model, dim_ens):
        self.dim_state_p = np.prod(model.nx_p)
        self.dim_state=np.prod(model.nx)
        self.dim_ens = dim_ens


class DAS:
    def __init__(self, pe, model, obs):
        self.pe = pe
        self.model = model
        self.obs = obs

    def init(self):
        # init model
        self.model.init_field('inputs_online/true_initial.txt', 
                                                   self.pe.mype_model)

        # init observations
        PDAFomi.init(len(self.obs))

        # Initialize PDAF
        self.assim_dim = AssimilationDimensions(self.model,
                                     dim_ens=self.pe.n_modeltasks)
        self.filter_options = FilterOptions(filtertype=6, subtype=0)
        self.filter_options.setTransformTypes(type_trans=0, 
                                              type_sqrt=0, 
                                              incremental=0, 
                                              covartype=1, 
                                              rank_analysis_enkf=0)
        self.infl = Inflation(type_forget=0, forget=1.)
        self.localization = Localization(loc_weight=0, local_range=0, 
                                            srange=0)
        PDAF_caller.init_pdaf(self.assim_dim, self.infl, 
                                        self.filter_options, 
                                        self.localization,
                                        self.model, self.pe, 
                                        self.obs, 0)

    def forward(self, istep):
        self.model.step(self.pe, istep)
        PDAF_caller.assimilate_pdaf(self.model, self.obs, self.pe, 
                        self.assim_dim, self.localization, 
                        self.filter_options.filtertype)

