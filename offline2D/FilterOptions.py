"""This file is part of pyPDAF

Copyright (C) 2022 University of Reading and
National Centre for Earth Observation

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""


class FilterOptions:

    """Summary

    Attributes
    ----------
    covartype : int
        definition of factor in covar. matrix
    filtertype : int
        filter type
    incremental : int
        (1) to perform incremental updating (only in SEIK/LSEIK!)
    rank_analysis_enkf : int
        rank to be considered for inversion of HPH in analysis of EnKF
    subtype : int
        sub type of a filter
    type_sqrt : int
        type of transform matrix square-root
    type_trans : int
        type of ensemble transformation
    """

    def __init__(self, filtertype, subtype):
        """Summary

        Parameters
        ----------
        filtertype : int
            filter type
        subtype : int
            subtype of the filter
        """
        # Select filter algorithm
        self.filtertype = filtertype
        # Subtype of filter algorithm
        self.subtype = subtype

    def setTransformTypes(self, type_trans, type_sqrt, incremental,
                          covartype, rank_analysis_enkf):
        """Summary

        Parameters
        ----------
        type_trans : int
            type of ensemble transformation
        type_sqrt : int
            Type of transform matrix square-root
            (0) symmetric square root, (1) Cholesky decomposition
        incremental : int
            (1) to perform incremental updating (only in SEIK/LSEIK!)
        covartype : int
            Definition of factor in covar. matrix used in SEIK
            - (0) for dim_ens^-1 (old SEIK)
            - (1) for (dim_ens-1)^-1 (real ensemble covariance matrix)

            This parameter has also to be set internally in PDAF_init.
        rank_analysis_enkf : int
            rank to be considered for inversion of HPH in analysis of EnKF

            (0) for analysis w/o eigendecomposition
            if set to >=ensemble size, it is reset to ensemble size - 1
        """
        self.type_trans = type_trans
        self.type_sqrt = type_sqrt
        self.incremental = incremental
        self.covartype = covartype

        self.rank_analysis_enkf = rank_analysis_enkf
