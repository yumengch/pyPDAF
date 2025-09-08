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
import config


class FilterOptions:
    """Here, we provide all filter options

    Attributes
    ----------
    filtertype : int
        the type of filter used
        1=SEIK, 2=EnKF, 3=LSEIK, 4=ETKF, 5=LETKF, 6=ESTKF, 7=LESTKF
        8=LEnKF, 9=NETF, 10=LNETF, 11=LKNETF, 12=PF, 100=GENOBS,
        200=3DVar, 0=SEEK
        For a simplified documentation, see:
        https://pdaf.awi.de/trac/wiki/AvailableOptionsforInitPDAF
        Different DA scheme requires different user-supplied functions
        More information can be found in the PDAF documentation
    subtype : int
        Variants of each DA scheme check https://pdaf.awi.de/trac/wiki/AvailableOptionsforInitPDAF
    type_forget : int
        type of forgetting factor
        - (0) fixed
        - (1) global adaptive
        - (2) local adaptive for LSEIK/LETKF/LESTKF
    forget : float
        forgetting factor
    type_trans : int
        type of ensemble transformation
    type_sqrt : int
        type of transform matrix square-root
    incremental : int
        (1) to perform incremental updating (only in SEIK/LSEIK!)
    covartype : int
        definition of factor in covar. matrix
    rank_analysis_enkf : int
        rank to be considered for inversion of HPH in analysis of EnKF
    """

    def __init__(self) -> None:
        # Select filter algorithm
        self.filtertype:int = config.filtertype
        # Subtype of filter algorithm
        self.subtype:int = config.subtype

        self.type_forget:int = config.type_forget
        # forgeting factor
        self.forget:float = config.forget

        # Set other parameters to default values
        self.type_trans:int = config.type_trans
        self.type_sqrt:int = config.type_sqrt
        self.incremental:int = config.incremental
        self.covartype:int = config.covartype

        self.rank_analysis_enkf:int = config.rank_analysis_enkf
