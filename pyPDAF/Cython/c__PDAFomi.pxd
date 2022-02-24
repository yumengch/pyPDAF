"""This file is part of pyPDAF

Copyright (C) 2022 University of Reading and National Centre for Earth Observation

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

Declaration of the user-defined PDAFomi subroutines
in Cython syntax. 
"""
cdef void c__init_dim_obs_PDAFomi(int* step, int* dim_obs);
cdef void c__obs_op_PDAFomi(int* step, int* dim_p, int* dim_obs, 
                            double* state_p, double* ostate);
cdef void c__init_dim_obs_l_PDAFomi(int* domain_p, int* step,
                                    int* dim_obs, int* dim_obs_l);
cdef void c__localize_covar_PDAFomi(int* dim_p, int* dim_obs, 
                                    double* hp_p, double* hph);