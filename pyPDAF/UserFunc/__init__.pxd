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

Declaration of the user-defined PDAF subroutines
in Cython syntax. 
"""
cdef void c__init_ens_pdaf(int* filtertype, int* dim_p, int* dim_ens, 
	                       double* state_p, double* uinv, 
                           double* ens_p, int* flag)

cdef void c__distribute_state_pdaf(int* dim_p, double* state_p);

cdef void c__collect_state_pdaf(int* dim_p, double* state_p);

cdef void c__next_observation_pdaf(int* stepnow, int* nsteps, 
                                   int* doexit, double* time);

cdef void c__prepoststep_ens_pdaf(int* step, int* dim_p, int* dim_ens,
                                  int* dim_ens_p, int* dim_obs_p, 
                                  double* state_p, double* uinv, 
                                  double* ens_p, int* flag);