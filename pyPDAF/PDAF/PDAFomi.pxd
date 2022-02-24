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

Cython declarations of PDAFomi Fortran iso_C_binding wrappers.
"""
cdef extern void c__init_pdafomi(int* n_obs);
cdef extern void c__set_pdafomi_doassim(int* i_obs, int* doassim);
cdef extern void c__set_pdafomi_disttype(int* i_obs, int* disttype);
cdef extern void c__set_pdafomi_ncoord(int* i_obs, int* ncoord);
cdef extern void c__set_pdafomi_id_obs_p(int* i_obs, int* nrows, 
                                         int* dim_obs_p, 
                                         int* id_obs_p);
cdef extern void c__set_pdafomi_icoeff_p(int* i_obs, int* nrows, 
                                         int* dim_obs_p, 
                                         double* icoeff_p);
cdef extern void c__set_pdafomi_domainsize(int* i_obs, int* ncoord, 
                                           double* domainsize);
cdef extern void c__set_pdafomi_obs_err_type(int* i_obs, 
                                             int* obs_err_type);
cdef extern void c__set_pdafomi_use_global_obs(int* i_obs, 
                                               int* use_global_obs);
cdef extern void c__pdafomi_gather_obs(int* i_obs, int* nrows, 
                                       int* dim_obs_p, double* obs_p, 
                                       double* ivar_obs_p, 
                                       double* ocoord_p, 
                                       double* local_range, 
                                       int* dim_obs);
cdef extern void c__pdafomi_set_domain_limits(double* lim_coords);
cdef extern void c__pdafomi_obs_op_gridpoint(int* i_obs, int* dim_p, 
                                             int* dim_obs, 
                                             double* state_p, 
                                             double* ostate);
cdef extern void c__pdafomi_init_dim_obs_l(int* i_obs, 
                                           double* coords_l, 
                                           int* locweight, 
                                           double* local_range, 
                                           double* srange,
                                           int* dim_obs_l);
cdef extern void c__pdafomi_localize_covar(int* i_obs, int* dim_p,
                                           int* dim_obs,
                                           int* dim_coords,
                                           int* locweight, 
                                           double* local_range, 
                                           double* srange, 
                                           double* coords_p, 
                                           double* hp_p, double* hph);
cdef extern void c__pdafomi_deallocate_obs(int* i_obs, int* step);
