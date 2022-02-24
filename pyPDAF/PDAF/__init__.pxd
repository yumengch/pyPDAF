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

Cython declarations of PDAF Fortran iso_C_binding wrappers.
"""
cdef extern void c__pdaf_init(int* filtertype, int* subtype,
                              int filter_param_i[], int* dim_pint,
                              double filter_param_r[], int* dim_preal,
                              int* comm_model, int* comm_filter,
                              int* comm_couple, int* task_id,
                              int* n_modeltasks, bint* filterpe,
                              void (*c__init_ens_pdaf) (int*, int*, int*,
                                                        double*, double*,
                                                        double*, int*),
                              int* screen, int* status_pdaf);


cdef extern void c__pdaf_get_state(int* steps, double* timenow,
                                   int* doexit,
                                   void (*c__next_observation_pdaf)(int*,
                                                                    int*,
                                                                    int*,
                                                                    double*),
                                   void (*c__distribute_state_pdaf)(int*,
                                                                    double*),
                                   void (*c__prepoststep_ens_pdaf)(int*, int*,
                                                                   int*, int*,
                                                                   int*,
                                                                   double*,
                                                                   double*,
                                                                   double*,
                                                                   int*),
                                   int* status_pdaf);


cdef extern void c__pdaf_get_localfilter(int* localfilter);


cdef extern void c__pdafomi_assimilate_local(
                 void (*c__collect_state_pdaf)(int*, double*),
                 void (*c__distribute_state_pdaf)(int*, double*),
                 void (*c__init_dim_obs_PDAFomi)(int*, int*),
                 void (*c__obs_op_PDAFomi)(int*, int*, int*,
                                           double*, double*),
                 void (*c__prepoststep_ens_pdaf)(int*, int*, int*, int*, int*,
                                                 double*, double*, double*, 
                                                 int*),
                 void (*c__init_n_domains_pdaf)(int*, int*),
                 void (*c__init_dim_l_pdaf)(int*, int*, int*),
                 void (*c__init_dim_obs_l_PDAFomi)(int*, int*, int*, int*),
                 void (*c__g2l_state_pdaf)(int*, int*, int*, double*, int*, 
                                           double*),
                 void (*c__l2g_state_pdaf)(int*, int*, int*, double*,
                                           int*, double*),
                 void (*c__next_observation_pdaf)(int*, int*, int*, double*),
                 int* status_pdaf);


cdef extern void c__pdafomi_assimilate_global(
                 void (*c__collect_state_pdaf)(int*, double*),
                 void (*c__distribute_state_pdaf)(int*, double*),
                 void (*c__init_dim_obs_PDAFomi)(int*, int*),
                 void (*c__obs_op_PDAFomi)(int*, int*, int*,
                                           double*, double*),
                 void (*c__prepoststep_ens_pdaf)(int*, int*, int*, int*, int*,
                                                 double*, double*, double*,
                                                 int*),
                 void (*c__next_observation_pdaf)(int*, int*, int*, double*),
                 int* status_pdaf);


cdef extern void c__pdafomi_assimilate_lenkf(
                 void (*c__collect_state_pdaf)(int*, double*),
                 void (*c__distribute_state_pdaf)(int*, double*),
                 void (*c__init_dim_obs_PDAFomi)(int*, int*),
                 void (*c__obs_op_PDAFomi)(int*, int*, int*,
                                           double*, double*),
                 void (*c__prepoststep_ens_pdaf)(int*, int*, int*, int*, int*,
                                                 double*, double*, double*,
                                                 int*),
                 void (*c__localize_covar_PDAFomi)(int*, int*, double*,
                                                   double*),
                 void (*c__next_observation_pdaf)(int*, int*, int*, double*),
                 int* status_pdaf);