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
for localization in Cython syntax. 
"""
cdef void c__init_dim_l_pdaf(int* step, int* domain_p, int* dim_l);

cdef void c__init_n_domains_pdaf(int* step, int* n_domains_p);

cdef void c__g2l_state_pdaf(int* step, int* domain_p, int* dim_p, 
                            double* state_p, int* dim_l, 
                            double* state_l);

cdef void c__l2g_state_pdaf(int* step, int* domain_p, int* dim_l,
                            double* state_l, int* dim_p, 
                            double* state_p);