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
import logging

import config
import model
import model_integrator
from parallelisation import parallelisation
from PDAF_system import PDAF_system

def main():
    pe = parallelisation(dim_ens=4, n_modeltasks=4)

    # Initial Screen output
    if pe.mype_ens == 0:
        logging.info('+++++ PDAF online mode +++++')
        logging.info('2D model with parallelization')

    # Initialise model
    # throughout this example and PDAF, we must assume that each ensemble member
    # uses the same domain decomposition.
    model_ens = [model.model(pe=pe) for i in range(pe.dim_ens_l)]
    model_ens[0].print_info(pe)

    if not config.USE_PDAF:
        model_ens.init_field(pe.mype_model)

    # Initialise model integrator
    integrator = model_integrator.model_integrator(model_ens)

    # Initialise PDAF system
    das = PDAF_system(pe, model_ens)
    if config.USE_PDAF:
        das.init_pdaf(screen=2)

    integrator.forward(config.nsteps, das)

    pe.finalize_parallel()

if __name__ == '__main__':
    main()
