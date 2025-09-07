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
import log

import config
import model
import model_integrator
from parallelisation import Parallelisation
from pdaf_system import PDAFsystem

def main():
    """main function for the example of online data assimilation"""
    pe = Parallelisation(dim_ens=config.dim_ens, n_modeltasks=config.n_modeltasks)

    # Initial Screen output
    if pe.mype_ens == 0:
        log.logger.info('+++++ PDAF online mode +++++')
        log.logger.info('2D model with parallelization')

    # Initialise model
    # throughout this example and PDAF, we must assume that each ensemble member
    # uses the same domain decomposition.
    model_t = model.Model(pe=pe)
    model_t.print_info(pe)

    # Initialise model integrator
    integrator = model_integrator.ModelIntegrator(model_t)

    # Initialise PDAF system
    das = PDAFsystem(pe, model_t)
    das.init_pdaf(screen=config.screen)

    integrator.forward(config.nsteps, das)

    das.finalise()

    pe.finalize_parallel()

if __name__ == '__main__':
    main()
