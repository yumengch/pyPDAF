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
from parallelisation import parallelisation
from PDAF_system import PDAF_system

def main():
    pe = parallelisation(dim_ens=config.dim_ens)

    # Initial Screen output
    if pe.mype_ens == 0:
        log.logger.info('+++++ PDAF offline mode +++++')
        log.logger.info('2D model with parallelization')

    # Initialise model grid for localisation/observation operator
    model_grid = model.model_grid(pe=pe)
    model_grid.print_info(pe)

    # Initialise PDAF system
    das = PDAF_system(pe, model_grid)

    das.init_pdaf(screen=config.screen)

    das.assimilate()

    das.finalise()

    pe.finalize_parallel()

if __name__ == '__main__':
    main()
