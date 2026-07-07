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
import configparser

import numpy as np
import pyPDAF

import collector
import filter_options
import io_nemo
import local
import model
import obs_factory
import parallel
import prepost
import sfields


class DAsystem:
    """Wire together model IO, state fields, observations, and PDAF callbacks.

    This class is the high-level orchestration layer. To adapt the workflow to
    another complex model, replace the NEMO-specific reader/writer/domain
    pieces while keeping the PDAF callback shape used here.
    """
    def __init__(self, pe:parallel.Parallel, config:configparser.ConfigParser) -> None:
        self.pe = pe
        self.reader = io_nemo.Reader(pe, config['io'])
        self.writer = io_nemo.Writer(pe, config['io'])
        self.model = model.NemoDomain.from_input(config['model'].getint('numcat', 0),
                                                 self.reader, pe)

        self.filter_options = filter_options.FilterOptions(config['filter_options'])
        self.sfields = sfields.Sfields(config, self.model.wet_grid, self.pe)
        self.local = local.Local(self.model, self.sfields.fields)
        self.obs = obs_factory.ObsFactory(config, self.model, self.sfields.fields)

    def init_pdaf(self, screen:int) -> None:
        """Initialise PDAF, PDAF-OMI, filter settings, and local-domain support.

        Parameters
        ----------
        screen : int
            verbosity of PDAF screen output
        """
        filter_param_i = np.array([max(self.sfields.dim_state_p, 1),
                                   self.pe.dim_ens], dtype=np.intc)
        filter_param_r = np.array([self.filter_options.forget, ])
        status = 0
        cltor = collector.Collector(self.reader, self.model.wet_grid,
                                    self.sfields.fields)
        # initialise PDAF filters, communicators, ensemble
        _, _, status = pyPDAF.init(self.filter_options.filtertype,
                                   self.filter_options.subtype,
                                   0, filter_param_i,
                                   len(filter_param_i),
                                   filter_param_r, 1,
                                   cltor.init_ens_pdaf, screen)
        assert status == 0, f'ERROR {status} \
                in initialization of PDAF - stopping! \
                (PE f{self.pe.mype_filter})'

        pyPDAF.PDAFomi.init(self.obs.nobs)

        lfilter = pyPDAF.PDAF.get_localfilter()
        self.local.local_filter = lfilter == 1
        if lfilter == 1:
            self.local.set_domain_limits()


    def assimilate(self) -> None:
        """offline assimilation function.
        """
        status = 0
        pre_post = prepost.Prepost(self.writer, self.model, self.sfields.fields)
        pyPDAF.assim_offline(self.obs.init_dim_obs_pdafomi,
                             self.obs.obs_op_pdafomi,
                             self.local.init_n_domains_pdaf,
                             self.local.init_dim_l_pdaf,
                             self.obs.init_dim_obs_l_pdafomi,
                             pre_post.prepostprocess, status)


        assert status == 0, f'ERROR {status} in pyPDAF.assimiate' \
            f' - stopping! (PE {self.pe.mype_filter})'

    def finalise(self) -> None:
        """finalise PDAF system"""
        pyPDAF.PDAF.print_info(11)
        if self.pe.mype_filter == 0:
            pyPDAF.PDAF.print_info(3)
        pyPDAF.PDAF.deallocate()
