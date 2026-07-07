"""Rank-aware logging helpers for the offline pyPDAF/NEMO workflow.

This module wraps Python logging so messages can be limited to selected MPI
ranks and PDAF screen levels.

This file is part of pyPDAF

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

_base_logger = logging.getLogger()
_base_logger.setLevel(logging.DEBUG)
if not _base_logger.handlers:
    ch = logging.StreamHandler()
    formatter = logging.Formatter('%(levelname)s: %(message)s')
    ch.setFormatter(formatter)
    _base_logger.addHandler(ch)


class PDAFLogger:
    """Logger wrapper with PDAF prefixing and rank-gated output."""

    def __init__(self, base_logger, prefix="NEMO-PDAF", pe=None, ranks=(0,), screen=None):
        self.base_logger = base_logger
        self.prefix = prefix
        self.pe = pe
        self.ranks = ranks
        self.screen = screen

    def configure(self, pe=None, ranks=None, prefix=None, screen=None):
        if pe is not None:
            self.pe = pe
        if ranks is not None:
            self.ranks = ranks
        if prefix is not None:
            self.prefix = prefix
        if screen is not None:
            self.screen = screen

    def should_log(self, ranks=None, screen=None):
        if self.screen > 3:
            return True
        if screen is not None and self.screen <= screen:
            return False
        rank = self.pe.mype_filter
        ranks_l = self.ranks if ranks is None else ranks
        if isinstance(ranks_l, int):
            return rank == ranks_l
        return rank in ranks_l

    def format(self, message, prefix=True):
        message = str(message)
        if not prefix or not self.prefix or message.startswith(self.prefix):
            return message
        return f"{self.prefix} {message}"

    def info(self, message, ranks=None, screen=None):
        if self.should_log(ranks=ranks, screen=screen):
            self.base_logger.info(self.format(message, prefix=self.prefix))

    def warning(self, message, ranks=None, screen=None):
        if self.should_log(ranks=ranks, screen=screen):
            self.base_logger.warning(self.format(message, prefix=self.prefix))

    def error(self, message, ranks=None, screen=None):
        if self.should_log(ranks=ranks, screen=screen):
            self.base_logger.error(self.format(message, prefix=self.prefix))

    def debug(self, message, ranks=None, screen=None):
        if self.should_log(ranks=ranks, screen=screen):
            self.base_logger.debug(self.format(message, prefix=self.prefix))

    def __getattr__(self, name):
        return getattr(self.base_logger, name)


logger = PDAFLogger(_base_logger)


def configure(pe=None, ranks=None, prefix=None, screen=None):
    logger.configure(pe=pe, ranks=ranks, prefix=prefix, screen=screen)


def set_printing_process(pe, ranks=0, screen=0):
    """Choose which MPI process rank(s) emit log messages.

    ``ranks`` can be an integer, an iterable of integers.
    """
    configure(pe=pe, ranks=ranks, screen=screen)


def should_log(ranks=None, screen=None):
    return logger.should_log(ranks=ranks, screen=screen)


def info(message, ranks=None, screen=None):
    logger.info(message, ranks=ranks, screen=screen)


def warning(message, ranks=None, screen=None):
    logger.warning(message, ranks=ranks, screen=screen)


def error(message, ranks=None, screen=None):
    logger.error(message, ranks=ranks, screen=screen)


def debug(message, ranks=None, screen=None):
    logger.debug(message, ranks=ranks, screen=screen)
