pyPDAF
======

pyPDAF is a Python interface to the `Parallel Data Assimilation Framwork (PDAF) <http://pdaf.awi.de/trac/wiki>`_ library written in Fortran. The latest pyPDAF supports PDAF-V2.1.

With a variety of packages in Python, it allows a simpler coding style for user-supplied functions, such as I/O of observations and post-processing. This is helpful for prototyping data assimilation systems, offline data assimilation systems. It can also benefit many Python-based numerical models, or models that can be interfaced with Python, with parallel and efficient data assimilation capability.

The core DA algorithm is as efficient as Fortran implementation in the interface. The efficiency of the Python-based user supplied functions can be improved if sufficient optimisations are used.