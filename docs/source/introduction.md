pyPDAF is a Python interface to the `Parallel Data Assimilation Framwork (PDAF) <http://pdaf.awi.de/trac/wiki>`_ written in Fortran.
The latest pyPDAF supports PDAF-V3.0.

As an interface to PDAF, pyPDAF supports all PDAF functionalities. You can use pyPDAF to construct
a parallel ensemble data assimilation system purely in Python. The pyPDAF is designed as a framework
that defines the workflow of given DA algorithms. Considering the versatility of the software,
information on the model and observations are passed to the DA algorithms through user-supplied
functions. With pyPDAF, all user-supplied functions can be implemented in Python. We expect that
the coding of the user-supplied functions will be easier and more flexible than in Fortran due to
the rich Python ecosystem.

pyPDAF can be used with two modes:
  - online mode: DA is performed without interrupting the model program. Here, the model code is
    extended by calling PDAF functions to generate a single program. In online mode,
    the filtering gets the model state and distributes the analysis to model by in-memory exchange.
    This is the recommended mode for better efficiency.
  - offline mode: DA is performed after the model program is finished. Here, a separate program
    is generated to perform DA. In offline mode,
    the filtering reads the model state from disk and writes the analysis to disk.


The potential applications of pyPDAF include:
  - online DA systems with Python models, e.g., machine-learning models
  - offline DA systems

This is a great tool for researchers who want to test and develop new DA systems.
Compared to Fortran systems, the efficiency is decreased mainly from user-supplied functions
and overhead for array conversions between Fortran and Python. The core DA algorithms are
as efficient as in PDAF. Note that, for computational
intensive user-supplied functions, the efficiency can be improved by using just-in-time compilation
tools such as `numba`.

To get started, we highly recommend to start from the Jupyter notebook
example for
`a serial ensemble DA system using a simple wave model <https://github.com/yumengch/pyPDAF/blob/main/tutorials/tutorial1_serial.ipynb>`_.

We also provide more structured offline and online examples. One can adapt these examples based on their needs
  - `A parallel online ensemble DA system using a simple wave model <https://github.com/yumengch/pyPDAF/tree/main/example/online>`_
  - `A parallel offline ensemble DA system using a simple wave model <https://github.com/yumengch/pyPDAF/tree/main/example/offline>`_
