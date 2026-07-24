# Data assimilation workflow

This page gives a conceptual map for readers who are new to data assimilation,
PDAF, or pyPDAF. It is intentionally practical: the goal is to explain what the
main objects mean and where user code fits into a pyPDAF program.

## What data assimilation does

Data assimilation combines three sources of information:

- a **model forecast**, which predicts the system state;
- **observations**, which measure part of the real system;
- estimates of **uncertainty**, which decide how strongly the analysis should
  trust the forecast or the observations.

The output is an **analysis**: an improved estimate of the state that can be
used for diagnostics or as the initial condition for the next forecast.

In ensemble methods, pyPDAF does this for an ensemble of model states. The spread
of the ensemble represents forecast uncertainty. Observations then adjust the
ensemble mean, and often the ensemble perturbations, according to the selected
PDAF filter.

## The pyPDAF division of responsibilities

PDAF provides the data assimilation algorithms. pyPDAF exposes those algorithms
to Python and lets user code provide model-specific information through callback
functions.

Your application normally supplies:

- how to initialise the model state and ensemble;
- how to convert model fields to a one-dimensional state vector;
- how to distribute an analysis state vector back into the model;
- what observations are available at the current analysis time;
- how the observation operator maps state variables to observation space;
- observation error information and, for local filters, localisation metadata.

pyPDAF handles:

- MPI communicator setup and ensemble organisation;
- calls into the selected PDAF algorithm;
- the timing of callback calls;
- observation handling through PDAFomi;
- analysis update and diagnostic helper routines.

## Online and offline modes

In **online mode**, the model and the assimilation code are part of the same
program. The model advances the ensemble, pyPDAF collects the forecast state in
memory, computes the analysis, and returns the updated state to the model. This
is usually the most efficient option when the model is written in Python or can
be controlled directly from Python.

In **offline mode**, the model run and the assimilation program are separate.
The model writes forecast states to files, pyPDAF reads those states, computes
the analysis, and writes updated states back to disk. This is useful when a
model is difficult to couple directly or when pyPDAF is being added to an
existing workflow.

## A typical online cycle

The high-level pattern is:

1. Configure MPI communicators with `pyPDAF.set_parallel` or
   `pyPDAF.init_parallel`.
2. Initialise PDAF with `pyPDAF.init`, passing the filter choice, dimensions,
   parameters, and callback functions.
3. Initialise PDAFomi with `pyPDAF.PDAFomi.init` for observation handling.
4. Call `pyPDAF.init_forecast` so PDAF can distribute initial ensemble states to the model tasks.
5. Advance each ensemble member with the model until the next observation time.
6. Call an assimilation driver such as `pyPDAF.assimilate`.
7. Repeat forecast and analysis cycles.
8. Call `pyPDAF.deallocate` before the program exits.

The assimilation driver calls the user-supplied functions at the right time. For
example, the driver may ask for observations, apply the observation operator,
compute local observation dimensions, or distribute the analysis state.

## A typical offline cycle

The offline workflow uses the same concepts, but the model state is exchanged
through files instead of in-memory model variables:

1. Run the model ensemble and write forecast restart files.
2. Start the pyPDAF offline assimilation program.
3. Read model states in the callback that collects the forecast state.
4. Provide observations and observation operators through PDAFomi callbacks.
5. Call an offline driver such as `pyPDAF.assim_offline`.
6. Write analysis states or restart files from the callback that distributes the
   analysis state.
7. Run the model ensemble again from the new analysis.

## State vectors and observations

PDAF works with one-dimensional state vectors. A model, however, usually stores
fields as structured arrays, for example temperature on a grid or variables on
different vertical levels. The user-supplied functions define the conversion
between model fields and PDAF state vectors.

Observations live in observation space. The observation operator maps a model
state to the quantities that would be observed. For example, it may select a
grid point, interpolate a gridded field to a station location, or compute a
derived quantity. PDAFomi provides helper functions for common observation
operators, interpolation coefficients, and localisation.

## Choosing where to start

For a first implementation, keep the system small:

- use a serial or lightly parallel example first;
- choose one observation type;
- start with an established ensemble filter such as an ESTKF or LETKF-style
  setup;
- implement the callback functions with explicit shape checks;
- add localisation only after the global or simple local case is working;
- compare diagnostics before and after assimilation to verify that the analysis
  changes in the expected direction.

After the basic cycle works, the API reference and user-function pages become
the main guides for adding observation types, diagnostics, localisation, and
parallel decomposition.

## Inspecting PDAF filter options

PDAF uses integer identifiers for filter types and subtypes. pyPDAF exposes the
PDAF information routines directly:

```python
import pyPDAF

pyPDAF.print_version()
pyPDAF.print_filter_types(verbose=0)
pyPDAF.options_filters(filtertype)
```

Use `pyPDAF.print_version()` to check the linked PDAF version,
`pyPDAF.print_filter_types(verbose=0)` to list available filter type numbers,
and `pyPDAF.options_filters(filtertype)` to inspect the options for one filter.

These routines are also installed as command-line tools:

```bash
pypdaf-version
pypdaf-filter-types --verbose 0
pypdaf-filter-options 7
```
