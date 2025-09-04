
# Variable naming conventions

The suffix of variables in (py)PDAF follows a few naming conventions.

`_p` typically means process-local variables. In weather and climate models, to
use multiple CPUs, the computational domain is decomposed into many sub-domains.
Each sub-domain is simulated by one processor. This is called domain decomposition.
The `_p` variables typically means that they only contains information of the
sub-domain. This is relevant for implementing observations. For example, `obs_p`
means observations reside in the region of corresponding sub-domains. If the model
is not parallel, `_p` is simply the global domain.

`_l` suffix is related to the local domain of domain localisation of ensemble filters.
This is different from the domain decomposition used for model parallelisation.
In domain localisation, each local domain does data assimilation independently.
One local domain typically only assimilates observations within given localisation radius.
Therefore, `_l` suffix corresponds to each local analysis domain.

`_f` denotes full observations. However, this does not necessarily mean all
observations or global observations. Instead, the meaning of these
variables depends on the `use_global_obs` set by :func:`pyPDAF.PDAFomi.set_use_global_obs`.
  - `use_global_obs=0`,
      - for filters using domain localisation: `_f` means observations within localisation radius.
      - for filters without domain localisation: `_f` means observations in each processor
        making it the same as `_p`.
  - `use_global_obs=1`, `_f` means all observations globally.
