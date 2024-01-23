# Userguide

Here, only a few conventions and caveats are provided. For more detailed explanation and documentation of each routines, we recommend to [PDAF documentation](http://pdaf.awi.de/trac/wiki).

## Naming conventions
Since we expect Python users use:
```Python
import pyPDAF.PDAF as PDAF
```
For all PDAF subroutines using name starting with `PDAF_` or `PDAFomi_`, the `PDAF` part is removed. For example, we expect `PDAF_get_state` is used as `PDAF.get_state` and `PDAFomi_gather_obs` is used as `PDAF.omi_gather_obs`. Hence, these prefixes in the subroutines are removed in `pyPDAF`.

## Arguments conventions
The arguments may be slightly different between `PDAF` and `pyPDAF`. We recommend a check at the API reference for more details.