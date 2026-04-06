# swash 2.0.0

## Breaking changes (Non-backwards compatible)
- Analyses are conducted via `infpan` objects rather than `sbm` objects
- Former method `plot_regions()` for class `sbm` is replaced by generic method `plot()` for `infpan` objects
- Former methods `growth()` and `growth_initial()` for `sbm` objects are are now methods for the `infpan` class
- Former function `plot_breakpoints()` is replaced by function `breaks_growth()` and the corresponding `plot()` method

## New features 
- Importing infections panel data via `load_infections_paneldata()`, which creates and instance of the new class `infpan`
- Calculation of spread indicators for `infpan` objects such as incidence and effective reproduction number $R_t$
- Function `hawkes_growth()` for parametrization of a Hawkes process equation for infections
- Method `growth_hawkes()` for parametrization of Hawkes process equations based on `infpan` objects
- Additional NLS estimation in function `exponential_growth()`
- Option to add a constant (if values of y equal to zero occur) in functions `logistic_growth()` and `exponential_growth()`, as well as in methods `growth(infpan)` and `growth_initial(infpan)` 

## Bugfixes
- Functions `metrics()`, `binary_metrics()` and `binary_metrics_glm()` now return results invisible
- Check for equal length of input vectors in `logistic_growth()`