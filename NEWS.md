# exdex 1.1.1

## New features

* The functions `kgaps()`, `kgaps_imt()` and `choose_uk()` can now accept a `data` argument that
    - is a matrix of independent subsets of data, such as monthly or seasonal time series from different years
    - contains missing values, that is, `NA`s 
* A new dataset `cheeseboro` is included, which is a matrix containing some missing values.
* In addition to `kgaps()`, the functions `kgaps_imt()` and `choose_uk()` now have an extra argument `inc_cens`, which allows contributions from censored K-gaps to be included in the log-likelihood for the extremal index.
* The default value of `inc_cens` in `kgaps()` (and in `kgaps_imt()` and `choose_uk()`) and is now `inc_cens = TRUE`.

## Bug fixes and minor improvements

* Plot and print methods have been added for objects of class `"confint_gaps"` returned from `confint.kgaps()`.
* In `confint.spm()` and `confint.kgaps()` the input confidence `level` is included in the output object.

# exdex 1.0.1

## Bug fixes and minor improvements

* An overloading ambiguity has been corrected to ensure installation on Solaris.

