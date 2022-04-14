# exdex 1.2.1

## New features

* A new estimator has been implemented, based on what we will call the D-gaps model of Holesovsky, J. and Fusek, M. Estimation of the extremal index using censored distributions. Extremes 23, 197–213 (2020). doi: 10.1007/s10687-020-00374-3

## Bug fixes and minor improvements

* The value returned by `nobs.kgaps()` was incorrect in cases where there are censored K-gaps that are equal to zero.  These K-gaps should not contribute to the number of observations. This has been corrected.

* In cases where the data used in `kgaps` are split into separate sequences, the threshold exceedance probability is estimated using all the data rather than locally within each sequence.

* A `logLik` method for objects inheriting from class `"kgaps"` has been added. 

* In the (unexported, internal) function `kgaps_conf_int()` the limits of the confidence intervals for the extremal index based on the K-gaps model are constrained manually to (0, 1) to avoid problems in calculating likelihood-based confidence intervals in cases where the the log-likelihood is greater than the interval cutoff when theta = 1.

* In the documentation of the argument `k` to `kgaps()` it is noted that in practice `k` should be no smaller than 1.

* The function `kgaps()` also return standard errors based on the expected information.

* In the package manual related functions have been arranged in sections for easier reading.

* Activated 3rd edition of the `testthat` package

# exdex 1.1.1

## New features

* The functions `kgaps()`, `kgaps_imt()` and `choose_uk()` can now accept a `data` argument that
    - is a matrix of independent subsets of data, such as monthly or seasonal time series from different years,
    - contains missing values, that is, `NA`s. 
* A new dataset `cheeseboro` is included, which is a matrix containing some missing values.
* In addition to `kgaps()`, the functions `kgaps_imt()` and `choose_uk()` now have an extra argument `inc_cens`, which allows contributions from censored K-gaps to be included in the log-likelihood for the extremal index.
* The default value of `inc_cens` in `kgaps()` (and in `kgaps_imt()` and `choose_uk()`) and is now `inc_cens = TRUE`.

## Bug fixes and minor improvements

* Plot and print methods have been added for objects of class `"confint_gaps"` returned from `confint.kgaps()`.
* In `confint.spm()` and `confint.kgaps()` the input confidence `level` is included in the output object.

# exdex 1.0.1

## Bug fixes and minor improvements

* An overloading ambiguity has been corrected to ensure installation on Solaris.

