
<!-- README.md is generated from README.Rmd. Please edit that file -->
exdex
=====

Estimation of the Extremal Index
--------------------------------

This package is under development: it is not ready for use.

### What does exdex do?

The extremal index is a measure of the degree of local dependence in the extremes of a stationary process. The `exdex` package performs frequentist inference about using two types of methodology. One type ([Northrop, 2015](https://doi.org/10.1007/s10687-015-0221-5)) is based on a model that relates the distribution of block maxima to the marginal distribution of the data. The other type uses a model for the distribution of threshold inter-exceedance times ([Ferro and Segers, 2003](https://doi.org/10.1111/1467-9868.00401)). Two versions of the the latter type of approach are provided, following [Suveges (2007)](http://dx.doi.org/10.1007/s10687-007-0034-2) and [Suveges and Davison (2010)](http://dx.doi.org/10.1214/09-AOAS292).

### A simple example

### Installation

You can install exdex from github with:

``` r
# install.packages("devtools") 
devtools::install_github("paulnorthrop/exdex")
```

### Vignette

See `vignette("exdex-vignette", package = "exdex")` for an overview of the package.
