
<!-- README.md is generated from README.Rmd. Please edit that file -->

# exdex

[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/paulnorthrop/exdex?branch=master&svg=true)](https://ci.appveyor.com/project/paulnorthrop/exdex)
[![Coverage
Status](https://codecov.io/github/paulnorthrop/exdex/coverage.svg?branch=master)](https://codecov.io/github/paulnorthrop/exdex?branch=master)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/exdex)](https://cran.r-project.org/package=exdex)
[![Downloads
(monthly)](https://cranlogs.r-pkg.org/badges/exdex?color=brightgreen)](https://cran.r-project.org/package=exdex)
[![Downloads
(total)](https://cranlogs.r-pkg.org/badges/grand-total/exdex?color=brightgreen)](https://cran.r-project.org/package=exdex)

## Estimation of the Extremal Index

### What does exdex do?

The extremal index *θ* is a measure of the degree of local dependence in
the extremes of a stationary process. The `exdex` package performs
frequentist inference about *θ* using two types of methodology.

One type ([Northrop, 2015](https://doi.org/10.1007/s10687-015-0221-5))
is based on a model that relates the distribution of block maxima to the
marginal distribution of the data, leading to a semiparametric maxima
estimator. Two versions of this type of estimator are provided,
following [Northrop, 2015](https://doi.org/10.1007/s10687-015-0221-5)
and [Berghaus and Bücher, 2018](https://doi.org/10.1214/17-AOS1621). A
slightly modified version of the latter is also provided. Estimates are
produced using both disjoint and sliding block maxima, the latter
providing greater precision of estimation. A graphical block size
diagnostic is provided.

The other type of methodology uses a model for the distribution of
threshold inter-exceedance times ([Ferro and Segers,
2003](https://doi.org/10.1111/1467-9868.00401)). Two versions of this
type of approach are provided: the iterated weight least squares
approach of [Süveges (2007)](https://doi.org/10.1007/s10687-007-0034-2)
and the *K*-gaps model of [Süveges and Davison
(2010)](https://doi.org/10.1214/09-AOAS292). For the *K*-gaps model the
`exdex` package allows missing values in the data, can accommodate
independent subsets of data, such as monthly or seasonal time series
from different years, and can incorporate information from censored
interexceedance times. A graphical diagnostic for the threshold level
and the runs parameter *K* is provided.

### A simple example

The following code estimates the extremal index using the semiparametric
maxima estimators, for an example dataset containing a time series of
sea surges measured at Newlyn, Cornwall, UK over the period 1971-1976.
The block size of 20 was chosen using a graphical diagnostic provided by
`choose_b()`.

``` r
library(exdex)
theta <- spm(newlyn, 20)
theta
#> 
#> Call:
#> spm(data = newlyn, b = 20)
#> 
#> Estimates of the extremal index theta:
#>           N2015   BB2018  BB2018b
#> sliding   0.2392  0.3078  0.2578 
#> disjoint  0.2350  0.3042  0.2542
summary(theta)
#> 
#> Call:
#> spm(data = newlyn, b = 20)
#> 
#>                   Estimate Std. Error Bias adj.
#> N2015, sliding      0.2392    0.01990  0.003317
#> BB2018, sliding     0.3078    0.01642  0.003026
#> BB2018b, sliding    0.2578    0.01642  0.053030
#> N2015, disjoint     0.2350    0.02222  0.003726
#> BB2018, disjoint    0.3042    0.02101  0.003571
#> BB2018b, disjoint   0.2542    0.02101  0.053570
```

Now we estimate *θ* using the *K*-gaps model. The threshold *u* and runs
parameter *K* were chosen using the graphical diagnostic provided by
`choose_uk()`.

``` r
u <- quantile(newlyn, probs = 0.60)
theta <- kgaps(newlyn, u, k = 2)
theta
#> 
#> Call:
#> kgaps(data = newlyn, u = u, k = 2)
#> 
#> Estimate of the extremal index theta:
#> [1]  0.1758
summary(theta)
#> 
#> Call:
#> kgaps(data = newlyn, u = u, k = 2)
#> 
#>       Estimate Std. Error
#> theta   0.1758   0.009211
```

### Installation

To get the current released version from CRAN:

``` r
install.packages("exdex")
```

### Vignette

See `vignette("exdex-vignette", package = "exdex")` for an overview of
the package.
