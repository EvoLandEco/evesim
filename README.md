
<!-- README.md is generated from README.Rmd. Please edit that file -->

# eve <a href='https://github.com/EvoLandEco/eve/'><img src='man/eve-logos/eve-logos_transparent.png' align="right" height="139" /></a>

[![CRAN
status](https://www.r-pkg.org/badges/version/eve)](https://cran.r-project.org/package=eve)
[![R build
status](https://github.com/DEvoLandEco/eve/workflows/R-CMD-check/badge.svg)](https://github.com/EvoLandEco/eve/actions)
[![DOI](https://zenodo.org/badge/382378337.svg)](https://zenodo.org/badge/latestdoi/382378337)

## Overview

The package `eve` is an evolution emulator which provides pipelines to
do phylogenetic-diversity-dependent simulation, analyse outputs and
generate publication-ready plots and statistics conveniently.

`eve` supports mainly three different scenarios, the species
diversification process is regulated by ***phylogenetic diversity
(PD)*** or ***evolutionary distinctiveness (ED)*** or ***nearest
neighbor distance (NND)*** respectively.

`eve` supports parallel computing , this feature is implemented with
[furrr](https://furrr.futureverse.org/), read through its documentation
before using this feature. Note that parralel computing may use a huge
amount of memory and may have large overhead.

## Installation

You can install the developmental version from
[GitHub](https://github.com/) by running the following commands in R
console:

``` r
# install.packages("remotes")
remotes::install_github("EvoLandEco/eve")
```

## Example

``` r
library(eve)

# make a combo of parameter sets
combo <- edd_combo_maker(
  la = c(0.5, 0.8),
  mu = c(0.1, 0.2),
  beta_n = c(-0.001, 0),
  beta_phi = c(-0.001, 0.001),
  gamma_n = c(-0.001, 0.001),
  gamma_phi = c(-0.001, 0.001),
  age = 5,
  model = "dsde2",
  metric = c("ed"),
  offset = "none"
)

# have a look at the combo
combo[1:3]
#> $`1`
#>   age model metric offset                                         pars
#> 1   5 dsde2     ed   none 0.500, 0.100, -0.001, -0.001, -0.001, -0.001
#> 
#> $`2`
#>   age model metric offset                                         pars
#> 2   5 dsde2     ed   none 0.800, 0.100, -0.001, -0.001, -0.001, -0.001
#> 
#> $`3`
#>   age model metric offset                                         pars
#> 3   5 dsde2     ed   none 0.500, 0.200, -0.001, -0.001, -0.001, -0.001
```

All the unique and possible combinations of parameters will be created
automatically.

`eve` supports parallel simulation, the following example shows how to
do an 8-session parallel EDD simulation given the combo, with 3
replications for each parameter set . The result will be saved to
`/result/example`. If `name` is not specified, a folder will be created
according to time and date at the moment.

``` r
# result will be save as .RData file
# edd_go(
#   combo = combo,
#   nrep = 3,
#   name = "example",
#   strategy = future::multisession,
#   workers = 8
# )
```

`edd_go` will not return any object unless you set `name = "no_save"`.

``` r
# set name = "no_save" if you want to assign the output to a variable
# output <- edd_go(
#   combo = combo,
#   nrep = 3,
#   name = "no_save",
#   strategy = future::multisession,
#   workers = 8
# )

# examine the result
# output[[1]]
```

Sequential simulation is also possible:

``` r
# edd <- edd_go(
#     combo = combo,
#     nrep = 3,
#     name = "example2",
#     strategy = future::sequential,
#     workers = 8
#   )
```
