<!-- badges: start -->
[![R-CMD-check](https://github.com/EvoLandEco/evesim/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EvoLandEco/evesim/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/eve)](https://CRAN.R-project.org/package=eve)
<!-- badges: end -->

# evesim

## How to install
```R
require(remotes)
remotes::install_github("EvoLandEco/evesim")
```

## How to use

```R
library(evesim)
library(ape)

# Optionally set the number of threads.
# Defaults to number of logical cores.
RcppParallel::setThreadOptions(numThreads = 1)

# Simulation
pars = c(0.5, 0.1, -0.001, -0.001, 0.0, 0.0)
sim <- edd_sim(pars = pars, age = 10, metric = "nnd", offset = "none")

# convert results to an ape::phylo object
phy <- SimTable.phylo(sim$sim, drop_extinct = TRUE)
plot(phy)

# convert results to an L table
L <- SimTable.ltable(sim$sim)
L
```
