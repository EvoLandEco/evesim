<!-- badges: start -->
[![R-CMD-check](https://github.com/EvoLandEco/evesim/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EvoLandEco/evesim/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/evesim)](https://CRAN.R-project.org/package=evesim)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

# Our publication of the 'evesim' package
Tianjian Qin, Luis Valente, Rampal Etienne. Impact of evolutionary relatedness on species diversification and tree shape. Journal of Theoretical Biology. DOI: https://doi.org/10.1016/j.jtbi.2024.111992

Special thanks to Hanno Hildenbrandt and Thijs Janzen

## How to install
```r
require(remotes)
remotes::install_github("EvoLandEco/evesim")
```

## How to use
```r
library(evesim)
library(ape)

# Optionally set the number of threads.
# Defaults to number of logical cores.
RcppParallel::setThreadOptions(numThreads = 1)

# Simulation
# In the below example, we set speciation rate to 0.5, extinction rate to 0.1
# We also set the effect sizes of species richness and evolutionary relatedness on speciation process to -0.001 and -0.001
# In the example, the effects on extinction are set to zeros
pars = c(0.5, 0.1, -0.001, -0.001, 0.0, 0.0)

# With pre-set parameters, we now run a simulation to generate trees with crown age fixed at 10 time units
# We use the NND (nearest neighbor distance) scenario
sim <- edd_sim(pars = pars, age = 10, metric = "nnd", offset = "none")

# Simulation outputs are embedded in the 'sim' object
# Plot extant tree
plot(sim$tes)

# Plot complete tree
plot(sim$tas)

# L table (to view the evolutionary history)
sim$L
```
