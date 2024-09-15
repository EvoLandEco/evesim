<!-- badges: start -->
[![R-CMD-check](https://github.com/HHildenbrandt/evesim/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/HHildenbrandt/evesim/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# evesim

## Example

```R
library(evesim)
library(ggplot2)
library(ape)

# optionally set the number of threads.
# defaults to number of logical cores.
RcppParallel::setThreadOptions(numThreads = 1)

pars = c(0.5, 0.1, -0.001, -0.001, 0.0, 0.0)
sim <- edd_sim(pars = pars, age = 50, metric = "nnd", offset = "none")

# convert to an ape::phylo object
phy <- SimTable.phylo(sim$sim, drop_extinct = TRUE)
plot(phy)

# extract some info
age <- SimTable.age(sim$sim)
nspecie <- SimTable.nspecie(sim$sim)                # tips
nclade_specie <- SimTable.nclade_specie(sim$sim)    # tips per crown lineage
size <- SimTable.size(sim$sim)                      # number of events in simulation

# convert to (full) ascending ltable
ltable <- SimTable.ltable(sim$sim)
prunned_ltable <- Ltable.prune(ltable, age = 10)
prunned_phy <- Ltable.phylo(prunned_ltable, drop_extinct = TRUE, age = 10)
plot(prunned_phy)

# tree metrices are also available for Ltables
ed <- Ltable.ed(prunned_ltable, drop_extince = TRUE, age = 10)

# and more...
```
