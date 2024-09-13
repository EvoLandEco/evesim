# From the example

library(evesim)
#detach("package:evesim", unload=TRUE)
#library(evesim)
library(ggplot2)
library(ape)

RcppParallel::setThreadOptions(numThreads = 1)

pars = c(0.5, 0.1, -0.001, -0.001, 0.0, 0.0)
sim <- edd_sim(pars = pars, age = 50, metric = "nnd", offset = "none")
phy <- SimTable.phylo(sim$sim, drop_extinct = TRUE)
plot(phy)
