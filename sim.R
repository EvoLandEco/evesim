# From the example

library(evesim)
#detach("package:evesim", unload=TRUE)
#library(evesim)
library(ggplot2)
library(ape)
library(DDD)

options("width" = 200)
options("digits" = 6)
RcppParallel::setThreadOptions(numThreads = 1)

pars = c(0.5, 0.1, -0.001, -0.001, 0.0, 0.0)
#pars = c(0.5, 0.1, 0.0, 0.0, 0.0, 0.0)
age = 10
model = "dsde2"
metric = "ed"
offset = "none"

sim_res <- evesim::edd_sim(pars = pars, age = age, metric = metric, offset = offset)
phy <- evesim::SimTable.phylo(sim_res$sim, drop_extinct = TRUE)
plot(phy)
