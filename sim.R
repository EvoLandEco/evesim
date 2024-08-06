# From the example

library(evesim)
#detach("package:evesim", unload=TRUE)
#library(evesim)
library(ggplot2)
library(ape)
library(DDD)

set.seed(42)
options("width" = 200)
options("digits" = 6)
RcppParallel::setThreadOptions(numThreads = 1)

pars = c(0.5, 0.1, -0.001, -0.001, 0.0, 0.0)
#pars = c(0.5, 0.1, 0.0, 0.0, 0.0, 0.0)
age = 10
model = "dsde2"
metric = "pd"
offset = "none"

set.seed(42)
#X <- eve::edd_sim(pars = pars, age = age, model = model, metric = metric, offset = offset, history = FALSE)
set.seed(42)
Y <- evesim::edd.sim(pars = pars, age = age, metric = metric, offset = offset)

dummy <- 0
