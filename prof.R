# Profiling eve

library(evesim)
#detach("package:evesim", unload=TRUE)
#library(evesim)
library(microbenchmark)
library(ggplot2)
library(DDDsim)

set.seed(42)
options("width"=200)

#pars = c(0.5, 0.1, -0.001, -0.001, 0.001, 0.001)
#age <- 100

#if (file.exists("./Ltable_prof.Rdata")) {
#  print("loading Ltable_prof.Rdata")
#  load("./Ltable_prof.Rdata")
#} else {
#  print("creating Ltable_prof.R")
#  phy <- evesim::edd_sim(
#    pars = pars,
#    age = age,
#    model = "dsde2",
#    metric = "ed",
#    offset = "none",
#    history = FALSE,
#    verbose = FALSE,
#    converter = "cpp")
#  L <- phy$l_table
#  save(L, file = "./Ltable_prof.Rdata")
#}

print("start profiling")

Phy2dist_ape <- function(phy) {
  d <- ape::cophenetic.phylo(phy)
  d[order(rownames(d)), order(colnames(d))]
}

RcppParallel::setThreadOptions(numThreads = 16)


L <- DDDsim::dd_sim_cpp(pars = c(0.5, 0.3, 10000), age = 100, ddmodel = 1.0)$L

mbms <- list()
for (age in c(100)) { #2,5,10,20,50,100)) {
  LP <- evesim::Ltable.prune(L, age)
  PHY <- evesim::Ltable.phylo(L, TRUE, age)
  print(c(age, dim(LP)[1], length(PHY$tip.label)))
  X <- Ltable.tree(LP, TRUE, age)
  mbm = microbenchmark(
    ddd.L2phylo = DDDsim::L2phylo(L),
    geiger.drop.extinct = geiger::drop.extinct(PHY),
    eve.Ltable.phylo = evesim::Ltable.phylo(LP, TRUE, age),
    eve.Xtree.phylo = evesim::Xtree.phylo(X),
    ape.cophenetic.phylo = Phy2dist_ape(PHY),
    eve.Xtree.dist = evesim::Xtree.dist(X),
    eve.Xtree.drop.extinct = evesim::Xtree.drop.extinct(X),
    mean_new = evesim::Ltable.mean(LP, TRUE, age),
    mean_tree = evesim::Xtree.mean(X),
    nnd_new = evesim::Ltable.nnd(LP, TRUE, age),
    nnd_tree = evesim::Xtree.nnd(X),
    pd_new = evesim::Ltable.pd(LP, TRUE, age),
    pd_tree = evesim::Xtree.pd(X),
    mpd_new = evesim::Ltable.mpd(LP, TRUE, age),
    mpd_tree = evesim::Xtree.mpd(X),
    times = 2
  )
  mbms <- append(mbms, list(mbm))
  print(mbm)
}
