# Profiling eve

library(evetiny)
#detach("package:evetiny", unload=TRUE)
#library(evetiny)
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
#  phy <- evetiny::edd_sim(
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
  LP <- evetiny::Ltable.prune(L, age)
  PHY <- evetiny::Ltable.phylo(L, TRUE, age)
  print(c(age, dim(LP)[1], length(PHY$tip.label)))
  X <- Ltable.tree(LP, TRUE, age)
  mbm = microbenchmark(
    ddd.L2phylo = DDDsim::L2phylo(L),
    geiger.drop.extinct = geiger::drop.extinct(PHY),
    eve.Ltable.phylo = evetiny::Ltable.phylo(LP, TRUE, age),
    eve.Xtree.phylo = evetiny::Xtree.phylo(X),
    ape.cophenetic.phylo = Phy2dist_ape(PHY),
    eve.Xtree.dist = evetiny::Xtree.dist(X),
    eve.Xtree.drop.extinct = evetiny::Xtree.drop.extinct(X),
    mean_new = evetiny::Ltable.mean(LP, TRUE, age),
    mean_tree = evetiny::Xtree.mean(X),
    nnd_new = evetiny::Ltable.nnd(LP, TRUE, age),
    nnd_tree = evetiny::Xtree.nnd(X),
    pd_new = evetiny::Ltable.pd(LP, TRUE, age),
    pd_tree = evetiny::Xtree.pd(X),
    mpd_new = evetiny::Ltable.mpd(LP, TRUE, age),
    mpd_tree = evetiny::Xtree.mpd(X),
    times = 2
  )
  mbms <- append(mbms, list(mbm))
  print(mbm)
}
