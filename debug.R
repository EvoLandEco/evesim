library(evesim)
#detach("package:evesim", unload = TRUE)
#library(evesim)
library(ape)
library(microbenchmark)
library(ggplot2)


set.seed(42)
options("width" = 200)
RcppParallel::setThreadOptions(numThreads = 1)


L2ED_cpp <- function(L, t, drop_extinct) {
  phy <- eve::l_to_phylo_ed(L, t, drop_extinct)
  dist_tips <-
    ape::cophenetic.phylo(phy)
  dist_means <- rowSums(dist_tips) / (dim(dist_tips)[1] - 1)
  dist_means_sorted <-
    dist_means[gtools::mixedorder(names(dist_means))]
  return(dist_means_sorted)
}

la <- 0.5
mu <- 0.3
K <- 10000
age <- 10
pars <- c(la, mu, K)

print("Running simulation")
sim <- evesim::dd_sim(pars = pars, age = age, ddmodel = 1)
phy <- SimTable.phylo(sim, drop_extinct = TRUE)
L <- Ltable.legacy.ascending(SimTable.ltable(sim), age)

LL <- evesim::SimTable.ltable(sim)
s <- evesim::SimTable(LL)
evesim::SimTable.ed(s)

D <- SimTable.cophenetic(sim)
phy <- SimTable.phylo(sim, drop_extinct = TRUE)

X <- SimTable.tree(sim, drop_extinct = FALSE)
phy1 <- Xtree.phylo(X)

return()

print("Starting benchmark")
mbm = microbenchmark(
  cophenetic_drop = SimTable.cophenetic(sim),
  ape_cophenetic_drop = ape::cophenetic.phylo(phy),
  ED = SimTable.ed(sim),
  L2ed_cpp_drop = L2ED_cpp(L, age, TRUE),
  cophenetic = Xtree.cophenetic(X),
  ape_cophenetic = ape::cophenetic.phylo(phy1),
  L2ed_cpp = L2ED_cpp(L, age, FALSE),
  NND = SimTable.nnd(sim),
  XED = Xtree.ed(X),
  XPD = Xtree.pd(X),
  XNND = Xtree.nnd(X),
  MPD = Xtree.mpd(X),
  ENND = ape::dist.nodes(phy1),
  times = 10,
  unit = "milliseconds"
)
print(mbm)
dummy <- 0
