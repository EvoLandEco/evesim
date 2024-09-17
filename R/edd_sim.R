edd_pars_check <- function(pars, age, metric, offset) {
  # pars range check
  if (pars[1] <= 0) {
    stop("speciation rate should be positive")
  }
  if (pars[2] < 0) {
    stop("extinction rate should be none-negative")
  }
  # pars and model match check
  if (length(pars) != 6) {
    stop("six parameters required")
  }
  # metric and offset match check
  if (metric != "pd" && offset != "none") {
    stop("only pd metric has offset methods")
  }
  if (age <= 0) {
    stop("age should be positive")
  }
  pars
}


edd_metric <- function(metric, offset) {
  if (metric == "ed") {
    fn <- SimTable.ed
  } else if (metric == "nnd") {
    fn <- SimTable.nnd
  } else if (metric == "pd") {
    if (offset == "none") {
      fn <- SimTable.pd
    } else if (offset == "simtime") {
      fn <- function(sim, t) {
        SimTable.pd(sim, t) - t
      }
    } else {
      stop("unknown offset")
    }
  } else {
    stop("unknown metric")
  }
  return(fn)
}


edd_lamu <- function(pars, metric, offset) {
  la0 <- pars[1]
  mu0 <- pars[2]
  beta_num <- pars[3]
  beta_phi <- pars[4]
  gamma_num <- pars[5]
  gamma_phi <- pars[6]
  ed_fn <- edd_metric(metric, offset)

  if ((metric == "nnd") || (metric == "ed")) {
    # metric returns vector
    fn <- function(sim, t) {
      ed <- ed_fn(sim, t)
      num <- SimTable.nspecie(sim)
      lamu <- pmax(0, c(
        la0 + beta_num * num + beta_phi * ed,
        mu0 + gamma_num * num + gamma_phi * ed
      ))
      list(
        lamu = lamu,
        prob = sum(lamu),
        num = num
      )
    }
  } else {
    # metric returns value
    fn <- function(sim, t) {
      ed <- ed_fn(sim, t)
      num <- SimTable.nspecie(sim)
      lamu <- pmax(0, c(
        la0 + beta_num * num + beta_phi * ed,
        mu0 + gamma_num * num + gamma_phi * ed
      ))
      list(
        lamu = lamu,
        prob = num * sum(lamu),
        num = num
      )
    }
  }
  return(fn)
}


edd_sampler <- function(metric) {
  if ((metric == "nnd") || (metric == "ed")) {
    # metric returns vector
    fn <- function(lamu) {
      sn <- sample.int(2 * lamu$num, 1, prob = lamu$lamu)
      list(
        specie = ((sn - 1) %% lamu$num) + 1,
        speciation = (sn <= lamu$num)
      )
    }
  } else {
    # metric returns value
    fn <- function(lamu) {
      list(
        specie = sample.int(lamu$num, 1),
        speciation = (1 == sample.int(2, 1, prob = lamu$lamu))
      )
    }
  }
  return(fn)
}


# returns valid SimTable or error message
edd_sim_single <- function(pars, age, metric, offset, size_limit) {
  sim <- SimTable(age)
  lamu_fn <- edd_lamu(pars, metric, offset)
  lamu <- lamu_fn(sim, 0.0)
  sampler <- edd_sampler(metric)
  t <- stats::rexp(1, lamu$prob)
  while (t <= age) {
    lamu_real <- lamu_fn(sim, t)
    if (lamu_real$prob == 0.0) {
      return("Simulation got stuck (before event happens)")
    }
    prob_fake <- max(0, lamu$prob - lamu_real$prob)
    if (1 == sample.int(2, 1, prob = c(lamu_real$prob, prob_fake))) {
      # real event
      sample <- sampler(lamu_real)
      if (sample$speciation) {
        SimTable.speciation(sim, sample$specie, t)
      } else {
        nc <- SimTable.extinction(sim, sample$specie, t)
        if (0 == (nc[1] * nc[2])) {
          return("crown lineage went extinct")
        }
      }
    }
    lamu <- lamu_fn(sim, t)
    if (lamu$prob == 0.0) {
      return("Simulation got stuck (after event happens)")
    }
    t <- t + stats::rexp(1, lamu$prob)
    if (SimTable.size(sim) > size_limit) {
      return("size limit exceeded")
    }
  }
  return(sim)
}


#' Simulate a phylogenetic tree using the eve model
#'
#' The `edd_sim` function simulates evolutionary relatedness dependent phylogenies based on the specified parameters, given a fixed crown age.
#' It provides functionality to retry the simulation multiple times in case of errors, with an optional limit on the size of the result.
#'
#' @param pars A numeric vector of 6 parameters for the simulation:
#' \describe{
#'   \item{\code{lambda_0}}{Intrinsic speciation rate (must be positive and larger than \code{mu_0}).}
#'   \item{\code{mu_0}}{Intrinsic extinction rate.}
#'   \item{\code{beta_N}}{Effect size of species richness on the speciation rate (can be any sign or zero).}
#'   \item{\code{beta_phi}}{Effect size of evolutionary relatedness on the speciation rate (can be any sign or zero).}
#'   \item{\code{gamma_N}}{Effect size of species richness on the extinction rate (can be any sign or zero).}
#'   \item{\code{gamma_phi}}{Effect size of evolutionary relatedness on the extinction rate (can be any sign or zero).}
#' }
#' @param age A numeric value representing the maximum crown age at which the simulation will stop.
#' @param metric A character string indicating the evolutionary relatedness measure to be used in the simulation.
#' Options are:
#' \describe{
#'   \item{\code{"pd"}}{Phylogenetic diversity as the evolutionary relatedness measure (a community-wise constraint on diversification).}
#'   \item{\code{"ed"}}{Evolutionary distinctiveness as a lineage-specific measure.}
#'   \item{\code{"nnd"}}{Phylogenetic distance to the nearest neighbor as a lineage-specific measure.}
#' }
#' @param offset A character string specifying the method for applying an offset. Currently, only `"simtime"` is available,
#' and it can only be used when \code{metric} is set to `"pd"`. Defaults to `"none"`.
#' @param size_limit An integer specifying the maximum size of the simulation result. Defaults to \code{10000}.
#' @param retry An integer specifying the number of retry attempts in case of failure. Defaults to \code{100}.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{tes}{A phylogeny with only extant lineages, or \code{NULL} if the simulation failed after all retry attempts.}
#'   \item{tas}{A phylogeny with all lineages, or \code{NULL} if the simulation failed after all retry attempts.}
#'   \item{L}{An L table recording the historical events, or \code{NULL} if the simulation failed after all retry attempts.}
#'   \item{msg}{A character string containing error messages (if any) accumulated over the retry attempts.}
#' }
#'
#' @details
#' The function simulates evolutionary relatedness dependent diversification based on a specified set of parameters (\code{pars}). The intrinsic
#' speciation rate (\code{lambda_0}) must be positive and greater than the intrinsic extinction rate (\code{mu_0}). The
#' remaining parameters, which account for the effects of species richness and evolutionary relatedness on diversification
#' rates, can be of any sign or zero. The simulation may be retried multiple times (up to the specified \code{retry} limit)
#' if failures occur during execution.
#'
#' The \code{metric} argument allows for different ways of measuring evolutionary relatedness, with three options: phylogenetic
#' diversity (\code{"pd"}), evolutionary distinctiveness (\code{"ed"}), and nearest neighbor distance (\code{"nnd"}). When using
#' the phylogenetic diversity metric, an optional offset method (\code{"simtime"}) is available. The simulation will terminate
#' once the specified \code{age} is reached or the \code{size_limit} is exceeded.
#'
#' @examples
#' # Example of simulation and plotting
#' # optionally set the number of threads.
#' # defaults to number of logical cores.
#' # RcppParallel::setThreadOptions(numThreads = 1)
#'
#' pars = c(0.5, 0.1, -0.001, -0.001, 0.0, 0.0)
#' sim <- edd_sim(pars = pars, age = 10, metric = "nnd", offset = "none")
#'
#' @references
#' Impact of Evolutionary Relatedness on Species Diversification: A New Birth-Death Model
#' Tianjian Qin, Luis Valente, Rampal S. Etienne
#' bioRxiv 2023.11.09.566365; doi: https://doi.org/10.1101/2023.11.09.566365
#'
#' @export
edd_sim <- function(pars,
                    age,
                    metric = "ed",
                    offset = "none",
                    size_limit = 10000,
                    retry = 100) {
  pars <- edd_pars_check(pars, age, metric, offset)
  age <- as.double(age)
  msg <- ""
  while (retry > 0) {
    res <- edd_sim_single(pars, age, metric, offset, size_limit)
    if (is.character(res)) {
      msg <- paste(msg, sep = "\n", res)
      retry <- retry - 1
    } else {
      # success
      tes <- SimTable.phylo(res, drop_extinct = TRUE)
      tas <- SimTable.phylo(res, drop_extinct = FALSE)
      L <- SimTable.ltable(res)
      return(list(tes = tes, tas = tas, L = L, msg = msg))
    }
  }
  list(tes = NULL, tas = NULL, L = NULL, msg = msg)
}
