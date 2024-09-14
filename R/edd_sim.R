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
    fn <- evesim::SimTable.ed
  } else if (metric == "nnd") {
    fn <- evesim::SimTable.nnd
  } else if (metric == "pd") {
    if (offset == "none") {
      fn <- evesim::SimTable.pd
    } else if (offset == "simtime") {
      fn <- function(sim, t) {
        evesim::SimTable.pd(sim, t) - t
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
      num <- evesim::SimTable.nspecie(sim)
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
      num <- evesim::SimTable.nspecie(sim)
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
  sim <- evesim::SimTable(age)
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
        evesim::SimTable.speciation(sim, sample$specie, t)
      } else {
        nc <- evesim::SimTable.extinction(sim, sample$specie, t)
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
    if (evesim::SimTable.size(sim) > size_limit) {
      return("size limit exceeded")
    }
  }
  return(sim)
}


#' @export edd_sim
edd_sim <- function(pars,
                    age,
                    metric = "ed",
                    offset = "none",
                    size_limit = 10000,
                    retry = 100) {
  pars <- edd_pars_check(pars, age, metric, offset)
  msg <- ""
  while (retry > 0) {
    res <- edd_sim_single(pars, age, metric, offset, size_limit)
    if (is.character(res)) {
      msg <- paste(msg, sep = "\n", res)
      retry <- retry - 1
    } else {
      # success
      return(list(sim = res, msg = msg))
    }
  }
  list(sim = NULL, msg = msg)
}
