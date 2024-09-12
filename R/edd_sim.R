edd_pars_check <- function(pars, age, model, metric, offset) {
  # pars range check
  if (pars[1] <= 0) {
    stop("speciation rate should be positive")
  }
  if (pars[2] < 0) {
    stop("extinction rate should be none-negative")
  }
  # pars and model match check
  if (model == "dsce2" && length(pars) != 6) {
    pars <- c(pars, 0.0, 0.0)   # zero out gamme_num & gamma_phi
  }
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


edd_lamu <- function(params) {
  la0 <- params[1]
  mu0 <- params[2]
  beta_num <- params[3]
  beta_phi <- params[4]
  gamma_num <- params[5]
  gamma_phi <- params[6]
  function(num, ed) {
    list(
      newlas = pmax(0, la0 + beta_num * num + beta_phi * ed),
      newmus = pmax(0, mu0 + gamma_num * num + gamma_phi * ed)
    )
  }
}


# retunrs valid SimTable or error message
edd_sim_single <- function(pars, age, model, offset, size_limit) {
  sim <- evesim::SimTable(age)
  return(sim)
}


#' @export edd_sim
edd_sim <- function(pars,
                    age,
                    model = "dsce2",
                    metric = "ed",
                    offset = "none",
                    size_limit = 1e6,
                    retry = 100) {
  pars <- edd_pars_check(pars, age, model, metric, offset)
  msg <- ""
  while (retry > 0) {
    res <- edd_sim_single(pars, age, model, offset, size_limit)
    if (is.character(res)) {
      msg <- paste(msg, sep="\n", res)
      retry <- retry - 1
    } else {
      return(list(sim = res, msg = msg))
    }
  }
  list(sim = NULL, msg = msg)
}
