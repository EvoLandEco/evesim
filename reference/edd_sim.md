# Simulate a phylogenetic tree using the eve model

The `edd_sim` function simulates evolutionary relatedness dependent
phylogenies based on the specified parameters, given a fixed crown age.
It provides functionality to retry the simulation multiple times in case
of errors, with an optional limit on the size of the result.

## Usage

``` r
edd_sim(
  pars,
  age,
  metric = "ed",
  offset = "none",
  size_limit = 10000,
  retry = 100
)
```

## Arguments

- pars:

  A numeric vector of 6 parameters for the simulation:

  `lambda_0`

  :   Intrinsic speciation rate (must be positive and larger than
      `mu_0`).

  `mu_0`

  :   Intrinsic extinction rate.

  `beta_N`

  :   Effect size of species richness on the speciation rate (can be any
      sign or zero).

  `beta_phi`

  :   Effect size of evolutionary relatedness on the speciation rate
      (can be any sign or zero).

  `gamma_N`

  :   Effect size of species richness on the extinction rate (can be any
      sign or zero).

  `gamma_phi`

  :   Effect size of evolutionary relatedness on the extinction rate
      (can be any sign or zero).

- age:

  A numeric value representing the maximum crown age at which the
  simulation will stop.

- metric:

  A character string indicating the evolutionary relatedness measure to
  be used in the simulation. Options are:

  `"pd"`

  :   Phylogenetic diversity as the evolutionary relatedness measure (a
      community-wise constraint on diversification).

  `"ed"`

  :   Evolutionary distinctiveness as a lineage-specific measure.

  `"nnd"`

  :   Phylogenetic distance to the nearest neighbor as a
      lineage-specific measure.

- offset:

  A character string specifying the method for applying an offset.
  Currently, only `"simtime"` is available, and it can only be used when
  `metric` is set to `"pd"`. Defaults to `"none"`.

- size_limit:

  An integer specifying the maximum size of the simulation result.
  Defaults to `10000`.

- retry:

  An integer specifying the number of retry attempts in case of failure.
  Defaults to `100`.

## Value

A list containing the following components:

- tes:

  A phylogeny with only extant lineages, or `NULL` if the simulation
  failed after all retry attempts.

- tas:

  A phylogeny with all lineages, or `NULL` if the simulation failed
  after all retry attempts.

- L:

  An L table recording the historical events, or `NULL` if the
  simulation failed after all retry attempts.

- msg:

  A character string containing error messages (if any) accumulated over
  the retry attempts.

## Details

The function simulates evolutionary relatedness dependent
diversification based on a specified set of parameters (`pars`). The
intrinsic speciation rate (`lambda_0`) must be positive and greater than
the intrinsic extinction rate (`mu_0`). The remaining parameters, which
account for the effects of species richness and evolutionary relatedness
on diversification rates, can be of any sign or zero. The simulation may
be retried multiple times (up to the specified `retry` limit) if
failures occur during execution.

The `metric` argument allows for different ways of measuring
evolutionary relatedness, with three options: phylogenetic diversity
(`"pd"`), evolutionary distinctiveness (`"ed"`), and nearest neighbor
distance (`"nnd"`). When using the phylogenetic diversity metric, an
optional offset method (`"simtime"`) is available. The simulation will
terminate once the specified `age` is reached or the `size_limit` is
exceeded.

## References

Impact of Evolutionary Relatedness on Species Diversification: A New
Birth-Death Model Tianjian Qin, Luis Valente, Rampal S. Etienne Journal
of Theoretical Biology; DOI: https://doi.org/10.1016/j.jtbi.2024.111992

## Examples

``` r
# Example of simulation and plotting
# optionally set the number of threads.
# defaults to number of logical cores.
# RcppParallel::setThreadOptions(numThreads = 1)

pars = c(0.5, 0.1, -0.001, -0.001, 0.0, 0.0)
sim <- edd_sim(pars = pars, age = 10, metric = "nnd", offset = "none")
```
