# Some notes

## [edd_sim](/R/edd_sim.R)

* Race condition in file output.
* Returns *ascending* LTable (in contrast to DDD, secsee et.al).
* Same old sample value, lookup index inefficiency (O(n) vs O(1)).
* Default `size-limit` too hight. The sim will crash before the threshold (78GB distance matrix just for tips).

```R
# edd_sim stripped by history and verbose messages
# slightly edited

edd_sim <- function(pars,
                    age,
                    model = "dsce2",
                    metric = "ed",
                    offset = "none",
                    history = TRUE,
                    verbose = FALSE,
                    converter = "cpp",
                    size_limit = 1e6) {
  edd_pars_check(pars, age, model, metric, offset)
  done <- 0
  while (done == 0) {
    i <- 1
    t <- rep(0, 1)
    l_table <- matrix(0, 2, 4)
    t[1] <- 0
    num <- 2
    l_table[1, 1:4] <- c(0, 0, -1, -1)
    l_table[2, 1:4] <- c(0, -1, 2, -1)
    linlist <- c(-1, 2)
    new_lin <- 2
    params <- c(num, pars)
    ed <- c(0, 0)
    ed_max <- edd_get_edmax(num, l_table, age, metric, offset, converter)
    lamu <- edd_update_lamu(ed, ed_max, params, model)
    t[i + 1] <- t[i] + stats::rexp(1, edd_sum_rates(lamu$newlas, lamu$newmus))
    while (t[i + 1] <= age & num[i] < size_limit) {
      i <- i + 1
      ed <- edd_get_ed(num[i - 1], l_table, t[i], metric, offset, converter)
      lamu_real <- edd_update_lamu(ed, ed_max, params, model)
      prob_real <- sum(lamu_real$newlas + lamu_real$newmus)
      prob_diff_la <- max(0, sum(lamu$newlas - lamu_real$newlas))
      prob_diff_mu <- max(0, sum(lamu$newmus - lamu_real$newmus))
      prob_fake <- sum(prob_diff_la + prob_diff_mu)
      event_type <- sample(c("real", "fake"), 1, prob = c(prob_real, prob_fake))
      if (event_type == "real") {
        event <- edd_sample_event(lamu_real$newlas, lamu_real$newmus, linlist)
        ran_lin <- c(linlist, linlist)[event]
        if (event <= length(linlist)) {
          num[i] <- num[i - 1] + 1
          new_lin <- new_lin + 1
          l_table <- rbind(l_table, c(t[i], ran_lin, sign(ran_lin) * new_lin, -1))
          linlist <- c(linlist, sign(ran_lin) * new_lin)
        } else {
          num[i] <- num[i - 1] - 1
          l_table[abs(ran_lin), 4] <- t[i]
          w <- which(linlist == ran_lin)
          linlist <- linlist[-w]
        }
      } else {
        num[i] <- num[i - 1]    # nothing to do
      }
      if (sum(linlist < 0) == 0 | sum(linlist > 0) == 0) {
        t[i + 1] <- Inf     # One of the crown lineage goes extinct, signal termination
      } else {
        ed <- edd_get_ed(num[i], l_table, t[i], metric, offset, converter)  
        ed_max <- edd_get_edmax(num[i], l_table, age, metric, offset, converter)
        params[1] <- num[i]
        lamu <- edd_update_lamu(ed, ed_max, params, model)
        if (edd_sum_rates(lamu$newlas, lamu$newmus) == 0) {
          t[i + 1] <- Inf   # lamu zero, signal termination
        } else {
          t[i + 1] <- t[i] + stats::rexp(1, edd_sum_rates(lamu$newlas, lamu$newmus))
          }
        }
      }
    }
    if (sum(linlist < 0) == 0 | sum(linlist > 0) == 0) {
      done <- 0     # both lineages extant
    } else {
      done <- 1
    }
  }
  ...
}
```

* Fail test `if (edd_sum_rates(lamu$newlas, lamu$newmus) == 0)` too generous.
* Event type is based on `edd_get_ed(num[i-1], l_table, t[i], ...)`, with `num[i-1] == current num`.<br>
next event time (Gillespie) is based on `edd_get_ed(num[i], l_table, t[i], ...)`, with `num[i] == num after event.`<br>
That is *instantaneous* after the current event. A bit wired but more importantly,<br>
this will create trees with *zero length edges*.<br>
**That has the potential to break a lot of algorithms downstream** - read: it will.<br>
First victim seems to be `L2phylo` that needed to be rewritten to handle `drop_extinct` 

## `L2phylo2`, `l_to_phylo_ed`

The old, tested conversion function `DDD::L2phylo` would not drop the following 
extinction event at t=4.5:

```R
L at t = 4.5
[1,]: num [1:4(1d)] 0 0 -1 4.5   // == t
[2,]: num [1:4(1d)] 0 -1 2 -1
[3,]: num [1:4(1d)] 4.16 -1 -3 -1
```

because `geiger::drop.extinct` et.al wouldn't.<br>
Note: in this example, `L2phylo2` creates a valid `phylo` object that
represents an invalid `Ltable` if converted back :(


```R
set.seed(42)
age <- 5
pars <- c(0.5, 0.1, -0.001, -0.001, 0.001, 0.001)
```

![](small_tree.png)

```
L:
[0] t: 5 ancestor: -1 death: 0.43719714520563535 flag: 0
[1] t: 5 ancestor: -1 death: -1 flag: 1
[2] t: 0.83669295891942053 ancestor: 0 death: -1 flag: 0
[3] t: 0.067530950604618667 ancestor: 1 death: -1 flag: 1
```

![](small_tree_ext.png)

```

{...}
[0] t: 6 ancestor: -1 death: 1.4371971452056354 flag: 0
[1] t: 6 ancestor: -1 death: -1 flag: 1
[2] t: 1.8366929589194205 ancestor: 0 death: -1 flag: 0
[3] t: 1.0675309506046187 ancestor: 1 death: -1 flag: 1
[4] t: 0.46745307649728396 ancestor: 2 death: -1 flag: 0
```

## Bug in `L2newick.h::index_of_parent`:

```C++
int index_of_parent(const std::vector< std::array< double, 4>>& ltable,
                    int parent) {
  int index = 0;
  bool found = false;
  for (; index < ltable.size(); ++index) {
   if (std::abs(ltable[index][2] - parent) < 0.0000001) {   // <-- Strange and buggy
      found = true;
      break;
    }
  }
  if (!found) index = -1;
  return index;
}
```

The function fails to find the index of 'parent == -2` -- flipped sign from DDD...<br>
Corrected:

```C++
   if (std::abs(ltable[index][2]) == std::abs(parent)) {
```

# edd_sim

* Unusual sampling:<br>
Time to next event (dt) based on lamu(t[i -1]) - current time.<br>
Event based on lamu(t[i]) - time in future?
* Why are "dsds2" and "dsce2" different models?
* Is "mpd" a good name?
* calling the pd-offset mode "simtime" is ambiguous. The "offset" is "simtime" for "ed" but "age" for "ed_max".
* Inconsistent uniform `ed_max = 2 * age` for "ed" and "nnd" compared to `ed(age)`.

## edd_update_lamu

```R
# edd_update_lamu
if (any(ed > ed_max)) {
  dummy <- 1  # breakpoint here never triggers (and it shouldn't based an namings)
}
if (beta_phi < 0) {
  newlas <- pmax(0, la0 + beta_num * num + beta_phi * ed)
} else {
  # pmin(1..., 2...) must return 1...
  newlas <- pmin(la0 + beta_num * num + beta_phi * ed, la0 + beta_num * num + beta_phi * ed_max)
  newlas <- pmax(0, newlas)
}
# same for gamma_phi
```

