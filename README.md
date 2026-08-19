
Replication material for *Protected Areas reduce but do not prevent cattle encroachment in the Brazilian Amazon* by Leonie Hodel, Marin Skidmore, Joseph App, Rodrigo Balbueno, Amintas Brandão, Cainã Couto-Silva, Julia Woods, Lisa Rausch, Holly K. Gibbs.
[under review]

---

## Overview

The analysis estimates the effect of protected area (PA) status on pasture cover,
deforestation and cattle stocking rates using a spatial regression discontinuity
design at PA boundaries.

Main results:

- **Figure 2** — pooled RDD plots for all four outcomes
- **Figure 3** — RDD plots by pasture-development cohort
- **Figure 4** — RDD plots Indigenous Lands vs. Strictly Protected Areas
- **Table S3** — pooled RDD estimates
- **Tables S4–S10** — robustness and subgroup estimates

---

## Requirements

R 4.4.2 or later. The RDD script requires:

```r
install.packages(c("dplyr", "tidyr", "purrr", "tibble", "readr", "ggplot2",
                   "rdrobust", "gt", "arrow","MatchIt","patchwork","stringr"))
```

---

## Reproducing the results

Download the data from Zenodo [link will be added soon] and place into the data/ folder.
Then, run the R scripts from the repository root:

```r
source("rdd_protectedareas_figures.R")
source("rdd_protectedareas_tables.R")

```
Output is written to `results/`.

---

## Method notes

Estimates come from `rdrobust` with a local polynomial of degree 1 (linear) and
0 (constant), MSE-optimal bandwidth (`bwselect = "mserd"`), and standard errors
clustered at the protected area. Tables report both the MSE-optimal bandwidth
(`h`) and twice that width (`2h`).

---

## Licence

All code in this repository is released under the MIT Licence (see `LICENSE`).

