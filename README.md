# HIV_Reconciliation

This is the code repository for the study:

> Antal T, Atkins KE, Hue S, Lythgoe KA, Regoes RR, Thompson RN, Villabona-Arenas CJ. "Reconciling HIV epidemiology by accounting for within-host dynamics." [Unpublished]

## Overview

We present a suite of seven mechanistic population-level models (M1--M7) that predict HIV acquisition probability and the probability of transmitting multiple variants. This work extends the population-level transmission model of [Thompson et al. (2019)](https://doi.org/10.1093/ve/vey038) by incorporating target cell limitation and stage-dependent viral infectivity. The models combine three mechanisms:

1. **Intermittent susceptibility**: infection only occurs during permissive exposures (probability $f$, `f_perm`)
2. **Stage-dependent virus infectivity**: the per-virion transmission probability ($p_E$, $p_A$, $p_L$) or the per-infected-cell establishment probability varies by infection stage ($r_E$, $r_A$, $r_L$)
3. **Target-cell-limited bottleneck**: there is a finite number of target cells $c$ (`nCells`) at the site of infection

The seven models arise from different combinations of these mechanisms:

| Model | Intermittent susceptibility | Stage-dependent infectivity | Target-cell-limited |
|-------|:--:|:--:|:--:|
| M1 | x |   |   |
| M2 |   | x |   |
| M3 | x | x |   |
| M4 |   |   | x |
| M5 | x |   | x |
| M6 |   | x | x |
| M7 | x | x | x |

## Repository structure

```
HIV_reconciliation/
├── README.md
├── models/
│   ├── populationmodel_acrossVL.R              # M1:M3 (non-cell-limited)
│   └── populationmodel_acrossVL_cellLimited.R  # M4:M7 (cell-limited)
├── src/
│   └── generate_tables.R                       # Generates h(x,t) tables
└── data/
    └── tables/
        ├── 7.csv                               # Pre-computed variant distribution table
        └── 100.csv                             # Pre-computed variant distribution table
```

- **`models/`** -- Population-level model functions. Each file defines `populationmodel_acrossVL()`, which computes acquisition probability and variant distributions given model parameters. The specific model (M1--M7) is selected by fixing parameters at the call site (see [Usage](#usage)).
- **`src/`** -- Scripts for generating pre-computed genetic variant distribution tables (see [Generating tables](#generating-variant-diversity-tables)).
- **`data/tables/`** -- Pre-computed genetic variant distribution tables as a function of established infections and time since infection. Each CSV is named by the maximum number of established infections it supports. In the cell-limited model, `nCells` is both a biological parameter (target cells available at the site of infection) and the hard cap on established infections. In the non-cell-limited model there is no biological cap, so `maxVirions` is a computational truncation of the Poisson tail (default 1000). This repository ships tables  `7.csv` and `100.csv`; additional tables can be generated using `src/generate_tables.R`.

## Requirements

### R

R >= 4.0 is recommended. Install the required packages:

```r
install.packages(c("doFuture", "dplyr", "fitdistrplus", "purrr", "furrr", "data.table", "parallel"))
```

| Package | Purpose |
|---------|---------|
| `doFuture` | Parallel processing backend |
| `dplyr` | Data manipulation |
| `fitdistrplus` | Distribution fitting (VL distributions) |
| `purrr` | Functional programming (`map`, `map2`) |
| `furrr` | Parallel map (`future_pmap`) |
| `data.table` | Fast CSV reading (`fread`) |
| `parallel` | Multicore table generation (base R) |


## Generating variant diversity tables

The models require pre-computed tables of the genetic variant distribution $h(x, t)$, which gives the probability that a transmitted virion is of genetic variant type $x$ given time since infection $t$. This distribution follows a gamma density:

$$h(x, t) = \frac{x^{-\gamma} \, e^{-\delta x / t}}{C(t)}$$

with $\gamma = 0.583$ and $\delta = 0.563$ (see manuscript for derivation).

Each table is a CSV with columns `nparticles`, `tvals`, `V1`, ..., `V9`, `V10plus`. For a given number of established virions/cells (`nparticles`) and time since infection (`tvals`), the V columns give the probability of observing exactly 1, 2, ..., 9, or 10+ distinct genetic variants. These probabilities are estimated via Monte Carlo simulation (100 replicates per condition).

The table file is named by the `nCells` value it supports. For example, `100.csv` is loaded when running the cell-limited model with `nCells = 100`, or the non-cell-limited model when `maxVirionsConsidered` evaluates to 100.

Tables only need to be generated once for each `nCells` value. Once generated, they are reused across all subsequent model runs.

```r
source("src/generate_tables.R")
create_table(100)            # generate only table 100
lapply(2:999, create_table)  # generate a range of tables

```

Tables are written to `data/tables/`.

## Usage

All commands should be run from the repository root directory.

### Running the cell-limited model (M4--M7)

```r
source("models/populationmodel_acrossVL_cellLimited.R")

# M7: full model (all three mechanisms)
result <- populationmodel_acrossVL(
  p_E    = 4.715e-8,   # per-virion probability, early stage
  p_A    = 4.715e-8,   # per-virion probability, asymptomatic stage
  p_L    = 4.715e-8,   # per-virion probability, late stage
  f_perm = 0.029,      # probability of permissive exposure
  nCells = 100,         # number of target cells
  r_E    = 0.5,         # establishment probability, early stage
  r_A    = 0.5,         # establishment probability, asymptomatic stage
  r_L    = 0.5,         # establishment probability, late stage
  tau_E  = 0.243         # duration of early stage (years)
)
```

The specific model is selected by fixing parameters:

| Model | How to call |
|-------|-------------|
| M4 | `f_perm = 1`, `p_E = p_A = p_L = p`, `r_E = r_A = r_L = r` |
| M5 | `p_E = p_A = p_L = p`, `r_E = r_A = r_L = r` |
| M6 | `f_perm = 1` |
| M7 | All parameters free |

### Running the non-cell-limited model (M1--M3)

```r
source("models/populationmodel_acrossVL.R")

# M3: intermittent susceptibility + stage-dependent infectivity
result <- populationmodel_acrossVL(
  p_E        = 4.715e-8,
  p_A        = 2e-8,
  p_L        = 6e-8,
  f_perm     = 0.029,
  tau_E      = 0.243,
  maxVirions = 1000        # Poisson tail truncation (requires matching table)
)
```

| Model | How to call |
|-------|-------------|
| M1 | `p_E = p_A = p_L = p` |
| M2 | `f_perm = 1` |
| M3 | All parameters free |

### Model output

Both models return a list with:

| Key | Description |
|-----|-------------|
| `probAcquisitionPerSexAct` | Overall per-act acquisition probability |
| `probAcquisitionPerSexAct_E` | Per-act acquisition probability from early-stage transmitters |
| `probAcquisitionPerSexAct_A` | Per-act acquisition probability from asymptomatic-stage transmitters |
| `probAcquisitionPerSexAct_L` | Per-act acquisition probability from late-stage transmitters |
| `multiple_variant_prob` | Overall probability of multiple variants (given transmission) |
| `multiple_variant_prob_E` | Multiple variant probability from early-stage transmitters |
| `multiple_variant_prob_A` | Multiple variant probability from asymptomatic-stage transmitters |
| `multiple_variant_prob_L` | Multiple variant probability from late-stage transmitters |

## Citation

```bibtex
@article{antal2025reconciling,
  title   = {Reconciling {HIV} epidemiology by accounting for within-host dynamics},
  author  = {Antal, Tibor and Atkins, Katherine E. and Hue, Stephane and
             Lythgoe, Katrina and Regoes, Roland and Thompson, Robin N. and
             Villabona-Arenas, Ch. Julian},
  year    = {2025},
  note    = {Unpublished}
}
```

## Acknowledgements

The original model formulation is available at [robin-thompson/MultiplicityOfInfection](https://github.com/robin-thompson/MultiplicityOfInfection).


