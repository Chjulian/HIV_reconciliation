#!/usr/bin/env Rscript

# Generate founder strain distribution tables h(x,t)
#
# For each value of nCells (= max virions or target cells considered),
# this script simulates the distribution of the number of distinct
# genetic variants Y, given N established virions/cells and time since
# infection t.
#
# The variant type distribution h(x,t) follows a gamma density:
#   h(x,t) ~ dgamma(x, shape = 0.417, scale = t / 0.563)
# which is equivalent to h(x,t) = x^{-gamma} exp(-delta x / t) / C(t)
# with gamma = 0.583, delta = 0.563 (see manuscript eq. for h).
#
# For each (nparticles, t) pair, nSims Monte Carlo draws determine
# the probability of Y = 1, 2, ..., 9, or 10+ distinct variants.
#
# Output: one CSV per nCells value, with columns:
#   nparticles, tvals, V1, V2, ..., V9, V10plus
# Trailing all-zero variant columns are trimmed to reduce file size.

library(parallel)

# ---- Parameters ----

# Stage durations (used only to determine the time grid)
tau_E <- 0.65          # max primary stage duration (years)
tau_L <- 10 / 12       # late stage duration (years)

# Hill function for asymptomatic duration: tau_A(v) = Dmax * D50^Dk / (v^Dk + D50^Dk)
Dmax   <- 25.4
Dfifty <- 3058
Dk     <- 0.41

# VL grid — needed to find the maximum possible infection duration
naVals <- round(c(10^2, 10^2.5, 10^3, 10^3.5, 10^4, 10^4.5, 10^5, 10^5.5, 10^6, 10^6.5, 10^7))
tau_A   <- Dmax * (Dfifty^Dk) / (naVals^Dk + Dfifty^Dk)
maximumTime <- max(tau_E + tau_A + tau_L)

# Time grid: 300 midpoints spanning [0, maximumTime]
nTimeSteps    <- 300
timeWindowEdges <- seq(0, maximumTime, by = maximumTime / nTimeSteps)
timeVals      <- timeWindowEdges[-length(timeWindowEdges)] + diff(timeWindowEdges) / 2

# Simulation settings
nSims      <- 100   # Monte Carlo replicates per (nparticles, t)
max_variants <- 10  # track V1-V9 individually; V10+ collapsed into V10plus

# Output
out_dir      <- "data/tables"
nCells_range <- 2:999   # range of nCells values to generate
nWorkers     <- 4       # parallel workers

# ---- Table generation ----

create_table <- function(nCells) {
  cat(sprintf("nCells = %d\n", nCells))

  n_time_vals <- length(timeVals)

  # Pre-allocate: nCells x nTimeSteps x max_variants
  result_densities <- array(0, dim = c(nCells, n_time_vals, max_variants))

  # Pre-calculate normalised gamma probabilities for each time step
  gamma_probs_list <- lapply(timeVals, function(t) {
    probs <- dgamma(1:nCells, shape = 0.417, scale = t / 0.563)
    probs / sum(probs)
  })

  # Simulate
  for (t_idx in 1:n_time_vals) {
    gamma_probs <- gamma_probs_list[[t_idx]]

    for (n_particles in 1:nCells) {
      if (n_particles == 1) {
        result_densities[n_particles, t_idx, 1] <- 1.0
      } else {
        strain_counts <- replicate(nSims, {
          ids <- sample(1:nCells, size = n_particles,
                        prob = gamma_probs, replace = TRUE)
          min(length(unique(ids)), max_variants)
        })
        hist_result <- hist(strain_counts, breaks = 0:max_variants, plot = FALSE)
        result_densities[n_particles, t_idx, ] <- hist_result$density
      }
    }
  }

  # Build dataframe
  result_rows <- vector("list", nCells * n_time_vals)
  idx <- 1
  for (n_particles in 1:nCells) {
    for (t_idx in 1:n_time_vals) {
      result_rows[[idx]] <- c(n_particles, timeVals[t_idx],
                              result_densities[n_particles, t_idx, ])
      idx <- idx + 1
    }
  }

  result_df <- as.data.frame(do.call(rbind, result_rows))
  colnames(result_df) <- c("nparticles", "tvals",
                            paste0("V", 1:(max_variants - 1)), "V10plus")

  # Trim trailing all-zero variant columns
  v_cols    <- 3:ncol(result_df)
  col_sums  <- colSums(result_df[, v_cols])
  last_nonzero <- max(which(col_sums > 0))
  if (last_nonzero < length(v_cols)) {
    keep_through <- min(last_nonzero + 1, length(v_cols))  # one zero buffer column
    result_df <- result_df[, 1:(keep_through + 2)]         # +2 for nparticles, tvals
  }

  write.csv(result_df,
            file = file.path(out_dir, paste0(nCells, ".csv")),
            row.names = FALSE)
}

# ---- Run ----
# Only execute batch generation when this script is run directly
# (not when sourced to load create_table).

if (sys.nframe() == 0L) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  cl <- makeCluster(nWorkers)
  clusterExport(cl, c("timeVals", "nSims", "max_variants", "out_dir"))

  cat(sprintf("Generating %d tables (%d to %d) with %d workers...\n",
              length(nCells_range), min(nCells_range), max(nCells_range), nWorkers))

  system.time({
    results <- parLapply(cl, nCells_range, create_table)
  })

  stopCluster(cl)
  cat("Done.\n")
}
