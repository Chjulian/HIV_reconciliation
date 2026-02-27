# Models M4--M7 (cell-limited)
# For given per-virion transmission probabilities (early, asymptomatic, late)
# and the probability of permissive conditions for infection (f_perm),
# this function calculates i) the probability of acquisition, and
# ii) the probability of multiple genetic variants
# by stage of infection.
# This code pre-loads the genetic variant distribution tables h(x,t),
# which give the probability of Y distinct variants given N_2
# established infections and time since infection t.

populationmodel_acrossVL <-
  function(p_E = 4.715e-8,
           p_A = 4.715e-8,
           p_L = 4.715e-8,
           f_perm = 0.029,
           nCells = 100,
           r_E = 1,
           r_A = 1,
           r_L = 1,
           tau_E = 0.243) {
    require(doFuture)
    require(dplyr)
    require(fitdistrplus)
    
    ##################################################
    # Viral load parameters
    ##################################################
    nVals <- round(10 ^ seq(2, 7, length.out = 10))
    
    ##################################################
    # Set durations for each stage (main text eq 24)
    ##################################################
    # tau_A(v_A) = tau_max * v_50^tau_k / (v_A^tau_k + v_50^tau_k)
    # (Thompson et al 2019)
    alpha = -3.55
    sigma = 0.78 / (sqrt(1 - ((2 * alpha ^ 2) / (pi * (
      1 + alpha ^ 2
    )))))
    mu = 4.74 - (2 * sigma * alpha) / (sqrt(2 * pi * (1 + alpha ^ 2)))
    Dmax = 25.4
    Dfifty = 3058
    Dk = 0.41
    tau_A = array(0, length(nVals))
    tau_A = Dmax * (Dfifty ^ Dk) / (nVals ^ Dk + (Dfifty ^ Dk))
    
    # tau_L (late stage), fixed value (Bellan et al 2015)
    tau_L = 10 / 12
    
    ##################################################
    # Set weights (g_E, g_A, g_L) for each VL (main text eq 1)
    ##################################################
    
    # Early stage: normal distribution fitted to Fiebig-stage viral load quantiles (Fiebig et al 2003)
    g_E <- dnorm(log10(nVals), 4.4588631, 0.7915453)
    g_E <- g_E / sum(g_E)
    
    # Asymptomatic stage: skew-normal distribution (Thompson et al 2019)
    g_A <- (2 / sigma) * dnorm((log10(nVals) - mu) / sigma) *
      pnorm(alpha * (log10(nVals) - mu) / sigma)
    # Save seroconversion weights before prevalence rescaling
    # (used in multiplicity calculation where the 1/T_i factor
    # from the uniform infection-age distribution is not explicit)
    g_A_sero <- g_A / sum(g_A)
    # total time of infectiousness across VL types
    totalOfTime <- length(g_A) * (tau_E + tau_L) + sum(tau_A)
    # Rescale g_A from seroconversion to prevalence (main text eq 25)
    # g_A(v_A) proportional to g_A^sero(v_A) * (tau_E + tau_A(v_A) + tau_L)
    # (Fraser et al 2007)
    g_A <- g_A * ((tau_E + tau_A + tau_L) / totalOfTime)
    # normalise g_A
    g_A <- g_A / sum(g_A)
    
    # Late stage: log-normal distribution fitted to viral load data
    # from Mellors et al (1995)
    L_fit <-
      fitdist(log(c(20190, 62530, 103900, 188000, 1235000), 10), distr = "norm", method = "mle")
    g_L <-
      dnorm(log(c(nVals), 10), L_fit$estimate[[1]], L_fit$estimate[[2]])
    g_L <- g_L / sum(g_L)
    
    ##################################################
    # Probability of acquisition
    ##################################################
    
    # The number of transmitted virions ~ Poisson(vp) (main text eq 13).
    # Virions per cell ~ Poisson(vp/c) (eq 14),
    # so kappa = 1 - exp(-vp/c) (eq 15)
    # and the number of infected cells I ~ Poisson(c * kappa) (eq 16).
    pois_fulldist_function <-
      function(nVals, p) {
        return(purrr::map(.x = nVals,
                          ~ dpois(0:nCells,
                                  nCells *
                                    (1 - exp(
                                      -(.x * p / nCells)
                                    )))))
      }
    
    #check that pois_fulldists sum to one (truncate if not)
    sumcheck <- function(mydist, threshold = 1e-9) {
      for (i in 1:length(mydist)) {
        current_sum <- sum(mydist[[i]])
        missing_prob <- 1 - current_sum
        
        if (abs(missing_prob) > threshold) {
          # distribution truncation - add missing mass to last bin
          mydist[[i]][length(mydist[[i]])] <-
            mydist[[i]][length(mydist[[i]])] + missing_prob
        } else if (abs(missing_prob) > 0) {
          # Minor rounding error - normalize
          mydist[[i]] <- mydist[[i]] / current_sum
        }
        # If missing_prob is exactly 0, do nothing
      }
      return(mydist)
    }
    
    E_prob_fulldist <-
      sumcheck(pois_fulldist_function(nVals, p_E))
    A_prob_fulldist <-
      sumcheck(pois_fulldist_function(nVals, p_A))
    L_prob_fulldist <-
      sumcheck(pois_fulldist_function(nVals, p_L))
    
    ##################################################
    # Prob. dist. of N_2 (established infections)
    # after permissive exposure
    ##################################################
    
    # Permissive conditions (main text eq 20) and
    # binomial thinning N_2 | I ~ Binom(I, r_s) (eq 17)
    myfun1 <- function(distr0, r) {
      distr0 <- distr0 * f_perm
      distr0[1] <-  distr0[1] + (1 - f_perm)
      return(Reduce(`+`, purrr::map2(
        .x = distr0,
        .y = lapply(0:nCells, function(x)
          dbinom(0:nCells, x, r)),
        .f = ~ (.x * .y)
      )))
    }
    
    # Triple summation (main text eq 21)
    compute_triple_summation <-
      function(nVals,
               g_A,
               g_E,
               g_L,
               tau_A,
               A_prob_fulldist,
               E_prob_fulldist,
               L_prob_fulldist,
               tau_E,
               tau_L,
               f_perm,
               nCells,
               r_E,
               r_A,
               r_L) {
        # PRE-COMPUTATION: Calculate all myfun1 results once
        env_distr_A_all <-
          purrr::map(A_prob_fulldist, ~ myfun1(.x, r_A))
        env_distr_E_all <-
          purrr::map(E_prob_fulldist, ~ myfun1(.x, r_E))
        env_distr_L_all <-
          purrr::map(L_prob_fulldist, ~ myfun1(.x, r_L))
        
        # PRE-COMPUTATION: Calculate all tau denominators
        tau_denominators <-
          tau_E + tau_A + tau_L  # Vector of length nVals
        tau_ratios_E <- tau_E / tau_denominators
        tau_ratios_A <- tau_A / tau_denominators
        tau_ratios_L <- tau_L / tau_denominators
        
        # Initialize result vector
        total_prob <- rep(0, nCells + 1)
        
        # Triple summation: nA, nE, nL
        for (i in 1:length(nVals)) {
          # nA loop (asymptomatic)
          gA <- g_A[i]
          env_distr_A <- env_distr_A_all[[i]]
          ratio_E <- tau_ratios_E[i]
          ratio_A <- tau_ratios_A[i]
          ratio_L <- tau_ratios_L[i]
          contrib_A <- env_distr_A * ratio_A
          for (j in 1:length(nVals)) {
            # nE loop (early/primary)
            gE <- g_E[j]
            env_distr_E <- env_distr_E_all[[j]]
            contrib_E <- env_distr_E * ratio_E
            for (k in 1:length(nVals)) {
              # nL loop (late)
              gL <- g_L[k]
              env_distr_L <- env_distr_L_all[[k]]
              contrib_L <- env_distr_L * ratio_L
              combined_contrib <- contrib_E + contrib_A + contrib_L
              total_prob <-
                total_prob + gA * gE * gL * combined_contrib
            }
          }
        }
        
        return(total_prob)
      }
    
    # Main calculation using triple summation
    prob_mLineagesTransmittedPerSexAct <- compute_triple_summation(
      nVals = nVals,
      g_A = g_A,
      g_E = g_E,
      g_L = g_L,
      tau_A,
      A_prob_fulldist = A_prob_fulldist,
      E_prob_fulldist = E_prob_fulldist,
      L_prob_fulldist = L_prob_fulldist,
      tau_E = tau_E,
      tau_L = tau_L,
      f_perm = f_perm,
      nCells = nCells,
      r_E = r_E,
      r_A = r_A,
      r_L = r_L
    )
    
    probAcquisitionPerSexAct = 1 - prob_mLineagesTransmittedPerSexAct[1]
    
    # Individual stage calculations - all follow the same pattern
    E_prob_fulldist_aggregated <- Reduce(`+`,
                                         purrr::map2(.x = g_E,
                                                     .y = E_prob_fulldist,
                                                     .f = ~ .x * .y))
    E_prob_fulldist_aggregated <-
      E_prob_fulldist_aggregated / sum(E_prob_fulldist_aggregated)
    prob_mLineagesTransmittedPerSexAct_E <-
      myfun1(E_prob_fulldist_aggregated, r_E)
    
    prob_mLineagesTransmittedPerSexAct_A_perVL <-
      purrr::map(A_prob_fulldist, ~ myfun1(.x, r_A))
    prob_mLineagesTransmittedPerSexAct_A <-
      Reduce(`+`,
             purrr::map2(
               .x = g_A,
               .y = prob_mLineagesTransmittedPerSexAct_A_perVL,
               .f = ~ (.x * .y)
             ))
    
    L_prob_fulldist_aggregated <- Reduce(`+`,
                                         purrr::map2(.x = g_L,
                                                     .y = L_prob_fulldist,
                                                     .f = ~ .x * .y))
    L_prob_fulldist_aggregated <-
      L_prob_fulldist_aggregated / sum(L_prob_fulldist_aggregated)
    prob_mLineagesTransmittedPerSexAct_L <-
      myfun1(L_prob_fulldist_aggregated, r_L)
    
    probAcquisitionPerSexAct_E = 1 - prob_mLineagesTransmittedPerSexAct_E[1]
    probAcquisitionPerSexAct_A = 1 - prob_mLineagesTransmittedPerSexAct_A[1]
    probAcquisitionPerSexAct_L = 1 - prob_mLineagesTransmittedPerSexAct_L[1]
    
    ##################################################
    # Probability of multiple genetic variants
    # (main text eq 22)
    ##################################################
    
    # Genetic variant distribution h(x,t) (main text eqs 11-12)
    # (Thompson et al 2019)
    
    if (nCells == 1) {
      multiple_variant_proportion =
        multiple_variant_proportion_E =
        multiple_variant_proportion_A =
        multiple_variant_proportion_L <- 0
    } else {
      numberFounderStrainDistribution <-
        data.table::fread(paste0('data/tables/', nCells, '.csv'))
      
      # Repeat variant distribution table for each viral load value
      listofnumberFounderStrainDistribution = rep(list(numberFounderStrainDistribution),
                                                  length(nVals))
      
      # Weight variant distribution by P(N_2) at each stage and time since infection
      weight_numberFounderStrainDistribution <-
        furrr::future_pmap(
          list(
            listofnumberFounderStrainDistribution,
            as.list(tau_A),
            prob_mLineagesTransmittedPerSexAct_A_perVL
          ),
          .f = function(tbl, tau_a, asym_prob) {
            mutate(
              tbl,
              prob_nparticles =
                case_when(
                  tvals <= tau_E ~ prob_mLineagesTransmittedPerSexAct_E[nparticles + 1],
                  (tvals > tau_E) &
                    (tvals <= (tau_E + tau_a)) ~ asym_prob[nparticles + 1],
                  tvals > (tau_E + tau_a) &
                    (tvals <= (tau_E + tau_a + tau_L)) ~ prob_mLineagesTransmittedPerSexAct_L[nparticles + 1],
                  (tvals > (tau_E + tau_a + tau_L)) ~ 0
                )
            ) %>%
              mutate(across(starts_with("V"), ~ (. * prob_nparticles)))
          }
        )
      rm(listofnumberFounderStrainDistribution)
      
      myfun0 <- function(tau0, w_nfd0, output) {
        V <- (grep('V', names(w_nfd0)))
        
        if (output == 'E' | output == 'ALL') {
          s <- colSums(w_nfd0[w_nfd0$tvals <= tau_E, ..V])
          toOut <- s
        }
        if (output == 'A' | output == 'ALL') {
          s <-
            colSums(w_nfd0[(w_nfd0$tvals > tau_E) &
                             (w_nfd0$tvals <= (tau_E + tau0)), ..V])
          toOut <- s
        }
        if (output == 'L' | output == 'ALL') {
          s <-
            colSums(w_nfd0[(w_nfd0$tvals > (tau_E + tau0)) &
                             (w_nfd0$tvals <= (tau_E + tau0 + tau_L)), ..V])
          toOut <- s
        }
        if (output == 'ALL') {
          ALL <- colSums(w_nfd0[w_nfd0$tvals <= (tau_E + tau0 + tau_L), ..V])
          toOut <- ALL
        }
        return(toOut)
        
      }
      
      variant_dist <- Reduce(`+`, purrr::map2(
        .x = g_A_sero,
        .y = purrr::map2(
          .x = tau_A,
          .y = weight_numberFounderStrainDistribution,
          .f = ~ myfun0(tau0 =
                          .x, w_nfd0 = .y, 'ALL')
        ),
        .f = ~ .x * .y
      ))
      variant_dist <- variant_dist / sum(variant_dist)
      multiple_variant_proportion <- 1 - as.numeric(variant_dist[1])
      
      
      variant_dist_E <- Reduce(`+`, purrr::map2(
        .x = g_A_sero,
        .y = purrr::map2(
          .x = tau_A,
          .y = weight_numberFounderStrainDistribution,
          .f = ~
            myfun0(tau0 = .x, w_nfd0 = .y, 'E')
        ),
        .f = ~ .x * .y
      ))
      variant_dist_E <- variant_dist_E / sum(variant_dist_E)
      multiple_variant_proportion_E <-
        1 - as.numeric(variant_dist_E[1])
      
      variant_dist_A <- Reduce(`+`, purrr::map2(
        .x = g_A_sero,
        .y = purrr::map2(
          .x = tau_A,
          .y = weight_numberFounderStrainDistribution,
          .f = ~
            myfun0(tau0 = .x, w_nfd0 = .y, 'A')
        ),
        .f = ~ .x * .y
      ))
      variant_dist_A <- variant_dist_A / sum(variant_dist_A)
      multiple_variant_proportion_A <-
        1 - as.numeric(variant_dist_A[1])
      
      variant_dist_L <- Reduce(`+`, purrr::map2(
        .x = g_A_sero,
        .y = purrr::map2(
          .x = tau_A,
          .y = weight_numberFounderStrainDistribution,
          .f = ~
            myfun0(tau0 = .x, w_nfd0 = .y, 'L')
        ),
        .f = ~ .x * .y
      ))
      variant_dist_L <- variant_dist_L / sum(variant_dist_L)
      multiple_variant_proportion_L <-
        1 - as.numeric(variant_dist_L[1])
      
      
      rm(weight_numberFounderStrainDistribution)
    }
    
    # Output: probability of acquisition and probability of multiple genetic variants
    output <-
      list(
        probAcquisitionPerSexAct = probAcquisitionPerSexAct,
        probAcquisitionPerSexAct_E = probAcquisitionPerSexAct_E,
        probAcquisitionPerSexAct_A = probAcquisitionPerSexAct_A,
        probAcquisitionPerSexAct_L = probAcquisitionPerSexAct_L,
        multiple_variant_prob = multiple_variant_proportion,
        multiple_variant_prob_E = multiple_variant_proportion_E,
        multiple_variant_prob_A = multiple_variant_proportion_A,
        multiple_variant_prob_L = multiple_variant_proportion_L
      )
    return(output)
    
    
  }

populationmodel_fixedVL_A <- function(SPVL = 10 ^ 5,
                                      p_A = 4.715e-8,
                                      f_perm = 0.029,
                                      nCells = 16,
                                      #c
                                      r_A = 1,
                                      #r
                                      nExpPerYear = 106.8) {
  SPVL <- round(SPVL)
  
  # I ~ Poisson(c * kappa) (main text eqs 15-16)
  A_prob_fulldist <- dpois(0:nCells,
                           nCells * (1 - exp(-(SPVL * p_A / nCells))))
  A_prob_fulldist <-
    A_prob_fulldist / sum(A_prob_fulldist)
  
  distInfectedCellAfterPermiExp <- f_perm * A_prob_fulldist
  distInfectedCellAfterPermiExp[1] <-
    distInfectedCellAfterPermiExp[1] + (1 - f_perm)
  
  # N_2 | I ~ Binomial(I, r_A) (main text eq 17)
  
  prob_mLineagesTransmittedPerSexAct <-
    Reduce(`+`,
           purrr::map2(
             .x = distInfectedCellAfterPermiExp,
             .y = lapply(0:nCells, function(x)
               dbinom(0:nCells, x, r_A)),
             .f = ~
               (.x * .y)
           ))
  
  probAcquisitionPerSexAct = 1 - prob_mLineagesTransmittedPerSexAct[1]
  
  # A(v_A) = 1 - (1 - f(1 - Pr(N_2=0|permissive)))^epsilon (main text eq 27)
  TransmissionRatePerYear = 1 - (1 - probAcquisitionPerSexAct) ^ nExpPerYear
  
  return(TransmissionRatePerYear = TransmissionRatePerYear)
}
