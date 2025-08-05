# This script loads the analytics-ready datasets created in script 002_
# and for each cohort:
#     estimates propensity scores using GBM (ps() in twang package)
#     estimates balancing weights via a variety of weighting methods

setwd('~/Desktop/EHR-mining/BacteremiaTreatmentEffects/')
library(dplyr)
library(twang)
library(cobalt)
library(MatchIt)
library(WeightIt)
library(ebal)
library(CBPS)
library(sbw)
source(file='~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/CreateNamedAbxAbbreviations.R') # this script can be found at: EHR-mining/UsefulDataForCleaning/antibiotic_names/...
Rdata_file_path <- '~/Desktop/EHR/EHR work/RdataFiles/BuildCohortDatasets/CleanedData/'

# Load cohort_info list.
load(file = 'CohortInfo/cohort_info.Rdata')

# Create directory for covariate balancing data.
dir.create('CohortInfo/balancing_weights_info/', showWarnings = FALSE)

# Loop over each cohort.
for (cohort in 7:8) { #seq_along(cohort_info)) {
   
   # Cohort information
   cohort_name <- cohort_info[[cohort]]$cohort_name
   cat(cohort_name, '\n')
   
   # Load clean data.
   load(file = paste0(Rdata_file_path, '/', cohort_name, '.Rdata'))
   
   # Minimal data preprocessing:
   if ('ABX_switch_time' %in% names(df))
      df <- df %>% dplyr::select(-ABX_switch_time)
   
   # Recode treatment as an integer (rather than a factor)
   col_vec <- c('ctrl' = '#3399FF', 'trt' = '#FFaa66')
   names(col_vec) <- names(sort(table(df$treatment), decreasing = TRUE))
   df$treatment <- as.integer(df$treatment == names(col_vec)[2]) # z = 0 (ctrl, more common, col_vec[1]) z = 1 (trt, less common, col_vec[2])
   
   # Remove outcome variables, recode logicals as integers (for modeling).
   data <- df %>%
      dplyr::select(-contains('OUT_')) %>%
      mutate(across(where(is.logical), as.integer)) %>%
      data.frame()
   
   # Function to calculate effective sample size
   ess <- function(w) sum(w)^2 / sum(w^2)
   
   # Treatment status vector
   z <- data$treatment
   
   # # Estimate "prognostic" score for each variable.
   # prog_data <- df %>%
   #    filter(treatment == 0L) %>%
   #    mutate(outcome = as.integer(!is.na(OUT_mortality_time) & OUT_mortality_time < 30L)) %>%
   #    dplyr::select(-contains('OUT_'), -treatment) %>%
   #    mutate(across(where(is.logical), as.integer))
   # 
   # prog_mod <- glm(formula = outcome ~ .,
   #                 family = 'binomial',
   #                 data = prog_data) %>%
   #    summary()
   # coefs <- prog_mod$coefficients[-1,]
   # 
   # # Compute covariate weights (binary (i.e., important or not) and weighted (weight by coefficient strength))
   # prog_DF <- data.frame(
   #    row.names = rownames(coefs),
   #    binary = as.integer(coefs[, 'Pr(>|z|)'] < 0.1),
   #    weights = abs(coefs[,'Estimate'])
   # )
   # 
   # # Normalize weights to the same as the binary case.
   # prog_DF$weights <- prog_DF$weights / sum(prog_DF$weights) * sum(prog_DF$binary)
   # 
   
   # Ensure control group covariate matrix is full rank.
   Treatment <- data$treatment
   X <- data %>% dplyr::select(-treatment)
   
   # Remove covariates that findLinearCombos thinks are too collinear.
   remove_vars <- caret::findLinearCombos(X)$remove
   if (length(remove_vars) > 0L)
      X <- X[-remove_vars]
   
   # Remove covariates that have > 0.9 correlation with another to avoid singularities during entropy balancing.
   r <- matrix(NA_real_, nrow=length(X), ncol=length(X), dimnames=list(names(X), names(X)))
   for (v1 in seq_along(X)[-length(X)]) {
      var1 <- X[[v1]]
      for (v2 in (v1+1):length(X)) {
         var2 <- X[[v2]]
         r[v1, v2] <- abs(cor(var1, var2))
      }
   }
   remove_vars <- rownames(which(r > 0.9, arr.ind=TRUE))
   if (length(remove_vars) > 0L)
      X <- X %>% dplyr::select(-!!remove_vars)
   
   
   # Using WeightIt package to get weights
   data_mod <- data.frame('treatment' = Treatment, X)
   
   cat('ATT ')
   w_ATT_energy <- weightit(treatment ~ .,
                            data = data_mod,
                            estimand = 'ATT',
                            method = 'energy')
   cat('ATC ')
   w_ATC_energy <- weightit(treatment ~ .,
                            data = data_mod,
                            estimand = 'ATC',
                            method = 'energy')
   cat('ATE ')
   w_ATE_energy <- weightit(treatment ~ .,
                            data = data_mod,
                            estimand = 'ATE',
                            method = 'energy')
   
   
   # Now, do entropy balancing (targets ATT) for both treatment and control, separately.
   get_eb_weights <- function(X, Treatment) {
      # Ensure "control" covariate matrix is full-rank.
      co.x <- X[Treatment == 0L, ]
      co.x <- cbind(Intercept = 1, co.x)
      qr_out <- qr(as.matrix(co.x))
      rank <- qr_out$rank
      co.x_fullrank <- co.x[, qr_out$pivot[1:rank]]
      X_fullrank <- X[, colnames(co.x_fullrank)[-1]]
      
      # Entropy balancing.
      eb.out <- ebalance(Treatment = Treatment,
                         X = scale(X_fullrank), 
                         method = 'AutoDiff')
      
      # Calculate weights.
      eb.w <- rep(NA_real_, length(Treatment))
      eb.w[Treatment == 1L] <- 1 / sum(Treatment)
      eb.w[Treatment == 0L] <- as.vector(eb.out$w)
      return(eb.w)
   }
   
   cat('Entropy balancing weights...\n')
   eb.w_list <- lapply(
      X = list(
         'att' = Treatment,
         'atc' = as.integer(!Treatment)
      ),
      FUN = function(treat_vec) get_eb_weights(X=X, Treatment=treat_vec)
   )
   # These weights did not work (all were 0, or near 0, 1 was near 1 - useless).
   if (grepl('faec', cohort_name) | cohort_name == 'S_aureus_MRSA_VAN_DAP_treatment') {
      eb.w_list$atc <- NULL
   }
   if (grepl('VRE|mirabilis', cohort_name)) {
      eb.w_list$att <- NULL
   }
   
   # Same thing for stabilizing weights.
   get_sb_weights <- function(X, Treatment, est='att') {
      # Ensure "control" covariate matrix is full-rank.
      if (est %in% c('att', 'atc')) {
         co.x <- X[Treatment == 0L, ]
         co.x <- cbind(Intercept = 1, co.x)
         qr_out <- qr(as.matrix(co.x))
         rank <- qr_out$rank
         co.x_fullrank <- co.x[, qr_out$pivot[1:rank]]
         X_fullrank <- X[, colnames(co.x_fullrank)[-1]]
      } else {
         X_fullrank <- X
      }
      
      # Stable balancing.
      sb.out <- sbw(dat = cbind(X_fullrank, 'treatment' = Treatment), 
                    ind = 'treatment',
                    bal = list(
                       bal_cov = names(X_fullrank),
                       bal_tol = 0.02
                    ),
                    par = list(par_est = est),
                    mes = FALSE)
      
      # Return weights.
      return(sb.out$dat_weights$sbw_weights)
   }
   
   
   # Stabilizing weights.
   cat('Stabilizing weights...\n')
   sb.w_list <- list()
   sb.w_list$att <- tryCatch(
      get_sb_weights(X=X, Treatment=Treatment),
      error = function(e) {
         message('sb_att')
         NULL
      }
   )
   sb.w_list$atc <- tryCatch(
      get_sb_weights(X=X, Treatment=as.integer(!Treatment)),
      error = function(e) {
         message('sb_atc')
         NULL
      }
   )
   sb.w_list$ate <- tryCatch(
      get_sb_weights(X=X, Treatment=Treatment, est='ate'),
      error = function(e) {
         message('sb_ate')
         NULL
      }
   )
   
   
   # CBPS - Covariate Balancing Propensity Score.
   cat('CBPS weights...\n')
   cbps.w_list <- lapply(
      X = c('ate' = 0,
            'att' = 1,
            'atc' = 2), # 
      FUN = function(att) {
         CBPS(formula = treatment ~ .,
              data = data,
              ATT = att)$weights
      }
   )
   
   
   # Propensity score estimation using twang - target ATE and minimize es.mean.
   cat('GBM propensity score estimation...\n')
   ps.gbm <- ps(formula = treatment ~ .,
                data = data,
                n.trees = 10000,
                keep.data = FALSE,
                stop.method = 'es.mean',
                estimand = 'ATE',
                verbose = FALSE,
                version = 'gbm')
   prop_scores <- ps.gbm$ps$es.mean.ATE
   
   
   # Using these propensity scores, construct a series of weights for estimating treatment effects:
   # ATE, ATO, ATS, ATT, ATC, matching at different ratios, EB, SB, CB, etc.
   weights_list <- list(
      'raw' = rep(1, nrow(data)),
      
      'ATT_energy' = w_ATT_energy$weights,
      'ATC_energy' = w_ATC_energy$weights,
      'ATE_energy' = w_ATE_energy$weights,
      
      'CBPS_ATE' = unname(cbps.w_list$ate),
      'CBPS_ATT' = unname(cbps.w_list$att),
      'CBPS_ATC' = unname(cbps.w_list$atc),
      'EB_ATT' = eb.w_list$att,
      'EB_ATC' = eb.w_list$atc,
      'SBW_ATT' = sb.w_list$att,
      'SBW_ATC' = sb.w_list$atc,
      'SBW_ATE' = sb.w_list$ate,
      
      'ATE' = (z / prop_scores) + ( (1 - z) / (1 - prop_scores) ),
      'ATO' = (z * (1 - prop_scores)) + ( (1 - z) * prop_scores ),
      'ATS' = ( z * ((sum(z == 1) / length(z)) / prop_scores) + ( (1 - z) * ((sum(z == 0) / length(z)) / (1 - prop_scores) )) ),
      'ATTw' = z + ( (1 - z) * prop_scores / (1 - prop_scores) ),
      'ATCw' = (1 - z) + ( z * (1 - prop_scores) / prop_scores )
   )
   weights_list <- weights_list[lengths(weights_list) > 0L]
   
   # Now, for each weighting scheme, assess and record balance across all variables.
   stats_list <- lapply(
      X = weights_list,
      FUN = function(w) {
         if (is.null(w)) return()
         vr <- col_w_vr(mat=data %>% dplyr::select(-treatment), treat=z, weights=w)
         vr[vr == 0 | is.nan(vr)] <- 1
         vr[vr < 1] <- 1 / vr[vr < 1]
         data.frame(
            'smd' = col_w_smd(mat=data %>% dplyr::select(-treatment), treat=z, weights=w),
            'ks' = col_w_ks(mat=data %>% dplyr::select(-treatment), treat=z, weights=w),
            'vr' = vr
         )
      }
   )
   
   # For each weighting scheme, get summary statistics and effective sample sizes.
   sum_balance_mat <- rbind(
      # Summary statistics.
      sapply(
         X = stats_list,
         FUN = function(x) {
            smd <- x$smd
            ks <- x$ks
            vr <- x$vr
            w_vr_inf <- which(is.infinite(vr))
            if (length(w_vr_inf) > 0L) {
               geomean_vr <- prod(vr[-w_vr_inf])^(1/length(vr[-w_vr_inf]))
            } else {
               geomean_vr <- prod(vr)^(1/length(vr))
            }
            c(
               'Mean Abs. Std. Mean Diffs.' = mean(abs(smd)),
               'Mean KS Test Statistic' = mean(ks),
               'Geom. Mean Var. Ratio' = geomean_vr
            )
         }
      ),
      # Effective sample sizes.
      sapply(
         X = weights_list,
         FUN = function(w) {
            c(
               'ess.ctrl' = ess(w[z == 0L]),
               'ess.trt' = ess(w[z == 1L]),
               'ess.total' = ess(w)
            )
         }
      )
   )
   
   # Combine all of this into one list per cohort.
   balance_list <- list(
      trt_df = data.frame(row.names=NULL, trt=c('ctrl','trt'), z = c(0,1), abx = names(col_vec), col_vec),
      treat_vec = z,
      prop_scores = prop_scores,
      summary_mat = sum_balance_mat,
      weights_list = weights_list,
      stats_list = stats_list,
      weightit_list = list(
         'ATT_energy' = w_ATT_energy,
         'ATC_energy' = w_ATC_energy,
         'ATE_energy' = w_ATE_energy
      )
   )

   # Remove patient-level data from weightit objects   
   for (i in grepv('energy', names(balance_list$weights_list))) {
      balance_list$weightit_list[[i]]$covs <- NULL
   }
   
   # Save it
   save(balance_list, file = paste0('CohortInfo/balancing_weights_info/', cohort_name, '.Rdata'))
   cat('\n')
   
   rm(balance_list, cbps.w_list, data, data_mod, df, eb.w_list, ps.gbm, r, sb.w_list, stats_list, sum_balance_mat, w_ATC_energy, 
      w_ATE_energy, w_ATT_energy, weights_list, X, alert, cohort_name, prop_scores, remove_vars, Treatment, v1, v2, var1, var2, z, col_vec)
}













