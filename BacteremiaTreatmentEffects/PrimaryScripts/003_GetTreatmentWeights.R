# This script loads the analytics-ready datasets created in script 002_
# and for each cohort:
#     estimates propensity scores using GBM (ps() in twang package),
#     estimates balancing weights via a variety of weighting methods, AND
#     assess covariate balance effective samples across a variety of weighting schemes.

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
base_path <- '~/Desktop/EHR-mining/BuildCohortDatasets/CohortInfo/'

# Load cohort_info list.
load(file = paste0(base_path, '/cohort_info.Rdata'))

# Create a directory for propensity score info for this treatment_var set (treatment or resistance).
prop_score_path <- paste0(base_path, '/propensity_score_info/')
if (!dir.exists(prop_score_path))
   dir.create(prop_score_path)

# Loop over each cohort.
for (cohort in seq_along(cohort_info)) {
   
   # Cohort information
   cohort_name <- cohort_info[[cohort]]$cohort_name
   cat(cohort_name, '\n')
   
   # Load clean data.
   load(file = paste0(Rdata_file_path, '/', cohort_name, '.Rdata'))
   df <- df %>% dplyr::select(-AMPICILLIN)
   
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
   # # These weights did not work (all were 0, or near 0, 1 was near 1 - useless).
   # if (cohort_name %in% c('S_aureus_MRSA_VAN_DAP_treatment', 'E_faecalis_VAN_AMP_treatment', 'E_faecium_VRE_DAP_LZD_treatment')) {
   #    eb.w_list$atc <- NULL
   # }
   
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
   
   # cat('Stabilizing weights...\n')
   # sb.w_list <- list()
   # sb.w_list$att <- sb.w_list$atc <- sb.w_list$ate <- NULL
   # # Sbw could not find optimal weights for some estimands for some cohorts.
   # if (!cohort_name %in% c('P_mirabilis_CRO_TZP_treatment', 'E_faecalis_VAN_AMP_treatment', 'E_faecium_VRE_DAP_LZD_treatment')) {
   #    sb.w_list$att <- get_sb_weights(X=X, Treatment=Treatment)
   # }
   # if (!cohort_name %in% c('K_pneumoniae_CRO_TZP_treatment', 'P_mirabilis_CRO_TZP_treatment', 'S_aureus_MRSA_switch_VAN_DAP_treatment',
   #                         'S_aureus_MRSA_VAN_DAP_treatment', 'E_faecalis_VAN_AMP_treatment', 'E_faecium_VRE_DAP_LZD_treatment')) {
   #    sb.w_list$atc <- get_sb_weights(X=X, Treatment=as.integer(!Treatment))
   # }
   # if (!cohort_name %in% c('K_pneumoniae_CRO_TZP_treatment', 'P_mirabilis_CRO_TZP_treatment', 'P_aeruginosa_FEP_TZP_treatment', 
   #                         'S_aureus_MRSA_switch_VAN_DAP_treatment', 'S_aureus_OXAorNFC_CFZtreatment', 'S_aureus_MRSA_VAN_DAP_treatment',
   #                         'E_faecalis_VAN_AMP_treatment', 'ampCprod_FEP_TZP_treatment', 'E_faecalis_switch_VAN_AMP_treatment',
   #                         'E_faecium_VRE_DAP_LZD_treatment')) {
   #    sb.w_list$ate <- get_sb_weights(X=X, Treatment=Treatment, est='ate')
   # }
   
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
   save(balance_list, file = paste0(prop_score_path, cohort_name, '.Rdata'))
   cat('\n')
   
   rm(balance_list, cbps.w_list, data, data_mod, df, eb.w_list, ps.gbm, r, sb.w_list, stats_list, sum_balance_mat, w_ATC_energy, 
      w_ATE_energy, w_ATT_energy, weights_list, X, alert, cohort_name, prop_scores, remove_vars, Treatment, v1, v2, var1, var2, z, col_vec)
}






# Mean balance and effective sample sizes plots
source(file='~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/CreateNamedAbxAbbreviations.R') 
base_path <- '~/Desktop/EHR-mining/BuildCohortDatasets/CohortInfo/treatment/'
load(file = paste0(base_path, '/cohort_info.Rdata'))
prop_score_path <- paste0(base_path, '/propensity_score_info/')
plot_path <- '~/Desktop/EHR-mining/BuildCohortDatasets/Plots/BalanceDiagnostics/summary_balance_ess_ps/'
if (!dir.exists(plot_path))
   dir.create(plot_path)

# Loop over cohorts
for (cohort in seq_along(cohort_info)) {
   
   # Cohort name.
   cohort_name <- cohort_info[[cohort]]$cohort_name
   clean_bug_name <- cohort_info[[cohort]]$clean_bug_name
   cat(clean_bug_name, '\n')
   
   # Load balance_list.
   load(file = paste0(prop_score_path, cohort_name, '.Rdata'))
   bal_mat <- balance_list$summary_mat
   bal_mat <- bal_mat[,grep('(AT(T|C)|COMB)m', colnames(bal_mat), invert=TRUE)]
   
   # cat(names(which(bal_mat[1,] < 0.01)), '\n\n')
   
   if (TRUE) {
      # Setup plot space.
      pdf(file = paste0(plot_path, cohort_name, '.pdf'), width=10, height=5.5)
      par(mfrow=c(2,3), tck=-0.015, mgp=c(1.5, 0.4, 0), mar=c(3, 6, 2, 2), cex.main=1)
      col_vec <- setNames(balance_list$trt_df$col_vec, 
                          balance_list$trt_df$abx)
      
      # Effective sample sizes.
      b <- barplot(bal_mat[c('ess.ctrl', 'ess.trt'),], beside=TRUE, horiz=TRUE, names.arg=rep(NA,ncol(bal_mat)), main='Effective sample sizes', 
                   xlab='n', col=col_vec, xlim=c(0, max(bal_mat[c('ess.ctrl', 'ess.trt'),]) * 1.02))
      axis(side=2, at=apply(b,2,mean), labels=colnames(bal_mat), las=1, tick=FALSE)
      abline(v=bal_mat[c('ess.ctrl', 'ess.trt'),'raw'], lty=2, xpd=FALSE)
      
      # Bug/cohort name.
      plot.new()
      title(main=clean_bug_name, line=-4, font.main=4, cex.main=1.5, xpd=NA)
      legend('center', legend=names(col_vec), pt.bg=col_vec, pch=22, pt.cex=3, bty='n', cex=1.5)
      
      # Propensity score distributions (pre- and post-weighting).
      plotPSdist <- function(e, z) {
         brks <- seq(0,1,0.02)
         h0 <- hist(e[z == 0L], breaks=brks, plot=FALSE)
         h1 <- hist(e[z == 1L],  breaks=brks, plot=FALSE)
         hist(e[z == 0L], col=paste0(col_vec[1], 'aa'), breaks=brks, xlim=c(0,1), ylim=c(0, max(c(h1$counts, h0$counts))), 
              ylab='', xlab='', xpd=NA, yaxt='n', main='Estimator: GBM')
         hist(e[z == 1L], col=paste0(col_vec[2], 'aa'), breaks=brks, add=TRUE)
         title(xlab='Propensity score', line=1.2)
         title(ylab='Frequency', line=2)
         axis(side=2, las=1)
      }
      plotPSdist(e=balance_list$prop_scores, z=balance_list$treat_vec)
      
      # Balance statistics.
      # xmaxs <- pmin(apply(bal_mat[1:3,], 1, max) * 1.1, c(0.5,0.5,2))
      xmaxs <- bal_mat[1:3,'raw'] * 1.1
      for (i in seq_len(nrow(bal_mat[1:3,]))) {
         v <- bal_mat[i,]
         b <- barplot(v, horiz=TRUE, names.arg=NA, main=rownames(bal_mat)[i], xlim=c(0, xmaxs[i]))
         axis(side=2, at=b, labels=names(v), las=1, tick=FALSE)
         abline(v=v['raw'], lty=2, xpd=FALSE)
         if (any(bal_mat[i,] > xmaxs[i])) {
            text(x=bal_mat[i,'raw']*1.025, y=b[bal_mat[i,] > xmaxs[i]], adj=0, xpd=NA, labels=round(bal_mat[i, bal_mat[i,] > xmaxs[i]], 2))
         }
      }
      
      dev.off()
   }
}

extra_weights <- list(
   'E. faecium (VRE)' = character(),
   'E. faecalis' = character(),
   'E. faecalis switch' = c('EB_ATT', 'EB_ATC'),
   'S. aureus (MRSA)' = c('EB_ATT', 'CBPS_ATT', 'SBW_ATT'),
   'S. aureus (MRSA) switch' = c('EB_ATT', 'CBPS_ATT', 'SBW_ATT'),
   'S. aureus' = c('EB_ATT', 'EB_ATC'),
   'ampC Producers' = c('EB_ATC', 'EB_ATT'),
   'P. aeruginosa' = c('EB_ATT', 'EB_ATC'),
   'P. mirabilis' = character(),
   'K. pneumoniae' = c('EB_ATT'),
   'E. coli' = c('EB_ATT', 'EB_ATC')
)
saveRDS(extra_weights, file='~/Desktop/EHR-mining/BuildCohortDatasets/CohortInfo/treatment/propensity_score_info/which_weighting_schemes_per_cohort.rds')















