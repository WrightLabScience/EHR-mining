# This script makes estimates of average treatment effects (ATE, ATT, ATC, ATO, ATS, etc.)
# For each bacteremia cohort
# Using propensity scores and weights calculated in script 003

setwd('~/Desktop/EHR-mining/BuildCohortDatasets/')
library(dplyr)
library(marginaleffects)
library(WeightIt)
source(file='~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/CreateNamedAbxAbbreviations.R')
Rdata_file_path <- '~/Desktop/EHR/EHR work/RdataFiles/BuildCohortDatasets/CleanedData/'

# Load cohort_info list.
load(file = paste0('CohortInfo/cohort_info.Rdata'))

# Create list object that will contain ATE estimates for each cohort.
ATE_estimates <- vector('list', length(cohort_info))

# Loop over each cohort (pathogen(s)/antibiotic(s)) and estimate average treatment effects.
for (cohort in seq_along(cohort_info)) {
   
   # Cohort information
   cohort_name <- cohort_info[[cohort]]$cohort_name
   clean_bug_name <- cohort_info[[cohort]]$clean_bug_name
   cat(cohort_name, '\n')
   
   # Use that clean cohort name for the ATE_estimates list object.
   names(ATE_estimates)[cohort] <- clean_bug_name
   
   # Load clean dataset.
   load(file = paste0(Rdata_file_path, '/', cohort_name, '.Rdata'))
   
   # Load previously estimated propensity scores and weights - balance_list.
   load(file = paste0('CohortInfo/propensity_score_info/', cohort_name, '.Rdata'))
   
   # Named color vector.
   col_vec <- setNames(balance_list$trt_df$col_vec, 
                       balance_list$trt_df$abx)
   
   # Create empty lists to fill with ATE estimates.
   outcome_vars <- c('mort', 'comp', 'LOS')
   ATE_estimates[[cohort]] <- setNames(vector('list', 3L), outcome_vars)
   
   # Next, get ATE estimates for each outcome using weighting from each ps method.
   threshold <- 30L # days
   
   for (o in c(1,3)) {# seq_along(outcome_vars)) {
      cat(outcome_vars[o], '--')
      
      # Data preparation - encode outcomes correctly.
      if (outcome_vars[o] == 'mort') {
         data <- df %>% 
            mutate(outcome = as.integer(!(is.na(OUT_mortality_time) | OUT_mortality_time > threshold)))
         
      } else if (outcome_vars[o] == 'comp') {
         data <- df %>%
            mutate(
               across(.cols = c(OUT_recur_time, OUT_mortality_time),
                      .fns = ~ ifelse(test = !is.na(.) & . <= threshold, 
                                      yes = ., 
                                      no = threshold + 1L)),
               outcome = OUT_mortality_time <= threshold | OUT_recur_time <= threshold
            )
         
      } else if (outcome_vars[o] == 'LOS') {
         data <- df %>%
            mutate(
               OUT_los_after_order = case_when(
                  OUT_los_after_order > threshold | (!is.na(OUT_mortality_time) & OUT_mortality_time < threshold) ~ threshold,
                  .default = OUT_los_after_order
               ),
               outcome = OUT_los_after_order > median(OUT_los_after_order)
            )
      }
      
      # Re-encode treatment as integer.
      data <- data %>% 
         mutate(treatment = as.integer(treatment == names(col_vec)[2])) %>%
         dplyr::select(-contains('OUT_'))
      
      # If the treatment involved switching, remove time to switch variable.
      if ('switch' %in% names(cohort_info[[cohort]]))
         data <- data %>% dplyr::select(-ABX_switch_time)
      
      
      # Function to parse the output of avg_comparisons
      getRiskDiff <- function(fit, weights=NULL) {
         if (is.null(fit)) return(NULL)
         if (is.null(weights)) {
            avg_comp <- avg_comparisons(fit, 
                                        variables = 'treatment',
                                        comparison = 'difference')
         } else {
            avg_comp <- avg_comparisons(fit, 
                                        variables = 'treatment',
                                        comparison = 'difference',
                                        wts = weights)
         }
         est <- avg_comp$estimate * 100
         marg_err <- qnorm(0.975) * avg_comp$std.error * 100
         return(c(
            'estimate' = est,
            'lowerCI' = est - marg_err,
            'upperCI' = est + marg_err,
            'pvalue' = avg_comp$p.value
         ))
      }
      
      # Create empty ATE estimates matrix.
      row_names <- c('unw uni', paste0('wtd ', rep(c('uni', 'mlt'), times=3), ' ', rep(c('ATE', 'ATT', 'ATC'), each=2L), ' - energy'))
      col_names <- c('estimate', 'lowerCI', 'upperCI', 'pvalue')
      coefs <- matrix(NA_real_, nrow=length(row_names), ncol=length(col_names), dimnames=list(row_names, col_names))
      
      # RAW - Unweight univariate.
      coefs['unw uni', ] <- getRiskDiff(fit = glm(formula = outcome ~ treatment,
                                                  data = data,
                                                  family = binomial()))
      
      for (estimand in c('ATE', 'ATT', 'ATC')) {
         weight_scheme <- paste0(estimand, '_energy')
         weightit_object <- balance_list$weightit_list[[weight_scheme]]
         
         # Weighted univariate.
         formula <- paste0('outcome ~ treatment')
         coefs[paste0('wtd uni ', estimand, ' - energy'), ] <- getRiskDiff(fit = glm_weightit(formula = formula,
                                                                                 data = data,
                                                                                 family = binomial(),
                                                                                 weightit = weightit_object))
         
         # ??? skip these?
         # if ( 
         #       (estimand == 'ATC' & cohort_name == 'E_faecalis_VAN_AMP_treatment') |
         #       (estimand == 'ATT' & cohort_name == 'E_faecium_VRE_DAP_LZD_treatment')
         #    ) {
         #    coefs <- coefs[rownames(coefs) != paste0('wtd mlt ', estimand, ' - energy'),]
         #    next
         # }
         
         # Weighted multivariate.
         add_covs <- rownames(balance_list$stats_list[[weight_scheme]])[abs(balance_list$stats_list[[weight_scheme]]$smd) > 0.05]
         if (length(add_covs) > 0L) {
            remove_vars <- caret::findLinearCombos(data %>% dplyr::select(!!add_covs))$remove
            if (length(remove_vars) > 0L)
               add_covs <- add_covs[-remove_vars]
            formula <- paste0(formula, ' + ', paste(add_covs, collapse=' + '))
         }
         coefs[paste0('wtd mlt ', estimand, ' - energy'), ] <- getRiskDiff(fit = glm_weightit(formula = formula,
                                                                                 data = data,
                                                                                 family = binomial(),
                                                                                 weightit = weightit_object))
      }
      
      
      # Run weighted average treatment estimation with several other weighting schemes.
      weight_schemes <- grepv('energy|raw', names(balance_list$weights_list), invert = TRUE)
      # weight_schemes <- extra_weights[[clean_bug_name]])
      for (weight_scheme in weight_schemes) {
         cat(weight_scheme, '')
         data$weights <- balance_list$weights_list[[weight_scheme]]
         
         formula <- paste0('outcome ~ treatment')
         new_row <- getRiskDiff(fit = survey::svyglm(formula = formula,
                                                     design = survey::svydesign(ids=~0, weights=~weights, data=data),
                                                     data = data,
                                                     family = quasibinomial()), 
                                weights = data$weights) %>% matrix() %>% t()
         rownames(new_row) <- paste0('wtd uni ', gsub('([A-Z]+)_([A-Z]+)', '\\2 - \\1', weight_scheme))
         coefs <- rbind(coefs, new_row)
         
         add_covs <- rownames(balance_list$stats_list[[weight_scheme]])[abs(balance_list$stats_list[[weight_scheme]]$smd) > 0.1]
         if (length(add_covs) > 0L) {
            remove_vars <- caret::findLinearCombos(data %>% dplyr::select(!!add_covs))$remove
            if (length(remove_vars) > 0L)
               add_covs <- add_covs[-remove_vars]
            formula <- paste0(formula, ' + ', paste(add_covs, collapse=' + '))
            
            # Fit the model, sometimes it fails...
            tryCatch(
               {
                  new_row <- getRiskDiff(fit = survey::svyglm(formula = formula,
                                                              design = survey::svydesign(ids=~0, weights=~weights, data=data),
                                                              data = data,
                                                              family = quasibinomial()), 
                                         weights = data$weights) %>% matrix() %>% t()
                  rownames(new_row) <- paste0('wtd mlt ', gsub('([A-Z]+)_([A-Z]+)', '\\2 - \\1', weight_scheme))
                  coefs <- rbind(coefs, new_row)
               },
               error = function(e) {
                  message('failed')
               }
            )
         }
      }
      
      
      # Add the treatment effect estimates, 95% CIs, and p-values for this outcome to the ATE estimates list.
      ATE_estimates[[cohort]][[o]] <- coefs
      cat('\n')
   }
   cat('\n')
}

save(ATE_estimates, file=paste0('CohortInfo/ATE_estimates_list.Rdata'))













