library(dplyr)
library(grf)

load(file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/MRSA_bacteremia/MRSA_BSI_DAP_VAN.Rdata')
df <- dfx %>%
   mutate(time = lubridate::decimal_date(ORDER_DAY) - 2017) %>%

load(file = '/Users/samblechman/Desktop/EHR-mining/UsefulDataForCleaning/UPMC_site_groups.Rdata')
Rdata_file_path <- '/Users/samblechman/Desktop/EHR/EHR work/RdataFiles/R01/OutcomePredictionModel/'
load(file = paste0(Rdata_file_path, 'data_feat_outcomes_processed_all.Rdata'))
target <- 'd30mortality'
data <- data %>%
   filter(MRSA == 1L) %>%
   mutate(time = lubridate::decimal_date(ORDER_DAY) - 2017) %>%
   select(
      time, 
      !!target,
      AGE, FEMALE, # demographics
      MRSA:Enterococci, # flag for bloodstream-infecting pathogen group
      RESPIRATORY_ISOLATE, # don't keep identity of respiratory isolate, just presence flag
      TIME_TO_CONC, # keep this instead of EmpDisc/DefDisc flags
      VAN_DEF, DAP_DEF,
      # matches('^[A-Z]{2,3}_(EMP|DEF)$'), # treatment variables - keep specific abx, not group
      matches('^(ICD|LAB|ENC)_') # keep ICD, LAB, and ENC variables
   ) %>% 
   filter(VAN_DEF != 0 | DAP_DEF != 0) %>%
   mutate(
      W = case_when(
         VAN_DEF == 1L & DAP_DEF == 1L ~ 'switch',
         VAN_DEF == 1L & DAP_DEF == 0L ~ 'van',
         VAN_DEF == 0L & DAP_DEF == 1L ~ 'dap'
      )
   ) %>%
   select(-VAN_DEF, -DAP_DEF) %>%
   relocate(W, .after=!!target) %>%
   mutate(
      across(.cols = c(ENC_FACILITY1, ENC_FACILITY2),
             .fns = ~ gsub('\\.| |-', '', as.character(.)))
   )

sites <- list(
   Presbyterian = 'Presbyterian',
   Mercy = 'Mercy',
   Altoona = 'Altoona',
   Shadyside = 'Shadyside',
   Hamot = 'Hamot',
   Regional = c('Jameson', 'Williamsport'),
   Community = c('East', 'StMargaret', 'McKeesport', 'Passavant', 'MageeW', 'SOL', 'Childrens'),
   Rural = c('Bedford', 'Northwest', 'Horizon', 'Chatauqua', 'LOC')
)

df <- data %>% filter(W != 'dap') %>% mutate(W = case_when(W == 'van' ~ 0L, .default=1L))
W <- df$W
Y <- as.integer(df$d30mortality) - 1L
X <- df %>% select(-W, -d30mortality)

# handle categorical variable encodings - one-hot:
for (s in seq_along(sites)) {
   if (names(sites)[s] == 'Rural')
      next
   X[[paste0(names(sites)[s], '_1')]] <- as.integer(X$ENC_FACILITY1 %in% sites[[s]])
   X[[paste0(names(sites)[s], '_2')]] <- as.integer(X$ENC_FACILITY2 %in% sites[[s]])
}
rm(s, sites)
X <- X %>% select(-ENC_FACILITY1, -ENC_FACILITY2)
X <- as.matrix(X)



# Train a causal forest.
tau.forest <- causal_forest(X, Y, W)

# Estimate treatment effects for the training data using out-of-bag prediction.
tau.hat.oob <- predict(tau.forest, estimate.variance = TRUE)
hist(tau.hat.oob$predictions)

plot(density(tau.hat.oob$predictions))

# Estimate the conditional average treatment effect on the full sample (CATE).
average_treatment_effect(tau.forest, target.sample = "all", method = 'TMLE')

# Estimate the conditional average treatment effect on the treated sample (CATT).
average_treatment_effect(tau.forest, target.sample = "treated")
average_treatment_effect(tau.forest, target.sample = "control")
average_treatment_effect(tau.forest, target.sample = "overlap")


rate.oob <- rank_average_treatment_effect(forest = tau.forest, 
                                          priorities = tau.hat.oob$predictions, 
                                          target = 'AUTOC')
plot(rate.oob)
t.stat.oob <- rate.oob$estimate / rate.oob$std.err
# Compute a two-sided p-value Pr(>|t|)
p.val = 2 * pnorm(-abs(t.stat.oob))
# Compute a one-sided p-value Pr(>t)
p.val.onesided = pnorm(t.stat.oob, lower.tail = FALSE)


# identifying heterogeneous treatment effects
train <- sample(seq_along(W), size=length(W) / 2)
train.forest <- causal_forest(X[train, ], Y[train], W[train])
eval.forest <- causal_forest(X[-train, ], Y[-train], W[-train])
rate <- rank_average_treatment_effect(eval.forest,
                                      predict(train.forest, X[-train, ])$predictions)
plot(rate)
paste("AUTOC:", round(rate$estimate, 3), "+/", round(qnorm(0.975) * rate$std.err, 3))



best_linear_projection(forest = tau.forest, 
                       A = X[,1:4])






