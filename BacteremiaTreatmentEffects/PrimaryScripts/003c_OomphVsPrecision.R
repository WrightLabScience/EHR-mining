


setwd('~/Desktop/EHR-mining/BuildCohortDatasets/')
load(file = 'ATE_estimates_list.Rdata')
load(file = 'cohort_info.Rdata')
cohort_name_dict <- setNames(sapply(cohort_info, '[[', 'clean_bug_name'), sapply(cohort_info, '[[', 'cohort_name'))

lf <- list.files('propensity_score_info/', pattern='treatment')
ATEs_all <- SMDs_all <- SIGs_all <- numeric()
for (i in seq_along(lf)) {
   load(file = paste0('propensity_score_info/', lf[i]))
   
   ATEs <- ATE_estimates[[unname(cohort_name_dict[gsub('.Rdata', '', lf[i], fixed=TRUE)])]]
   SIGs <- setNames(ATEs$mort[,'pvalue'], rownames(ATEs$mort))
   ATEs <- setNames(ATEs$mort[,'estimate'], rownames(ATEs$mort))
   SIGs <- SIGs[grep('uni', names(SIGs))]
   ATEs <- ATEs[grep('uni', names(ATEs))]
   names(ATEs)[1] <- 'raw'
   names(ATEs) <- gsub('wtd uni ', '', names(ATEs))
   w <- grep(' - (CBPS|EB|SBW)', names(ATEs))
   names(ATEs)[w] <- gsub('([A-Z]{3}) - ([A-z]+)', '\\2_\\1', names(ATEs)[w])
   names(ATEs) <- gsub(' - ', '_', names(ATEs))
   names(SIGs) <- names(ATEs)
   
   SMDs <- balance_list$summary_mat[1, match(names(ATEs), colnames(balance_list$summary_mat))]
   
   SMDs <- SMDs[-1]
   ATEs <- abs(ATEs)[-1]
   SIGs <- SIGs[-1]
   
   plot(x=SMDs, y=ATEs, xlab='mean |SMDs|', ylab='effect size', main=lf[i])
   SMDs_all <- c(SMDs_all, SMDs)
   ATEs_all <- c(ATEs_all, ATEs)
   SIGs_all <- c(SIGs_all, SIGs)
}

plot(x=SMDs_all, y=ATEs_all, xlab='mean |SMDs|', ylab='effect size')
l <- lm(ATEs_all ~ SMDs_all)
abline(l)
summary(l)

plot(x=log(SIGs_all), y=SMDs_all)
l <- lm(SMDs_all ~ log(SIGs_all))
abline(l)
summary(l)


