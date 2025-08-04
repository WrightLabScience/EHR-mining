setwd('~/Desktop/EHR-mining/BuildCohortDatasets/')
library(dplyr)
library(causalTree)
library(survey)
library(twang)
source(file='~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/CreateNamedAbxAbbreviations.R')
source(file='~/Desktop/EHR-mining/UsefulDataForCleaning/handle_UPMC_site_names.R')
rm(site_info, code2cat, code2cat_fxn, code2fac, code2fac_fxn, code2mixed, code2mixed_fxn, fac2mixed, fac2mixed_fxn)
Rdata_file_path <- '~/Desktop/EHR/EHR work/RdataFiles/BuildCohortDatasets/CleanedData/'
base_path <- '~/Desktop/EHR-mining/BuildCohortDatasets/'
treatment_vars <- c('treatment', 'resistance')


for (s in 1) { # seq_along(treatment_vars)) {
   treatment_var <- treatment_vars[s]
   
   # Load cohort_info list.
   load(file = paste0('CohortInfo/', treatment_var, '/cohort_info.Rdata'))
   
   CATE_list <- vector('list', length(cohort_info))
   
   # Loop over each cohort.
   # for (cohort in seq_along(cohort_info)) {
   for (cohort in 6:7) {
      
      # Cohort information
      cohort_name <- cohort_info[[cohort]]$cohort_name
      clean_bug_name <- cohort_info[[cohort]]$clean_bug_name
      cat(cohort_name, '\n')
      names(CATE_list)[[cohort]] <- cohort_name
      
      # Load clean data.
      load(file = paste0(Rdata_file_path, treatment_var, '/', cohort_name, '.Rdata'))
      
      # Load balancing weights and stats.
      load(file = paste0('CohortInfo/', treatment_var, '/propensity_score_info/', cohort_name, '.Rdata'))
      
      # Create named color vector.
      col_vec <- setNames(balance_list$trt_df$col_vec, 
                          balance_list$trt_df$abx)
      
      # Data preparation for modeling.
      if ('switch' %in% names(cohort_info[[cohort]])) {
         df %>% 
            mutate(SDAY = ifelse(treatment == names(col_vec)[1], -1, floor(ABX_switch_time))) %>% 
            summarise(n=n(), sum(OUT_mortality_time < 30L, na.rm=T) / n, .by=SDAY) %>% 
            arrange(SDAY)
         df <- df %>% select(-ABX_switch_time)
      }
      
      df <- df %>%
         mutate(
            outcome = as.integer(!(is.na(OUT_mortality_time) | OUT_mortality_time > 30)),
            treatment = as.integer(treatment == names(col_vec)[2])
         ) %>%
         select(-contains('OUT', ignore.case = FALSE))
      
      data <- df %>%
         select(
            # -DEM_year_decimal,
            -contains('ENC', ignore.case = FALSE)
         )
      
      Y <- data$outcome
      W <- data$treatment
      X <- data %>% 
         select(-outcome, -treatment) %>%
         data.frame()
      remove_vars <- caret::findLinearCombos(X)$remove
      if (length(remove_vars) > 0L)
         X <- X[-remove_vars]
      
      
      set.seed(124L)
      split_by <- c('random', 'median_year')[2]
      if (split_by == 'random') 
         idx <- sample(seq_along(Y), size=floor(length(Y)*0.5), replace=FALSE)
      if (split_by == 'median_year') 
         idx <- which(X$DEM_year_decimal < median(X$DEM_year_decimal))
      
      for (i in c(-1,1)) {
         train_data <- data.frame(Y, W, X, 'weights'=balance_list$weights_list$ATE_energy)[-i * idx,]
         est_data <- data.frame(Y, W, X, 'weights'=balance_list$weights_list$ATE_energy)[i * idx,]
         
         rpart.plot(honest.causalTree(formula = Y ~ .,
                                      data = train_data %>% select(-weights, -DEM_year_decimal),
                                      treatment = train_data$W,
                                      est_data = est_data %>% select(-weights, -DEM_year_decimal),
                                      est_treatment = est_data$W,
                                      weights = train_data$weights,
                                      est_weights = est_data$weights,
                                      split.Rule = "CT",
                                      split.Honest = TRUE,
                                      cv.option = "CT", 
                                      split.Bucket = TRUE,
                                      HonestSampleSize = nrow(est_data),
                                      minsize = nrow(est_data)*0.025))
      }
      
      
      train_data <- data.frame(Y, W, X, 'weights'=balance_list$weights_list$ATE_energy)[-i * idx,]
      est_data <- data.frame(Y, W, X, 'weights'=balance_list$weights_list$ATE_energy)[i * idx,]
      rpart.plot(causalTree(formula = Y ~ .,
                            data = train_data$X,
                            treatment = train_data$W,
                            weights = balance_list$weights_list$ATE_energy,
                            split.Rule = "CT",
                            split.Honest = TRUE,
                            cv.option = "CT", 
                            split.Bucket = TRUE,
                            minsize = nrow(X)*0.025))
      
      tree <- causalTree(formula = Y ~ .,
                            data = data.frame('Age'=X$DEM_age, Y, W),
                            treatment = W,
                            weights = balance_list$weights_list$ATE_energy,
                            split.Rule = "CT",
                            split.Honest = TRUE,
                            cv.option = "CT",
                            split.Bucket = TRUE,
                            minsize = nrow(X)*0.025)
      rpart.plot(tree)
      
      est_tree <- estimate.causalTree(tree, data.frame(X, Y, W), treatment=W)
      rpart.plot(est_tree)
      
      
      data_cf <- data.frame(Y, W, X %>% select(-DEM_year_decimal))
      # By default = 1 tree per row!!
      cf <- causalForest(formula = as.formula(paste0('Y ~ ', paste(names(X %>% select(-DEM_year_decimal)), collapse=' + '))),
                         data = data_cf,
                         treatment = data_cf$W,
                         ncov_sample = ncol(data_cf) - 2L,
                         ncolx = ncol(data_cf) - 2L,
                         num.trees = 1000L,
                         weights = balance_list$weights_list$ATE_energy,
                         minsize = nrow(X) * 0.025)
      
      rpart.plot(cf$trees[[1]], roundint=FALSE)
      
      var1 <- sapply(
         X = cf$trees,
         FUN = function(x) {
            x$frame$var[rownames(x$frame) %in% 1:3 & x$frame$var != '<leaf>']
         }
      )
      sort(table(unlist(var1)))
      
      
      data_mod <- data
      w <- balance_list$weights_list$ATE_energy
      z <- data_mod$treatment
      o <- data_mod$outcome
      
      age <- seq(40, 90, 5)
      num_yth0 <- num_yth1 <- num_eld0 <- num_eld1 <- num_eld <- num_yth <- eld1 <- eld0 <- yth1 <- yth0 <- numeric(length(age))
      for (i in seq_along(age)) {
         e <- data$DEM_age > age[i]
         if (length(unique(e)) == 1L)
            next
         
         num_eld[i] <- sum(e)
         num_yth[i] <- sum(!e)
         
         num_eld1[i] <- sum(e & z)
         num_yth1[i] <- sum(!e & z)
         num_eld0[i] <- sum(e & !z)
         num_yth0[i] <- sum(!e & !z)
         
         # Elderly: treated - control
         eld1[i] <- weighted.mean(o[z & e], w[z & e]) * 100
         eld0[i] <- weighted.mean(o[!z & e], w[!z & e]) * 100
         
         # Youth: treated - control
         yth1[i] <- weighted.mean(o[z & !e], w[z & !e]) * 100
         yth0[i] <- weighted.mean(o[!z & !e], w[!z & !e]) * 100
      }
      age_df <- data.frame(age, 
                           num_eld, num_yth,
                           num_eld1, num_eld0, num_yth1, num_yth0,
                           eld = eld1 - eld0, yth = yth1 - yth0, 
                           eld1, eld0, yth1, yth0)

      # difference in % difference between old and young
      {
         par(mfrow=c(2,1), tck=-0.015, mgp=c(2, 0.4, 0), mar=c(3.5, 4, 2, 4))
         age_cols <- c('yth' = '#000000cc', 'eld' = '#FF0000cc')
         
         plot(NA, xlim=range(age_df$age), ylim=range(c(age_df$eld, age_df$yth), na.rm=T),
              xlab = 'Age (years)', ylab = '% difference 30-day mortality (DAP - VAN)',
              yaxt = 'n')
         legend('right', inset=-0.16, xpd=NA, legend=paste0(c('≤', '>'), ' age'), col=age_cols, pch=16, pt.cex=2)
         axis(side=2, las=1)
         abline(h=0, lty=2)
         pnts_lns <- function(x=age_df$age, y, pch=16, col, cex=NULL) {
            if (is.null(cex)) { 
               cex <- 1
            } else {
               omm <- range(age_df %>% select(contains('num_')))
               or <- diff(omm)
               nmm <- c(1,3)
               nr <- diff(nmm)
               cex <- (((cex - omm[1]) * nr) / or) + nmm[1]
            }
            points(x=x, y=y, col=col, pch=pch, cex=cex)
            lines(x=x, y=y, col=col)
         }
         pnts_lns(y=age_df$yth, cex=age_df$num_yth, col=age_cols['yth'])
         pnts_lns(y=age_df$eld, cex=age_df$num_eld, col=age_cols['eld'])
         
         # Show the four lines separately
         plot(NA, xlim=range(age_df$age), ylim=c(0, max(age_df %>% select(eld1, eld0, yth1, yth0))),
              xlab = 'Age (years)', ylab = '30-day mortality %', yaxt='n')
         legend('right', inset=-0.145, xpd=NA, legend=c('DAP', 'VAN'), pch=c(16,1), col=age_cols['yth'], pt.cex=2)
         axis(side=2, las=1)
         pnts_lns(y=age_df$yth1, cex=age_df$num_yth1, col=age_cols['yth'])
         pnts_lns(y=age_df$yth0, cex=age_df$num_yth0, col=age_cols['yth'], pch=1)
         pnts_lns(y=age_df$eld1, cex=age_df$num_eld1, col=age_cols['eld'])
         pnts_lns(y=age_df$eld0, cex=age_df$num_eld0, col=age_cols['eld'], pch=1)
      }
      # 
      

      if (length(grep('leaf', tree$frame$var)) <= 4L) {
         bush <- tree
      } else {
         w <- unname(tree$cptable[,'nsplit'])
         w <- w[which.min(abs(w - 3))]
         bush <- prune(tree, tree$cptable[w,1])
      }
      rpart.plot(bush)
      title(main=clean_bug_name)

      
      # Assess balance in leaves and estimate ATE.
      
      
      # Add leaves and then per-leaf estimate ATE using PS.
      data <- df %>%
         mutate(leaf = bush$where) %>%
         rename(W = treatment,
                Y = outcome) %>%
         relocate(leaf, .after=Y) #%>% filter(n() >= min_leaves, .by=leaf)
      
      
      uniq_leaves <- unique(data$leaf)
      cnames <- c('OR', 'lowerOR', 'upperOR', 'pval', 'num.trt', 'num.ctrl', 'num.trt.ess', 'num.ctrl.ess')
      coefs <- matrix(NA, nrow=length(uniq_leaves), ncol=length(cnames), dimnames=list(uniq_leaves, cnames))
      
      for (l in seq_along(uniq_leaves)) {
         cat(l, '\t')
         data_leaf <- data %>% filter(leaf == uniq_leaves[l])
         
         # Check number of treated and control units in this leaf
         num_trt <- sum(data_leaf$W == 1L)
         num_ctrl <- sum(data_leaf$W == 0L)
         if (num_trt < 50L | num_ctrl < 50L) {
            cat('Too few treated or control units.\n')
            cat('Num trt:', num_trt, '\n')
            cat('Num ctrl:', num_ctrl, '\n')
            coefs[l, c('num.trt', 'num.ctrl')] <- c(num_trt, num_ctrl)
            next
         }
         
         data_leaf <- data_leaf %>% 
            select(!c(!!names(which(sapply(X = data_leaf, 
                                           FUN = function(x) length(unique(x)) == 1L)))))
         
         # ps estimation
         ps_model <- ps(formula = W ~ .,
                        data = data_leaf %>%
                           select(-Y) %>%
                           mutate(across(where(is.logical), as.integer)) %>%
                           data.frame(),
                        n.trees = 10000,
                        keep.data = FALSE,
                        stop.method = 'es.mean',
                        estimand = estimand,
                        verbose = FALSE,
                        version = 'gbm')
         
         # find still imbalanced covariates
         tab <- bal.table(ps_model)
         tab <- tab$es.mean.ATE
         vars <- rownames(tab)[which(abs(tab[['std.eff.sz']]) > 0.1)]
         ess.sizes <- summary(ps_model)[2, c('ess.treat', 'ess.ctrl')]
         names(ess.sizes) <- c('num.trt', 'num.ctrl')
         ess.size <- sum(ess.sizes)
         np_ratio <- 0.1
         if (length(vars) > np_ratio * ess.size) {
            vars <- rownames(tail(tab[order(abs(tab$std.eff.sz)), ], n=floor(np_ratio * ess.size)))
         }
         cat(length(vars), 'will be included in outcome regression.\n')
         
         data_leaf_mod <- data_leaf %>% select(W, Y, !!vars)
         data_leaf_mod$weights <- ps_model$w[[1]]
         formula <- as.formula(paste0('Y ~ W + ', paste(vars, collapse = ' + ')))
         
         mod <- data_leaf_mod %>%
            svyglm(
               formula = formula,
               design = svydesign(ids=~1, data=data_leaf_mod, weights=~weights),
               family = quasibinomial(),
               data = .
            ) %>%
            summary()
         
         # Extract and store model coefficients on treatment
         mod <- mod$coefficients['W',]
         coefs[l,] <- c(exp(mod['Estimate'] + c(0,-1,1) * qnorm(0.975) * mod['Std. Error']), mod['Pr(>|t|)'], num_trt, num_ctrl, ess.sizes)
      }
      cat('\n\n')
      
      CATE_list[[cohort]] <- list('ATEs' = coefs, 'bush' = bush)
   }
   
   # sapply(CATE_list, '[', 'ATEs')
   
   names(CATE_list) <- gsub('_treatment', '', names(CATE_list))
   names(CATE_list)[1:5] <- gsub('([A-Z])_([a-z]+)_([A-Z]+)_([A-Z]+)', '*\\1. \\2* - \\3 vs. \\4', names(CATE_list)[1:5])
   names(CATE_list)[6] <- '*S. aureus* (MRSA) - VAN vs. switch to DAP'
   names(CATE_list)[7] <- '*E. faecalis* - VAN vs. switch to AMP'
   names(CATE_list)[8] <- '*E. faecium* (VRE) - DAP vs. LZD'
   
   
   save(CATE_list, file = paste0('~/Desktop/EHR-mining/BuildCohortDatasets/CohortInfo/', treatment_var, ifelse(only_sepsis_flags[s], '_sep', ''), '/CATE_estimates_list.Rdata'))
}






# Loads CATEs and pruned trees (bushes).
# Plot everything.
library(dplyr)
library(rpart.plot)
# load(file = '~/Desktop/EHR-mining/BuildCohortDatasets/CohortInfo/treatment/CATE_estimates_list.Rdata')
plot_path <- '~/Desktop/EHR-mining/BuildCohortDatasets/Plots/Trees/'
options(scipen=3)
treatment_vars <- c('treatment', 'treatment', 'resistance')
only_sepsis_flags <- c(FALSE, TRUE, FALSE)


## Figure 5 - main text - focal cohorts.
for (s in 1:2) {
   treatment_var <- treatment_vars[1]
   load(file = paste0('~/Desktop/EHR-mining/BuildCohortDatasets/CohortInfo/', treatment_var, ifelse(only_sepsis_flags[s], '_sep', ''), '/CATE_estimates_list.Rdata'))
   
   pdf(file=paste0(plot_path, 'CATE_trees_', ifelse(only_sepsis_flags[s], 'sep', ''), '.pdf'), height=6.75, width=7.5)
   # {
   layout(mat=matrix(c(1:2,3,3), nrow=2, byrow=T), widths=c(0.45,0.55), height=c(1,0.85))
   par(oma=c(0,0,0,1))
   for (i in 1:2) {
      if (i == 1) par(mar=c(0,4,2,0))
      if (i == 2) par(mar=c(1.5,3.5,3,1))
      coefs <- CATE_list[[i+5]]$ATEs
      bush <- CATE_list[[i+5]]$bush
      rpart.plot(bush, roundint=FALSE, digits=2, do.par=FALSE, tweak=0.95)
      title(main=gsub('*', '', names(CATE_list)[i+5], fixed=TRUE), font.main=4, xpd=NA, line=ifelse(i == 1, 0, 1), cex.main=1)
      
      if (s == 1) {
         if (i == 1) {
            xpos <- c(0.085, 0.455, 0.84)
            ypos <- -0.05
         } else if (i == 2) {
            xpos <- c(0.06, 0.35314962, 0.64010020, 0.92705079)
            ypos <- -0.04
         }  
      } else {
         xpos <- c(0.06, 0.35314962, 0.64010020, 0.92705079)
         ypos <- -0.04
      }
      text(x=xpos, y=ypos, labels=paste0('leaf: ', seq_along(xpos), letters[i]), xpd=NA, font=3)
   }
   
   par(mar=c(3,3.5,2,0.5), tck=-0.015, mgp=c(2.2, 0.4, 0))
   ylim <- 200
   plot(NA, xlim=c(1, 10), yaxs='i', ylim=c(1/ylim, ylim), log='y', xaxt='n', ylab='Odds ratio', xlab='', font.main=4, yaxt='n')
   sq <- c(2, 4, 8, 16, 32)
   axis(side=2, at=c(1/rev(sq), 1, sq), labels=paste0(rep(c('1/', ''), times=c(length(sq), length(sq)+1L)), c(rev(sq), 1, sq)), las=1)
   abline(h=1, lty=2, lwd=0.5)
   xpos <- list()
   if (s == 1) {
      xpos[[1]] <- c(1.186, 2.56, 3.95)
   } else {
      xpos[[1]] <- c(1.33, 2.29, 3.25, 4.2)
   }
   xpos[[2]] <- c(5.93, 7.18, 8.46, 9.721)
   # abline(v = unlist(xpos), xpd=NA)
   axis(side=1, at=unlist(xpos), labels=paste0('leaf: ', c(1:(s+2), 1:4), rep(letters[1:2], times=c(s+2,4))))
   # }
   for (i in 1:2) {
      coefs <- CATE_list[[i+5]]$ATEs
      coefs <- coefs[order(rownames(coefs)),]
      points(x=xpos[[i]], y=coefs[,'OR'], pch=16, cex=1.4)
      arrows(x0=xpos[[i]], y0=coefs[,'lowerOR'], y1=coefs[,'upperOR'], code=3, angle=90, length=0.075, lwd=1.25, xpd=NA)
      text(x=xpos[[i]], y=coefs[,'upperOR']*1.5, labels=paste0('p=', unname(formatC(coefs[,'pval'], format='g', digits=2, flag='#'))), 
           xpd=NA, font=ifelse(coefs[,'pval'] < 0.05, 2, 1), cex=0.9)
      text(x=rep(xpos[[i]], each=2), y=ylim*rep(c(0.4,0.7), times=nrow(coefs)), labels=paste0('n=', round(as.vector(t(coefs[,c('num.trt.ess', 'num.ctrl.ess')])), 2)), 
           col=rep(c('#FFaa66', '#3399FF'), times=nrow(coefs)), cex=0.9)
      text(x=xpos[[i]]+0.075, y=coefs[,'OR'], labels=paste0('aOR=', round(coefs[,'OR'], 2)), adj=0, xpd=NA, cex=0.9)
   }
   
   dev.off()
}






library(dplyr)
library(rpart.plot)
plot_path <- '~/Desktop/EHR-mining/BuildCohortDatasets/Plots/Trees/'
options(scipen=3)
treatment_vars <- c('treatment', 'treatment', 'resistance')
only_sepsis_flags <- c(FALSE, TRUE, FALSE)

s <- 2
treatment_var <- treatment_vars[s]

load(file = paste0('~/Desktop/EHR-mining/BuildCohortDatasets/CohortInfo/', treatment_var, ifelse(only_sepsis_flags[s], '_sep', ''), '/CATE_estimates_list.Rdata'))
load(file = paste0('~/Desktop/EHR-mining/BuildCohortDatasets/CohortInfo/', treatment_var, ifelse(only_sepsis_flags[s], '_sep', ''), '/ATE_estimates_list.Rdata'))

## Appx figures
# for (i in seq_along(CATE_list)) {
for (i in 1:2) {
   col_vec <- ATE_estimates[[i+5]]$info$ABX
   coefs <- CATE_list[[i+5]]$ATEs
   coefs <- coefs[order(rownames(coefs)),]
   bush <- CATE_list[[i+5]]$bush
   cohort_name <- gsub('*', '', names(CATE_list)[i+5], fixed=TRUE)
   num_leaves <- nrow(coefs)
   
   pdf(file=paste0(plot_path, cohort_name, ifelse(only_sepsis_flags[s], 'sep', ''), '.pdf'), height=6, width=7)
   layout(mat=matrix(1:2,nrow=2), height=c(0.6,0.5))
   par(oma=c(0,4,2,3))
   
   rpart.plot(bush, roundint=FALSE)
   # title(main=c('MRSA', 'E. faecalis')[i], font.main=c(2,4)[i], line=)
   # title(main = cohort_name, line=4, xpd=NA, font.main=4)
   
   ylim <- 60
   par(mar=c(4,0,1,0), tck=-0.015, mgp=c(2.5, 0.4, 0))
   plot(NA, xlim=c(1, 4), yaxs='i', ylim=c(1/ylim, ylim), log='y', xaxt='n', ylab='Odds ratio', xlab='', font.main=4, yaxt='n', xpd=NA)
   sq <- c(2, 4)
   # axis(side=2, at=c(1/rev(sq), 1, sq), labels=paste0(rep(c('1/', ''), times=c(length(sq), length(sq)+1L)), c(rev(sq), 1, sq)), las=1)
   axis(side=2, at=c(1/rev(sq), 1, sq), labels=c(0.25, 0.5, 1, 2, 4), las=1)
   abline(h=1, lty=2, lwd=0.5)
   
   yscale <- 15
   axis(side=4, at=1/(ylim/yscale), labels=paste0('favors\n', names(col_vec)[1]), las=1, col.axis=col_vec[1], tick=F)
   axis(side=4, at=ylim/yscale, labels=paste0('favors\n', names(col_vec)[2]), las=1, col.axis=col_vec[2], tick=F)
   
   
   if (num_leaves == 2) xpos <- c(1.96, 3.05)
   if (num_leaves == 3) xpos <- c(1.25, 2.51, 3.75)
   if (num_leaves == 4) xpos <- c(1.25, 2.09, 2.93, 3.77)
   
   if (i == 7) xpos[2:4] <- xpos[2:4] * 0.99
   if (i == 3 | i == 4) xpos[1] <- xpos[1] + 0.03
   if (i == 3) xpos[2:3] <- xpos[2:3] * 1.005
   
   axis(side=1, at=xpos, labels=paste0('leaf: ', seq_len(num_leaves)), xpd=NA, font=3)
   axis(side=3, at=xpos, labels=paste0('leaf: ', seq_len(num_leaves)), xpd=NA, font=3, tick=F)
   #axis(side=1, at=xpos, labels=paste0('ess=', round(as.vector(t(coefs[,'num.trt.ess'])), 2)), col.axis=col_vec[1], line=2, tick=F)
   #axis(side=1, at=xpos, labels=paste0('ess=', round(as.vector(t(coefs[,'num.ctrl.ess'])), 2)), col.axis=col_vec[2], line=1, tick=F)
   
   points(x=xpos, y=coefs[,'OR'], pch=16, cex=1.2)
   arrows(x0=xpos, y0=coefs[,'lowerOR'], y1=coefs[,'upperOR'], code=3, angle=90, length=0.075, lwd=1.25, xpd=NA)
   text(x=xpos, y=coefs[,'upperOR']*1.2, labels=paste0('p=', unname(formatC(coefs[,'pval'], format='g', digits=2, flag='#'))), xpd=NA, font=ifelse(coefs[,'pval'] < 0.05, 2, 1), adj=c(0.5,0))
   text(x=xpos, y=coefs[,'lowerOR']/1.2, labels=paste0('aOR=', round(coefs[,'OR'], 2)), adj=c(0.5,1), xpd=NA)
   
   dev.off()
}












# cf <- causal_forest(
#    X = X,
#    Y = Y,
#    W = W,
#    Y.hat = NULL,
#    W.hat = NULL #ps[[estimand]]$ps[[1]]
# )
# 
# data.frame(sort(setNames(variable_importance(cf)[,1], colnames(X))))
# 
# # cfpred <- predict(cf)
# cfpred <- predict(cf, estimate.variance = TRUE)
# hist(cfpred$predictions)
# mean(cfpred$predictions)
# 
# n <- nrow(df)
# train <- sample(1:n, n / 2)
# train.forest <- causal_forest(X[train, ], Y[train], W[train])
# eval.forest <- causal_forest(X[-train, ], Y[-train], W[-train])
# rate <- rank_average_treatment_effect(eval.forest,
#                                       predict(train.forest, X[-train, ])$predictions)
# plot(rate)
# paste("AUTOC:", round(rate$estimate, 2), "+/", round(1.96 * rate$std.err, 2))
# 
# tau.hat <- cfpred$predictions
# high.effect <- tau.hat > median(tau.hat)
# ate.high <- average_treatment_effect(cf, subset = high.effect)
# ate.low <- average_treatment_effect(cf, subset = !high.effect)
# 
# ate.high[["estimate"]] - ate.low[["estimate"]] +
#    c(-1, 1) * qnorm(0.975) * sqrt(ate.high[["std.err"]]^2 + ate.low[["std.err"]]^2)
# 
# 
# e.hat <- cf$W.hat
# p <- mean(W)
# Y.hat.0 <- cf$Y.hat - e.hat * tau.hat
# Y.hat.1 <- cf$Y.hat + (1 - e.hat) * tau.hat
# 
# bias <- (e.hat - p) * (p * (Y.hat.0 - mean(Y.hat.0))  + (1 - p) * (Y.hat.1 - mean(Y.hat.1)))
# 
# hist(bias / sd(Y))
# 
# # p <- cfpred$predictions
# # pu <- cfpred$predictions + qnorm(0.975) * cfpred$variance.estimates
# # pl <- cfpred$predictions - qnorm(0.975) * cfpred$variance.estimates
# # o <- order(p)
# # p <- p[o]
# # pu <- pu[o]
# # pl <- pl[o]
# # plot(p, pch=16, cex=0.1)
# # points(pu)
# # points(pl)
# 
# 
# library(rpart)
# tree <- rpart(tau.hat ~ ., data=data.frame(X))
# rpart.plot(tree)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# library(dplyr)
# library(h2o)
# h2o.init()
# 
# data <- df %>%
#    mutate(
#       outcome = as.factor(!(is.na(OUT_mortality_time) | OUT_mortality_time > 30))
#    ) %>%
#    relocate(outcome, .after=treatment) %>%
#    select(!contains('OUT', ignore.case=FALSE)) %>%
#    select(-treatment)
# 
# data <- as.h2o(data)
# 
# y <- 'outcome'
# X <- setdiff(names(data), y)
# 
# splits <- h2o.splitFrame(data, ratios = 0.75, seed = 1234)
# train <- splits[[1]]
# valid <- splits[[2]]
# 
# # Train XGBoost model
# xgb_model <- h2o.xgboost(
#    x = X,
#    y = y,
#    training_frame = train,
#    validation_frame = valid,
#    ntrees = 100,
#    max_depth = 4,
#    learn_rate = 0.1,
#    seed = 1234
# )
# 
# # View model performance
# h2o.performance(xgb_model, valid = TRUE)
# 
# 
# 
# df_og <- df
# 
# 
# df <- df_og %>%
#    mutate(
#       outcome = as.factor(!(is.na(OUT_mortality_time) | OUT_mortality_time > 30)),
#       treatment = as.integer(treatment == names(col_vec)[1])
#    ) %>%
#    relocate(outcome, .after=treatment) %>%
#    select(!contains('OUT', ignore.case=FALSE))
# 
# # Define variables
# treatment_col <- 'treatment'
# outcome_col <- 'outcome'
# feature_cols <- setdiff(colnames(df), c(treatment_col, outcome_col))
# 
# # Split data into treated and control groups
# treated <- df[df[[treatment_col]] == 1, ]
# control <- df[df[[treatment_col]] == 0, ]
# 
# # ---- Train XGBoost on each group with cross-validation ----
# 
# model_all <- h2o.xgboost(
#    x = feature_cols,
#    y = outcome_col,
#    training_frame = as.h2o(df),
#    nfolds = 5,
#    keep_cross_validation_predictions = TRUE,
#    ntrees = 100,
#    max_depth = 4,
#    learn_rate = 0.1,
#    seed = 1
# )
# 
# # Treated model
# model_treated <- h2o.xgboost(
#    x = feature_cols,
#    y = outcome_col,
#    training_frame = as.h2o(treated),
#    nfolds = 5,
#    keep_cross_validation_predictions = TRUE,
#    ntrees = 100,
#    max_depth = 4,
#    learn_rate = 0.1,
#    seed = 1
# )
# 
# # Control model
# model_control <- h2o.xgboost(
#    x = feature_cols,
#    y = outcome_col,
#    training_frame = as.h2o(control),
#    nfolds = 5,
#    keep_cross_validation_predictions = TRUE,
#    ntrees = 100,
#    max_depth = 4,
#    learn_rate = 0.1,
#    seed = 1
# )
# 
# # ---- Get predictions for all units using both models ----
# 
# mu_hat <- h2o.cross_validation_holdout_predictions(model_all)
# mu_hat <- as.vector(mu_hat[['TRUE']])
# 
# # Predict potential outcome under treatment (T=1)
# mu1_hat <- h2o.predict(model_treated, as.h2o(df))
# #mu1_hat <- h2o.cross_validation_holdout_predictions(model_treated)
# 
# # Predict potential outcome under control (T=0)
# mu0_hat <- h2o.predict(model_control, as.h2o(df))
# #mu0_hat <- h2o.cross_validation_holdout_predictions(model_control)
# 
# mu1_hat <- as.vector(mu1_hat[['TRUE']])
# mu0_hat <- as.vector(mu0_hat[['TRUE']])
# 
# plot(x=mu0_hat, y=mu1_hat - mu0_hat, xlim=c(0,1), ylim=c(-1, 1), pch=16, col='#00000022')
# 
# 
# # Estimate ITEs
# ite <- mu1_hat - mu0_hat
# ite <- as.data.frame(ite)
# ite <- 
#    
#    df_with_ite <- h2o.cbind(df, mu1_hat, mu0_hat, ite)
# colnames(df_with_ite)[(ncol(df)+1):(ncol(df)+3)] <- c("mu1_hat", "mu0_hat", "ITE_hat")
# 
# # View results
# head(as.data.frame(df_with_ite))
# 

