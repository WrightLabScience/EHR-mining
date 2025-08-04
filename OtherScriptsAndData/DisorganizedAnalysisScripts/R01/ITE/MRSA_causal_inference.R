library(dplyr)
cohort <- 'gramneg'
years <- 2017:2023
Rdata_file_path <- '~/Desktop/EHR/EHR work/RdataFiles/BuildCohortDatasets/'
df <- tibble(read.csv(paste0(Rdata_file_path, 'CleanedData/', cohort, '_modded_', paste(range(years), collapse='_'), '_bacteremia_withCATEs.csv'), sep='\t'))

vars <- read.csv('~/Desktop/EHR-mining/Scripts/AnalysisScripts/BuildCohortDatasets/mrsa/features.csv', header=F)[[1]]


# cohort <- 'mrsa'
plot_path <- paste0('~/Desktop/EHR-mining/Scripts/AnalysisScripts/BuildCohortDatasets/', cohort, '/plots/')
col_vec <- c('#3399FF', '#FFaa66')
if (cohort == 'mrsa') {
   trt <- c('VAN', 'DAP')
} else if (cohort == 'gramneg') {
   trt <- c('TZP', 'CRO')
}
names(col_vec) <- trt
ptcex <- 0.9
lwd <- 3
bw <- 0.05

cate_range <- c(-0.6, 0.6)

{
   pdf(file = paste0(plot_path, '006c_CATE_dist_by_outcome_trt.pdf'), height=5.5, width=6)
   par(mfrow=c(1,1), mgp=c(1.5, 0.3, 0), tck=-0.01, mar=c(3, 3, 2, 0.5))
   plot(NA, xlim=cate_range, ylim=c(0,4), xlab='', ylab='Density')
   title(xlab=expression(hat(tau)), cex.lab=1.2)
   lines(density(df$ite[df$outcome == 1L & df$treatment == 1L]), col='black', lwd=lwd)
   lines(density(df$ite[df$outcome == 0L & df$treatment == 1L]), col='red', lwd=lwd)
   lines(density(df$ite[df$outcome == 1L & df$treatment == 0L]), col='blue', lwd=lwd)
   lines(density(df$ite[df$outcome == 0L & df$treatment == 0L]), col='gray', lwd=lwd)
   legend('topleft', legend=c('expired-DAP', 'survived-DAP', 'expired-VAN', 'survived-VAN'), col=c('black', 'red', 'blue', 'gray'), lwd=lwd)
   dev.off()
}


{
   pdf(file = paste0(plot_path, '006a_CATE_dists.pdf'), height=8, width=8)
   par(mfrow=c(2,2), mgp=c(1.5, 0.3, 0), tck=-0.01, mar=c(3, 3, 1, 0.5))
   plot(NA, xlim=cate_range, ylim=c(0,4), xlab='', ylab='Density')
   title(xlab=expression(hat(tau)), cex.lab=1.2)
   # lines(ecdf(df$ite[df$treatment == 1L]), col=col_vec[trt[2]], lwd=lwd)
   # lines(ecdf(df$ite[df$treatment == 0L]), col=col_vec[trt[1]], lwd=lwd)
   lines(density(df$ite[df$treatment == 1L], bw=bw), col=col_vec[trt[2]], lwd=lwd)
   lines(density(df$ite[df$treatment == 0L], bw=bw), col=col_vec[trt[1]], lwd=lwd)
   legend('topright', legend=names(col_vec), col=col_vec, lwd=lwd)
   
   plot(NA, xlim=cate_range, ylim=c(0,4), xlab='', ylab='Density')
   title(xlab=expression(hat(tau)), cex.lab=1.2)
   lines(density(df$ite[df$outcome == 1L]), col='black', lwd=lwd)
   lines(density(df$ite[df$outcome == 0L]), col='red', lwd=lwd)
   legend('topright', legend=c('expired', 'survived'), col=c('black', 'red'), lwd=lwd)
#    dev.off()
# }
# 
# 
# 
# {
   #pdf(file = paste0(plot_path, '006_CATE_by_age.pdf'), height=4.5, width=10)
   #par(mfrow=c(1,2), mgp=c(1.5, 0.3, 0), tck=-0.01, mar=c(3, 3, 0.5, 0.5))
   plotCont <- function(var='age', strat='trt') {
      vec <- df[[var]]
      if (strat == 'trt') {
         w <- df$treatment == 1L
         cols <- col_vec
      } else if (strat == 'out') {
         w <- df$outcome == 1L
         cols <- c('survived' = '#009900', 'expired' = '#000000')
      }
      
      plot(x=vec, y=df$ite, pch=ifelse(w, 15, 16), col=paste0(cols[w + 1], 'aa'), xlim=c(0, 100), ylim=cate_range,
           xlab=stringr::str_to_sentence(var), ylab='', cex=ptcex, yaxt='n')
           #main=paste0(cohort, ' bacteremia - ', trt[1], ' vs. ', trt[2]))
      axis(side=2, las=1)
      mtext(text=expression(hat(tau)), side=2, line=2.5, at=mean(cate_range), las=1, adj=0.5, cex=1.2)
      legend('topleft', legend=names(cols), col=paste0(cols, 'aa'), pch=c(16, 15))
      # 'E[Y(1) - Y(0) | X]\n(CATE)'
      abline(h=0, lty=2)
      lines(smooth.spline(x=vec, y=df$ite, spar = 1))
      xpos <- 10
      yshift <- 0.05
      arrows(x0=xpos, y0=-yshift, y1=yshift, code=3, angle=45, length=0.15, lwd=2)
      text(x=xpos, y=yshift * 1.5 * c(-1,1), labels=paste0(c(trt[2], trt[1]), ' better'), font=2)
   }
   plotCont(var='age', strat='trt')
   plotCont(var='age', strat='out')
   dev.off()
}


plotStratScat <- function(var, xval_var) {
   vec <- df[[var]]
   xval <- df[[xval_var]]
   plot(x=xval, y=df$ite, pch=16, col=paste0(col_vec[vec + 1L], '22'), ylim=cate_range,
        xlab=stringr::str_to_sentence(xval_var), ylab='', main=var, yaxt='n')
   title(ylab=expression(hat(tau)), line=1.4, cex.lab=1.2)
   axis(side=2, las=1)
   abline(h=0, lty=3)
   legend('topleft', col=col_vec, legend=c(F, T), pch=16, bty='n')
}
plotCDFs <- function(var) {
   vec <- df[[var]]
   plot(ecdf(df$ite[vec == 1L]), xlim=cate_range, main=var, lwd=2, col=col_vec[1L + 1L], yaxt='n', xlab='', yaxs='i')
   title(xlab=expression(hat(tau)), line=1.4, cex.lab=1.2)
   axis(side=2, las=1)
   lines(ecdf(df$ite[vec == 0L]), lwd=2, col=col_vec[1L + 0L])
   abline(v=0, lty=3)
   legend('topleft', col=col_vec, legend=c(F, T), lwd=2, bty='n')
}


par(mfrow=c(2,2), mgp=c(1.5, 0.25, 0), tck=-0.01, mar=c(3, 3.5, 2, 1))
s <- which(vars %in% c('enc_nursing_home_hospice', 'icd_sepsis_1w', 'icd_aki_1w', 'icd_endocarditis_1w'))
s <- which(vars %in% c('enc_nursing_home_hospice', 'icd_sepsis_1w', 'enc_emergency_dept', 'icd_cellulitis_1w'))
sapply(s, function(x) plotCDFs(vars[x]))
sapply(s, function(x) plotStratScat(vars[x], 'age'))
s <- which(vars %in% c('icd_osteomyelitis_1w', 'enc_transfer', 'enc_emergency_dept', 'female'))
sapply(s, function(x) plotCDFs(vars[x]))
sapply(s, function(x) plotStratScat(vars[x], 'age'))



x <- df %>% filter(enc_transfer == 0L, icd_osteomyelitis_1w == 1L, icd_sepsis_1w == 0L) %>% pull(ite)

plot(mean(x), pch=16, ylim=c(-1, 1))
arrows()



{
   par(mfrow=c(1,2), tck=-0.008, mgp=c(1.75, 0.25, 0), mar=c(3, 3.5, 2, 1))
   
   plot(1.5, ylim=c(1/2, 3), log='y', pch=16, xaxt='n', ylab='Odds ratio', xlab='', yaxt='n',
        main='Average treatment effect')
   axis(side=2,las=1)
   abline(h=1, lty=2)
   arrows(x0=1, y0=1.2, y1=1.8, code=3, angle=90, length=0.1)
   text(x=0.58, y=c(0.92, 1.08), adj=0, labels=c('Treatment B more effective', 'Treatment A more effective'), col='darkgray')
   
   x <- seq(0, 1, 0.01)
   y <- rev(sqrt(x)) + 0.7
   plot(x=x, y=y, type='l', ylim=c(1/2, 3), log='y', yaxt='n', xlab='feature', ylab='Odds ratio', lwd=3, lty=3, xaxs='i',
        main = 'Conditional average treatment effect')
   axis(side=2,las=1)
   abline(h=1, lty=2)
   polygon(x=c(x, rev(x)), y=c(y+0.2, rev(y)-0.2), border=NA, col='#00000055')
}
# 





bvars <- df %>% select(-d30mortality, -MRSA, -Staphylococci) %>% select(where(is.integer)) %>% names
bvars <- setNames(vector('list', length(bvars)),
                  bvars)

for (v in seq_along(bvars)) {
   cat(v)
   var <- names(bvars)[v]
   
   vec1_r <- df$ITE_RAW[df[[var]] == 1L]
   vec0_r <- df$ITE_RAW[df[[var]] == 0L]
   
   vec1_p <- df$ITE_PROB[df[[var]] == 1L]
   vec0_p <- df$ITE_PROB[df[[var]] == 0L]
   
   ks <- ks.test(vec1_r, vec0_r)
   wx <- wilcox.test(vec1_r, vec0_r)
   
   bvars[[v]] <- list(
      vec1_r = vec1_r,
      vec0_r = vec0_r,
      vec1_p = vec1_p,
      vec0_p = vec0_p,
      pval = ks$p.value,
      Dstat = ks$statistic,
      pvalw = wx$p.value
   )
}
rm(vec1_r, vec0_r, vec1_p, vec0_p, ks, v, var, wx)

# 
{
   pdf(file = '~/Desktop/EHR-mining/Scripts/AnalysisScripts/R01/ITE/MRSA_ITE_plots/KS_test_pvals.pdf', width=6, height=5.2)
   par(mfcol=c(1,1), mar=c(4, 4, 3, 1), mgp=c(3, 0.5, 0), tck=-0.015)
   
   pval_threshold <- 0.01 / length(bvars)
   pvals <- sapply(bvars, '[[', 'pval')
   Dstats <- sapply(bvars, '[[', 'Dstat')
   
   plot(x = pvals, y = Dstats, log='xy', ylab='D statistic', xlab='', yaxt='n', main='KS test for ITE (strat. by features)', 
        #ylim=c(min(Dstats), max(Dstats)*1.1),
        col = ifelse(pvals < pval_threshold, 'black', 'gray'),
        pch = ifelse(Dstats > 0.1, 16, 1))
   title(xlab='p-value', line=1.5, xpd=NA)
   axis(side=2, las=1)
   axis(side=2, at=seq(0.05, 0.5, 0.05), labels=rep('', length(seq(0.05, 0.5, 0.05))))
   abline(v=pval_threshold, h=0.1, lty=3)
   # legend('topleft', pch=rep(c(16,1), each=2), col=rep(c('black', 'gray'), times=2), 
   #        legend=c('sig. p-val; D > 0.1', 'n.s. p-val; D > 0.1', 'sig. p-val; D < 0.1', 'n.s. p-val; D < 0.1'))
   dev.off()
   
   sig_vars <- which(pvals < pval_threshold & Dstats > 0.1)
   rm(pvals, pval_threshold, Dstats)
}



# PLOT CDFs - stratified ITEs by binary variables
# For example: CPT_DEF refers to receipt of ceftaroline
# patients who receiveed cefaroline (CPT_DEF = 1; black line in CDF plot)
# had lower ITE on average
{
   plot_type <- 'density'
   plot_type <- 'ecdf'
   lwd <- 1
   pdf(file = '~/Desktop/EHR-mining/Scripts/AnalysisScripts/R01/ITE/MRSA_ITE_plots/CDFs.pdf', width=10)
   par(mfrow=c(5,5), mar=c(0,0,0,0))
   
   # plot CDFs
   for (v in seq_along(sig_vars)) {
      var_data <- bvars[[sig_vars[v]]]
      vec1 <- var_data$vec1_p
      vec0 <- var_data$vec0_p
      
      if (plot_type == 'ecdf') {
         p1 <- ecdf(vec1)
         p0 <- ecdf(vec0)
      } else if (plot_type == 'density') {
         p1 <- density(vec1)
         p0 <- density(vec0)
      } else {
         cat('valid plot type not specified')
         break
      }
      plot(p1, xlim=range(c(vec1, vec0)), lwd=lwd, main='', xaxt='n', yaxt='n', xlab='', ylab='')
      lines(p0, col='blue', lwd=lwd)
      
      # add variable name
      mtext(text=names(sig_vars)[v], line=-2, at=0.01, adj=0, cex=0.8)
   }
   
   # add legend
   plot(NA, xlim=c(0,1), ylim=c(0,1), axes=F, ann=F)
   legend('center', legend=1:0, col=c('black', 'blue'), lty=1, lwd=2, bty='n')
  
   dev.off()
   rm(v, lwd, var_data, vec1, vec0, p1, p0, plot_type)
}



# Do sicker patients have higher ITEs??
{
   pdf(file = '~/Desktop/EHR-mining/Scripts/AnalysisScripts/R01/ITE/MRSA_ITE_plots/num_comorb_ITEs.pdf', width=6, height=8)
   par(mfrow=c(2,1), mgp=c(1.5,0.5,0), tck=-0.015, mar=c(3, 3, 2, 1))
   
   num_comorb <- apply(X = df %>% select(contains('ICD_')), 
                       MARGIN = 1, 
                       FUN = sum)
   # hist(num_comorb, breaks=diff(range(num_comorb)), right = FALSE,
   #      main='', xlab='Number of comorbidities (single/grouped ICD codes)')
   
   plot(x=num_comorb, y=df$ITE_PROB, pch=16, col='#00000022')
   
   spline_fit <- smooth.spline(x=num_comorb, y=df$ITE_PROB, spar = 0, all.knots = FALSE, nknots = 5)
   lines(spline_fit, col="blue")
   
   num_comorb_g <- round(num_comorb / 5) * 5
   boxplot(split(x=df$ITE_PROB, f=num_comorb_g), varwidth=TRUE, pch=16, cex=0.5, 
           main='ITEs grouped by number comorbidities', xlab='Num comorbidity (binned)', ylab='ITE')
   abline(h = 0.5, lty=3)
   
   dev.off()
}









cvars <- df %>% select(-W, -ITE_RAW, -ITE_PROB, -d30mortality, -MRSA, -Staphylococci) %>% select(!where(is.integer)) %>% names

for (var in cvars) {
   plot(x=df[[var]], y=df$ITE_PROB, xlab=var, pch=16, cex=0.4)
   abline(lm(df$ITE_PROB ~ df[[var]]))
}


hist(df$ITE_PROB, breaks=100)
abline(v = 0.5, lty=2, lwd=3)

mean(df$ITE_PROB)
qnorm(0.975) * sd(df$ITE_PROB)


mod <- df %>%
   select(-W, -ITE_RAW, -d30mortality, -MRSA, -Staphylococci) %>% 
   glm(formula = ITE_PROB ~ ., family = 'gaussian', data = .)

ite_pred <- predict(mod, newdata=df %>%
           select(-W, -ITE_RAW, -ITE_PROB, -d30mortality, -MRSA, -Staphylococci))


plot(ite_pred, df$ITE_PROB, )




















df <- data %>% filter(W != 'dap') %>% mutate(W = case_when(W == 'van' ~ 0L, .default=1L))


W <- df$W
Y <- as.integer(df$d30mortality) - 1L
X <- df %>% select(-W, -d30mortality)


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








pathDF %>%
   filter(ESBL == 1L) %>%
   select(PERSON_ID, ORDER_DAY) %>%
   left_join(
      x = .,
      y = abxtDF,
      by = join_by(PERSON_ID, ORDER_DAY)
   ) %>%
   select(matches('.+_DEF', ignore.case=FALSE)) %>%
   mutate(across(everything(), ~ as.integer(. > 0L))) %>%
   count(cephalosporin_3g_DEF | cephalosporin_4g_DEF, carbapenem_DEF)

















