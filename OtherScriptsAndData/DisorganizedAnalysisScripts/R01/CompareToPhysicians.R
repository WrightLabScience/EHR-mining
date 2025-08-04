# physicians "predict" ESBL by prescribing a carbapenem empirically
# our model does the same
library(dplyr)
setwd('~/Desktop/EHR/EHR work/RdataFiles/R01/')
load(file = 'model_metrics/evalDF.Rdata')

# fit <- evalDF$RxC_carbapenem
fit <- evalDF$ESBL

load(file = '~/Desktop/EHR/EHR work/RdataFiles/R01/bact_outcomes_variables.Rdata')
phys_confmat <- table(outcomesDF$ESBL, outcomesDF$RxC_carbapenem)
fisher.test(phys_confmat)

calcMetrics <- function(t) {
   sens <- t[1,1] / sum(t[1,])
   spec <- t[2,2] / sum(t[2,])
   prec <- t[1,1] / sum(t[,1])
   ngpv <- t[2,2] / sum(t[,2])
   # r <- function(x) round(x * 100, 1)
   # cat('Sens:', r(sens), '\n')
   # cat('Prec:', r(prec), '\n')
   # cat('Spec:', r(spec), '\n')
   # cat('Ngpv:', r(ngpv), '\n')
   return(c('sens' = sens, 'spec' = spec, 'prec' = prec))
}
phys_metrics <- calcMetrics(phys_confmat)
phys_metrics


sens <- fit$table$sens
fpr <- fit$table$fpr
prec <- fit$table$prec

par(mfrow=c(2,1), mar=c(3,3,2,1), mgp=c(2, 0.5, 0), tck=-0.015)
plot(x=fpr, y=sens, xlim=c(0,1), ylim=c(0,1), type='l', xaxs='i', yaxs='i', xlab='FPR', ylab='TPR', yaxt='n', main='ROC')
axis(side=2, las=1)
text(x=0.01, y=c(0.95, 0.85), adj=0, labels=c('ESBL prediction', round(fit$data['AUROC'], 2)))
abline(a=0, b=1, lty=3)

points(x=1-phys_metrics['spec'], y=phys_metrics['sens'])

plot(x=sens, y=prec, xlim=c(0,1), ylim=c(0,1), type='l', xaxs='i', yaxs='i', xlab='Recall', ylab='Precision', yaxt='n', main='Precision-recall')
axis(side=2, las=1)
text(x=0.01, y=c(0.95, 0.85), adj=0, labels=c('ESBL prediction', round(fit$data['AUROC'], 2)))
abline(h=fit$data['Prevalence'], lty=3)

points(x=phys_metrics['sens'], y=phys_metrics['prec'])















