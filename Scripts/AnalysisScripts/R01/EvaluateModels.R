library(dplyr)
setwd('~/Desktop/EHR/EHR work/RdataFiles/R01/')

model_files <- list.files('glmnet_fits/', pattern='_fit\\.Rdata')

evalDF <- vector(mode = 'list', length = length(model_files))
names(evalDF) <- gsub('^(.+)_fit.Rdata$', '\\1', model_files)

for (i in seq_along(model_files)) {
   cat(i, '\n')
   load(file = paste0('glmnet_fits/', model_files[i]))
   target_var <- names(evalDF)[i]
   
   labels <- cv_fit$labels
   prevalence <- sum(labels == 1L) / length(labels)
   probs <- cv_fit$fit.preval[, which(cv_fit$lambda == cv_fit$lambda.min)]
   
   # sort the prediction probabilities and labels so decreasing
   o <- order(probs, decreasing=TRUE)
   probs <- probs[o]
   labels <- labels[o]
   
   # calculate metrics across every possible threshold
   tp_cumsum <- cumsum(labels == 1)
   fp_cumsum <- cumsum(labels == 0)
   total_pos <- sum(labels == 1)
   total_neg <- sum(labels == 0)
   
   sens <- tp_cumsum / total_pos  # Sensitivity
   spec <- (total_neg - fp_cumsum) / total_neg  # Specificity
   prec <- tp_cumsum / (tp_cumsum + fp_cumsum)  # Precision
   prec[is.nan(prec)] <- 0
   
   fpr <- 1 - spec
   auroc <- sum(( (fpr[-1] - fpr[-length(fpr)]) * (sens[-1] + sens[-length(sens)]) ) / 2)
   auprc <- sum(( (sens[-1] - sens[-length(sens)]) * (prec[-1] + prec[-length(prec)]) ) / 2)
   
   # add metrics data to list
   evalDF[[i]] <- list(
      table = data.frame(sens, fpr, prec),
      data = c(
         'Prevalence' = prevalence,
         'AUROC' = auroc,
         'AUPRC' = auprc,
         'AUPRC_adj' = auprc - prevalence 
      )
   )
   
   # remove all variables just in case
   rm(cv_fit, target_var, labels, prevalence, probs, o, tp_cumsum, fp_cumsum, total_pos, total_neg, sens, spec, prec, fpr, auroc, auprc)
}

save(evalDF, file = 'model_metrics/evalDF.Rdata')


o <- order(sapply(evalDF, '[[', 2)['AUROC',])
evalDF <- evalDF[o]
{
   pdf(file = 'model_metrics/ROC.pdf', width=11)
   par(mfrow=c(11,11), mar=c(0,0,0,0), mgp=c(2, 0.35, 0), tck=-0.015, oma=c(0,0,0,0))
   for (i in seq_along(evalDF)) {
      sens <- evalDF[[i]]$table$sens
      fpr <- evalDF[[i]]$table$fpr
      
      plot(x=fpr, y=sens, xlim=c(0,1), ylim=c(0,1), type='l', xaxs='i', yaxs='i', xlab='FPR', ylab='TPR', yaxt='n', xaxt='n')
      text(x=0.01, y=c(0.95, 0.85), adj=0, labels=c(names(evalDF)[i], round(evalDF[[i]]$data['AUROC'], 2)))
      abline(a=0, b=1, lty=3)
   }
   dev.off()
}

o <- order(sapply(evalDF, '[[', 2)['AUPRC_adj',])
evalDF <- evalDF[o]
{
   pdf(file = 'model_metrics/PRC.pdf', width=11)
   par(mfrow=c(11,11), mar=c(0,0,0,0), mgp=c(2, 0.35, 0), tck=-0.015, oma=c(0,0,0,0))
   for (i in seq_along(evalDF)) {
      sens <- evalDF[[i]]$table$sens
      prec <- evalDF[[i]]$table$prec
      
      plot(x=sens, y=prec, xlim=c(0,1), ylim=c(0,1), type='l', xaxs='i', yaxs='i', xlab='Recall', ylab='Precision', yaxt='n', xaxt='n')
      text(x=0.01, y=c(0.95, 0.85), adj=0, labels=c(names(evalDF)[i], round(evalDF[[i]]$data['AUPRC_adj'], 2)))
      abline(h=evalDF[[i]]$data['Prevalence'], lty=3)
   }
   dev.off()
}













