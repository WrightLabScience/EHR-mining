
evalModel_h2o_all <- function(model, datasets, labels, test=TRUE, xval=TRUE, xval_probs=NULL, pred_col='p1', mode='classification', plot_data=TRUE) {
   l <- list(train = evalModel_h2o(model=model, dataset=datasets$train, labels=labels$train, pred_col=pred_col, mode=mode, plot_data=plot_data))
   if (test) {
      l$test <- evalModel_h2o(model=model, dataset=datasets$test, labels=labels$test, pred_col=pred_col, mode=mode, plot_data=plot_data)  
   }
   if (xval) {
      if (is.null(xval_probs)) {
         xval_probs <- as.vector(h2o.cross_validation_predictions(model)[[pred_col]])
      }
      l$cross_val <- evalModel(probs = xval_probs, labels = labels$train, mode='classification', plot_data=plot_data)
   }
   return(l)
}

evalModel_h2o <- function(model, dataset, labels, pred_col='p1', mode='classification', plot_data=TRUE) {
   pred <- h2o.predict(object = model,
                       newdata = as.h2o(dataset))
   pred <- as.data.frame(pred)
   return(evalModel(probs = pred[[pred_col]],
                    labels = labels,
                    mode = mode,
                    plot_data = plot_data))
}

calc_auroc <- function(fpr, sens) {
   return(sum(( (fpr[-1] - fpr[-length(fpr)]) * (sens[-1] + sens[-length(sens)]) ) / 2))
}

calc_auprc <- function(sens, prec) {
   return(sum(( (sens[-1] - sens[-length(sens)]) * (prec[-1] + prec[-length(prec)]) ) / 2))
}

evalModel <- function(probs, labels, mode='classification', plot_data=TRUE) {
   if (mode == 'classification') {
      prevalence <- sum(labels == 1L) / length(labels)
      
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
      auroc <- calc_auroc(fpr, sens)
      auprc <- calc_auprc(sens, prec)
      
      l <- list(
         data = c(
            'Prevalence' = prevalence,
            'AUROC' = auroc,
            'AUPRC' = auprc,
            'AUPRC_adj' = auprc - prevalence 
         )
      )
      
      if (plot_data) {
         l$table <- data.frame(sens, fpr, prec)
      }
      
      # add metrics data to list
      return(l)  
   }
   
   if (mode == 'regression') {
      
   }
}

plotROC_PRC <- function(evalDF, evalDF_others=NULL) {
   col_vec <- c('blue', 'red', 'darkgreen', 'gray', 'purple')
   xpos <- 0.99
   ypos <- seq(0.04, 1, by=0.06)
   
   fpr <- evalDF$table$fpr
   sens <- evalDF$table$sens
   prec <- evalDF$table$prec
   
   ### RECEIVER OPERATING CHARACTERISTIC curve
   par(mfrow=c(2,1), mar=c(3,3,2,1), mgp=c(2, 0.5, 0), tck=-0.015)
   plot(x=fpr, y=sens, xlim=c(0,1), ylim=c(0,1), type='l', xaxs='i', yaxs='i', xlab='FPR', ylab='TPR', main='ROC curve', yaxt='n')
   axis(side=2, las=1)
   text(x=xpos, y=ypos[1], adj=1, labels=paste0('AUC (train): ', round(evalDF$data['AUROC'], 3)))
   abline(a=0, b=1, lty=3)
   
   if (!is.null(evalDF_others)) {
      # there were additional evalDF objects given: validation, testing, cross-validation, other models
      # let's loop through them and add them to the above plot
      for (i in seq_along(evalDF_others)) {
         evalDF_ <- evalDF_others[[i]]
         name <- names(evalDF_others)[i]
         fpr_ <- evalDF_$table$fpr
         sens_ <- evalDF_$table$sens
         lines(x=fpr_, y=sens_, col=col_vec[i])
         text(x=xpos, y=ypos[i+1], adj=1, 
              labels=paste0('AUC (', name, '): ', round(evalDF_$data['AUROC'], 3)), col=col_vec[i])  
      }
   }
   
   ### PRECISION - RECALL curve
   plot(x=sens, y=prec, xlim=c(0,1), ylim=c(0,1), type='l', xaxs='i', yaxs='i', xlab='Recall', ylab='Precision', main='PR curve', yaxt='n')
   axis(side=2, las=1)
   text(x=xpos, y=ypos[1], adj=1, labels=paste0('AUC (train):', round(evalDF$data['AUPRC'], 3)))
   abline(h=evalDF$data['Prevalence'], lty=3)  
   
   if (!is.null(evalDF_others)) {
      # there were additional evalDF objects given: validation, testing, cross-validation, other models
      # let's loop through them and add them to the above plot
      for (i in seq_along(evalDF_others)) {
         evalDF_ <- evalDF_others[[i]]
         name <- names(evalDF_others)[i]
         sens_ <- evalDF_$table$sens
         prec_ <- evalDF_$table$prec
         lines(x=sens_, y=prec_, col=col_vec[i])
         text(x=xpos, y=ypos[i+1], adj=1, 
              labels=paste0('AUC (', name, '): ', round(evalDF_$data['AUPRC'], 3)), col=col_vec[i])  
      }
   }
}

plotCDFs <- function() {
   
}