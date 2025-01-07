library(dplyr)
load(file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/MRSA_bacteremia/DF_before_AbxAdmin.Rdata')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_CLEANED_2017_2023_AbxAdmin.Rdata')
df <- ast; rm(ast)

abxDF <- abxDF %>%
   filter(PERSON_ID %in% unique(ast$PERSON_ID)) %>%
   filter(ABX %in% c('VANCOMYCIN', 'DAPTOMYCIN', 'LINEZOLID', 'SYNERICD', 'TRIMETHOPRIM/SULFAMETHOXAZOLE', 'CLINDAMYCIN', 'TIGECYCLINE')) %>%
   select(-END_DATE) %>%
   distinct()

# join ASTs + AbxAdmin
ast_abx <- df %>%
   mutate(JOIN_START = ORDER_DAY - 1,
          JOIN_END = ORDER_DAY + 7) %>%
   left_join(y = abxDF %>% filter(ABX %in% c('VANCOMYCIN', 'DAPTOMYCIN')),
             by = join_by(
                PERSON_ID,
                JOIN_START <= START_DAY,
                JOIN_END >= START_DAY
             )) %>%
   select(-JOIN_START, -JOIN_END) %>%
   mutate(X = as.integer(START_DAY - ORDER_DAY),
          XT = as.numeric(lubridate::as.duration(START_DATE - ORDER_DATE)) / 86400) %>% 
   filter(sum(unique(X) %in% 0:3) >= 2L, .by=c(PERSON_ID, ORDER_DAY)) %>%
   arrange(PERSON_ID, ORDER_DAY, XT)

# take only patients started on vancomycin, then switched (or didn't) to daptomycin
ast_abx <- ast_abx %>%
   group_by(PERSON_ID, ORDER_DAY) %>%
   filter(any(ABX == 'VANCOMYCIN'), first(ABX) == 'VANCOMYCIN')


# when did the switch to DAP happen?
ast_abx <- ast_abx %>%
   mutate(ALL_DAP_TIMES = list(XT[X %in% -1:7 & ABX == 'DAPTOMYCIN']),
          ALL_VAN_TIMES = list(XT[X %in% -1:7 & ABX == 'VANCOMYCIN']),
          FIRST_DAP_DAY = list(X[X %in% -1:7 & ABX == 'DAPTOMYCIN']),
          FIRST_VAN_DAY = list(X[X %in% -1:7 & ABX == 'VANCOMYCIN'])) %>%
   ungroup() %>%
   summarise(FIRST_ABX_TIME = min(XT[X %in% -1:3]),
             switched_status = any(ABX[X %in% -1:7] == 'DAPTOMYCIN'),
             .by = c(PERSON_ID, ORDER_DATE, RESULT_DATE, ORDER_DAY, RESULT_DAY, FIRST_DAP_DAY, FIRST_VAN_DAY,
                     ALL_DAP_TIMES, ALL_VAN_TIMES, GENDER, time, AGE, status, time_censored,
                     FACILITY, DAYS_SINCE_PRV, HospAcq, EARLY_CULTURE, NURSING_HOME, EMERGENCY_DEPT, TIME_BW_ADMIT_ORDER))
ast_abx$FIRST_DAP_DAY  <- sapply(ast_abx$FIRST_DAP_DAY, function(l) ifelse(length(l) == 0L, NA, min(l)))
ast_abx$FIRST_VAN_DAY  <- sapply(ast_abx$FIRST_VAN_DAY, function(l) ifelse(length(l) == 0L, NA, min(l)))
ast_abx$FIRST_DAP_TIME <- sapply(ast_abx$ALL_DAP_TIMES, function(l) ifelse(length(l) == 0L, NA, min(l)))
ast_abx$FIRST_VAN_TIME <- sapply(ast_abx$ALL_VAN_TIMES, function(l) ifelse(length(l) == 0L, NA, min(l)))

# how much overlap in time spent on VAN and DAP concurrently? and is the final VAN day LATER than the DAP day?
df <- ast_abx %>% filter(switched_status) %>% select(ALL_VAN_TIMES, ALL_DAP_TIMES)
keep <- logical(nrow(df))
dap_switch_time_r <- dap_switch_time_a <- rep(NA, nrow(df))
for (i in seq_len(nrow(df))) {
   van <- df$ALL_VAN_TIMES[[i]]
   dap <- df$ALL_DAP_TIMES[[i]]
   overlap <- min(c(max(van), max(dap))) - max(c(min(van), min(dap)))
   if (max(dap) >= max(van) & overlap <= 1) {
      keep[i] <- TRUE
      dap_switch_time_r[i] <- min(dap) - min(van)
      dap_switch_time_a[i] <- min(dap)
   }
}
w <- which(ast_abx$switched_status) # 693
ast_abx$DAP_SWITCH_TIME_A <- ast_abx$DAP_SWITCH_TIME_R <- NA
ast_abx$DAP_SWITCH_TIME_A[w] <- dap_switch_time_a
ast_abx$DAP_SWITCH_TIME_R[w] <- dap_switch_time_r
ast_abx <- ast_abx[-w[!keep], ]
rm(df, keep, i, van, dap, overlap, w, dap_switch_time_a, dap_switch_time_r)


df <- ast_abx %>%
   mutate(YEAR = lubridate::year(ORDER_DAY))
df <- df %>%
   mutate(switched_status = as.factor(switched_status))

# Feature engineering: prescriptions in the last 30 - 90 - 180 days??
# TODO






############## MODELING ###############
# start with a base model that uses just age, gender, facility, year, etc. to predict whether swtich or not
library(randomForest)
labels <- df$switched_status

# function to calculate metrics at every threshold
evaluatePredictions <- function(preds, labels) {
   o <- order(preds)
   labels_sorted <- labels[o]
   preds_sorted <- preds[o]
   
   sens <- spec <- prec <- numeric(length(preds_sorted))
   for (i in seq_along(preds_sorted)) {
      if (i %% 100 == 0L) print(i)
      # get predicted labels with given threshold
      # calculated sens and 1-spec for pred_labs vs. labs
      threshold <- preds_sorted[i]
      pred_labels_sorted <- ifelse(preds_sorted < threshold, 0, 1)
      t <- table(labels_sorted, pred_labels_sorted)
      if (ncol(t) == 1L) {
         if ('1' %in% colnames(t)) t <- cbind(t, 0)
         if ('0' %in% colnames(t)) t <- cbind(0, t)
      } else {
         t <- t[2:1, 2:1]
      }
      sens[i] <- t[1,1] / sum(t[1,])
      spec[i] <- t[2,2] / sum(t[2,])
      prec[i] <- t[1,1] / sum(t[,1])
   }
   rm(i, t, threshold, pred_labels_sorted)
   
   c_stat <- unname(wilcox.test(preds[labels], preds[!labels], paired=FALSE)$statistic / prod(table(labels)))
   
   return(list(sens=sens, spec=spec, prec=prec, auc=c_stat))
}

# TODO: implement 5-fold cross-validation...
set.seed(123L)
kfolds_idx <- sample(1:5, nrow(df), replace=TRUE)
train_metrics_lr <- train_metrics_rf <- test_metrics_lr <- test_metrics_rf <- list()
for (k in 1:5) {
   fold_id <- which(kfolds_idx == k)
   
   test_labels <- labels[fold_id]
   train_labels <- labels[-fold_id]
   
   test_df <- df[fold_id, ]
   train_df <- df[-fold_id, ]
   
   # LOGISTIC REGRESSION
   lr <- glm(switched_status ~ AGE + GENDER + FACILITY + YEAR, family='binomial', data=train_df)
   train_preds_lr <- predict(lr, newdata=train_df, type='response')
   test_preds_lr <- predict(lr, newdata=test_df, type='response')
   
   # RANDOM FOREST
   rf <- randomForest(formula = switched_status ~ AGE + GENDER + FACILITY + YEAR, data=train_df, nodesize=10)
   train_preds_rf <- predict(rf, newdata=train_df, type='prob')
   train_preds_rf <- train_preds_rf[, 'TRUE']
   test_preds_rf <- predict(rf, newdata=test_df, type='prob')
   test_preds_rf <- test_preds_rf[, 'TRUE']
   
   train_metrics_lr[[k]] <- evaluatePredictions(train_preds_lr, train_labels)
   train_metrics_rf[[k]] <- evaluatePredictions(train_preds_rf, train_labels)
      
   test_metrics_lr[[k]] <- evaluatePredictions(test_preds_lr, test_labels)
   test_metrics_rf[[k]] <- evaluatePredictions(test_preds_rf, test_labels)
}
rm(fold_id, k, test_df, train_df, lr, rf, test_labels, train_labels, train_preds_lr, train_preds_rf, test_preds_lr, test_preds_rf)


plotCurves <- function(train_metrics, test_metrics, model_name='Model unknown') {
   pdf(file = paste0('~/Desktop/EHR/EHR-mining/Scripts/AnalysisScripts/PredictTreatment/VAN_DAP_switch/', model_name, '.pdf'), width=11, height=6)
   par(mfrow=c(1,2), mgp=c(2, 0.5, 0), tck=-0.015)
   
   plot(NA, xlim=c(0,1), ylim=c(0,1), xaxs='i', yaxs='i', xlab='1 - Specificity', ylab='Sensitivity', main='ROC curve')
   sapply(train_metrics, function(x) lines(x = 1 - x$spec, y = x$sens, lwd=2, col='gray'))
   sapply(test_metrics,  function(x) lines(x = 1 - x$spec, y = x$sens, lwd=2))
   abline(a=0, b=1, lty=2, lwd=1.5)
   text(x=0.02, y=0.97, labels=round(mean(sapply(test_metrics, function(x) x$auc)), 3), adj=0)
   legend('bottomright', legend=c('train', 'test'), col=c('gray', 'black'), lwd=2)
   
   plot(NA, xlim=c(0,1), ylim=c(0,1), xaxs='i', yaxs='i', xlab='Recall', ylab='Precision', main='PR curve')
   sapply(train_metrics, function(x) lines(x = x$sens, y = x$prec, lwd=2, col='gray'))
   sapply(test_metrics,  function(x) lines(x = x$sens, y = x$prec, lwd=2))
   legend('bottomright', legend=c('train', 'test'), col=c('gray', 'black'), lwd=2)
   
   mtext(text=model_name, outer=TRUE, line=-1, at=0.5, font=2, cex=1.2)
   
   dev.off()
}

# plot ROC and PR curves
plotCurves(train_metrics_lr, test_metrics_lr, 'Logistic_regression')
plotCurves(train_metrics_rf, test_metrics_rf, 'Random_forest')






















