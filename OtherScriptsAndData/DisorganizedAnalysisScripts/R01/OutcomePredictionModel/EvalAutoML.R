library(dplyr)
library(h2o)
h2o.init()
# h2o.no_progress()
source('~/Desktop/EHR-mining/Scripts/CleaningScripts/evaluate_model_outputs_and_plot.R')
Rdata_file_path <- '/Users/samblechman/Desktop/EHR/EHR work/RdataFiles/R01/OutcomePredictionModel/'
load(file = paste0(Rdata_file_path, 'data_feat_outcomes_processed_all.Rdata'))
target <- 'd30mortality'
data <- data %>%
   mutate(year = as.factor(lubridate::year(ORDER_DAY))) %>% 
   select(
      year, 
      !!target,
      AGE, FEMALE, # demographics
      MRSA:Enterococci, # flag for bloodstream-infecting pathogen group
      RESPIRATORY_ISOLATE, # don't keep identity of respiratory isolate, just presence flag
      TIME_TO_CONC, # keep this instead of EmpDisc/DefDisc flags
      matches('^[A-Z]{2,3}_(EMP|DEF)$'), # treatment variables - keep specific abx, not group
      matches('^(ICD|LAB|ENC)_') # keep ICD, LAB, and ENC variables
   )

data <- data
labels <- as.integer(data[[target]]) - 1L


# TRAIN MODEL
xgb <- h2o.xgboost(
   y = target,
   training_frame = as.h2o(data[w,]),
   fold_column = 'year',
   nfolds = 0,
   keep_cross_validation_models = TRUE,
   keep_cross_validation_predictions = TRUE,
   keep_cross_validation_fold_assignment = TRUE,
   colsample_bylevel = 0.7,
   colsample_bytree = 1,
   learn_rate = 0.3,
   max_depth = 15, # (2 did barely worse)
   ntrees = 20, # (not huge difference from 200)
   reg_alpha = 10,
   reg_lambda = 1,
   sample_rate = 0.7,
   seed = 123L
)



## collapse cross validation predictions into one dataframe
xval_probs <- rep(NA, nrow(data))
idxs_list <- tapply(X = w, # list, each element contains the indices for the corresponding folds
                    INDEX = data$year, 
                    FUN = list)
xval_pred_list <- h2o.cross_validation_predictions(xgb)
for (i in seq_along(idxs_list)) {
   xval_probs[idxs_list[[i]]] <- as.vector(xval_pred_list[[i]]$p1[idxs_list[[i]]])
}
rm(idxs_list, xval_pred_list, i)


# PLOT ROC and PR curves to assess performance on train set, cross-validation folds, and (optionally) test set
evalDF <- evalModel_h2o_all(
   model = xgb,
   datasets = list(
      train = data
   ),
   labels = list(
      train = labels
   ),
   xval = TRUE,
   xval_probs = xval_probs,
   test = FALSE
)
plotROC_PRC(evalDF = evalDF$train,
            evalDF_others = list(
               'cross-val' = evalDF$cross_val
            ))


# I have each xval model here:
# loop through and make interventional predictions, then combine

# interventional distributions
# make a function that takes as input:
#     filter by bug
#     two treatments to "swap"
# an outputs:
#     six pairwise comparisons: (0,0)->(1,0), (0,0)->(0,1), (1,0)->(0,1), (0,1)->(1,1), (1,0)->(1,1)
#     for each comparison:
#        - football plot
#        - statistical summary (how many with response above or below X%)
#     dataframe
#        6 columns - one for each pairwise difference (effect size)
#        4 columns - raw predicted probability of outcome

# get the right data...
bug <- 'MRSA'
df <- data[data[[bug]] == 1L, ]

# define treatments to intervene on
trts <- c('VAN_DEF', 'DAP_DEF')
trt_groups <- apply(
   X = expand.grid(
      paste0(gsub('_DEF', '', trts[1]), 0:1), 
      paste0(gsub('_DEF', '', trts[2]), 0:1)
   ),
   MARGIN = 1, 
   FUN = paste, collapse='_'
)
trt_assignments <- setNames(expand.grid(0:1, 0:1), trts)

# create empty dataframes to hold: predicted probs (4), pairwise diffs (6) data
x <- rep(NA, nrow(data))
pred_prob_DF <- data.frame(x,x,x,x)
names(pred_prob_DF) <- trt_groups
pair_diffs_DF <- data.frame(
   add1 = x,
   add2 = x,
   swap = x,
   add1_2 = x,
   add2_1 = x,
   add12 = x
)
rm(x)

# loop over each fold:
#     get model
#     loop over each intervention:
#           make predictions - add to pred_prob_DF

folds <- 2017:2023
for (f in seq_along(folds)) {
   wrows <- which(df$year == folds[f])
   model <- h2o.cross_validation_models(xgb)[[f]]

   for (i in seq_len(nrow(trt_assignments))) { # i for intervention
      # intervene on the data
      df[[trts[1]]] <- trt_assignments[[trts[1]]][i]
      df[[trts[2]]] <- trt_assignments[[trts[2]]][i]
      
      # make predictions
      pred_prob_DF[[i]][wrows] <- as.vector(h2o.predict(object=model, newdata=as.h2o(df[wrows,] %>% select(-year)))$p1)
   }
}



trtDF <- data.frame(row.names=seq_len(nrow(df)))

df_mod <- df
for (trt1_status in 0:1) {
   df_mod[[trts[1]]] <- trt1_status
   for (trt2_status in 0:1) {
      df_mod[[trts[2]]] <- trt2_status
      pred <- h2o.predict(object=xgb_cal, newdata=as.h2o(df_mod), )
      col_name <- paste0(trts[1], '_', trt1_status, '__', trts[2], '_', trt2_status)
      cat(col_name, table(df_mod[[trts[1]]], df_mod[[trts[2]]]), '\n')
      trtDF[col_name] <- as.vector(pred$p1)
   }
}


diffs <- vector('list', (length(trtDF) * (length(trtDF) - 1)) / 2)
idx <- 1
for (i in seq_along(trtDF)[-length(trtDF)]) {
   for (j in (i+1):length(trtDF)) {
      diffs[[idx]] <- trtDF[[i]] - trtDF[[j]] # P(d30mort | i) - P(d30mort | j)
      names(diffs)[idx] <- paste(names(trtDF)[c(i,j)], collapse='____')
      idx <- idx + 1
      #hist(diffs, xlim=c(-0.26, 0.26), breaks=seq(-0.26, 0.26, 0.02),
      #     main = paste(names(trtDF)[c(i,j)], collapse='\n'))
      #title(sub=paste0(round(mean(diffs), 4) * 100, '%'))
   }
}

# negative ~ dap - van < 0 ~ dap better
# d <- diffs$MEM_DEF_0__ETP_DEF_1____MEM_DEF_1__ETP_DEF_0
d <- diffs$VAN_DEF_0__DAP_DEF_1____VAN_DEF_1__DAP_DEF_0 # swapping DAP for VAN
#d <- diffs$VAN_DEF_1__DAP_DEF_0____VAN_DEF_1__DAP_DEF_1 # add DAP to VAN
#d <- diffs$VAN_DEF_0__DAP_DEF_1____VAN_DEF_1__DAP_DEF_1 # add VAN to DAP
sum(d == 0) / length(d) * 100 # 69% no difference
sum(d < 0.02 & d > -0.02) / length(d) * 100 # 89.8% within 2% no difference
sum(d < 0.05 & d > -0.05) / length(d) * 100 # 96% within 5% no difference
sum(d <= -0.05) / length(d) * 100 # ~2.93% dap better
sum(d >= 0.05) / length(d) * 100 # ~0.91% van better

# is there consistency to who responds according to the model?
diffsDF <- data.frame(do.call(cbind, diffs))

# yes! 950 had 4 pairs different, 13 had 6 different
table(apply(X = diffsDF, 
            MARGIN = 1,
            FUN = function(x) {
               sum(x != 0)
            }))

# some comparisons are non-zero more often:
data.frame(
   non0 = sapply(X = diffsDF, 
                 FUN = function(x) {
                    sum(x != 0)
                 })
)


# plot X (P(outcome | trt1)) vs. Y (P(outcome | trt2) - P(outcome | trt1))
# ~ P(outcome | reference treatment) vs. predicted effect size of alternative treatment
par(mfrow=c(2,1))
plot(x =                              trtDF$VAN_DEF_1__DAP_DEF_0, # P(outcome | reference treatment)
     y = trtDF$VAN_DEF_0__DAP_DEF_1 - trtDF$VAN_DEF_1__DAP_DEF_0, # P(outcome | alternative) - P(outcome | reference)
     xlab = 'P(outcome | VAN=1, DAP=0)',
     ylab = 'See title...',
     main = 'y-axis: P(outcome | VAN=0, DAP=1) - P(outcome | VAN=1, DAP=0)',
     pch = 16, ylim=c(-0.32, 0.32), col='#00000022')
abline(h = 0)
abline(h = -0.03, lty=2)

plot(x =                              trtDF$MEM_DEF_1__ETP_DEF_0, # P(outcome | reference treatment)
     y = trtDF$MEM_DEF_0__ETP_DEF_1 - trtDF$MEM_DEF_1__ETP_DEF_0, # P(outcome | alternative) - P(outcome | reference)
     xlab = 'P(outcome | VAN=1, DAP=0)',
     ylab = 'See title...',
     main = 'y-axis: P(outcome | VAN=0, DAP=1) - P(outcome | VAN=1, DAP=0)',
     pch = 16, ylim=c(-0.2, 0.2), col='#00000022')


# some patients expected to respond highly to abx swapping (effect sizes ~ 5-10%)
# who are they??
d <- trtDF$VAN_DEF_0__DAP_DEF_1 - trtDF$VAN_DEF_1__DAP_DEF_0
hist(d, main='Predicted effect size of DAP vs. VAN', xlim=c(-0.15, 0.15))
sum(d < -0.03) # 181 with 3% or greater effect size
sum(d < -0.05) # 90 with 5% or greater effect size

w <- which(d < -0.03)

df$DAP_RESPONDER <- FALSE
df$DAP_RESPONDER[w] <- TRUE


df %>%
   summarise(
      n = n(),
      mean(AGE),
      sum(FEMALE) / n,
      sum(ICD_SepticShock_1w) / n,
      sum(ICD_Sepsis_1w) / n,
      .by = DAP_RESPONDER
   )


plot(x=df$AGE, y=d)
fisher.test(table(df$ICD_SepticShock_1w, df$DAP_RESPONDER))














# model explainability
par(mfrow=c(1,1))
h2o.varimp_plot(model = xgb, num_of_features = 20L)

# where are the antibiotic prescription features in the rankings of feature importance?
impDF <- as.data.frame(h2o.varimp(object = xgb))

w_emp <- grep('.+_EMP', impDF$variable)
w_def <- grep('.+_DEF', impDF$variable)

mean(w_emp) / nrow(impDF) # 69%
mean(w_def) / nrow(impDF) # 64%

sum(impDF$percentage[w_emp]) * 100 # 0.67%
sum(impDF$percentage[w_def]) * 100 # 1.8%
sum(impDF$percentage[impDF$variable == 'AGE']) * 100 # 4.9%
sum(impDF$percentage[grep('FACILITY', impDF$variable)]) * 100 # 1.05%

# using h2o.explain() functions
h2o.explain(object = xgb, newdata = as.h2o(data %>% select(!c(year))))
h2o.explain_row(object = xgb, newdata = as.h2o(data %>% select(!c(year))), row_index = 1)



h2o.explain(object = xgb, 
            newdata = as.h2o(
               data %>% 
                  select(!c(year)) %>%
                  filter(MRSA == 1L)
            ),
            top_n_features = 25L,
            include_explanations = c('confusion_matrix', 'varimp', 'shap_summary'),
            plot_overrides = list(
               shap_summary_plot = list(top_n_features = 50L)
            ))

h2o.shap_explain_row_plot(
   model = xgb, 
   newdata = as.h2o(data %>% select(-year) %>% filter(ESBL == 1L)), 
   row_index = 1,
   columns = grep('EMP|DEF', names(data), value=TRUE)
)

h2o.shap_summary_plot(
   model = xgb, 
   newdata = as.h2o(data %>% select(-year) %>% filter(MRSA == 1L)),
   columns = apply(expand.grid(c('VAN', 'LZD', 'CPT'), c('EMP', 'DEF')), 1, paste, collapse='_')
)

h2o.partialPlot(
   object = xgb,
   newdata = as.h2o(data %>% select(-year) %>% filter(MRSA == 1L)),
   cols = apply(expand.grid(c('VAN', 'LZD', 'CPT'), c('EMP', 'DEF')), 1, paste, collapse='_')  
)



generateCF <- function(model, bug, filter_var, switch_vars) {
   df <- data %>% select(-year)
   for (i in seq_along(filter_vars)) {
      df <- df[df[[names(filter_vars)[i]]] == filter_vars[[i]],]
      cat(nrow(df), '')
   }
   pred_og <- h2o.predict(object = model, as.h2o(df))
   
   for (i in switch_vars) {
      switch_var <- filter_vars[i]
      new_value <- as.integer(1 & !unlist(switch_var))
      df[names(switch_var)]
   }
}







# # save model
# model_path <- h2o.saveModel(
#    object = xgb,
#    path = '/Users/samblechman/Desktop/EHR-mining/Scripts/AnalysisScripts/R01/OutcomePredictionModel/models/',
#    export_cross_validation_predictions = TRUE, 
#    force = TRUE,
#    filename = paste0('xgb_all_', target)
# )
# 
# # load model and evaluation functions
# xgb <- h2o.loadModel(path = '/Users/samblechman/Desktop/EHR-mining/Scripts/AnalysisScripts/R01/OutcomePredictionModel/models/xgb_all_d30mortality')















