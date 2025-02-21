library(dplyr)
Rdata_file_path <- '/Users/samblechman/Desktop/EHR/EHR work/RdataFiles/R01/OutcomePredictionModel/'
load(file = paste0(Rdata_file_path, 'data_feat_outcomes_processed_all.Rdata'))
data_og <- data
target <- 'd30mortality'
data <- data_og %>%
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

set.seed(123L)
w <- sample(seq_len(nrow(data)), size = 0.8 * nrow(data))
train <- data[w,]
test <- data[-w,] %>% select(-year)
labels_train <- as.integer(train[[target]]) - 1L
labels_test <- as.integer(test[[target]]) - 1L



library(h2o)
h2o.init()
h2o.no_progress()


# Perform grid search
grid <- h2o.grid(
   algorithm = 'xgboost',
   y = target,
   training_frame = as.h2o(train),
   seed = 123L,
   fold_column = 'year',
   nfolds = 0,
   search_criteria = list(
      strategy = 'RandomDiscrete',
      max_runtime_secs = 3600
   ),
   hyper_params = list(
      ntrees = c(20, 50, 200),            # default: 50
      max_depth = c(2, 6, 15, 0),          # default: 6 - 0 is unlimited
      learn_rate = c(0.05, 0.3),       # default: 0.3 - 0.0 to 1.0
      sample_rate = c(0.7, 1),        # default: 1 - 0.0 to 1.0 - also called subsample
      colsample_bytree = c(0.7, 1),   # default: 1 - 0.0 to 1.0 - also called col_sample_rate_per_tree
      colsample_bylevel = c(0.7, 1),  # default: 1 - 0.0 to 1.0 - also called col_sample_rate
      reg_lambda = c(1, 10),   # L2 regularization on weights
      reg_alpha = c(0, 10) 
   )
)



df <- as.data.frame(grid@summary_table)

for (i in 1:8) {
   print(names(df)[i])
   print(tapply(X = df$logloss,
                INDEX = df[[i]],
                FUN = function(x) mean(x)))
}

colsample_bylevel <- 0.7
colsample_bytree <- 1
learn_rate <- 0.3
max_depth <- 15 # (2 did barely worse)
ntrees <- 20 # (not huge difference)
reg_alpha <- 10
reg_lambda <- 1
sample_rate <- 0.7


model_ids <- df$model_ids

res <- sapply(X = model_ids,
              FUN = function(id) {
                 model <- h2o.getModel(id)
                 return(c(
                    'auroc' = h2o.auc(object = model, train=TRUE, xval=TRUE),
                    'auprc' = h2o.aucpr(object = model, train=TRUE, xval=TRUE)
                 ))
              })

df$train_auroc <- res['auroc.train',]
df$xval_auroc <- res['auroc.xval',]
df$train_auprc <- res['auprc.train',]
df$xval_auprc <- res['auprc.xval',]



test_res <- sapply(X = model_ids,
                   FUN = function(id) {
                      model <- h2o.getModel(id)
                      return(evalModel_h2o(model=model, dataset=test, labels=labels_test, plot_data=FALSE)$data)
                   })
df$test_auroc <- test_res['AUROC',]
df$test_auprc <- test_res['AUPRC',]



par(mfrow=c(2,1))
plot(NA, xlim=c(0.5, 1), ylim=c(0.5, 1), xlab='Train AUROC', ylab='Cross-val AUROC')
text(x=df$train_auroc, y=df$xval_auroc, labels=gsub('.+_([0-9]+)', '\\1', df$model_ids))

plot(NA, xlim=c(0, 1), ylim=c(0, 1), xlab='Train AUPRC', ylab='Cross-val AUPRC')
text(x=df$train_auprc, y=df$xval_auprc, labels=gsub('.+_([0-9]+)', '\\1', df$model_ids))


plot(NA, xlim=c(0.5, 1), ylim=c(0.5, 1), xlab='Train AUROC', ylab='Test AUROC')
text(x=df$train_auroc, y=df$test_auroc, labels=gsub('.+_([0-9]+)', '\\1', df$model_ids))

plot(NA, xlim=c(0, 1), ylim=c(0, 1), xlab='Train AUPRC', ylab='Test AUPRC')
text(x=df$train_auprc, y=df$test_auprc, labels=gsub('.+_([0-9]+)', '\\1', df$model_ids))


plot(NA, xlim=c(0.5, 1), ylim=c(0.5, 1), xlab='Cross-val AUROC', ylab='Test AUROC')
abline(a=0, b=1, lty=3)
text(x=df$xval_auroc, y=df$test_auroc, labels=gsub('.+_([0-9]+)', '\\1', df$model_ids))

plot(NA, xlim=c(0, 1), ylim=c(0, 1), xlab='Cross-val AUPRC', ylab='Test AUPRC')
abline(a=0, b=1, lty=3)
text(x=df$xval_auprc, y=df$test_auprc, labels=gsub('.+_([0-9]+)', '\\1', df$model_ids))


# using year as cross-validation folds means that performance in cross validation
# and test set are nearly identical





source('~/Desktop/EHR-mining/Scripts/CleaningScripts/evaluate_model_outputs_and_plot.R')


evalDF <- evalModel_h2o_all(
   model = h2o.getModel(model_ids[1]),
   datasets = list(
      train = train,
      test = test
   ),
   labels = list(
      train = labels_train,
      test = labels_test
   ),
   xval = FALSE
)

plotROC_PRC(evalDF)





