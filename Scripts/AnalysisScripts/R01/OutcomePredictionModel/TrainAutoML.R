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
labels <- as.integer(data[[target]]) - 1


# install.packages("h2o", type="source", repos=(c("http://h2o-release.s3.amazonaws.com/h2o/latest_stable_R")))
library(h2o)

h2o.init()

aml <- h2o.automl(
   y = target,
   training_frame = as.h2o(data),
   max_runtime_secs = (60/60) * 3600,
   fold_column = 'year',
   nfolds = 0,
   keep_cross_validation_predictions = TRUE,
   seed = 123L
)
# cross-validation? Default explained in "fold_column" argument description:
#     fold_column:   Specifies a column with cross-validation fold index assignment per observation. 
#                    This is used to override the default, randomized, 5-fold cross-validation scheme 
#                    for individual models in the AutoML run.

h2o.saveModel(aml, export_cross_validation_predictions = TRUE, force = TRUE,
              path = 'Users/samblechman/Desktop/EHR-mining/Scripts/AnalysisScripts/R01/OutcomePredictionModel/models/',
              filename = 'automl_all')

save(aml, file = '~/Desktop/EHR-mining/Scripts/AnalysisScripts/R01/OutcomePredictionModel/AutoML_model_object_d30mortality_allmodels1hr_foldyear.Rdata')
# save(aml, file = '~/Desktop/EHR-mining/Scripts/AnalysisScripts/R01/OutcomePredictionModel/AutoML_model_object_d90mortality_allmodels1hr_foldyear.Rdata')
# save(aml, file = '~/Desktop/EHR-mining/Scripts/AnalysisScripts/R01/OutcomePredictionModel/AutoML_model_object_d90event_allmodels1hr_foldyear.Rdata')


lb <- as.data.frame(h2o.get_leaderboard(aml, extra_columns="all"))
model_id_1st_non_ensemble <- lb$model_id[which(lb$algo != 'StackedEnsemble')[1]]
model_1st_non_ensemble <- h2o.getModel(model_id_1st_non_ensemble)


source('~/Desktop/EHR-mining/Scripts/CleaningScripts/evaluate_model_outputs_and_plot.R')


pred_ens <- h2o.predict(object = aml@leader, 
                        newdata = as.h2o(data))
evalDF_ens <- evalModel(probs = as.vector(pred_ens$p1),
                        labels = as.integer(data[[target]]) - 1)

pred <- h2o.predict(object = model_1st_non_ensemble, 
                    newdata = as.h2o(data))
evalDF <- evalModel(probs = as.vector(pred$p1),
                    labels = as.integer(data[[target]]) - 1)
# plotROC_PRC(evalDF_train = evalDF)

pred_cv <- h2o.cross_validation_holdout_predictions(object = model_1st_non_ensemble)
evalDF_cv <- evalModel(probs = as.vector(pred_cv$p1),
                       labels = as.integer(data[[target]]) - 1)

plotROC_PRC(evalDF_train = evalDF,
            evalDF_test = evalDF_cv)






# regression comparison
plot(x = trainDF[[target]],
     y = as.vector(pred_train$predict),
     xlim = c(0, 100),
     ylim = c(0, 100),
     pch = 16,
     cex = 0.5)

plot(x = testDF[[target]],
     y = as.vector(pred_test$predict),
     xlim = c(0, 100),
     ylim = c(0, 100),
     pch = 16,
     cex = 0.5)


# variable importance
impDF <- h2o.varimp(aml@leader)
h2o.varimp_plot(aml@leader, num_of_features = 20L)
barplot(impDF$percentage)






vanDF <- trainDF %>% 
   filter(MRSA == 1L, VAN_EMP == 1L, DAP_EMP == 0L)

pred_og <- h2o.predict(aml, as.h2o(vanDF))

vanDF$VAN_EMP <- 0L
vanDF$DAP_EMP <- 1L

# vanDF$MRSA <- 1L

pred_new <- h2o.predict(aml, as.h2o(vanDF))

diffs <- as.vector(pred_new$p1) - as.vector(pred_og$p1)
mean(diffs)
hist(diffs, breaks=20, xlim=c(-.1, .1))










# get and combine cross valiation results
cv_res <- h2o.cross_validation_predictions(aml@leader)
folds <- as.vector(h2o.cross_validation_fold_assignment(aml@leader))
cvDF <- data.frame(
   predict = rep(NA_real_, nrow(trainDF)),
   p0 = rep(NA_real_, nrow(trainDF)),
   p1 = rep(NA_real_, nrow(trainDF))
)
for (i in seq_along(cv_res)) {
   w <- which(folds == (i - 1))
   cvDF[w,] <- as.numeric(cv_res[[i]][w,])
}
rm(i, w, cv_res, folds)



