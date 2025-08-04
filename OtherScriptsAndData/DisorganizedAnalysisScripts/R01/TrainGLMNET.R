library(dplyr)
library(glmnet)
source('~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/CreateNamedAbxAbbreviations.R')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/R01/bact_outcomes_variables.Rdata')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/R01/bact_features_variables.Rdata')
featuresDF <- full_join(
   x = outcomesDF %>% select(PERSON_ID, ORDER_DAY),
   y = featuresDF,
   by = join_by(PERSON_ID, ORDER_DAY)
) %>%
   mutate(
      across(.cols = !c(PERSON_ID, ORDER_DAY),
             .fns = ~ ifelse(is.na(.), 0L, .)),
   )


for (target_var in outcomesDF %>% select(-PERSON_ID, -ORDER_DAY) %>% names) {
   cat(target_var, '\n')
   df <- featuresDF
   df[target_var] <- outcomesDF[target_var]
   X <- df %>% select(-PERSON_ID, -ORDER_DAY, -!!target_var) %>% as.matrix()
   labels <- df[[target_var]]
   if (length(unique(labels)) != 2L) next
   cv_fit <- cv.glmnet(
      x = X,
      y = labels,
      family = 'binomial',
      nfolds = 5,
      keep = TRUE
   )
   cv_fit$labels <- labels
   save(cv_fit, file = paste0('~/Desktop/EHR/EHR work/RdataFiles/R01/glmnet_fits/', target_var, '_fit.Rdata'))
}


















