library(dplyr)
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

id_df <- featuresDF %>% select(PERSON_ID, ORDER_DAY) %>% arrange(PERSON_ID, ORDER_DAY) %>% mutate(ID = row_number())

save(id_df, file = '~/Desktop/EHR/EHR work/RdataFiles/R01/identification_numbers.Rdata')
