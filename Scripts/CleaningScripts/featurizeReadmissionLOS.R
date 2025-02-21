featurizeReadmissionLOS <- function(df, encsDF) {
   # length of stay
   source('~/Desktop/EHR-mining/Scripts/CleaningScripts/CombineOverlappingEncountersFxn.R')
   encsDFc <- combineOverlappingEncounters(encsDF)
   losDF <- df %>%
      select(PERSON_ID, ORDER_DAY) %>%
      distinct() %>%
      mutate(WEEK_AFTER_CULTURE = ORDER_DAY + 7L) %>%
      left_join(
         x = .,
         y = encsDFc,
         by = join_by(
            PERSON_ID,
            overlaps(x_lower = ORDER_DAY, x_upper = WEEK_AFTER_CULTURE,
                     y_lower = ADMIT_DAY, y_upper = DISCHARGE_DAY)
         )
      ) %>%
      group_by(PERSON_ID, ORDER_DAY) %>%
      mutate(TOTAL_LOS = sum(LENGTH_OF_STAY)) %>%
      ungroup() %>%
      select(PERSON_ID, ORDER_DAY, TOTAL_LOS) %>%
      distinct()
   
   # time to readmission (after 14 days?)
   readDF <- df %>%
      select(PERSON_ID, ORDER_DAY) %>%
      distinct() %>%
      mutate(AFTER_CULTURE = ORDER_DAY + 14L) %>%
      left_join(
         x = .,
         y = encsDF,
         by = join_by(
            PERSON_ID,
            closest(y$ADMIT_DAY >= x$AFTER_CULTURE)
         )
      ) %>%
      mutate(readmit_time = as.integer(ADMIT_DAY - ORDER_DAY)) %>%
      select(PERSON_ID, ORDER_DAY, readmit_time) %>%
      distinct()
   readDF <- XdayOutcome(readDF, col_name='readmit_time')
   
   df <- inner_join(
      x = losDF,
      y = readDF,
      by = join_by(PERSON_ID, ORDER_DAY)
   )
   
   return(df)
}
