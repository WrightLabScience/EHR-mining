featurizeReadmissionLOS <- function(df, encsDF, index_var='ORDER', var_location='treatment') {
   index_date <- paste0(index_var, '_DATE')
   index_day <- paste0(index_var, '_DAY')
   cat(index_date, index_day, 'as indexes.\n')
   
   source('~/Desktop/EHR-mining/Scripts/CleaningScripts/CombineOverlappingEncountersFxn.R')
   encsDFc <- combineOverlappingEncounters(encsDF)
   
   # Get total length of stay in (mostly) overlapping encounters
   df <- df %>%
      mutate(WEEK_AFTER_INDEX = get(index_day) + 7L) %>%
      left_join(
         x = .,
         y = encsDFc,
         by = join_by(
            PERSON_ID,
            overlaps(x_lower = !!index_day, x_upper = WEEK_AFTER_INDEX,
                     y_lower = ADMIT_DAY, y_upper = DISCHARGE_DAY)
         )
      ) %>%
      group_by(PERSON_ID, !!index_day) %>%
      mutate(
         OUT_length_of_stay = sum(LENGTH_OF_STAY),
         FINAL_DISCHARGE_DATE = max(DISCHARGE_DATE),
         OUT_los_after_index = as.numeric(lubridate::as.duration(FINAL_DISCHARGE_DATE - get(index_date))) / 86400
      ) %>%
      ungroup() %>%
      relocate(OUT_length_of_stay, OUT_los_after_index, .after=!!var_location) %>%
      select(!c(WEEK_AFTER_INDEX, ADMIT_DATE, DISCHARGE_DATE, ADMIT_DAY, DISCHARGE_DAY, FINAL_DISCHARGE_DATE,
                COMBINED_ADJ, ENCOUNTER_TYPE, FACILITY, ADMIT_SOURCE, ADMIT_TYPE, LENGTH_OF_STAY)) %>%
      distinct()
   
   cat(round(sum(is.na(df$OUT_los_after_index)) / nrow(df) * 100 , 2), '% of units did not match to a visit to calculate LOS. Removing them.\n', sep='')
   
   
   # get time to readmission (since time of blood culture...)
   # time to readmission (after 14 days?)
   # readDF <- df %>%
   #    select(PERSON_ID, ORDER_DAY) %>%
   #    distinct() %>%
   #    mutate(AFTER_CULTURE = ORDER_DAY + 14L) %>%
   #    left_join(
   #       x = .,
   #       y = encsDF,
   #       by = join_by(
   #          PERSON_ID,
   #          closest(y$ADMIT_DAY >= x$AFTER_CULTURE)
   #       )
   #    ) %>%
   #    mutate(readmit_time = as.integer(ADMIT_DAY - ORDER_DAY)) %>%
   #    select(PERSON_ID, ORDER_DAY, readmit_time) %>%
   #    distinct()
   # readDF <- XdayOutcome(readDF, col_name='readmit_time')
   # 
   # df <- inner_join(
   #    x = losDF,
   #    y = readDF,
   #    by = join_by(PERSON_ID, ORDER_DAY)
   # )
   
   return(df)
}
