featurizeLabs <- function(df, labsDF, time_frame=c(-4L, 4L)) {
   df <- df %>%
      select(PERSON_ID, ORDER_DAY) %>%
      distinct() %>%
      rename(CULTURE_DAY = ORDER_DAY) %>%
      mutate(
         BEFORE_CULTURE_DAY = CULTURE_DAY + time_frame[1],
         AFTER_CULTURE_DAY = CULTURE_DAY + time_frame[2]
      ) %>%
      left_join(
         x = .,
         y = labsDF,
         by = join_by(
            PERSON_ID,
            between(y$ORDER_DAY, x$BEFORE_CULTURE_DAY, x$AFTER_CULTURE_DAY)
         )
      ) %>%
      mutate(DAYS = as.integer(ORDER_DAY - CULTURE_DAY)) %>%
      select(!c(BEFORE_CULTURE_DAY, AFTER_CULTURE_DAY, ORDER_DATE, RESULT_DATE, ORDER_DAY, CLEAN_REF_UNIT, RESULT_LAB_NAME)) %>%
      rename(ORDER_DAY = CULTURE_DAY)
   
   df <- df %>%
      summarise(
         RESULT_VALUE = mean(RESULT_VALUE),
         .by = c(PERSON_ID, ORDER_DAY, LAB)
      ) %>%
      tidyr::pivot_wider(
         values_from = RESULT_VALUE,
         names_from = LAB
      ) %>%
      select(-`NA`)
   
   names(df)[!names(df) %in% c('PERSON_ID', 'ORDER_DAY')] <- paste0('LAB_', names(df)[!names(df) %in% c('PERSON_ID', 'ORDER_DAY')])
   
   return(df)
}