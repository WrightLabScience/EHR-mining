featurizeVasopressors <- function(df, vasDF, preprocess=TRUE, time_frame) {
   start_time_frame <- time_frame[1]
   end_time_frame <- time_frame[2]
   
   if (preprocess) {
      vasDF <- vasDF %>%
         filter(PERSON_ID %in% unique(df$PERSON_ID)) %>%
         rename(START_DATE = ADMIN_START_DATE) %>%
         mutate(PERSON_ID = as.character(PERSON_ID)) %>%
         select(PERSON_ID, START_DATE) %>%
         distinct()
   }
   
   df_og <- df
   df <- df_og %>% select(PERSON_ID, ORDER_DATE)
   
   df <- df %>%
      mutate(
         JOIN_START = ORDER_DATE + start_time_frame * 86400,
         JOIN_END = ORDER_DATE + end_time_frame * 86400
      ) %>%
      left_join(
         x = .,
         y = vasDF,
         by = join_by(
            PERSON_ID,
            between(y$START_DATE, x$JOIN_START, x$JOIN_END)
         )
      ) %>%
      summarise(
         ENC_vasopressors = sum(!is.na(START_DATE)),
         .by = c(PERSON_ID, ORDER_DATE)
      )
   
   df <- df %>%
      left_join(
         x = .,
         y = df_og,
         by = join_by(PERSON_ID, ORDER_DATE)
      )
   
   return(df)
}