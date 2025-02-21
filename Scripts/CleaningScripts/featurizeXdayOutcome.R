XdayOutcome <- function(df, col_name='mortality_time', threshold=c(14, 21, 30, 60, 90), includes_pid=TRUE) {
   if (includes_pid)
      df <- df %>%
         select(PERSON_ID, ORDER_DAY, !!col_name) %>%
         distinct()
   for (x in threshold) {
      df[[paste0('d', x, gsub('^(.+)_time$', '\\1', col_name))]] <- ifelse(test = is.na(df[[col_name]]) | df[[col_name]] > x, 
                                                                           yes = 0L,
                                                                           no = 1L)
   }
   return(df)
}
