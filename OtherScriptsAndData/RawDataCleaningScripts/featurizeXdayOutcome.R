XdayOutcome <- function(df, col_name='OUT_mortality_time', threshold=c(14L, 21L, 30L, 60L, 90L), index_day='ORDER_DAY') {
   stopifnot('provided time column must be named as `OUT_()_time`' = grepl('^OUT_.+_time$', col_name))
   
   og_vars <- names(df)
   for (x in threshold) {
      df[[paste0('OUT_', gsub('^OUT_(.+)_time$', '\\1', col_name), '_', x, 'd')]] <- ifelse(test = is.na(df[[col_name]]) | df[[col_name]] > x, 
                                                                                            yes = 0L,
                                                                                            no = 1L)
   }
   new_vars <- names(df)[!names(df) %in% og_vars]
   df <- df %>% relocate(!!new_vars, .after=!!col_name)
   return(df)
}
