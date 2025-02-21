getBSIrecurrence <- function(df) {
   load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_clean_ASTs.Rdata')
   load(file = '~/Desktop/EHR-mining/UsefulDataForCleaning/bug_groups.Rdata')
   astDF <- astDF %>%
      filter(
         BLOOD,
         PERSON_ID %in% unique(df$PERSON_ID),
         BUG %in% unname(unlist(bug_groups))
      )
   w <- which(astDF$BUG == 'Staphylococcus aureus' & is.na(astDF$OXACILLIN))
   if (length(w) > 0L)
      astDF <- astDF[-w,]
   
   # join with index cases
   df <- df %>%
      select(PERSON_ID, ORDER_DAY) %>%
      distinct() %>%
      mutate(AFTER_CULTURE = ORDER_DAY + 7L) %>%
      left_join(
         x = .,
         y = astDF %>% select(PERSON_ID, ORDER_DAY) %>% distinct() %>% rename(NEXT_CULTURE_DAY = ORDER_DAY),
         by = join_by(
            PERSON_ID,
            closest(y$NEXT_CULTURE_DAY >= x$AFTER_CULTURE)
         )
      ) %>%
      mutate(BSIrecur_time = as.integer(NEXT_CULTURE_DAY - ORDER_DAY)) %>%
      select(PERSON_ID, ORDER_DAY, BSIrecur_time) %>%
      distinct()
   
   source('~/Desktop/EHR-mining/Scripts/CleaningScripts/featurizeDemographics.R')
   df <- XdayOutcome(df, col_name='BSIrecur_time')
   
   return(df)
}
