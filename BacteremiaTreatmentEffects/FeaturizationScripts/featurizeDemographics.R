featurizeDemographics <- function(df, demoDF, earliest_death_date=4L) {
   # Pre-process demographic data
   demoDF <- demoDF %>%
      mutate(
         across(
            .cols = c(DOB, DEATH_DATE),
            .fns = ~ as.Date(., format='%m/%d/%Y')
         )
      )
   
   # Execute join and create features.
   df <- df %>%
      left_join(
         x = .,
         y = demoDF,
         by = join_by(PERSON_ID)
      ) %>%
      mutate(
         DEM_age = as.integer(ORDER_DAY - DOB) / 365,
         mortality_time = as.integer(DEATH_DATE - ORDER_DAY),
         DEM_female = GENDER == 'FEMALE'
      ) %>%
      select(-DOB, -GENDER, -PATIENT_STATUS, -DEATH_DATE)
   
   df <- df %>% 
      filter(is.na(mortality_time) | mortality_time >= earliest_death_date) %>%
      rename(OUT_mortality_time = mortality_time)
   cat('Removed patients that died within ', earliest_death_date, ' days: ', nrow(df), '\n', sep='')
   
   # Report and remove children.
   df <- df %>% filter(DEM_age >= 18)
   cat('Removed patients younger than 18 years of age:', nrow(df), '\n')
   
   return(df)
}
