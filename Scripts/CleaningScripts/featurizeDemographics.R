getDemographics <- function(df) {
   load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_DEMO.Rdata')
   df <- df %>%
      select(PERSON_ID, ORDER_DAY) %>%
      distinct() %>%
      left_join(
         x = .,
         y = dth,
         by = join_by(PERSON_ID)
      ) %>%
      mutate(
         AGE = as.integer(ORDER_DAY - DOB) / 365,
         mortality_time = as.integer(DEATH_DATE - ORDER_DAY),
         FEMALE = as.integer(GENDER == 'FEMALE')
      ) %>%
      select(PERSON_ID, ORDER_DAY, FEMALE, AGE, mortality_time) %>%
   return(df)
}
