source('~/Desktop/EHR/EHR work/config_file.R')


meds <- c('PHENYLEPHRINE', 'NOREPINEPHRINE', 'EPINEPHRINE', 'VASOPRESSIN', 'DOPAMINE', 
          'SEPSIS POWERPLAN ACTIVE 1ST DOSE ABX ARE STAT', 'SEPSIS POWERPLAN ACTIVE 1ST DOSE ABX ARE')

medsDF <- tibble()
for (view in c('', '_IV')) {
   for (med in meds) {
      query <- paste0("select * from AMB_ETL.SENS_MED_ADMIN", view, "_VW where MEDICATION = '", med, "'")
      cat(query, '\n')
      
      medsDF <- rbind(medsDF,
                      dbGetQuery(conn, query))
   }
}

medsDF <- tibble(medsDF)

save(medsDF, file = '~/Desktop/EHR/EHR work/RdataFiles/SENS_raw_vasopressor_meds.Rdata')

medsDF <- medsDF %>%
   mutate(PERSON_ID = as.character(PERSON_ID)) %>%
   rename(START_DATE = ADMIN_START_DATE) %>%
   select(PERSON_ID, START_DATE, MEDICATION) %>%
   distinct() %>%
   mutate(START_DAY = lubridate::as_date(START_DATE))

save(medsDF, file = '~/Desktop/EHR/EHR work/RdataFiles/SENS_clean_vasopressor_meds.Rdata')
