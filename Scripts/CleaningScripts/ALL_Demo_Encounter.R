source(file = '~/Desktop/EHR/EHR work/config_file.R')

dth <- tbl(conn, in_schema('AMB_ETL', 'SENS_PATIENT_DEMO_VW')) %>%
   collect()
length(unique(dth$PERSON_ID)) # 374,305
dth <- dth %>%
   arrange(PERSON_ID) %>%
   mutate(across(c(DOB, DEATH_DATE), ~ as.Date(., format='%m/%d/%Y')))

save(dth, file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_DEMO.Rdata')

enc <- tbl(conn, in_schema('AMB_ETL', 'SENS_ENCOUNTER_VW')) %>%
   collect()

enc %>% count(substr(ADMIT_DATE,7,10))
enc %>% count(substr(DISCHARGE_DATE,7,10))
length(unique(enc$PERSON_ID)) # 210,153
enc %>% count(ADMIT_SOURCE, sort=TRUE)
enc %>% count(ENCOUNTER_TYPE)
enc %>% count(substr(ADMIT_DATE,12,19) == '00:00:00')
enc %>% count(substr(DISCHARGE_DATE,12,19) == '00:00:00')




































