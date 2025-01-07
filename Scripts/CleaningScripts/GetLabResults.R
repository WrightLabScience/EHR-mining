source('~/Desktop/EHR/EHR work/config_file.R')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/MRSA_bacteremia/MRSA_BSI_DAP_VAN.Rdata')
person_IDs <- unique(dfx$PERSON_ID)
chunks <- list(1:832, 833:1664, 1665:2496)
labs <- tibble()
for (chunk in chunks) {
   print(chunk[1])
   for (year in 2017:2023) {
      print(year)
      labs <- rbind(labs,
                    tbl(conn, in_schema('AMB_ETL', paste0('SENS_LAB_RESULT_', year, '_VW'))) %>% 
                       filter(PERSON_ID %in% local(person_IDs[chunk])) %>%
                       collect())
   }
}
save(labs, file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/MRSA_bacteremia/LabResults_raw.Rdata')
