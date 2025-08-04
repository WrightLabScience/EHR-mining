# This script gets the raw data for a cohort of patients and timeframe
# As of now, it grabs patients with any positive blood culture between 2017-2023
# Gets raw labs and cleans them
# Gets raw encounters
# Gets raw ICDs

# TODO: should this exclude non-priority/contaminant/fungal pathogens?


# Get person_ids for 2017-2023 blood cultures
source(file = '~/Desktop/EHR/EHR work/config_file.R')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_clean_ASTs.Rdata')
Rdata_file_path <- '~/Desktop/EHR/EHR work/RdataFiles/BuildCohortDatasets/RawData/'
years <- 2017:2023
bacteremia_person_ids <- astDF %>%
   filter(BLOOD) %>%
   filter(lubridate::year(ORDER_DAY) %in% years) %>%
   select(PERSON_ID) %>%
   distinct() %>%
   pull(PERSON_ID) # 45,362
chunks <- mapply(seq(1, floor(length(bacteremia_person_ids) / 1000) * 1000 + 1, 1000),
                 c(seq(1000, floor(length(bacteremia_person_ids) / 1000) * 1000, 1000), length(bacteremia_person_ids)),
                 FUN = ':') # length = 46


#### START LABS ####
# Get raw lab results from Neptune
labsDF <- tibble()
for (year in years) {
   cat(year, '- ')
   for (chunk in chunks) {
      cat(chunk[1], '')
      labsDF <- rbind(
         labsDF,
         tbl(conn, in_schema('AMB_ETL', paste0('SENS_LAB_RESULT_', year, '_VW'))) %>%
            filter(PERSON_ID %in% local(bacteremia_person_ids[chunk])) %>%
            collect()
      )
   }
   cat('\n')
}

# Clean date-time columns and arrange
labsDF <- labsDF %>%
   mutate(
      across(.cols = c(ORDER_DATE, RESULT_DATE),
             .fns = ~ strptime(., format='%m/%d/%Y %T'))
   ) %>%
   arrange(PERSON_ID, ORDER_DATE)

# Save raw labs dataframe
save(labsDF, file = paste0(Rdata_file_path, 'labsDF_raw_2017_2023_bacteremia.Rdata'))

# Clean whole labs dataframe
source('~/Desktop/EHR-mining/Scripts/CleaningScripts/CleanLabResults.R')
labsDF <- cleanLabResults(labsDF)

# Save cleaned labs dataframe
save(labsDF, file = paste0(Rdata_file_path, 'labsDF_clean_2017_2023_bacteremia.Rdata'))

# remove unnecessary variables
rm(getRawLabsData, cleanLabResults, labsDF)
#### END LABS ####


#### START ENCOUNTERS ####
# Get data from Neptune
encsDF <- tibble()
for (chunk in chunks) {
   cat(chunk[1], '')
   encsDF <- rbind(
      encsDF,
      tbl(conn, in_schema('AMB_ETL', 'SENS_ENCOUNTER_VW')) %>%
         filter(PERSON_ID %in% local(bacteremia_person_ids[chunk]),
                substr(ADMIT_DATE, 7, 10) %in% years) %>%
         collect()
   )
}

# Clean up date-time columns and arrange
encsDF <- encsDF %>%
   mutate(
      across(.cols = c(ADMIT_DATE, DISCHARGE_DATE),
             .fns = ~ strptime(., format='%m/%d/%Y %T'))
   ) %>%
   mutate(
      ADMIT_DAY = lubridate::as_date(ADMIT_DATE),
      DISCHARGE_DAY = lubridate::as_date(DISCHARGE_DATE)
   ) %>%
   arrange(PERSON_ID, ADMIT_DATE, DISCHARGE_DATE)

# Save raw encounters data
save(encsDF, file = paste0(Rdata_file_path, 'encsDF_raw_2017_2023_bacteremia.Rdata'))


site_groups <- list(
   academic = c("Mercy", "Presbyterian", "Shadyside", "MageeW"), 
   regional = c("Hamot", "Williamsport", "Jameson", "Altoona"), 
   community = c("Passavant", "East", "McKeesport", "StMargaret", 'SOL', 'LOC'), 
   rural = c("Bedford", "Northwest", "Horizon", "Chatauqua")
)
save(site_groups, '~/Desktop/EHR-mining/UsefulDataForCleaning/UPMC_site_groups.Rdata')
#### END ENCOUNTERS ####


#### START ICDS ####
icdsDF <- tibble()
for (chunk in chunks) {
   cat(chunk[1], '')
   icdsDF <- rbind(
      icdsDF,
      tbl(conn, in_schema('AMB_ETL', 'LAB_SENS_DX_VW')) %>%
         filter(PERSON_ID %in% local(bacteremia_person_ids[chunk]),
                lubridate::year(DX_FROM_DATE) %in% (min(years)-2L):max(years)) %>%
         collect()
   )
}

# Clean and arrange data
icdsDF <- icdsDF %>%
   mutate(
      across(.cols = c(CODE_TYPE, CODE_DESCRIPTION),
             .fns = ~ ifelse(is.na(.), 'unknown', .)), # replace NA code descriptions with character string (so doesn't behave as NA)
      DX_DATE = lubridate::as_date(DX_FROM_DATE)
   ) %>% 
   arrange(PERSON_ID, DX_FROM_DATE, PRIMARY_DX_IND, CODE_DESCRIPTION)

# Save raw ICDs
save(icdsDF, file=paste0(Rdata_file_path, 'icdsDF_raw_2017_2023_bacteremia'))
#### END ICDS ####













