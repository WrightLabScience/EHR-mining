getRawEncountersData <- function(person_ids, years=2023) {
   chunks <- mapply(seq(1, floor(length(person_ids) / 1000) * 1000 + 1, 1000),
                    c(seq(1000, floor(length(person_ids) / 1000) * 1000, 1000), length(person_ids)),
                    FUN = ':')
   years <- (min(years) - 1):max(years)
   source(file = '~/Desktop/EHR/EHR work/config_file.R')
   
   cat('Getting raw data from database.\n')
   cat('For', length(chunks), 'chunks of 1000 patients each.\n')
   
   encsDF <- tibble()
   for (chunk in chunks) {
      cat(chunk[1], '')
      encsDF <- rbind(
         encsDF,
         tbl(conn, in_schema('AMB_ETL', 'SENS_ENCOUNTER_VW')) %>%
            filter(PERSON_ID %in% local(person_ids[chunk]),
                   substr(ADMIT_DATE, 7, 10) %in% years) %>%
            collect()
      )
   }
   
   cat('Processing data.\n')
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
   
   return(encsDF)
}
