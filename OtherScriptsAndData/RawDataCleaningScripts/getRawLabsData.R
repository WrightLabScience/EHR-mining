getRawLabsData <- function(person_ids, years=2023L) {
   chunks <- mapply(seq(1, floor(length(person_ids) / 1000) * 1000 + 1, 1000),
                    c(seq(1000, floor(length(person_ids) / 1000) * 1000, 1000), length(person_ids)),
                    FUN = ':')
   source(file = '~/Desktop/EHR/EHR work/config_file.R')
   years <- intersect(2017:2023, years)
   
   cat('Getting raw labs data from database.\n')
   cat('Getting years: ', years, '\n')
   cat('For', length(chunks), 'chunks of 1000 patients each.\n')
   labsDF <- tibble()
   for (year in years) {
      cat(year, '- ')
      for (chunk in chunks) {
         cat(chunk[1], '')
         labsDF <- rbind(
            labsDF,
            tbl(conn, in_schema('AMB_ETL', paste0('SENS_LAB_RESULT_', year, '_VW'))) %>%
               filter(PERSON_ID %in% local(person_ids[chunk])) %>%
               collect()
         )
      }
      cat('\n')
   }
   
   cat('Final processing.')
   labsDF <- labsDF %>%
      mutate(
         across(.cols = c(ORDER_DATE, RESULT_DATE),
                .fns = ~ strptime(., format='%m/%d/%Y %T'))
      ) %>%
      arrange(PERSON_ID, ORDER_DATE)
   
   return(labsDF)
}
