getRawICDsData <- function(person_ids, years=2023L) {
   years <- (min(years) - 2):max(years)
   chunks <- mapply(seq(1, floor(length(person_ids) / 1000) * 1000 + 1, 1000),
                    c(seq(1000, floor(length(person_ids) / 1000) * 1000, 1000), length(person_ids)),
                    FUN = ':')
   source(file = '~/Desktop/EHR/EHR work/config_file.R')
   
   cat('Getting raw ICDs data from database.\n')
   cat('Getting years: ', years, '\n')
   cat('For', length(chunks), 'chunks of 1000 patients each.\n')
   icdsDF <- tibble()
   for (chunk in chunks) {
      cat(chunk[1], '')
      icdsDF <- rbind(
         icdsDF,
         tbl(conn, in_schema('AMB_ETL', 'LAB_SENS_DX_VW')) %>%
            filter(PERSON_ID %in% local(person_ids[chunk]),
                   lubridate::year(DX_FROM_DATE) %in% years) %>%
            collect()
      )
   }
   
   cat('Final processing.\n')
   icdsDF <- icdsDF %>% 
      arrange(PERSON_ID, DX_FROM_DATE, PRIMARY_DX_IND, CODE_DESCRIPTION) %>%
      mutate(
         across(.cols = c(CODE_TYPE, CODE_DESCRIPTION),
                .fns = ~ ifelse(is.na(.), 'unknown', .)),
         DX_DATE = lubridate::as_date(DX_FROM_DATE)
      )
   
   return(icdsDF)
}
