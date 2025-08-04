featurizeICUstays <- function(df, icusDF, preprocess=TRUE, time_frame=c(-7L, 2L)) {
   
   # hs <- tbl(conn, in_schema('AMB_ETL', 'SENS_ENCOUNTER_ICU_LOCS_VW')) %>%
   #    count(HOSP_SERVICE, sort=TRUE) %>%
   #    collect()
   icu_service_groups <- list(
      'surgery' = c(
         'General Surgery', 'Cardiothoracic Surgery', 'Thoracic Surgery', 'Cardiac Surgery',
         'Vascular Surgery', 'Colon/Rectal Surgery', 'Surgical Oncology', 'Orthopedic Surgery',
         'Orthopaedic Surgery', 'Pediatric Surgery', 'Plastic Surgery', 'Breast Surgery',
         'Gynecologic Oncology', 'Gynocologic Oncology', 'BARIATRIC/SURG (MIS)',
         'Oral Surgery', 'Oral and Maxillofacial Surg', 'Urogynecology', 'Transplant Medicine'
      ),
      
      'critical' = c(
         'Critical Care Med', 'Critical Care Medicine', 'Emergency Medicine', 'Emergency Care',
         'Intensive Care', 'Trauma', 'Anesthesiology', 'Anesthesia'
      ),
      
      'oncology' = c(
         'Hematology/Oncology', 'Surgical Oncology', 'Hematology Oncology / BMT', 
         'Radiation Oncology', 'Gynecologic Oncology', 'Gynocologic Oncology'
      ),
      
      'ob_gyn' = c(
         'Obstetrics & Gynecology', 'Gynecology', 'Maternal-Fetal Medicine',
         'Maternal & Fetal Medicine-Obs & Gynecolo', 'Midwife', 'Obstetrics',
         'COMMUNITY MED (OB)', 'Gynecologic Oncology'
      )
   )
   
   if (preprocess) {
      icusDF <- icusDF %>%
         mutate(
            across(
               .cols = c(ADMIT_DATE, DISCHARGE_DATE, ICU_EFFECTIVE_DATE),
               .fns = ~ strptime(., format = '%m/%d/%Y %T')
            ),
         )
   }
   
   # Join and make sure that the ICU_EFFECTIVE_DATE for a given row in df
   # is within 48-72 hours of blood culture, rather than some other time during the visit.
   df <- df %>%
      mutate(
         JOIN_START = ORDER_DATE + time_frame[1] * 86400,
         JOIN_END = ORDER_DATE + time_frame[2] * 86400
      ) %>%
      left_join(
         x = .,
         y = icusDF %>% select(PERSON_ID, HOSP_SERVICE, ICU_EFFECTIVE_DATE) %>% distinct(),
         by = join_by(
            PERSON_ID,
            between(y$ICU_EFFECTIVE_DATE, x$JOIN_START, x$JOIN_END)
         )
      ) %>%
      mutate(ENC_icu_stay = !is.na(ICU_EFFECTIVE_DATE)) %>%
      select(-JOIN_START, -JOIN_END, -ICU_EFFECTIVE_DATE) %>%
      distinct()
   
   # Featurize the icu service type, then reduce back to one row per visit.
   for (g in seq_along(icu_service_groups)) {
      df[[paste0('ENC_icu_', names(icu_service_groups)[g])]] <- df$HOSP_SERVICE %in% icu_service_groups[[g]]
   }
   df <- df %>%
      select(-HOSP_SERVICE) %>%
      group_by(PERSON_ID, ORDER_DAY) %>%
      mutate(
         across(
            .cols = contains('ENC_icu_'),
            .fns = any
         )
      ) %>%
      ungroup() %>%
      distinct()
   
   # Remove an icu service type column if there aren't at least 2% of rows that have it.
   remove_icu_cols <- names(which(sapply(df %>% select(contains('ENC_icu_')), sum) / nrow(df) < 0.02))
   if (length(remove_icu_cols) > 0L)
      df <- df %>% select(-!!remove_icu_cols)
   
   return(df)
}