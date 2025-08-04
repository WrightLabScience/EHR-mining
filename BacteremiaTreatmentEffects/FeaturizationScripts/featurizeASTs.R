featurizeASTs <- function(df, astDF, cohort_name, time_frame) {
   start_time_frame <- time_frame[1]
   end_time_frame <- time_frame[2]
   
   # Prep columns we will keep for later.
   keep_cols <- c('PERSON_ID', 'ORDER_DATE', 'BUG', 'SITE')
   if (grepl('E_faecalis', cohort_name)) 
      keep_cols <- c(keep_cols, 'AMPICILLIN')
   
   # Heavily modify astDF data for analysis.
   other_astDF <- astDF %>%
      filter(
         PERSON_ID %in% unique(df$PERSON_ID),
         lubridate::year(ORDER_DAY) %in% years,
         BUG %in% unname(unlist(bug_groups)),
         !SITE %in% c('urine', 'other')
      ) %>%
      mutate(
         BUG = gsub('^([A-Z])[a-z]+ ([a-z]+)$', '\\1_\\2', BUG),
         BUG = case_when(
            ESBL == 1L & BUG %in% gsub('^([A-Z])[a-z]+ ([a-z]+)$', '\\1_\\2', bug_groups$Enterobacterales) ~ paste0(BUG, '_ESBL'),
            OXACILLIN == 1L & BUG == 'S_aureus' ~ paste0(BUG, '_MRSA'),
            VANCOMYCIN == 1L & BUG %in% c('E_faecalis', 'E_faecium') ~ paste0(BUG, '_VRE'),
            .default = BUG
         )
      ) %>%
      filter(!BUG %in% cohort_bug) %>%
      mutate(
         BUG = case_when(
            BUG == 'S_lugdunensis' & OXACILLIN == 1L ~ 'MRSL',
            BUG == 'S_aureus' & OXACILLIN == 1L ~ 'MRSA',
            BUG == 'S_aureus' & OXACILLIN == 0L ~ 'MSSA',
            BUG == 'S_lugdunensis' & OXACILLIN == 0L ~ 'MSSL',
            ESBL == 1L & BUG %in% gsub('^([A-Z])[a-z]+ ([a-z]+)$', '\\1_\\2', bug_groups$Enterobacterales) ~ 'ESBL',
            ESBL == 0L & BUG %in% gsub('^([A-Z])[a-z]+ ([a-z]+)$', '\\1_\\2', bug_groups$Enterobacterales) ~ 'Enterobacterales',
            BUG %in% paste0('E_faec', c('ium', 'alis')) & VANCOMYCIN == 1L ~ 'VRE',
            BUG %in% gsub('^([A-Z])[a-z]+ ([a-z]+)$', '\\1_\\2', bug_groups$Streptococci) ~ 'Streptococci',
            BUG %in% gsub('^([A-Z])[a-z]+ ([a-z]+)$', '\\1_\\2', bug_groups$NonFermGN) ~ 'NonFermGN',
            .default = BUG
         )
      ) %>%
      select(!!keep_cols) %>%
      distinct() %>%
      rename(
         ORDER_DATE2 = ORDER_DATE,
         other_BUG = BUG
      )
   
   # Join with df data to keep relevant events
   other_astDF <- other_astDF %>%
      right_join(
         x = .,
         y = df %>%
            select(PERSON_ID, ORDER_DATE) %>%
            mutate(
               JOIN_START = ORDER_DATE + start_time_frame * 86400,
               JOIN_END = ORDER_DATE + end_time_frame * 86400
            ),
         by = join_by(
            PERSON_ID,
            between(ORDER_DATE2, JOIN_START, JOIN_END)
         )
      ) %>%
      select(-JOIN_START, -JOIN_END, -ORDER_DATE2) %>%
      distinct()
   
   # Remove contraindicating concomitant infections.
   if (grepl('S_aureus_MRSA', cohort_name)) {
      other_astDF <- other_astDF %>%
         group_by(PERSON_ID, ORDER_DATE) %>%
         filter(is.na(other_BUG) | !any(SITE == 'respiratory' & other_BUG == 'MRSA')) %>%
         ungroup()
      
   } else if (grepl('E_faecium', cohort_name)) {
      other_astDF <- other_astDF %>%
         group_by(PERSON_ID, ORDER_DATE) %>%
         filter(is.na(other_BUG) | !any(SITE == 'respiratory' & other_BUG == 'VRE')) %>%
         ungroup()
      
   } else if (grepl('E_faecalis', cohort_name)) {
      other_astDF <- other_astDF %>%
         mutate(AMPICILLIN = ifelse(is.na(AMPICILLIN), 0L, AMPICILLIN)) %>%
         group_by(PERSON_ID, ORDER_DATE) %>%
         filter(!(any(other_BUG %in% c('MRSA', 'MRSL', 'MSSA') | (other_BUG == 'E_faecium_VRE' & AMPICILLIN == 1L) | grepl('ESBL', other_BUG) | other_BUG %in% gsub('^([A-Z])[a-z]+ ([a-z]+)$', '\\1_\\2', bug_groups$NonFermGN)))) %>%
         ungroup() %>%
         select(-AMPICILLIN)
   }
   
   # Pivot to wide format.
   other_astDF <- other_astDF %>%
      select(-SITE) %>%
      distinct() %>%
      mutate(
         X = TRUE,
         other_BUG = ifelse(is.na(other_BUG), 'no_other_bug', paste0('BUG_', other_BUG))
      ) %>%
      tidyr::pivot_wider(
         names_from = other_BUG,
         values_from = X,
         values_fill = FALSE
      )
   
   # Remove 'no other bug' column, if it exists.
   if ('no_other_bug' %in% names(other_astDF))
      other_astDF <- other_astDF %>% select(-no_other_bug)
   
   # Keep only other bug columns if they exceed 1% prevalence.
   keep_bug_cols <- c('PERSON_ID', 'ORDER_DATE', names(which(sapply(other_astDF %>% select(-PERSON_ID, -ORDER_DATE), sum) / nrow(other_astDF) > 0.01)))
   
   # Do a right join because some rows of other_astDF may have been removed (due to contraindications).
   df <- df %>%
      right_join(
         x = .,
         y = other_astDF %>% select(PERSON_ID, ORDER_DATE, !!keep_bug_cols),
         by = join_by(PERSON_ID, ORDER_DATE)
      )
   
   
   # Add flag for if MRSA isolate in last year. Also add vancomycin MIC values.
   if (cohort_name == 'S_aureus_MRSA_switch_VAN_DAP_treatment' | cohort_name == 'S_aureus_MRSA_VAN_DAP_treatment') {
      df <- left_join(
         x = df %>% 
            mutate(
               JOIN_START = ORDER_DAY - 365L,
               JOIN_END = ORDER_DAY - 1L
            ),
         y = astDF_og %>%
            filter(
               PERSON_ID %in% unique(df$PERSON_ID),
               lubridate::year(ORDER_DAY) %in% years,
               BUG == 'Staphylococcus aureus', OXACILLIN == 1L,
               !SITE %in% c('urine', 'other')
            ) %>%
            select(PERSON_ID, ORDER_DAY) %>%
            distinct() %>%
            rename(LAST_MRSA_DAY = ORDER_DAY),
         by = join_by(
            PERSON_ID,
            between(y$LAST_MRSA_DAY, x$JOIN_START, x$JOIN_END)
         )
      ) %>%
         mutate(DEM_MRSA_1yr = !is.na(LAST_MRSA_DAY)) %>%
         select(-LAST_MRSA_DAY, -JOIN_START, -JOIN_END) %>%
         distinct()
      
      # Add vancomycin MIC values.
      load(file = '~/Desktop/EHR/EHR work/RdataFiles/BuildCohortDatasets/RawData/mrsa_vancomycin_MIC_values_2015_2024.Rdata')
      df <- df %>%
         left_join(
            x = .,
            y = micDF,
            by = join_by(PERSON_ID, ORDER_DAY)
         ) %>%
         mutate(
            VAN_MIC_LOG2 = as.integer(round(log2(VAN_MIC))),
            VAN_MIC_LOG2 = ifelse(is.na(VAN_MIC_LOG2), median(VAN_MIC_LOG2, na.rm=T), VAN_MIC_LOG2)
         ) %>%
         select(-VAN_MIC)
   }
   
   cat('Removed units with concomitant infections that contraindicate either treatment:', nrow(df), '\n')
   
   # Add BSI recurrence days.
   recur_astDF <- astDF %>%
      filter(
         PERSON_ID %in% unique(df$PERSON_ID),
         SITE == 'blood',
         lubridate::year(ORDER_DAY) %in% min(years):(max(years)+1),
         BUG %in% unname(unlist(bug_groups))
      ) %>%
      select(PERSON_ID, ORDER_DAY) %>% 
      distinct() %>% 
      rename(NEXT_CULTURE_DAY = ORDER_DAY)
   
   # Execute join.
   df <- df %>%
      mutate(AFTER_CULTURE = ORDER_DAY + start_time_frame) %>%
      left_join(
         x = .,
         y = recur_astDF,
         by = join_by(
            PERSON_ID,
            closest(y$NEXT_CULTURE_DAY >= x$AFTER_CULTURE)
         )
      ) %>%
      mutate(OUT_recur_time = as.integer(NEXT_CULTURE_DAY - ORDER_DAY)) %>%
      select(-NEXT_CULTURE_DAY, -AFTER_CULTURE) %>%
      distinct()
   
   return(df)
}
