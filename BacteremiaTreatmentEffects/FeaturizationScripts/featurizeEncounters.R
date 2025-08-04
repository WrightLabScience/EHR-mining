featurizeEncounters <- function(df, encsDF, time_frame, preprocess=TRUE) {
   og_vars <- names(df)
   source('~/Desktop/EHR-mining/UsefulDataForCleaning/UPMC_facilities/handle_UPMC_site_names.R')
   
   if (preprocess) {
      encsDF <- encsDF %>%
         mutate(
            across(
               .cols = c(ADMIT_DATE, DISCHARGE_DATE),
               .fns = ~ strptime(., format = '%m/%d/%Y %T')
            ),
            ADMIT_DAY = lubridate::as_date(ADMIT_DATE),
            DISCHARGE_DAY = lubridate::as_date(DISCHARGE_DATE)
         )
   }
   
   start_time_frame <- time_frame[1]
   end_time_frame <- time_frame[2]
   
   df <- df %>%
      mutate(
         JOIN_START = ORDER_DATE + start_time_frame * 86400,
         JOIN_END = ORDER_DATE + end_time_frame * 86400
      ) %>%
      left_join(
         x = .,
         y = encsDF %>%
            filter(PERSON_ID %in% unique(df$PERSON_ID)) %>%
            select(PERSON_ID, ADMIT_DATE, DISCHARGE_DATE, ADMIT_DAY, DISCHARGE_DAY) %>% 
            distinct(),
         by = join_by(
            PERSON_ID,
            overlaps(x_lower = JOIN_START,
                     x_upper = JOIN_END,
                     y_lower = ADMIT_DAY,
                     y_upper = DISCHARGE_DAY)
         )
      ) %>%
      select(-JOIN_START, -JOIN_END)
   if (nrow(df %>% filter(n() > 1L, .by=c(PERSON_ID, ORDER_DAY))) > 0L) {
      df <- rbind(
         df %>% filter(n() == 1L, .by=c(PERSON_ID, ORDER_DAY)),
         df %>% 
            filter(n() > 1L, .by=c(PERSON_ID, ORDER_DAY)) %>%
            summarise(
               ADMIT_DATE = min(ADMIT_DATE),
               ADMIT_DAY = min(ADMIT_DAY),
               DISCHARGE_DATE = max(DISCHARGE_DATE),
               DISCHARGE_DAY = max(DISCHARGE_DAY),
               .by = c(!!og_vars)
            )
      ) %>%
         arrange(PERSON_ID, ORDER_DATE)
   }
   
   df <- df %>% filter(!is.na(ADMIT_DAY))
   df <- df %>% filter(DISCHARGE_DATE > ORDER_DATE)
   cat('Removed units without any encounters data:', nrow(df), '\n')
   
   df <- df %>% 
      mutate(
         ENC_time_admit_order = as.numeric(lubridate::as.duration(ORDER_DATE - ADMIT_DATE)) / 86400,
         OUT_los_after_order = as.numeric(lubridate::as.duration(DISCHARGE_DATE - ORDER_DATE)) / 86400,
      ) %>%
      select(!!og_vars, OUT_los_after_order, ENC_time_admit_order)
   
   df <- df %>% filter(OUT_los_after_order > 3)
   
   
   df <- df %>%
      left_join(
         x = .,
         y = encsDF %>% select(-FACILITY) %>% distinct(),
         by = join_by(
            PERSON_ID,
            closest(y$DISCHARGE_DAY > x$ORDER_DAY)
         )
      )
   df <- df %>% group_by(PERSON_ID, ORDER_DAY)
   df <- rbind(
      df %>% 
         filter(n() == 1L) %>% 
         ungroup(),
      df %>%
         filter(n() > 1L) %>%
         slice_min(ADMIT_DATE, with_ties = FALSE) %>%
         ungroup()
   ) %>%
      arrange(PERSON_ID, ORDER_DATE)
   
   if (any(is.na(df$ADMIT_DAY))) {
      df <- df %>% filter(!is.na(ADMIT_DAY))
   }
   cat('Removed units discharged before end of treatment period:', nrow(df), '\n')
   
   # Admit type, source, ICU
   df <- df %>%
      rename(
         ADMIT_TYPE_OG = ADMIT_TYPE,
         ADMIT_SOURCE_OG = ADMIT_SOURCE
      ) %>%
      mutate(
         
         ADMIT_TYPE = case_when(
            ADMIT_TYPE_OG %in% c('ED', 'Emergency', 'DEC', 'DEA') ~ 'ED',
            grepl('SDS', ADMIT_TYPE_OG) ~ 'Surgery',
            ADMIT_TYPE_OG == 'Trauma Center' ~ 'Trauma',
            ADMIT_TYPE_OG %in% c('Urgent', 'Elective') ~ ADMIT_TYPE_OG,
            .default = 'Other'
         ),
         ADMIT_TYPE = as.factor(ADMIT_TYPE),
         
         ADMIT_SOURCE = case_when(
            ADMIT_SOURCE_OG %in% c('Nursing Home', 'Trans Frm hospice') ~ 'NursingHome',
            grepl('Trans', ADMIT_SOURCE_OG) ~ 'Transfer',
            ADMIT_SOURCE_OG %in% c('Other Health Facility', 'Court/Law Enforcement', 'Long Term Acute Care Hospital', 'Clinic Referral') ~ 'Transfer',
            ADMIT_SOURCE_OG %in% c('Med Staff Referral', 'HMO Referral') ~ 'MedReferral',
            ADMIT_SOURCE_OG %in% c('Non-Staff Referral', 'Non-Healthcare Facility�') ~ 'Other',
            .default = 'Other'
         ),
         ADMIT_SOURCE = case_when(
            ENCOUNTER_TYPE == 'Inpt Maternity' ~ 'Maternity',
            .default = ADMIT_SOURCE
         ),
         ADMIT_SOURCE = as.factor(ADMIT_SOURCE),
         
         ICU_FLAG = recode(.x=ICU_FLAG, Y=TRUE, N=FALSE)
      ) %>%
      select(!c(ADMIT_TYPE_OG, ADMIT_SOURCE_OG, ENCOUNTER_TYPE, ADMIT_DATE, DISCHARGE_DATE, LENGTH_OF_STAY, ADMIT_DAY, DISCHARGE_DAY)) %>%
      rename(
         ENC_ADMIT_TYPE = ADMIT_TYPE,
         ENC_ADMIT_SOURCE = ADMIT_SOURCE,
         ENC_ICU_stay = ICU_FLAG
      )
   
   
   # what about transfers from regional/rural/community to quaternary??
   df <- df %>%
      mutate(
         JOIN_START = ORDER_DATE,
         JOIN_END = ORDER_DATE + 4 * 86400
      ) %>%
      left_join(
         x = .,
         y = encsDF %>% 
            filter(PERSON_ID %in% unique(df$PERSON_ID)) %>%
            select(PERSON_ID, ADMIT_DATE, DISCHARGE_DATE, ADMIT_DAY, DISCHARGE_DAY, FACILITY) %>% 
            distinct(),
         by = join_by(
            PERSON_ID,
            overlaps(x_lower = JOIN_START,
                     x_upper = JOIN_END,
                     y_lower = ADMIT_DATE,
                     y_upper = DISCHARGE_DATE)
         )
      )
   if (any(is.na(df$ADMIT_DAY))) {
      df <- df %>% filter(!is.na(ADMIT_DAY))
   }
   
   df <- df %>% group_by(PERSON_ID, ORDER_DAY)
   df <- rbind(
      df %>%
         filter(length(unique(FACILITY)) == 1L) %>%
         ungroup() %>%
         mutate(
            ENC_ADMIT_FACILITY = FACILITY,
            ENC_NEXT_FACILITY = FACILITY,
            ENC_TRANSFER_DURING_TRT = FALSE
         ),
      df %>%
         filter(length(unique(FACILITY)) > 1L) %>%
         mutate(
            ADMIT_FAC = DISCHARGE_DATE < JOIN_END,
            ENC_ADMIT_FACILITY = ifelse(test = any(ADMIT_FAC),
                                        yes = first(FACILITY[ADMIT_FAC]), 
                                        no = NA),
            TRANSFER_DATE = first(DISCHARGE_DATE),
            ENC_NEXT_FACILITY = last(FACILITY)
         ) %>%
         ungroup() %>%
         mutate(
            ENC_TRANSFER_DURING_TRT = as.numeric(lubridate::as.duration(TRANSFER_DATE - ORDER_DATE)) / 86400 < 4
         ) %>%
         select(-ADMIT_FAC, -TRANSFER_DATE)
   ) %>%
      select(!c(ADMIT_DATE, DISCHARGE_DATE, ADMIT_DAY, DISCHARGE_DAY, FACILITY, JOIN_START, JOIN_END)) %>%
      distinct() %>%
      mutate(
         ENC_ADMIT_FACILITY = as.factor(code2fac_fxn(ENC_ADMIT_FACILITY)),
         ENC_NEXT_FACILITY = as.factor(code2fac_fxn(ENC_NEXT_FACILITY))
      ) %>%
      arrange(PERSON_ID, ORDER_DATE)
   
   
   # mark as any transfer?
   df <- df %>%
      mutate(
         ENC_ADMIT_FAC_CAT = fac2cat_fxn(ENC_ADMIT_FACILITY),
         ENC_NEXT_FAC_CAT = fac2cat_fxn(ENC_NEXT_FACILITY)
      ) %>%
      mutate(
         ENC_TRANSFER_2_QUATERN = ENC_ADMIT_FAC_CAT != 'Academic' & ENC_NEXT_FAC_CAT == 'Academic'
      ) %>%
      select(-ENC_ADMIT_FAC_CAT, -ENC_NEXT_FAC_CAT)
   
   # how many days spent inpatient in last year?
   df <- left_join(
      x = df %>% 
         mutate(
            JOIN_START = ORDER_DAY - 365L,
            JOIN_END = ORDER_DAY
         ),
      y = encsDF %>% 
         filter(ENCOUNTER_TYPE %in% c('Inpatient', 'Inpt Maternity', 'Inpatient Nursery')) %>%
         select(PERSON_ID, ADMIT_DAY, DISCHARGE_DAY) %>%
         distinct(),
      by = join_by(
         PERSON_ID,
         overlaps(x_lower = JOIN_START,
                  x_upper = JOIN_END,
                  y_lower = ADMIT_DAY,
                  y_upper = DISCHARGE_DAY)
      )
   ) %>%
      mutate(
         X1 = ifelse(ADMIT_DAY < JOIN_START, as.integer(JOIN_START - ADMIT_DAY), 0L),
         X2 = ifelse(DISCHARGE_DAY > JOIN_END, as.integer(DISCHARGE_DAY - JOIN_END), 0L),
         ENC_days_inpt_1yr = as.integer(DISCHARGE_DAY - ADMIT_DAY) - X1 - X2,
         across(c(X1, X2, ENC_days_inpt_1yr), ~ ifelse(is.na(.), 0L, .)),
      ) %>%
      group_by(PERSON_ID, ORDER_DAY) %>%
      mutate(ENC_days_inpt_1yr = sum(ENC_days_inpt_1yr)) %>%
      ungroup() %>%
      select(!c(X1, X2, JOIN_START, JOIN_END, ADMIT_DAY, DISCHARGE_DAY)) %>%
      distinct()
   
   # rename to fix capitalization
   w <- grep('^ENC_.+', names(df))
   enc_vars <- names(df)[w]
   enc_vars <- paste0('ENC_', tolower(gsub('^ENC_(.+)$', '\\1', enc_vars)))
   names(df)[w] <- enc_vars
   
   # To simplify and reduce number of variables slightly, only keep admit facility.
   df <- df %>% 
      mutate(ENC_admit_facility = fac2cat_fxn(ENC_admit_facility)) %>%
      select(-ENC_next_facility)
   
   # One-hot encode as specific large academic facility or just facility category (regional, community, or rural).
   fct_vars <- c('ENC_admit_type', 'ENC_admit_source', 'ENC_admit_facility')
   for (var in fct_vars) {
      uniq_vals <- unique(df[[var]])
      for (val in uniq_vals) {
         df[[paste0(var, '_', val)]] <- df[[var]] == val
      }
      df <- df %>% select(-!!var)
   }
   
   # Remove if one category is <1% prevalence.
   for (var in fct_vars) {
      rmv_cols <- names(which(sapply(df %>% select(contains(!!var)), sum) / nrow(df) < 0.01))
      if (length(rmv_cols) > 0L)
         df <- df[!names(df) %in% rmv_cols]
   }
   
   return(df)
}
