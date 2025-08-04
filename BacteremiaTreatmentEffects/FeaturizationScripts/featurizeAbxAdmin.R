featurizeAbxAdmin <- function(df, abxDF, preprocess=TRUE) {
   #     Tasks:
   # preprocess abxDF
   # any ABX in last month
   # remove if no abx in 0-48 hours after blood culture
   # first abx time, number of empiric doses
   # deal with antibiotic switches
   # featurize early/empiric administrations
   
   
   # Transform raw antibiotic administration data.
   if (preprocess) {
      # Load antibiotic names.
      load(file = '~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/abx_names_for_abxDF.Rdata')
      
      # Rename, modify, and reduce raw abx admin data.
      abxDF <- abxDF %>%
         rename(START_DATE = ADMIN_START_DATE) %>%
         mutate(PERSON_ID = as.character(PERSON_ID)) %>%
         select(PERSON_ID, START_DATE, MEDICATION) %>%
         distinct()
      
      # Standardized abx names.
      meds <- abxDF %>% count(MEDICATION, sort=TRUE)
      meds$ABX <- ''
      for (abx in abx_all) {
         w <- grep(abx, meds$MEDICATION)
         if (length(w) > 0L) 
            meds$ABX[w] <- paste(meds$ABX[w], abx, sep=',')
      }
      rm(w, abx_all, abx)
      
      meds <- meds %>%
         mutate(ABX = gsub('^,', '', ABX)) %>%
         mutate(ABX = case_when(
            ABX == 'PIPERACILLIN,TAZOBACTAM' ~ 'PIPERACILLIN/TAZOBACTAM',
            ABX == 'AMPICILLIN,SULBACTAM' ~ 'AMPICILLIN/SULBACTAM',
            ABX == 'CIPROFLOXACIN,OFLOXACIN' ~ 'CIPROFLOXACIN',
            ABX == 'AMOXICILLIN,CLAVULANATE' ~ 'AMOXICILLIN/CLAVULANATE',
            ABX == 'LEVOFLOXACIN,OFLOXACIN' ~ 'LEVOFLOXACIN',
            ABX == 'SULFAMETHOXAZOLE,TRIMETHOPRIM' ~ 'TRIMETHOPRIM/SULFAMETHOXAZOLE',
            ABX == 'DOXYCYCLIN,DOXYCYCLINE' ~ 'DOXYCYCLINE',
            ABX == 'AVIBACTAM,CEFTAZIDIME' ~ 'CEFTAZIDIME/AVIBACTAM',
            ABX == 'CEFTOLOZANE,TAZOBACTAM' ~ 'CEFTOLOZANE/TAZOBACTAM',
            ABX == 'IMIPENEM,RELEBACTAM' ~ 'IMIPENEM/RELEBACTAM',
            ABX == 'CLOXACILLIN,DICLOXACILLIN,OXACILLIN' ~ 'DICLOXACILLIN',
            ABX == 'MEROPENEM,VABORBACTAM' ~ 'MEROPENEM/VABORBACTAM',
            ABX == 'DALFOPRISTIN,QUINUPRISTIN' ~ 'SYNERCID',
            ABX == 'PENICILLIN,PENICILLIN G' ~ 'PENICILLIN G',
            ABX == 'PENICILLIN,PENICILLIN V' ~ 'PENICILLIN V',
            ABX == 'CILASTATIN,IMIPENEM,RELEBACTAM' ~ 'IMIPENEM/RELEBACTAM',
            ABX == 'CILASTATIN,IMIPENEM' ~ 'IMIPENEM',
            ABX == 'CLAVULANATE,TICARCILLIN' ~ 'TICARCILLIN/CLAVULANATE',
            ABX == 'DURLOBACTAM,SULBACTAM' ~ 'SULBACTAM/DURLOBACTAM',
            .default = ABX
         ))# %>% filter(grepl(',', ABX)) %>% count(ABX, sort)
      meds <- setNames(meds$ABX, meds$MEDICATION)
      
      # Add additional column for standardized abx name
      abxDF <- abxDF %>% mutate(ABX = unname(meds[MEDICATION]))
      
      # Handle multi antibiotic administrations that should be considered separately
      separateAndExtend <- function(abxDF, abx_list) {
         w <- which(abxDF$ABX == abx_list)
         if (length(w) == 0L)
            return(abxDF)
         
         abx_list_sep <- strsplit(abx_list, split=',')[[1]]
         df1 <- df2 <- abxDF[w,]
         df1$ABX <- abx_list_sep[1]
         df2$ABX <- abx_list_sep[2]
         df <- rbind(df1, df2)
         if (length(abx_list_sep) == 3L) {
            df3 <- abxDF[w,]
            df3$ABX <- abx_list_sep[3]
            df <- rbind(df, df3)
         }
         return(rbind(abxDF[-w,], df))
      }
      abxDF <- separateAndExtend(abxDF, 'BACITRACIN,NEOMYCIN,POLYMYXIN B')
      abxDF <- separateAndExtend(abxDF, 'BACITRACIN,POLYMYXIN B')
      abxDF <- separateAndExtend(abxDF, 'POLYMYXIN B,TRIMETHOPRIM')
      abxDF <- separateAndExtend(abxDF, 'NEOMYCIN,POLYMYXIN B')
      abxDF <- separateAndExtend(abxDF, 'COLISTIN,GENTAMICIN')
      abxDF <- separateAndExtend(abxDF, 'AMPHOTERICIN,GENTAMICIN,VANCOMYCIN')
      abxDF <- separateAndExtend(abxDF, 'COLISTIN,NEOMYCIN')
      rm(meds, separateAndExtend)
      
      # Remove medication column, arrange.
      abxDF <- abxDF %>%
         select(-MEDICATION) %>%
         distinct() %>%
         arrange(PERSON_ID, START_DATE, ABX)
   }
   
   # First, do all of the abx admin work with just necessary cohort data.
   df_og <- df
   df <- df_og %>% select(PERSON_ID, ORDER_DATE)
   
   # Get number of days received any abx in last month (-2 through -30 days).
   df <- df %>%
      mutate(
         JOIN_START = ORDER_DATE - 30 * 86400,
         JOIN_END = ORDER_DATE - 2 * 86400
      ) %>%
      left_join(
         x = .,
         y = abxDF %>% 
            select(PERSON_ID, START_DATE) %>%
            distinct(),
         by = join_by(
            PERSON_ID, 
            between(y$START_DATE, x$JOIN_START, x$JOIN_END)
         )
      ) %>%
      mutate(START_DATE = lubridate::as_date(START_DATE)) %>%
      distinct() %>%
      mutate(ABX_last1m = !is.na(START_DATE)) %>%
      select(PERSON_ID, ORDER_DATE, ABX_last1m) %>%
      group_by(PERSON_ID, ORDER_DATE) %>%
      mutate(ABX_last1m = sum(ABX_last1m)) %>%
      ungroup() %>%
      distinct()
   
   # Remove if no antibiotics within 48 hours of blood culture.
   df <- df %>%
      mutate(
         JOIN_START = ORDER_DATE,
         JOIN_END = ORDER_DATE + 2 * 86400
      ) %>%
      semi_join(
         x = .,
         y = abxDF %>%
            select(PERSON_ID, START_DATE) %>%
            distinct(),
         by = join_by(
            PERSON_ID, 
            between(y$START_DATE, x$JOIN_START, x$JOIN_END)
         )
      ) %>% 
      select(-JOIN_START, -JOIN_END)
   cat('Removed units with no abx record:', nrow(df), '\n')
   
   # Get time of first antibiotic and number of empiric doses.
   df <- df %>%
      mutate(
         JOIN_START = ORDER_DATE - 2 * 86400,
         JOIN_END = ORDER_DATE + 2 * 86400
      ) %>%
      left_join(
         x = .,
         y = abxDF %>% 
            select(PERSON_ID, START_DATE) %>%
            distinct(),
         by = join_by(
            PERSON_ID, 
            between(y$START_DATE, x$JOIN_START, x$JOIN_END)
         )
      ) %>%
      mutate(ABX_TIME = as.numeric(lubridate::as.duration(START_DATE - ORDER_DATE)) / 86400) %>%
      summarise(
         ABX_numEmpDoses = length(unique(ABX_TIME[ABX_TIME > -2 & ABX_TIME < 0.5])),
         ABX_firstAbxTime = min(ABX_TIME) * 24,
         .by = c(PERSON_ID, ORDER_DATE, ABX_last1m)
      )
   
   # Define treatment window and cohort treatment definitions.
   df <- df %>%
      mutate(
         JOIN_START = ORDER_DATE - 2 * 86400,
         JOIN_END = ORDER_DATE + 4 * 86400
      ) %>%
      left_join(
         x = .,
         y = abxDF %>% 
            filter(ABX %in% unlist(strsplit(cohort_abx, split='+', fixed=TRUE))) %>%
            select(PERSON_ID, START_DATE, ABX) %>%
            distinct(),
         by = join_by(
            PERSON_ID,
            between(y$START_DATE, x$JOIN_START, x$JOIN_END)
         )
      ) %>%
      mutate(ABX_TIME = as.numeric(lubridate::as.duration(START_DATE - ORDER_DATE)) / 86400) %>%
      filter(!is.na(START_DATE)) %>%
      select(-JOIN_START, -JOIN_END, -START_DATE)
   cat('Removed patients that did not receive either trt [0-4]:', df %>% group_by(PERSON_ID, ORDER_DATE) %>% n_groups(), '\n')
   
   # Abbreviate antibiotic names (there might be >1 abx per "treatment").
   if (any(grepl('+', cohort_abx, fixed=TRUE))) {
      df <- df %>%
         mutate(
            ABX = case_when(
               ABX %in% strsplit(cohort_abx, split='+', fixed=TRUE)[[1]] ~ paste(abbr[strsplit(cohort_abx, split='+', fixed=TRUE)[[1]]], collapse='/'),
               ABX %in% strsplit(cohort_abx, split='+', fixed=TRUE)[[2]] ~ paste(abbr[strsplit(cohort_abx, split='+', fixed=TRUE)[[2]]], collapse='/')
            )
         )
   } else {
      df <- df %>% mutate(ABX = unname(abbr[ABX]))
   }
   
   # Classify "treatment" as what was received during this period.
   df <- df %>%
      group_by(PERSON_ID, ORDER_DATE) %>%
      mutate(treatment = paste(unique(ABX), collapse=',')) %>%
      ungroup()
   
   # Deal with overlap and order for ABX switch cases.
   if (!'switch' %in% names(cohort_info[[cohort]])) { # If no switch, quite simple.
      df <- df %>%
         select(-ABX, -ABX_TIME) %>%
         distinct()
      
   } else { # If switch, do this (correct order AND <24 hours overlap)
      cohort_abx_fixed <- unname(sapply(cohort_abx, function(str) paste(abbr[unlist(strsplit(str, split='+', fixed=TRUE))], collapse='/')))
      
      OTH <- cohort_abx_fixed[2]
      ONE <- cohort_abx_fixed[1]
      keep_combo <- paste0(ONE, ',', OTH)
      w <- which(df$treatment == keep_combo)
      
      df <- rbind(
         df[-w,] %>% 
            select(-ABX, -ABX_TIME) %>% 
            distinct() %>%
            mutate(ABX_switch_time = NA_real_),
         df[w,] %>%
            group_by(PERSON_ID, ORDER_DATE) %>%
            mutate(
               max_one = max(ABX_TIME[ABX == ONE]),
               min_one = min(ABX_TIME[ABX == ONE]),
               max_oth = max(ABX_TIME[ABX == OTH]),
               min_oth = min(ABX_TIME[ABX == OTH]),
               min1 = min(max_one, max_oth),
               max1 = max(min_one, min_oth),
               overlap = min1 - max1
            ) %>%
            ungroup() %>%
            filter(max_oth >= max_one, overlap <= 1) %>%
            rename(ABX_switch_time = min_oth) %>%
            select(!c(max_one, min_one, max_oth, min1, max1, overlap, ABX, ABX_TIME)) %>%
            distinct()
      ) %>%
         arrange(PERSON_ID, ORDER_DATE)
      
      df <- df %>%
         filter(treatment %in% c(keep_combo, ONE)) %>%
         mutate(treatment = ifelse(treatment == keep_combo, OTH, ONE))
   }
   
   # Finalize treatment variables.
   df <- df %>%
      filter(!grepl(',', treatment)) %>%
      mutate(treatment = as.factor(treatment))
   cat('Removed patients that did not match trt definition:', nrow(df), '\n')
   
   # Featurize early/empiric antibiotic treatments.
   df <- df %>%
      mutate(
         JOIN_START = ORDER_DATE - 2 * 86400,
         JOIN_END = ORDER_DATE + 1 * 86400
      ) %>%
      left_join(
         x = .,
         y = abxDF %>%
            filter(
               !ABX %in% unlist(strsplit(cohort_abx, split='+', fixed=TRUE)),
               PERSON_ID %in% unique(df$PERSON_ID)
            ),
         by = join_by(
            PERSON_ID,
            between(y$START_DATE, x$JOIN_START, x$JOIN_END)
         )
      ) %>%
      select(-JOIN_START, -JOIN_END, -START_DATE) %>%
      mutate(
         ABX = unname(abbr[ABX]),
         ABX = case_when(
            ABX %in% c('MEM', 'ETP', 'DOR', 'IPM', 'MVB') ~ 'CAR',
            ABX %in% c('BAC', 'PB', 'NEO', 'NIT', 'FOF', 'RIF', 'TMP', 'EMB', 'CST', 'LEX') ~ 'OTH',
            ABX %in% c('CFZ', 'ATM', 'SXT', 'CPT', 'CLI') ~ ABX,
            ABX %in% c('LVX', 'CIP', 'MXF', 'OFX', 'GAT') ~ 'FLQ',
            ABX %in% c('GEN', 'TOB', 'AMK') ~ 'AMG',
            ABX %in% c('AZM', 'ERY', 'CAM') ~ 'MAC',
            ABX %in% c('DOX', 'TGC', 'MIN', 'TET') ~ 'TET',
            ABX %in% c('TZP', 'SAM', 'AMC') ~ 'BLI',
            ABX %in% c('AMP', 'OXA', 'NFC', 'AMX') ~ 'BLA',
            ABX %in% c('LZD', 'DAP') ~ 'GRP',
            ABX %in% c('CXM', 'FOX', 'CTT') ~ 'CF2',
            ABX %in% c('CRO', 'CAZ', 'POD', 'CTX') ~ 'CF3',
            ABX %in% c('FEP', 'CZA', 'CFT') ~ 'CF4',
            ABX %in% c('VAN', 'TLV') ~ 'GLP',
            .default = NA
         )
      ) %>%
      distinct() %>%
      mutate(
         ABX = paste0('ABX_early_', ABX),
         X = TRUE
      ) %>%
      tidyr::pivot_wider(
         values_from = X,
         names_from = ABX,
         values_fill = FALSE
      )
   if ('ABX_early_NA' %in% names(df))
      df <- df %>% select(-ABX_early_NA)
   
   # Keep only those abx early columns with at least 5% prevalence.
   keep_abx_early_cols <- c(df %>% select(-contains('ABX_early_')) %>% names(),
                            names(which(sapply(df %>% select(contains('ABX_early_')), function(x) sum(x) / nrow(df) >= 0.05))))
   df <- df %>% select(!!keep_abx_early_cols)
   
   # Recombine with original data to get back other columns.
   df <- df %>%
      left_join(
         x = .,
         y = df_og,
         by = join_by(PERSON_ID, ORDER_DATE)
      )
   
   return(df)
}