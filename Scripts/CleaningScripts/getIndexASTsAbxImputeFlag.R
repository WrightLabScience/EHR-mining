getIndexASTsABXsImputeFlag <- function(year=2023L, index_culture_time_filter=42L, end_time_frame=4L) {
   cat('Load ASTs and AbxAdmin datasets.\n')
   load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_clean_ASTs.Rdata')
   load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_CLEANED_2015_2024_AbxAdmin.Rdata')
   load(file = '~/Desktop/EHR-mining/UsefulDataForCleaning/bug_groups.Rdata')
   
   cat('Preprocess ASTs: filter to blood, year, index, remove problems.\n')
   # get just blood cultures for common bugs from the requested year(s)
   astDF <- astDF %>%
      filter(
         BLOOD,
         lubridate::year(ORDER_DAY) %in% year,
         BUG %in% unname(unlist(bug_groups))
      ) %>%
      select(-BLOOD, -RESPIRATORY)
   
   # take only index cultures
   astDF <- astDF %>%
      group_by(PERSON_ID) %>%
      mutate(DAYS_SINCE_PREVIOUS_BLOOD = ORDER_DAY - lag(ORDER_DAY)) %>%
      ungroup() %>%
      filter(is.na(DAYS_SINCE_PREVIOUS_BLOOD) | DAYS_SINCE_PREVIOUS_BLOOD >= index_culture_time_filter)
   
   
   # remove problematic rows
   astDF <- astDF %>%
      mutate(RESULT_DELAY = as.numeric(lubridate::as.duration(RESULT_DATE - ORDER_DATE)) / 86400) %>%
      filter(RESULT_DELAY <= 8)
   w <- which(astDF$BUG == 'Staphylococcus aureus' & is.na(astDF$OXACILLIN)) # 5
   if (length(w) > 0L)
      astDF <- astDF[-w,]
   
   # impute missing AST results
   cat('Impute ASTs.\n')
   source('~/Desktop/EHR-mining/UsefulDataForCleaning/ASTimputation/ASTimputation.R')
   astDF <- imputeASTs(astDF)
   
   # filter abx admin events to this cohort/time-frame
   cat('Filter AbxAdmin data.\n')
   abxDF <- abxDF %>% 
      filter(PERSON_ID %in% unique(astDF$PERSON_ID)) %>%
      filter(lubridate::year(START_DAY) %in% year) %>%
      filter(!ABX %in% c('CASPOFUNGIN', 'FLUCONAZOLE', 'METRONIDAZOLE', 'VORICONAZOLE', 'POSACONAZOLE', 'CLOTRIMAZOLE', 
                         'MICAFUNGIN', 'KETOCONAZOLE', 'AMPHOTERICIN', 'TERBINAFINE', 'ITRACONAZOLE', 'FLUCYTOSINE', 'BACITRACIN')) %>%
      select(-START_DAY, -END_DATE) %>%
      distinct()
   w <- which(abxDF$ABX == 'PENICILLIN G')
   if (length(w) > 0L)
      abxDF$ABX[w] <- 'BENZYLPENICILLIN'
   
   
   
   cat('Join ASTs + AbxAdmin - label empiric vs. definitive.\n')
   # join abx
   ast_abx <- astDF %>%
      mutate( # empiric prescription time frame -48 hours to +12 hours
         JOIN_START = ORDER_DATE - 86400 * 2,
         JOIN_END = ORDER_DATE + 86400 * end_time_frame
      ) %>%
      left_join(
         x = .,
         y = abxDF,
         by = join_by(
            PERSON_ID,
            between(y$START_DATE, x$JOIN_START, x$JOIN_END)
         )
      ) %>%
      select(-JOIN_START, -JOIN_END) %>%
      relocate(ABX, .after=ORDER_DATE) %>%
      mutate(TIME_DIFF_DAYS = as.numeric(lubridate::as.duration(START_DATE - ORDER_DATE)) / 86400)
   
   common_25_abx <- ast_abx %>% 
      select(PERSON_ID, ORDER_DAY, ABX) %>% 
      filter(!is.na(ABX)) %>% 
      count(ABX, sort=TRUE) %>% 
      slice(1:25) %>% 
      pull(ABX)
   ast_abx <- ast_abx %>% filter(is.na(ABX) | ABX %in% common_25_abx)
   
   # for each person-infection-abx, calculate number of days received abx
   num_trt_days <- ast_abx %>%
      mutate(TIME_DIFF_WDAYS = round(TIME_DIFF_DAYS)) %>%
      select(PERSON_ID, ORDER_DAY, ABX, TIME_DIFF_WDAYS) %>%
      distinct() %>%
      count(PERSON_ID, ORDER_DAY, ABX) %>%
      rename(NUM_DAYS_TRT = n)
      
   df <- ast_abx %>%
      mutate(
         ABX_TIMING = case_when(
            TIME_DIFF_DAYS < 0.5 ~ 'EMP',
            TIME_DIFF_DAYS >= 0.5 ~ 'DEF'
         )
      ) %>%
      select(-START_DATE) %>%
      distinct()
   
   # concordance/discordance flags
   cat('Flag abx as concordant/discordant.\n')
   source('~/Desktop/EHR-mining/Scripts/CleaningScripts/getDiscordanceFlag.R')
   df <- getDiscordanceFlag(df)
   
   df <- df %>%
      group_by(PERSON_ID, ORDER_DAY) %>%
      mutate(TIME_TO_CONC = list(unique(TIME_DIFF_DAYS[FLAG == 'Concordant']))) %>%
      ungroup()
   
   # time to first concordant abx
   df <- df %>%
      select(PERSON_ID, ORDER_DAY, ORDER_DATE, ABX, FLAG, ABX_TIMING, TIME_TO_CONC) %>%
      distinct()
   
   discDF <- df %>%
      summarise(
         EmpDisc = !any(FLAG[ABX_TIMING == 'EMP'] == 'Concordant'),
         DefDisc = !any(FLAG[ABX_TIMING == 'DEF'] == 'Concordant'),
         .by = c(PERSON_ID, ORDER_DAY, ORDER_DATE, TIME_TO_CONC)
      )
   discDF$TIME_TO_CONC <- sapply(discDF$TIME_TO_CONC,
                                 FUN = function(vec) {
                                    if (all(is.na(vec)))
                                       return(NA_real_)
                                    return(min(vec))
                                 })
   w <- which(is.na(discDF$TIME_TO_CONC)) # 4,462
   # discDF[w,] %>% count(EmpDisc, DefDisc) # hmm, most missing didn't receive abx, some only ever received discordant abx (according to my ruleset...)
   if (length(w) > 0L)
      discDF$TIME_TO_CONC[w] <- max(discDF$TIME_TO_CONC[-w])
   
   df <- df %>% select(-TIME_TO_CONC)
   
   # get abx class
   cat('Get abx class - pivot to wide format.\n')
   source('~/Desktop/EHR-mining/UsefulDataForCleaning/getAbxBugClassFxn.R')
   df <- getAbxClass(df)
   
   # join with number of days on each abx dataframe
   df <- df %>%
      left_join(
         x = .,
         y = num_trt_days,
         by = join_by(PERSON_ID, ORDER_DAY, ABX)
      )
   
   
   source('~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/CreateNamedAbxAbbreviations.R')
   df$ABX_LABEL <- paste0(abbr[df$ABX], '_', df$ABX_TIMING)
   w <- which(df$ABX_LABEL %in% c('NA_NA', 'NA_EMP', 'NA_DEF'))
   if (length(w) > 0L)
      df$ABX_LABEL[w] <- NA
   
   df$ABXC_LABEL <- paste0(df$ABXC, '_', df$ABX_TIMING)
   w <- which(df$ABXC_LABEL %in% c('NA_NA', 'NA_EMP', 'NA_DEF'))
   if (length(w) > 0L)
      df$ABXC_LABEL[w] <- NA
   
   df <- df %>%
      select(PERSON_ID, ORDER_DAY, ORDER_DATE, ABX_LABEL, ABXC_LABEL, NUM_DAYS_TRT) %>%
      distinct()
   df <- rbind(
      df %>% select(-ABX_LABEL) %>% rename(TRT_LABEL = ABXC_LABEL),
      df %>% select(-ABXC_LABEL) %>% rename(TRT_LABEL = ABX_LABEL)
   )
   df <- df %>%
      #mutate(X = 1L) %>%
      distinct() %>%
      tidyr::pivot_wider(
         id_cols = c(PERSON_ID, ORDER_DAY, ORDER_DATE),
         names_from = TRT_LABEL,
         values_from = NUM_DAYS_TRT,
         values_fill = 0L
      ) %>%
      select(-`NA`)
   
   
   cat('Combine discordant and treatment labels.\n')
   df <- inner_join(
      x = discDF,
      y = df,
      by = join_by(PERSON_ID, ORDER_DAY, ORDER_DATE)
   )
   
   w <- which(is.na(df$EmpDisc))
   if (length(w) > 0L) {
      for (col in df %>% select(-PERSON_ID, -ORDER_DAY, -ORDER_DATE, -EmpDisc, -DefDisc) %>% names) {
         df[[col]][w] <- NA_integer_
      }
   }
   
   return(df)
}
