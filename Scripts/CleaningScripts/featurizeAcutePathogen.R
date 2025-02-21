featurizeAcutePathRes <- function(df, end_time_frame=4L) {
   cat('Preprocess input dataframe.\n')
   df <- df %>% 
      select(PERSON_ID, ORDER_DAY) %>% 
      distinct()
   
   cat('Loading, filtering, and preprocessing AST data.\n')
   load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_clean_ASTs.Rdata')
   astDF <- astDF %>% 
      filter(PERSON_ID %in% unique(df$PERSON_ID)) %>%
      select(PERSON_ID, ORDER_DAY, BUG, BLOOD, RESPIRATORY, CEFEPIME:ESBL) %>%
      distinct() %>%
      rename(CULTURE_DAY = ORDER_DAY)
   
   
   # impute missing ASTs
   source('~/Desktop/EHR-mining/UsefulDataForCleaning/ASTimputation/ASTimputation.R')
   astDF <- imputeASTs(astDF)
   
   
   cat('Joining ASTs with input dataframe.\n')
   df <- df %>%
      mutate(AFTER_ORDER_DAY = ORDER_DAY + end_time_frame) %>%
      left_join(
         x = .,
         y = astDF %>% filter(BLOOD | RESPIRATORY),
         by = join_by(
            PERSON_ID,
            between(y$CULTURE_DAY, x$ORDER_DAY, x$AFTER_ORDER_DAY)
         )
      ) %>%
      select(-AFTER_ORDER_DAY, -CULTURE_DAY) %>%
      mutate(
         APP_R = case_when(`PIPERACILLIN/TAZOBACTAM` == 1L ~ 1L, .default = 0L),
         CEF_R = case_when(CEFTAZIDIME == 1L | CEFEPIME == 1L ~ 1L, .default = 0L),
         FLQ_R = case_when(CIPROFLOXACIN == 1L | LEVOFLOXACIN == 1L ~ 1L, .default = 0L),
         AMI_R = case_when(TOBRAMYCIN == 1L | GENTAMICIN == 1L | AMIKACIN == 1L ~ 1L, .default = 0L),
         CAR_R = case_when(MEROPENEM == 1L | IMIPENEM == 1L ~ 1L, .default = 0L)
      ) %>%
      mutate(MDR_SCORE = APP_R + CEF_R + FLQ_R + AMI_R + CAR_R) %>%
      mutate(
         MRSA = BUG == 'Staphylococcus aureus' & OXACILLIN == 1L,
         VRE = grepl('Enterococcus', BUG) & VANCOMYCIN == 1L,
         ESBL = !is.na(ESBL) & ESBL == 1L,
         MDR_Pseud = BUG == 'Pseudomonas aeruginosa' & MDR_SCORE >= 3L,
      ) %>%
      select(PERSON_ID, ORDER_DAY, MRSA, VRE, ESBL, MDR_Pseud, BUGC, BLOOD, RESPIRATORY) %>%
      distinct()
   
   df_blood <- df %>%
      filter(BLOOD) %>%
      select(-BLOOD, -RESPIRATORY) %>%
      group_by(PERSON_ID, ORDER_DAY) %>%
      mutate(
         MRSA = as.integer(any(MRSA)),
         VRE = as.integer(any(VRE)),
         ESBL = as.integer(any(ESBL)),
         MDR_Pseud = as.integer(any(MDR_Pseud))
      ) %>%
      ungroup() %>%
      distinct() %>%
      mutate(X = 1L) %>%
      tidyr::pivot_wider(
         values_from = X,
         names_from = BUGC,
         values_fill = 0L
      )
   
   df_resp <- df %>%
      filter(RESPIRATORY) %>%
      select(-BLOOD, -RESPIRATORY) %>%
      group_by(PERSON_ID, ORDER_DAY) %>%
      mutate(
         MRSA = as.integer(any(MRSA)),
         VRE = as.integer(any(VRE)),
         ESBL = as.integer(any(ESBL)),
         MDR_Pseud = as.integer(any(MDR_Pseud))
      ) %>%
      ungroup() %>%
      distinct() %>%
      mutate(X = 1L) %>%
      tidyr::pivot_wider(
         values_from = X,
         names_from = BUGC,
         values_fill = 0L
      )
   names(df_resp)[!names(df_resp) %in% c('PERSON_ID', 'ORDER_DAY')] <- paste0('Respiratory_',
                                                                              names(df_resp)[!names(df_resp) %in% c('PERSON_ID', 'ORDER_DAY')])
   df_resp$RESPIRATORY_ISOLATE <- 1L
   
   df <- left_join(
      x = df_blood,
      y = df_resp,
      by = join_by(PERSON_ID, ORDER_DAY)
   ) %>%
      mutate(across(.cols = contains('Respiratory_'),
                    .fns = ~ ifelse(is.na(.), 0L, .)))
   
   return(df)
}