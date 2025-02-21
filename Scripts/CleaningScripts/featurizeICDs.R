featurizeICDs <- function(df, icdsDF, end_time_frame = 4L) {
   df <- df %>% 
      select(PERSON_ID, ORDER_DAY) %>%
      distinct() %>%
      mutate(
         TIME_AFTER_ORDER = ORDER_DAY + end_time_frame,
         WEEK_BEFORE_ORDER = ORDER_DAY - 7,
         MONTH_BEFORE_ORDER = ORDER_DAY - 30,
         TWO_YEARS_BEFORE_ORDER = ORDER_DAY - 730
      )
   icdsDF <- icdsDF %>%
      select(PERSON_ID, DX_CODE, CODE_DESCRIPTION, DX_DATE) %>%
      distinct()
   cat('Preprocessed data.\n')
   
   #### WITHIN WEEK BEFORE BLOOD CULTURE
   dfx <- df %>%
      left_join(
         x = .,
         y = icdsDF,
         by = join_by(
            PERSON_ID,
            WEEK_BEFORE_ORDER <= DX_DATE,
            TIME_AFTER_ORDER >= DX_DATE
         )
      ) %>%
      group_by(PERSON_ID, ORDER_DAY) %>%
      mutate(
         SepticShock_1w = any(grepl('with septic shock', CODE_DESCRIPTION)),
         Sepsis_1w = any(grepl('sepsis', CODE_DESCRIPTION)), 
         AKI_1w = any(grepl('^N17.9', DX_CODE)),
         Endocarditis_1w = any(grepl('^I33.0', DX_CODE)),
         Osteomyelitis_1w = any(grepl('^M86', DX_CODE) & !grepl('chronic', CODE_DESCRIPTION)),
         Cellulitis_1w = any(grepl('^L03.90', DX_CODE)),
         Peritonitis_1w = any(grepl('^K65.9', DX_CODE)),
         Respiratory_1w = any(grepl('^J', DX_CODE) & grepl('infect', CODE_DESCRIPTION))
      ) %>%
      ungroup() %>%
      select(-DX_DATE, -DX_CODE, -CODE_DESCRIPTION) %>%
      distinct()
   cat('Last week done.\n')
   
   
   #### WITHIN MONTH OF BLOOD CULTURE
   dfx <- dfx %>%
      left_join(
         x = .,
         y = icdsDF,
         by = join_by(
            PERSON_ID,
            MONTH_BEFORE_ORDER <= DX_DATE,
            TIME_AFTER_ORDER >= DX_DATE
         )
      ) %>%
      group_by(PERSON_ID, ORDER_DAY) %>%
      mutate(
         PulmCircDis_1m = any(grepl('^I26|^I27|^I28\\.[089]', DX_CODE)),
         CPD_1m = any(grepl('^I27\\.28|^I27\\.9|J4[0-7]|^J6[0-7]|^J68\\.4|^J70\\.[13]', DX_CODE)),
         MyocInfarc_1m = any(grepl('^I21|^I22|^I25', DX_CODE)),
         MetastSolidTumor_1m = any(grepl('^C7[7-9]|^C80', DX_CODE)),
         Malignancy_1m = any(grepl('^C[01][0-9]|^C2[0-6]|^C3[01234789]|^C4[013]|^C4[5-9]|^C5[0-8]|^C6[0-9]|^C7[0-6]|^C8[123458]|^C9[0-7]', DX_CODE)),
      ) %>%
      ungroup() %>%
      select(-DX_DATE, -DX_CODE, -CODE_DESCRIPTION) %>%
      distinct()
   cat('Last month done.\n')
   
   #### WITHIN 2 YEARS OF BLOOD CULTURE
   dfx <- dfx %>%
      left_join(
         x = .,
         y = icdsDF,
         by = join_by(
            PERSON_ID,
            TWO_YEARS_BEFORE_ORDER <= DX_DATE,
            TIME_AFTER_ORDER >= DX_DATE
         )
      ) %>%
      group_by(PERSON_ID, ORDER_DAY) %>%
      mutate(
         OsteoChronic = any(grepl('M86', DX_CODE) & grepl('chronic', CODE_DESCRIPTION)),
         OnDialysis = any(CODE_DESCRIPTION == 'dependence on renal dialysis'),
         Obesity = any(grepl('^E66', DX_CODE)),
         WeightLoss = any(grepl('^E4[0-6]|^R63\\.4|^R64', DX_CODE)),
         Anemia = any(grepl('^D50\\.[089]|^D5[1-3]', DX_CODE)),
         Hypothyroid = any(grepl('^E0[0-3]|^E89', DX_CODE) & !grepl('E000\\.0|E030|E015\\.2|E001\\.0', DX_CODE)),
         FluidElectroDis = any(grepl('^E22\\.2|^E86|^E87', DX_CODE)),
         Coagulopathy = any(grepl('^D6[5-8]|^D69\\.[16789]', DX_CODE)),
         Alcohol = any(grepl('^F10|^E52|^G62\\.1|^I42\\.6|^K29\\.2|^K70\\.0|^K70\\.3|^K70\\.9|^T51|^Z50\\.2|^Z71\\.4|^Z72\\.1', DX_CODE)),
         Drugs = any(grepl('^F1[12345689]|Z71\\.5|Z72\\.2', DX_CODE)),
         Psychoses = any(grepl('^F20|^F2[234589]|F30\\.2|F31\\.2|F31\\.5', DX_CODE)),
         Depression = any(grepl('^F20\\.4|^F31\\.[3-5]|^F32|^F33|^F34\\.1|^F41\\.2|^F43\\.2', DX_CODE)),
         NeuroDisease = any(grepl('^G1[0-3]|^G2[0-2]|^G25\\.4|^G25\\.5|^G31\\.[289]|^G32|^G3[5-7]|^G4[01]|^G93\\.[14]|^R47\\.0|^R56', DX_CODE)),
         CardiacArrythm = any(grepl('^I44\\.[1-3]|^I45\\.6|^I45\\.9|^I4[7-9]|^R00\\.[018]|^T82\\.1|^Z45\\.0|^Z95\\.0', DX_CODE)),
         MyocInfarc = any(grepl('^I21|^I22|^I25', DX_CODE)),
         CompHypertension = any(grepl('^I1[135]', DX_CODE)),
         UncompHypertension = any(grepl('^I10', DX_CODE)),
         CongHeartFailure = any(grepl('^I42|^I43|^I50', DX_CODE)),
         PeriphVasDis = any(grepl('^I70|^I71|^I73|^I79|^K55|^Z95', DX_CODE)),
         CereVasDis = any(grepl('^G45|^G46|^I6[0-9]', DX_CODE)),
         Dementia = any(grepl('^F0[01235]|^G30', DX_CODE)),
         CPD_Pneum = any(grepl('^I26|^I27|^I28\\.[089]|^J4[0-7]|^J6[0-8]|^J70', DX_CODE)), # this includes additional codes from Elixhauser
         RheumaticDis = any(grepl('^L94\\.[013]|^M05|^M06|^M08|^M12\\.[03]|^M3[0-6]|^M45|M46\\.[189]', DX_CODE)), # this includes additional codes from Elixhauser
         PepticUlcerDis = any(grepl('^K2[5-8]', DX_CODE)),
         MildLiverDis = any(grepl('^B18|^K7[1346]|^Z94', DX_CODE)),
         Diabetes = any(grepl('^E1[1-4]', DX_CODE)),
         HemiParaplegia = any(grepl('^G8[0-4]', DX_CODE)),
         RenalDis = any(grepl('^N03|^N05|^N18|^N19|^Z49', DX_CODE)),
         Malignancy = any(grepl('^C[01][0-9]|^C2[0-6]|^C3[01234789]|^C4[013]|^C4[5-9]|^C5[0-8]|^C6[0-9]|^C7[0-6]|^C8[123458]|^C9[0-7]', DX_CODE)),
         ModSevLivDis = any(grepl('I85|K72', DX_CODE)),
         MetastSolidTumor = any(grepl('^C7[7-9]|^C80', DX_CODE)),
         AIDS_HIV = any(grepl('^B2[0124]', DX_CODE)),
         Hyperlipid = any(grepl('(^| )hyperlipid', CODE_DESCRIPTION)),
         Smoking = any(grepl('nicotine', CODE_DESCRIPTION))
      ) %>%
      ungroup() %>%
      select(-DX_DATE, -DX_CODE, -CODE_DESCRIPTION) %>%
      distinct()
   cat('Last 2 years done.\n')
   
   
   dfx <- dfx %>%
      mutate(
         across(.cols = SepticShock_1w:Smoking,
                .fns = as.integer)
      ) %>%
      select(-TIME_AFTER_ORDER, -WEEK_BEFORE_ORDER, -MONTH_BEFORE_ORDER, -TWO_YEARS_BEFORE_ORDER)
   cat('Final processing.\n')
   
   names(dfx)[!names(dfx) %in% c('PERSON_ID', 'ORDER_DAY')] <- paste0('ICD_', names(dfx)[!names(dfx) %in% c('PERSON_ID', 'ORDER_DAY')])
   
   return(dfx)
}
