featurizeICDs <- function(df, icdsDF, end_time_frame = 2L, cohort_name, preprocess = TRUE,
                          week_vars=TRUE, month_vars=TRUE, twoyear_vars=TRUE, verbose=FALSE) {
   # Keep a copy of the original data and a reduced version - join later.
   df_og <- df
   df <- df_og %>% select(PERSON_ID, ORDER_DATE)
   
   # Clean up icdsDF raw data.
   if (preprocess) {
      icdsDF <- icdsDF %>%
         select(PERSON_ID, DX_CODE, CODE_DESCRIPTION, DX_FROM_DATE) %>%
         distinct() %>%
         na.omit() %>%
         mutate(CODE_DESCRIPTION = tolower(CODE_DESCRIPTION)) %>%
         rename(DX_DATE = DX_FROM_DATE)
   }
   
   # Add columns to df for making joins later.
   df <- df %>%
      mutate(
         TIME_AFTER_INDEX = ORDER_DATE + end_time_frame * 86400,
         WEEK_BEFORE_INDEX = ORDER_DATE - 7 * 86400,
         MONTH_BEFORE_INDEX = ORDER_DATE - 30 * 86400,
         TWO_YEARS_BEFORE_INDEX = ORDER_DATE - 730 * 86400
      )
   
   #### WITHIN WEEK BEFORE BLOOD CULTURE
   if (week_vars) {
      df <- df %>%
         left_join(
            x = .,
            y = icdsDF,
            by = join_by(
               PERSON_ID,
               between(y$DX_DATE, x$WEEK_BEFORE_INDEX, x$TIME_AFTER_INDEX)
            )
         ) %>%
         group_by(PERSON_ID, ORDER_DATE) %>%
         mutate(
            SepticShock_1w = any(grepl('with septic shock', CODE_DESCRIPTION)),
            Sepsis_1w = any(grepl('sepsis', CODE_DESCRIPTION)), 
            AKI_1w = any(grepl('^N17.9', DX_CODE)),
            Endocarditis_1w = any(grepl('^I33.0', DX_CODE)),
            Osteomyelitis_1w = any(grepl('^M86', DX_CODE) & !grepl('chronic', CODE_DESCRIPTION)),
            Cellulitis_1w = any(grepl('^L03.90', DX_CODE)),
            Peritonitis_1w = any(grepl('^K65.9', DX_CODE)),
            Respiratory_1w = any(DX_CODE %in% c('J06.9', 'J44.0', 'J22', 'J47.0')),
            Pneumonia_1w = any(DX_CODE %in% c('J16.8', 'J18.9', 'J18.1', 'J15.9', 'J15.1', 'J15.211', 'J13', 'J95.851', 'J18.0', 'J15.6', 'J15.0', 
                                              'J18.8', 'J85.1', 'J15.5', 'J15.8', 'J16.8', 'J14', 'J84.116', 'J15.7', 'J15.4', 'J15.20', 'J17')),
            PneumoniaMRSA_1w = any(DX_CODE == 'J15.212')
         ) %>%
         ungroup() %>%
         select(-DX_DATE, -DX_CODE, -CODE_DESCRIPTION) %>%
         distinct()
   }
   
   #### WITHIN MONTH OF BLOOD CULTURE
   if (month_vars) {
      df <- df %>%
         left_join(
            x = .,
            y = icdsDF,
            by = join_by(
               PERSON_ID,
               between(y$DX_DATE, x$MONTH_BEFORE_INDEX, x$TIME_AFTER_INDEX)
            )
         ) %>%
         group_by(PERSON_ID, ORDER_DATE) %>%
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
   }
   
   #### WITHIN 2 YEARS OF BLOOD CULTURE
   if (twoyear_vars) {
      df <- df %>%
         left_join(
            x = .,
            y = icdsDF,
            by = join_by(
               PERSON_ID,
               between(y$DX_DATE, x$TWO_YEARS_BEFORE_INDEX, x$TIME_AFTER_INDEX)
            )
         ) %>%
         group_by(PERSON_ID, ORDER_DATE) %>%
         mutate(
            OsteoChronic_2y = any(grepl('M86', DX_CODE) & grepl('chronic', CODE_DESCRIPTION)),
            OnDialysis_2y = any(CODE_DESCRIPTION == 'dependence on renal dialysis'),
            Obesity_2y = any(grepl('^E66', DX_CODE)),
            WeightLoss_2y = any(grepl('^E4[0-6]|^R63\\.4|^R64', DX_CODE)),
            Anemia_2y = any(grepl('^D50\\.[089]|^D5[1-3]', DX_CODE)),
            Hypothyroid_2y = any(grepl('^E0[0-3]|^E89', DX_CODE) & !grepl('E000\\.0|E030|E015\\.2|E001\\.0', DX_CODE)),
            FluidElectroDis_2y = any(grepl('^E22\\.2|^E86|^E87', DX_CODE)),
            Coagulopathy_2y = any(grepl('^D6[5-8]|^D69\\.[16789]', DX_CODE)),
            Alcohol_2y = any(grepl('^F10|^E52|^G62\\.1|^I42\\.6|^K29\\.2|^K70\\.0|^K70\\.3|^K70\\.9|^T51|^Z50\\.2|^Z71\\.4|^Z72\\.1', DX_CODE)),
            Drugs_2y = any(grepl('^F1[12345689]|Z71\\.5|Z72\\.2', DX_CODE)),
            Psychoses_2y = any(grepl('^F20|^F2[234589]|F30\\.2|F31\\.2|F31\\.5', DX_CODE)),
            Depression_2y = any(grepl('^F20\\.4|^F31\\.[3-5]|^F32|^F33|^F34\\.1|^F41\\.2|^F43\\.2', DX_CODE)),
            NeuroDisease_2y = any(grepl('^G1[0-3]|^G2[0-2]|^G25\\.4|^G25\\.5|^G31\\.[289]|^G32|^G3[5-7]|^G4[01]|^G93\\.[14]|^R47\\.0|^R56', DX_CODE)),
            CardiacArrythm_2y = any(grepl('^I44\\.[1-3]|^I45\\.6|^I45\\.9|^I4[7-9]|^R00\\.[018]|^T82\\.1|^Z45\\.0|^Z95\\.0', DX_CODE)),
            MyocInfarc_2y = any(grepl('^I21|^I22|^I25', DX_CODE)),
            CompHypertension_2y = any(grepl('^I1[135]', DX_CODE)),
            UncompHypertension_2y = any(grepl('^I10', DX_CODE)),
            CongHeartFailure_2y = any(grepl('^I42|^I43|^I50', DX_CODE)),
            PeriphVasDis_2y = any(grepl('^I70|^I71|^I73|^I79|^K55|^Z95', DX_CODE)),
            CereVasDis_2y = any(grepl('^G45|^G46|^I6[0-9]', DX_CODE)),
            Dementia_2y = any(grepl('^F0[01235]|^G30', DX_CODE)),
            CPD_Pneum_2y = any(grepl('^I26|^I27|^I28\\.[089]|^J4[0-7]|^J6[0-8]|^J70', DX_CODE)), # this includes additional codes from Elixhauser
            RheumaticDis_2y = any(grepl('^L94\\.[013]|^M05|^M06|^M08|^M12\\.[03]|^M3[0-6]|^M45|M46\\.[189]', DX_CODE)), # this includes additional codes from Elixhauser
            PepticUlcerDis_2y = any(grepl('^K2[5-8]', DX_CODE)),
            MildLiverDis_2y = any(grepl('^B18|^K7[1346]|^Z94', DX_CODE)),
            Diabetes_2y = any(grepl('^E1[1-4]', DX_CODE)),
            HemiParaplegia_2y = any(grepl('^G8[0-4]', DX_CODE)),
            RenalDis_2y = any(grepl('^N03|^N05|^N18|^N19|^Z49', DX_CODE)),
            Malignancy_2y = any(grepl('^C[01][0-9]|^C2[0-6]|^C3[01234789]|^C4[013]|^C4[5-9]|^C5[0-8]|^C6[0-9]|^C7[0-6]|^C8[123458]|^C9[0-7]', DX_CODE)),
            ModSevLivDis_2y = any(grepl('I85|K72', DX_CODE)),
            MetastSolidTumor_2y = any(grepl('^C7[7-9]|^C80', DX_CODE)),
            AIDS_HIV_2y = any(grepl('^B2[0124]', DX_CODE)),
            Hyperlipid_2y = any(grepl('(^| )hyperlipid', CODE_DESCRIPTION)),
            Smoking_2y = any(grepl('nicotine', CODE_DESCRIPTION))
         ) %>%
         ungroup() %>%
         select(-DX_DATE, -DX_CODE, -CODE_DESCRIPTION) %>%
         distinct()
   }
   
   # Clean-up df.
   df <- df %>% select(-TIME_AFTER_INDEX, -WEEK_BEFORE_INDEX, -MONTH_BEFORE_INDEX, -TWO_YEARS_BEFORE_INDEX)
   names(df)[!names(df) %in% c('PERSON_ID', 'ORDER_DATE')] <- paste0('ICD_', names(df)[!names(df) %in% c('PERSON_ID', 'ORDER_DATE')])
   
   # Keep only this with >1% prevalence.
   remove_icd_cols <- names(which(sapply(df %>% select(contains('ICD_')), sum) / nrow(df) < 0.01))
   if (length(remove_icd_cols) > 0L)
      df <- df %>% select(-!!remove_icd_cols)
   
   # remove MRSA patients with pneumonia
   if (any(names(df) == 'ICD_PneumoniaMRSA_1w')) {
      if (grepl('S_aureus_MRSA', cohort_name))
         df <- df %>% filter(!ICD_PneumoniaMRSA_1w)
      df <- df %>% select(-ICD_PneumoniaMRSA_1w)
   }
   
   # Re-join with df_og that has all variables.
   df <- df %>%
      left_join(
         x = .,
         y = df_og,
         by = join_by(PERSON_ID, ORDER_DATE)
      )
   
   return(df)
}
