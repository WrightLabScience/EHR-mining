featurizeLabs <- function(df, labsDF, preprocess=TRUE, time_frame, missingness_threshold=0.3) {
   # Do a massive cleaning job, standardizing lab names, units, etc.
   if (preprocess) {
      cml <- list(
         GLUCOSE = c('GLUCOSE(BEDSIDE TEST)', 'GLUCOSE IN SER/PLAS (MG/DL)', 'GLUCOSE, POC', 'GLUCOSE WHOLE BLOOD', 
                     'WHOLE BLOOD GLUCOSE (POCT)', 'GLUCOSE ISTAT', 'GLUCOSE, WHOLE BLOOD ISTAT'),
         HEMOGLOBIN = c('HGB', 'HEMOGLOBIN-ARTERIAL', 'HEMOGLOBIN ISTAT', 'TOTAL HGB(VENOUS)', 'HEMOGLOBIN (POCT)'),
         HEMOGLOBIN_A1C = c('HEMOGLOBIN A1C'),
         HEMATOCRIT = c('HEMATOCRIT(HCT)', 'HEMATOCRIT ISTAT', 'HEMATOCRIT (POCT)', 'HEMATOCRIT DERIVED', 'HEMATOCRIT(HCT) MANUAL PCV &&',
                        'HEMATOCRIT-BODY FLUID (HCT)', 'HEMATOCRIT VENOUS POC', 'HEMATOCRIT DERIVED - MIXED VEN'),
         POTASSIUM = c('POTASSIUM(K)', 'POTASSIUM ISTAT', 'POTASSIUM(K) WHOLE BLOOD', 'POTASSIUM(K) WHOLE BLOOD &&', 'POTASSIUM (K) (POCT)', 'POTASSIUM K WHOLE BLOOD'),
         CARBON_DIOXIDE = c('CARBON DIOXIDE(CO2)'),
         CHLORIDE = c('CHLORIDE(CL)', 'CHLORIDE, WHOLE BLOOD', 'CHLORIDE ISTAT'),
         SODIUM = c('SODIUM(NA)', 'SODIUM(NA) WHOLE BLOOD', 'SODIUM ISTAT', 'SODIUM (NA) (POCT)', 'SODIUM NA WHOLE BLOOD', 'SODIUM VENOUS POC'),
         UREA_NITROGEN = c('UREA NITROGEN', 'UREA NITROGEN ISTAT', 'UREA NITROGEN, WHOLE BLOOD IST'),
         CREATININE = c('CREATININE', 'CREATININE, WHOLE BLOOD', 'CREATININE ISTAT', 'CREATININE VENOUS ISTAT'),
         RBC = c('RBC'),
         MCV = c('MCV'),
         MCH = c('MCH'),
         MCHC = c('MCHC'),
         WBC = c('WBC'),
         ANION_GAP = c('ANION GAP', 'ANION GAP ISTAT'),
         PLATELETS = c('PLATELETS'),
         MEAN_PLATELET_VOLUME = c('MEAN PLATELET VOLUME'),
         CALCIUM = c('CALCIUM(CA)', 'CALCIUM I-STAT ARTERIAL'),
         EGFR = c('EGFR', 'EGFR (CKD-EPI 2021), SER/PLAS/WHOLE BLD'),
         RDW = c('RDW', 'RDW-CV'),
         ABS_NEUTROPHILS = c('ABS NEUTROPHILS'),
         PCT_NEUTROPHILS = c('NEUTROPHILS', 'NEUTROPHILS %'),
         ABS_LYMPHOCYTES = c('ABS LYMPHOCYTES'),
         PCT_LYMPHOCYTES = c('LYMPHOCYTES'),
         ATYP_PCT_LYMPHOCYTES = c('ATYPICAL LYMPHOCYTES', 'ATYPICAL LYMPHOCYTES @'),
         ATYP_ABS_LYMPHOCYTES = c('ABS ATYPICAL LYMPHOCYTES'),
         ABS_MONOCYTES = c('ABS MONOCYTES'),
         PCT_MONOCYTES = c('MONOCYTES'),
         ABS_BASOPHILS = c('ABS BASOPHILS', 'BASOPHILS, ABS'),
         PCT_BASOPHILS = c('BASOPHILS'),
         ABS_EOSINOPHILS = c('ABS EOSINOPHILS'),
         PCT_EOSINOPHILS = c('EOSINOPHILS'),
         PHOSPHORUS = c('PHOSPHORUS'),
         ALBUMIN = 'ALBUMIN',
         TOTAL_BILIRUBIN = 'TOTAL BILIRUBIN',
         AST = c('ASPARTATE AMINOT.(AST)', 'ASPARTATE AMINOTRANSF. (AST)'),
         ALT = c('ALANINE AMINOTRANS(ALT)'),
         ALP = c('ALKALINE PHOSPHATASE'),
         TOTAL_PROTEIN = 'TOTAL PROTEIN',
         INR = c('INR', 'INR (POC)', 'INR ISTAT'),
         PROTHROMBIN_TIME = c('PROTHROMBIN TIME'),
         LACTATE = c('LACTATE BLOOD', 'LACTATE WHOLE BLOOD', 'LACTATE ARTERIAL (POC)', 'LACTATE ISTAT', 'LACTATE VENOUS (POC)')
      )
      
      # Translate component names to sensible, consistent lab names
      cns <- labsDF %>% 
         count(COMPONENT_NAME, sort=TRUE) %>%
         mutate(LAB = NA)
      for (i in seq_along(cml)) {
         w <- which(cns$COMPONENT_NAME %in% cml[[i]])
         cns$LAB[w] <- names(cml)[i]
      }
      # cns %>% filter(is.na(LAB))
      cns <- cns %>% filter(!is.na(LAB))
      cns <- setNames(cns$LAB, cns$COMPONENT_NAME)
      labsDF <- labsDF %>%
         mutate(LAB = unname(cns[COMPONENT_NAME])) %>%
         relocate(LAB, .before=COMPONENT_NAME) %>%
         filter(!is.na(LAB)) %>%
         filter(!is.na(RESULT_VALUE))
      
      # Handle cutoffs.
      labsDF$CUTOFF <- 0
      w <- grep('<', labsDF$RESULT_VALUE)
      if (length(w) > 0L) {
         labsDF$RESULT_VALUE[w] <- gsub('<', '', labsDF$RESULT_VALUE[w])
         labsDF$CUTOFF[w] <- -1
      }
      w <- grep('>', labsDF$RESULT_VALUE)
      if (length(w) > 0L) {
         labsDF$RESULT_VALUE[w] <- gsub('>', '', labsDF$RESULT_VALUE[w])
         labsDF$CUTOFF[w] <- 1
      }
      
      # Standardize reference units.
      labsDF <- labsDF %>% mutate(REFERENCE_UNIT = tolower(REFERENCE_UNIT))
      units <- list(
         'pcnt' = c('%', '% of total hgb'),
         'thou' = c('x10e+09/l', '10e+9/l', '10*9/l',
                    'x10e*3', 'x10e3/ul', 'x10e+03/ul', '10s3/ul', '10s3ul', 'x10e3/cumm', '10^3', '10[3]',
                    'thousand/ul', 'th/cumm', 'th/ul', 'th/mm3', 'thous/mcl', '/cmm', 'thou',
                    '1000/mcl',
                    'k/ul'),
         'mill' = c('x10e+12/l', '10e+12/l', '10*12/l',
                    '/cmm',
                    '10^3/ul', '10x3 ul',
                    'x10e*6','x10e+06/ul', 'x10e6/ul', '10s6/ul', '10^6', '10[6]',
                    'mil/ul', 'mill/mcl', 'm/ul',  'mil/mm3', 'mill/ul', 'million/ul', 'mil'),
         'mmol' = c('mmol/l', 'meq/l', 'mm/l', 'mmoll'),
         'cells' = c('cells/mcl', 'cells/ul'),
         'iuL' = c('iu/l', 'u/l'),
         'rate' = c('ml/min/1.73**2', 'ml/min/1.73m**2', 'ml/min', 'ml/min/1.73m2', 'ml/min/1.73'),
         'g_dl' = c('g/dl', 'gm/dl', 'g%'),
         'mg_dl' = c('mg/dl'),
         'secs' = c('seconds', 'sec.', 'sec', 'secs.'),
         'fl' = c('fl', 'um3', 'ume3', 'um'),
         'pg' = c('pg', 'uug'),
         'none' = c('therapeutic', 'ratio')
      )
      
      # ABS BASOPHILS SHOULDNT HAVE PERCENTAGE, remove those
      w <- which(grepl('^ABS_', labsDF$LAB) & labsDF$REFERENCE_UNIT == '%')
      if (length(w) > 0L)
         labsDF <- labsDF[-w,]
      
      # PCT_NEUTROPHILS shouldnt have absolute counts, remove them
      w <- which(grepl('^PCT_', labsDF$LAB) & labsDF$REFERENCE_UNIT == 'absolute')
      if (length(w) > 0L)
         labsDF <- labsDF[-w,]
      
      # INR and x10 - not a thing
      w <- which(labsDF$LAB == 'INR' & labsDF$REFERENCE_UNIT == 'x10') # 19
      if (length(w) > 0L)
         labsDF <- labsDF[-w,]
      
      # CLEAN UP UNITS
      labsDF$CLEAN_REF_UNIT <- NA_character_
      for (i in seq_along(units)) {
         labsDF$CLEAN_REF_UNIT[labsDF$REFERENCE_UNIT %in% units[[i]]] <- names(units)[i]
      }
      labsDF <- labsDF %>%
         mutate(CLEAN_REF_UNIT = case_when(
            grepl('^PCT_', LAB) & is.na(REFERENCE_UNIT) ~ 'pcnt',
            LAB == 'INR' ~ 'none',
            LAB == 'EGFR' & is.na(REFERENCE_UNIT) ~ 'rate',
            .default = CLEAN_REF_UNIT
         ))
      labsDF <- labsDF %>% 
         filter(!is.na(CLEAN_REF_UNIT)) %>% 
         relocate(CLEAN_REF_UNIT, .before=REFERENCE_UNIT)
      
      
      # Handle numerical result values.
      w <- grep('%', labsDF$RESULT_VALUE) # 1243
      if (length(w) > 0L)
         labsDF$RESULT_VALUE[w] <- gsub('%', '', labsDF$RESULT_VALUE[w])
      
      # some words??
      labsDF <- labsDF %>% mutate(RESULT_VALUE = tolower(RESULT_VALUE))
      w <- which(labsDF$RESULT_VALUE == 'less than 1')
      if (length(w) > 0L)
         labsDF$RESULT_VALUE[w] <- 1
      # remove rows with words
      labsDF <- labsDF %>% filter(grepl('^[0-9.]+$', RESULT_VALUE))
      
      # Make result values into actual numbers!!!
      labsDF <- labsDF %>% mutate(RESULT_VALUE = as.numeric(RESULT_VALUE))
      
      # Scale some of the result values according to their units
      labsDF <- labsDF %>%
         mutate(RESULT_VALUE = case_when(
            CLEAN_REF_UNIT == 'cells' ~ RESULT_VALUE / 1000,
            .default = RESULT_VALUE
         ))
      
      # Final processing.
      labsDF <- labsDF %>%
         rename(LAB_DATE = ORDER_DATE) %>%
         select(PERSON_ID, LAB_DATE, LAB, RESULT_VALUE) %>%
         distinct() %>%
         mutate(LAB_DATE = strptime(LAB_DATE, format='%m/%d/%Y %T'))
   }
   
   start_time_frame <- time_frame[1]
   end_time_frame <- time_frame[2]
   
   df_og <- df
   df <- df_og %>%
      select(PERSON_ID, ORDER_DATE) %>%
      mutate(
         JOIN_START = ORDER_DATE + start_time_frame * 86400,
         JOIN_END = ORDER_DATE + end_time_frame * 86400
      ) %>%
      left_join(
         x = .,
         y = labsDF,
         by = join_by(
            PERSON_ID,
            between(y$LAB_DATE, x$JOIN_START, x$JOIN_END)
         )
      ) %>%
      select(PERSON_ID, ORDER_DATE, LAB, RESULT_VALUE)
   
   df <- df %>%
      mutate(LAB = paste0('LAB_', LAB)) %>%
      summarise(
         RESULT_VALUE = mean(RESULT_VALUE),
         .by = c(PERSON_ID, ORDER_DATE, LAB)
      ) %>%
      tidyr::pivot_wider(
         values_from = RESULT_VALUE,
         names_from = LAB
      )
   
   if (any(names(df) == 'LAB_NA'))
      df <- df %>% select(-LAB_NA)
   
   
   # Remove lab components if too many missing values across all rows.
   miss_vars <- sapply(X = df %>% select(contains('LAB_')),
                       FUN = function(x) sum(is.na(x)) / length(x))
   if (sum(miss_vars > missingness_threshold) > 0L) {
      remove_vars <- names(miss_vars)[miss_vars > missingness_threshold]
      df <- df[!names(df) %in% remove_vars]
   }
   
   # Remove rows if >50% missingness across all lab components.
   df$LAB_pcnt_missing <- apply(df %>% select(contains('LAB_')), 1, function(x) sum(is.na(x))) / sum(grepl('^LAB_', names(df)))
   df <- df %>% filter(LAB_pcnt_missing < 0.5)
   
   # Impute variables (with mean) if <30% missingness.
   miss_vars <- sapply(X = df %>% select(contains('LAB_')),
                       FUN = function(x) sum(is.na(x)) / length(x))
   if (sum(miss_vars <= missingness_threshold & miss_vars > 0) > 0L) {
      impt_vars <- names(miss_vars)[miss_vars <= missingness_threshold & miss_vars > 0]
      for (var in impt_vars) {
         w <- which(is.na(df[[var]]))
         if (length(w) > 0L) {
            df[[var]][w] <- mean(df[[var]][-w])  
         }
      }
   }
   
   # Rename to fix capitalization
   w <- grep('^LAB_.+', names(df))
   og_lab_names <- names(df)[w]
   og_lab_names <- gsub('^LAB_', '', og_lab_names)
   str <- rep('LAB_', length(w))
   pat <- '^(PCT|ABS|ATYP|ATYP_PCT|ATYP_ABS)_(.+)$'
   
   # first, the ones that need to stay capitalized
   w1 <- which(og_lab_names %in% c('EGFR', 'RBC', 'MCV', 'MCH', 'MCHC', 'WBC', 'RDW', 'AST', 'ALT', 'ALP', 'INR'))
   if (length(w1) > 0L)
      str[w1] <- paste0(str[w1], og_lab_names[w1])
   
   # now the ones that have abbreviations
   w2 <- grep(pat, og_lab_names)
   if (length(w2) > 0L)
      str[w2] <- paste0(str[w2], gsub(pat, '\\1', og_lab_names[w2]), '_', tolower(gsub(pat, '\\2', og_lab_names[w2])))
   
   # now the easy ones
   w3 <- which(str == 'LAB_')
   if (length(w3) > 0L)
      str[w3] <- paste0(str[w3], tolower(og_lab_names[w3]))
   
   names(df)[w] <- str
   
   # Re-join with full dataset with all variables.
   # Left join with labs-containing df because we may have removed rows.
   df <- df %>%
      left_join(
         x = .,
         y = df_og,
         by = join_by(PERSON_ID, ORDER_DATE)
      )
   cat('Removed rows with >50% missingness in laboratory components:', nrow(df), '\n')
   
   return(df)
}
