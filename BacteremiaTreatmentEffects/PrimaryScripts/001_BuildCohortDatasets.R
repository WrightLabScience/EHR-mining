# This script builds analytics-ready datasets from EHR data.
# It will build datasets that are pathogen(s)-specific (e.g., MRSA)
# which include EHR-derived features related to:
# antibiotic treatments, encounters, comorbidities, lab results, etc.
# Importantly, the data built in this script will not include any dates or person_ids,
# but STILL TREAT THE DATA AS SENSITIVE AS ALWAYS.
# The next step(s) beyond this script are to perform causal inference using
# propensity score-based approaches (see scripts: 003_, 004_, 005, and 006_)

# Running this script requires the ability to establish a connection to the Neptune database
# hosted by R3 at DBMI. The schema used for this project is AMB_ETL.

# For my MS thesis and the R01 grant proposal related to antibiotics superiority in bloodstream infections,
# there were 8 bacteremia cohorts included
# Since then, I have included a few additional bacteremia cohorts.

# As a brief side project, Erik was interested in applying the same causal inference
# approach to resistance (bacteremia) cohorts where the 'treatment' was specific antibiotic resistance
# phenotypes (e.g., MRSA vs. MSSA [oxacillin resistance in S. aureus]).
# To run that analysis, many changes will need to be made to this and subsequent scripts.

# Enjoy! Email sam.blechman@gmail.com with questions, I will try to be good about responding and helping.


# Setup environment, load relevant objects, etc.
START <- Sys.time()
setwd('~/Desktop/EHR-mining/BacteremiaTreatmentEffects/')
source(file = '~/Desktop/EHR/EHR work/config_file.R') # Establish connection to Neptune.
load(file = 'CohortInfo/cohort_info.Rdata') # Bug and treatment info for each cohort - this was created in script 001.
load(file = '~/Desktop/EHR-mining/UsefulDataForCleaning/CleanPathogenNames/bug_groups.Rdata') # Load priority pathogen names and groups.
source('~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/CreateNamedAbxAbbreviations.R') # Get a named vector of antibiotic abbreviations.

# Load all featurization functions.
sapply(c('Demographics', 'Encounters', 'ICUstays', 'ICDs', 'Labs', 'AbxAdmin', 'Vasopressors', 'ASTs'), 
       \(fname) source(paste0('FeaturizationScripts/featurize', fname, '.R')))

# Global variables used throughout. Could do sensitivity analyses varying these:
years <- 2016:2024
days_separating_bacteremia <- 60L
start_time_frame <- -5L
end_time_frame <- 2L
earliest_death_date <- 4L
encounters_time_frame <- c(-2L, 4L)
icus_start_time_frame <- -7L
vasopressors_time_frame <- c(start_time_frame, end_time_frame)
lab_start_time_frame <- -2L
lab_missingness_threshold <- 0.3

# Prepare AST data:
#     Keep only priority pathogens - ie, no fungus or possible contaminants
#     Rename bugs to be consistent with cohorts in bug_abx_list
#     Keep only 'index' blood cultures - first blood culture in record or first in 60 days
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_clean_ASTs.Rdata')
astDF_og <- astDF
astDF <- astDF_og %>%
   filter(
      SITE == 'blood',
      lubridate::year(ORDER_DAY) %in% years,
      BUG %in% unname(unlist(bug_groups))
   ) %>%
   mutate(
      BUG = gsub('^([A-Z])[a-z]+ ([a-z]+)$', '\\1_\\2', BUG)
   )
astDF$DAYS_SINCE_PRV <- astDF$ORDER_DAY - lag(astDF$ORDER_DAY)
astDF$DAYS_SINCE_PRV[c(1, which(astDF$PERSON_ID != lag(astDF$PERSON_ID)))] <- NA_integer_
astDF <- astDF %>%
   filter(is.na(DAYS_SINCE_PRV) | DAYS_SINCE_PRV > days_separating_bacteremia) %>%
   mutate(ENC_bacteremia_1yr = !is.na(DAYS_SINCE_PRV) & DAYS_SINCE_PRV <= 365L) %>%
   select(-DAYS_SINCE_PRV)

blood_astDF <- astDF %>%
   mutate(
      BUG = case_when(
         ESBL == 1L & BUG %in% gsub('^([A-Z])[a-z]+ ([a-z]+)$', '\\1_\\2', bug_groups$Enterobacterales) ~ paste0(BUG, '_ESBL'),
         OXACILLIN == 1L & BUG == 'S_aureus' ~ paste0(BUG, '_MRSA'),
         VANCOMYCIN == 1L & BUG %in% c('E_faecalis', 'E_faecium') ~ paste0(BUG, '_VRE'),
         .default = BUG
      )
   )

# Ensure the cleaned datasets have a place to go
clean_data_path <- '~/Desktop/EHR/EHR work/RdataFiles/BuildCohortDatasets/CleanedData/'
if (!dir.exists(clean_data_path))
   dir.create(clean_data_path)

# Directory for .txt files containing output of this script.
output_path <- 'CohortInfo/cat_output_df_building/'
if (!dir.exists(output_path))
   dir.create(output_path)

# Function to break up the unique PERSON_IDs into chunks of 1000 (limit for DBI queries).
get_id_chunks <- function(ids) {
   ids <- unique(ids)
   if (length(ids) < 1000)
      return(list(ids))
   idx_chunks <- mapply(
      seq(1, floor(length(ids) / 1000) * 1000 + 1, 1000),
      c(seq(1000, floor(length(ids) / 1000) * 1000, 1000), length(ids)),
      FUN = ':'
   )
   return(lapply(idx_chunks, \(x) ids[x]))
}


# Loop over each cohort (pathogen(s)/antibiotic(s)) and build 
# a dataset according to cohort_info list.
for (cohort in seq_along(cohort_info)) {
   start <- Sys.time()
   
   # Get bug and treatment variables
   cohort_bug <- cohort_info[[cohort]]$bug
   cohort_abx <- cohort_info[[cohort]]$abx
   cohort_name <- cohort_info[[cohort]]$cohort_name
   cat(cohort_name, '\n')
   
   # Open text file to dump all output to.
   sink(file = paste0(output_path, cohort_name, '.txt'))
   
   
   # Keep only requested pathogen(s).
   df <- blood_astDF %>% filter(BUG %in% cohort_bug)
   cat('Filtered to requested pathogens:', nrow(df), '\n')
   
   # Remove isolates if they were resistant to any of the treatments we're investigating.
   for (abx in unlist(strsplit(cohort_abx, split='+', fixed=TRUE))) 
      df <- df %>% filter(get(abx) %in% c(NA_integer_, 0L))
   cat('Removed isolates resistant to treatments:', nrow(df), '\n')
   rm(abx)
   
   # Remove cultures with AST result delay longer than 1 week.
   df <- df %>% filter(as.numeric(lubridate::as.duration(RESULT_DATE - ORDER_DATE) / 86400) < 7)
   cat('Removed long AST delay cultures:', nrow(df), '\n')
   
   # If multiple bugs, one-hot encode them.
   if (length(cohort_bug) > 1L) {
      for (bug in cohort_bug) {
         df[[paste0('BUG_BSI', bug)]] <- df$BUG == bug
      }
      df <- df %>% select(-BUG)
   }
   
   # Keep only certain variables. Add time as a variable.
   df <- df %>% 
      select(PERSON_ID, ORDER_DATE, ORDER_DAY, ENC_bacteremia_1yr) %>%
      mutate(DEM_year_decimal = lubridate::decimal_date(ORDER_DATE) - min(years))
   
   
   # Featurize comcomitant infections.
   cat('ASTs and infection data.\n')
   df <- featurizeASTs(df = df,
                       astDF = astDF_og,
                       cohort_name = cohort_name,
                       time_frame = c(start_time_frame, end_time_frame))
   
   # Get raw demographic data and featurize.
   id_chunks <- get_id_chunks(df$PERSON_ID)
   demoDF <- tibble()
   cat('Demographics data.\n')
   for (i in seq_along(id_chunks)) {
      cat(i, '')
      demoDF <- rbind(
         demoDF,
         tbl(conn, in_schema('AMB_ETL', 'SENS_PATIENT_DEMO_VW')) %>%
            filter(PERSON_ID %in% local(id_chunks[[i]])) %>%
            collect()
      )
   }; cat('\n')
   df <- featurizeDemographics(df = df, 
                               demoDF = demoDF,
                               earliest_death_date = earliest_death_date)
   
   # Get raw encounters data and featurize.
   id_chunks <- get_id_chunks(df$PERSON_ID)
   encsDF <- tibble()
   cat('Encounters data.\n')
   for (i in seq_along(id_chunks)) {
      cat(i, '')
      encsDF <- rbind(
         encsDF,
         tbl(conn, in_schema('AMB_ETL', 'SENS_ENCOUNTER_VW')) %>%
            filter(
               PERSON_ID %in% local(id_chunks[[i]]),
               substr(ADMIT_DATE,7,10) %in% years
            ) %>%
            collect()
      )
   }; cat('\n')
   df <- featurizeEncounters(df = df, 
                             encsDF = encsDF,
                             preprocess = TRUE,
                             time_frame = encounters_time_frame)
   rm(site_info, code2cat, code2cat_fxn, code2fac, code2fac_fxn, code2mixed, code2mixed_fxn, fac2cat, fac2cat_fxn, mixed2cat, mixed2cat_fxn)
   
   # Get raw ICU stay data and featurize.
   id_chunks <- get_id_chunks(df %>% filter(ENC_icu_stay) %>% pull(PERSON_ID) %>% unique())
   icusDF <- tibble()
   cat('ICU stay data.\n')
   for (i in seq_along(id_chunks)) {
      cat(i, '')
      icusDF <- rbind(
         icusDF,
         tbl(conn, in_schema('AMB_ETL', 'SENS_ENCOUNTER_ICU_LOCS_VW')) %>%
            filter(
               PERSON_ID %in% local(id_chunks[[i]]),
               substr(ADMIT_DATE,7,10) %in% years,
               !is.na(ICU_EFFECTIVE_DATE)
            ) %>%
            collect()
      )
   }; cat('\n')
   df <- featurizeICUstays(df = df,
                           icusDF = icusDF,
                           preprocess = TRUE,
                           time_frame = c(icus_start_time_frame, end_time_frame))
   
   
   # Get raw medication administration data.
   load(file = '~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/abx_names_in_med_admin_vw.Rdata')
   load(file = '~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/abx_names_in_med_admin_iv_vw.Rdata')
   load(file = '~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/vasopressor_names_in_med_admin_vw.Rdata')
   load(file = '~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/vasopressor_names_in_med_admin_iv_vw.Rdata')
   
   id_chunks <- get_id_chunks(df$PERSON_ID)
   medsDF <- tibble()
   cat('Antibiotic administration.\n')
   for (i in seq_along(id_chunks)) {
      cat(i, '')
      medsDF <- rbind(
         medsDF,
         tbl(conn, in_schema('AMB_ETL', 'SENS_MED_ADMIN_VW')) %>%
            filter(
               PERSON_ID %in% local(id_chunks[[i]]),
               lubridate::year(ADMIN_START_DATE) %in% years,
               MEDICATION %in% local(abx_admin)
            ) %>%
            collect()
      )
   }; cat('\n')
   cat('Vasopressor administration.\n')
   for (i in seq_along(id_chunks)) {
      cat(i, '')
      medsDF <- rbind(
         medsDF,
         tbl(conn, in_schema('AMB_ETL', 'SENS_MED_ADMIN_VW')) %>%
            filter(
               PERSON_ID %in% local(id_chunks[[i]]),
               lubridate::year(ADMIN_START_DATE) %in% years,
               MEDICATION %in% local(vas_admin)
            ) %>%
            collect()
      )
   }; cat('\n')
   cat('Medication administration IV.\n')
   for (i in seq_along(id_chunks)) {
      cat(i, '')
      medsDF <- rbind(
         medsDF,
         tbl(conn, in_schema('AMB_ETL', 'SENS_MED_ADMIN_IV_VW')) %>%
            filter(
               PERSON_ID %in% local(id_chunks[[i]]),
               lubridate::year(ADMIN_START_DATE) %in% years,
               MEDICATION %in% local(c(abx_admin_iv, vas_admin_iv))
            ) %>%
            collect()
      )
   }; cat('\n')
   
   # Separate and featurize antibiotics and vasopressors.
   df <- featurizeAbxAdmin(df = df,
                           abxDF = medsDF %>% filter(MEDICATION %in% c(abx_admin, abx_admin_iv)),
                           preprocess = TRUE)
   df <- featurizeVasopressors(df = df,
                               vasDF = medsDF %>% filter(MEDICATION %in% c(vas_admin, vas_admin_iv)),
                               preprocess = TRUE,
                               time_frame = vasopressors_time_frame)
   
   
   # Get ICD codes for comorbidities.
   id_chunks <- get_id_chunks(df$PERSON_ID)
   icdsDF <- tibble()
   cat('ICD-10 code data.\n')
   for (i in seq_along(id_chunks)) {
      cat(i, '')
      icdsDF <- rbind(
         icdsDF,
         tbl(conn, in_schema('AMB_ETL', 'LAB_SENS_DX_VW')) %>%
            filter(
               PERSON_ID %in% local(id_chunks[[i]]),
               lubridate::year(DX_FROM_DATE) %in% local((min(years)-2L):max(years))
            ) %>%
            collect()
      )
   }; cat('\n')
   df <- featurizeICDs(df = df,
                       icdsDF = icdsDF,
                       end_time_frame = end_time_frame,
                       cohort_name = cohort_name,
                       preprocess = TRUE)
   
   # Get raw laboratory result values.
   id_chunks <- get_id_chunks(df$PERSON_ID)
   labsDF <- tibble()
   for (i in seq_along(id_chunks)) {
      cat(i, '- ')
      for (year in years) {
         cat(year, '')
         labsDF <- rbind(
            labsDF,
            tbl(conn, in_schema('AMB_ETL', paste0('SENS_LAB_RESULT_', year, '_VW'))) %>%
               filter(PERSON_ID %in% local(id_chunks[[i]])) %>%
               select(PERSON_ID, ORDER_DATE, COMPONENT_NAME, RESULT_VALUE, REFERENCE_UNIT) %>%
               collect()
         )
      }
      cat('\n')
   }
   df <- featurizeLabs(df = df, 
                       labsDF = labsDF,
                       preprocess = TRUE,
                       time_frame = c(lab_start_time_frame, end_time_frame),
                       missingness_threshold = lab_missingness_threshold)
   
   # Final output about treatment group counts.
   t <- sort(table(df$treatment))
   cat(names(t)[1], ': ', t[1], '\n', sep='')
   cat(names(t)[2], ': ', t[2], '\n', sep='')
   cat(Sys.time() - start)
   sink()
   
   # Re-order all variables
   df <- df %>%
      relocate(treatment, contains('OUT'), .after=PERSON_ID) %>%
      relocate(contains('DEM_'), .after=last_col()) %>%
      relocate(contains('ABX_'), .after=last_col()) %>%
      relocate(contains('ENC_'), .after=last_col()) %>%
      relocate(contains('LAB_'), .after=last_col()) %>%
      relocate(contains('ICD_'), .after=last_col()) %>%
      select(-PERSON_ID, -ORDER_DATE, -ORDER_DAY)
   
   # Remove variables with no variation
   remove_vars <- names(which(sapply(df, function(x) length(unique(x))) == 1))
   if (length(remove_vars) > 0L)
      df <- df %>% select(-c(!!remove_vars))
   
   
   # Save data
   save(df, file = paste0(clean_data_path, cohort_name, '.Rdata'))
}

print(Sys.time() - START)












