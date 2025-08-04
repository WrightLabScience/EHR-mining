# This script will get the vancmoycin MIC values reported in the raw susceptibility data
# In this case, only for MRSA bacteremia cases
library(dplyr)

# Load raw susceptibility data.
load(file = '~/Desktop/EHR/EHR work/RdataFiles/lab_micro_sens_all_vw.Rdata')

# Get the MRSA bacteremia data.
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_clean_ASTs.Rdata')
astDF <- astDF %>%
   filter(
      BLOOD,
      BUG == 'Staphylococcus aureus',
      OXACILLIN == 1L
   ) %>%
   select(PERSON_ID, ORDER_DAY, RESULT_DAY) %>%
   distinct()

# Minimal filtering of raw data to facilitate join speed.
astrDF <- astrDF %>%
   filter(
      PERSON_ID %in% unique(astDF$PERSON_ID),
      ANITBIOTIC_NAME == 'VANCOMYCIN'
   )
astrDF <- astrDF %>%
   mutate(RESULT_DAY = as.Date(RESULT_DATE, format='%m/%d/%Y'))

astrDF <- astrDF_og %>% 
   mutate(ORGANISM_NAME = gsub(' ', '', ORGANISM_NAME)) %>%
   filter(
      grepl('MRSA|STAPH|AUREUS', ORGANISM_NAME)
   )

# Execute join.
micDF <- left_join(
   x = astDF,
   y = astrDF %>% select(PERSON_ID, RESULT_DAY, ORGANISM_NAME, SUSCEPTIBILITY_NAME, SENSITIVITY_VALUE) %>% distinct(),
   by = join_by(PERSON_ID, RESULT_DAY),
   relationship = 'many-to-many'
) %>%
   filter(lubridate::year(ORDER_DAY) %in% 2015:2024)

micDF$SENS_VAL <- gsub('I +| +Intermediate| +Sensitive|S ?|\\(|\\)|<=', '', micDF$SENSITIVITY_VALUE)
micDF$SENS_VAL[which(micDF$SENS_VAL == 'Not Done')] <- NA
micDF$VAN_MIC <- as.numeric(micDF$SENS_VAL)

micDF <- micDF %>%
   select(!c(RESULT_DAY, ORGANISM_NAME, SUSCEPTIBILITY_NAME, SENSITIVITY_VALUE, SENS_VAL)) %>%
   distinct() %>%
   group_by(PERSON_ID, ORDER_DAY) %>%
   slice_max(VAN_MIC) %>%
   ungroup()


save(micDF, file='~/Desktop/EHR/EHR work/RdataFiles/BuildCohortDatasets/RawData/mrsa_vancomycin_MIC_values_2015_2024.Rdata')

left_join(
   x = micDF_ %>% filter(is.na(VAN_MIC)) %>% select(-VAN_MIC),
   y = micDF,
   by = join_by(PERSON_ID, ORDER_DAY)
) %>%
   count(SUSCEPTIBILITY_NAME, SENSITIVITY_VALUE, SENS_VAL, VAN_MIC)





# Get the MRSA bacteremia data.
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_clean_ASTs.Rdata')
astDF <- astDF %>%
   filter(
      SITE == 'blood',
      BUG == 'Enterococcus faecalis'
   ) %>%
   select(PERSON_ID, ORDER_DAY, RESULT_DAY, VANCOMYCIN) %>%
   distinct()

astDF$DAYS_SINCE_PRV <- astDF$ORDER_DAY - lag(astDF$ORDER_DAY)
astDF$DAYS_SINCE_PRV[c(1, which(astDF$PERSON_ID != lag(astDF$PERSON_ID)))] <- NA_integer_
astDF <- astDF %>%
   filter(is.na(DAYS_SINCE_PRV) | DAYS_SINCE_PRV > 60L) %>%
   #mutate(ENC_bacteremia_1yr = !is.na(DAYS_SINCE_PRV) & DAYS_SINCE_PRV <= 365L) %>%
   select(-DAYS_SINCE_PRV)

# Minimal filtering of raw data to facilitate join speed.
astrDF <- astrDF %>%
   filter(
      PERSON_ID %in% unique(astDF$PERSON_ID),
      ANITBIOTIC_NAME == 'VANCOMYCIN'
   )
astrDF <- astrDF %>%
   mutate(RESULT_DAY = as.Date(RESULT_DATE, format='%m/%d/%Y'))
astrDF <- astrDF %>%
   mutate(ORGANISM_NAME = gsub(' ', '', ORGANISM_NAME)) %>%
   filter(grepl('ENTEROCOCCUSFAECALIS', ORGANISM_NAME))

mics <- astDF %>%
   filter(VANCOMYCIN == 0L) %>%
   left_join(
      x = .,
      y = astrDF %>% 
         select(PERSON_ID, RESULT_DAY, ORGANISM_NAME, SUSCEPTIBILITY_NAME, SENSITIVITY_VALUE, SENSITIVITY_COMMENT) %>% 
         distinct(),
      by = join_by(
         PERSON_ID, RESULT_DAY
      ),
      relationship = 'one-to-many'
   ) %>%
   mutate(VAN_MIC <- as.numeric(gsub('[^0-9.]', '', mics$SENSITIVITY_VALUE)))

mics <- mics %>%
   select(PERSON_ID, ORDER_DAY, ) %>%


mics$SENSITIVITY_VALUE <- as.numeric(gsub('[^0-9.]', '', mics$SENSITIVITY_VALUE))
table(mics$SENSITIVITY_VALUE)

















