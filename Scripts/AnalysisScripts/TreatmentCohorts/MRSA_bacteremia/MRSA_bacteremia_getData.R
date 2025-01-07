library(dplyr)
source('~/Desktop/WrightLab/UsefulRfuncs/newBarplotFxn.R')


##### GET MRSA BSI ASTS #####
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_clean_ASTs.Rdata')
ast <- astDF %>%
   filter(BUG == 'Staphylococcus aureus',
          OXACILLIN == 1L, 
          BLOOD,
          lubridate::year(ORDER_DAY) %in% 2017:2023) %>%
   summarise(
      ORDER_DATE = min(ORDER_DATE),
      RESULT_DATE = min(RESULT_DATE),
      .by = c(PERSON_ID, ORDER_DAY)
   ) %>%
   mutate(RESULT_DAY = lubridate::as_date(RESULT_DATE)) %>%
   filter(as.integer(RESULT_DAY - ORDER_DAY) <= 10L)

length(unique(ast$PERSON_ID)) # 3,468

ast <- ast %>%
   group_by(PERSON_ID) %>%
   mutate(DAYS_SINCE_PRV = as.integer(ORDER_DAY - lag(RESULT_DAY))) %>%
   ungroup()

sum(ast$DAYS_SINCE_PRV <= 0, na.rm=T) / nrow(ast) # 27%
sum(ast$DAYS_SINCE_PRV <= 30, na.rm=T) / nrow(ast) # 35.5%

plotBarplot(ast$DAYS_SINCE_PRV[ast$DAYS_SINCE_PRV < 90])

ast <- ast %>% filter(is.na(DAYS_SINCE_PRV) | DAYS_SINCE_PRV >= 30L)
ast %>% count(is.na(DAYS_SINCE_PRV)) # only 10% have any other MRSA BSI

resp <- astDF %>%
   filter(BUG == 'Staphylococcus aureus',
          OXACILLIN == 1L, 
          RESPIRATORY,
          lubridate::year(ORDER_DAY) %in% 2017:2023) %>%
   select(PERSON_ID, ORDER_DAY) %>%
   distinct()
ast <- ast %>%
   left_join(y=resp %>% 
                mutate(JOIN_START = ORDER_DAY - 10, 
                       JOIN_END = ORDER_DAY + 4) %>% 
                rename(resp_O_DAY = ORDER_DAY), 
             by=join_by(PERSON_ID, between(ORDER_DAY, JOIN_START, JOIN_END))) %>%
   mutate(RESP_CULT_DAY = as.integer(ORDER_DAY - resp_O_DAY)) %>%
   slice_min(resp_O_DAY, by = c(PERSON_ID, ORDER_DAY)) %>%
   select(-resp_O_DAY, -JOIN_START, -JOIN_END)
ast %>% pull(RESP_CULT_DAY) %>% plotBarplot()
ast %>% count(RESP_CULT_DAY)
rm(astDF, resp); gc()
##### END #####


##### JOIN DEMO #####
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_DEMO.Rdata')
censor_time <- 30
ast <- ast %>%
   left_join(y = dth, by = join_by(PERSON_ID)) %>%
   mutate(ORDER_DAY = as.Date(substr(ORDER_DATE,1,10))) %>%
   mutate(
      time = as.integer(DEATH_DATE - ORDER_DAY),
      AGE = as.integer(ORDER_DAY - DOB) / 365
   ) %>%
   mutate(status = case_when(
      is.na(time) | time > censor_time ~ 0,
      .default = 1
   )) %>%
   mutate(time_censored = case_when(
      is.na(time) | time > censor_time ~ censor_time,
      .default = time
   )) %>%
   select(!c(PATIENT_STATUS, DEATH_DATE, DOB))
ast %>% count(time)
rm(censor_time, dth); gc()
##### END #####


##### JOIN ENCOUNTERS #####
if (FALSE) {
   source('~/Desktop/EHR/EHR work/config_file.R')
   ids <- unique(ast$PERSON_ID)
   chunks <- mapply(FUN = ':',
                    seq(1, floor(length(ids) / 1000) * 1000 + 1, 1000),
                    c(seq(1000, floor(length(ids) / 1000) * 1000, 1000), length(ids)))
   encs_og <- tibble()
   for (chunk in chunks) {
      print(chunk[1])
      encs_og <- rbind(encs_og,
                        tbl(conn, in_schema('AMB_ETL', 'SENS_ENCOUNTER_VW')) %>%
                           filter(PERSON_ID %in% local(ids[chunk])) %>%
                           collect())
   }
   encs_og <- encs_og %>%
      mutate(across(.cols = c(ADMIT_DATE, DISCHARGE_DATE),
                    .fns = ~ strptime(., format='%m/%d/%Y %T'))) %>%
      mutate(ADMIT_DAY = lubridate::as_date(ADMIT_DATE), 
             DISCHARGE_DAY = lubridate::as_date(DISCHARGE_DATE))
   source('~/Desktop/EHR/EHR-mining/Scripts/CleaningScripts/CombineOverlappingEncountersFxn.R')
   encs <- combineOverlappingEncounters(encs_og)
   
   save(encs, file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/MRSA_bacteremia/EncountersCleaned.Rdata')
   save(encs_og, file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/MRSA_bacteremia/EncountersRaw.Rdata')
   
   rm(combineOverlappingEncounters, chunk, chunks, ids)
}

load(file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/MRSA_bacteremia/EncountersCleaned.Rdata')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/MRSA_bacteremia/EncountersRaw.Rdata')
###

# join with asts
df <- ast %>%
   left_join(y = encs %>% mutate(BEFORE_ADMIT = ADMIT_DAY - 1),
             by = join_by(
                PERSON_ID,
                between(ORDER_DAY, BEFORE_ADMIT, DISCHARGE_DAY)
             )) %>%
   filter(!is.na(ADMIT_DAY)) %>% # lost 423
   mutate(TIME_BW_ADMIT_ORDER = as.numeric(lubridate::as.duration(ORDER_DATE - ADMIT_DATE)) / 86400,
          DAYS_BW_ADMIT_ORDER = as.integer(ORDER_DAY - ADMIT_DAY)) %>% #pull(DAYS_BW_ADMIT_ORDER) %>% plotBarplot()
   mutate(HospAcq = TIME_BW_ADMIT_ORDER >= 2,
          EARLY_CULTURE = TIME_BW_ADMIT_ORDER < -0.5 & DAYS_BW_ADMIT_ORDER < 0) %>%
   select(-BEFORE_ADMIT, -DAYS_BW_ADMIT_ORDER, -ADMIT_DATE, -DISCHARGE_DATE, -ADMIT_DAY, -DISCHARGE_DAY)

# some overlapping (or nearly so) encounters were combined and may have multiple facilities
df <- df %>%
   rename(FACILITIES = FACILITY) %>%
   left_join(y = encs_og %>% 
                select(PERSON_ID, ADMIT_DAY, DISCHARGE_DAY, FACILITY) %>% 
                mutate(BEFORE_ADMIT = ADMIT_DAY - 1), 
             by = join_by(PERSON_ID, 
                          between(ORDER_DAY, BEFORE_ADMIT, DISCHARGE_DAY))) %>%
   distinct() %>%
   group_by(PERSON_ID, ORDER_DAY) %>%
   slice_max(DISCHARGE_DAY) %>%
   ungroup() %>%
   select(-ADMIT_DAY, -DISCHARGE_DAY, -BEFORE_ADMIT)
any(is.na(df$FACILITY))


# hospital-acquired MRSA bacteremia = death-sentence!
df %>% summarise(n=n(), sum(time < 30, na.rm=T) / n, .by=HospAcq)
df %>% summarise(n=n(), sum(time < 4, na.rm=T) / n, .by=HospAcq)
df %>% filter(is.na(time) | time > 3) %>% summarise(n=n(), sum(time < 30, na.rm=T) / n, .by=HospAcq)
sum(df$HospAcq, na.rm=T) / nrow(df) # 16% hospital-acquired

# admit to ED first?
df %>% count(ENCOUNTER_TYPE, sort=TRUE) # basically all inpatients, some from "maternity" or "nursery"
df %>% filter(grepl('maternity|nursery', ENCOUNTER_TYPE, ignore.case=TRUE)) %>% count(ENCOUNTER_TYPE)
df %>% count(ADMIT_TYPE, sort=TRUE)
data.frame(sort(table(unlist(strsplit(df$ADMIT_TYPE, ', '))), decreasing=TRUE))
data.frame(sort(table(unlist(strsplit(df$ADMIT_SOURCE, ', '))), decreasing=TRUE))

df$NURSING_HOME <- grepl('Nursing', df$ADMIT_SOURCE) # much higher mortality rate, nearly double 30-day
df %>% summarise(n = n(), sum(time < 30, na.rm=T) / n, .by=NURSING_HOME)
df$EMERGENCY_DEPT <- grepl('ED|Emergency|Trauma|Urgent', df$ADMIT_TYPE) # 4% higher 30-day mortality rate
df %>% summarise(n = n(), sum(time < 30, na.rm=T) / n, .by=EMERGENCY_DEPT)

df %>%
   summarise(n = n(), 
             ec = sum(EARLY_CULTURE) / n,
             ha = sum(HospAcq) / n,
             nh = sum(NURSING_HOME) / n,
             .by=FACILITY) %>%
   filter(n > 50L) %>%
   arrange(desc(ec)) # UPMCWIL is doube/triple most other facilities

# remove unnecessary variables
ast <- df %>% select(PERSON_ID, ORDER_DATE, RESULT_DATE, ORDER_DAY, RESULT_DAY, DAYS_SINCE_PRV, GENDER, time, time_censored, 
                    status, AGE, TIME_BW_ADMIT_ORDER, FACILITY, HospAcq, EARLY_CULTURE, NURSING_HOME, EMERGENCY_DEPT, RESP_CULT_DAY)
##### END #####



save(ast, file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/MRSA_bacteremia/DF_before_AbxAdmin.Rdata')
##########################################################################################################
library(dplyr)
source('~/Desktop/WrightLab/UsefulRfuncs/newBarplotFxn.R')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/MRSA_bacteremia/DF_before_AbxAdmin.Rdata')



##### JOIN ABX ADMIN #####
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_CLEANED_2017_2023_AbxAdmin.Rdata')
abxDF <- abxDF %>% 
   filter(PERSON_ID %in% unique(ast$PERSON_ID),
          ABX %in% c('VANCOMYCIN', 'DAPTOMYCIN')) %>% # 108K
   select(-END_DATE) %>%
   distinct()

# join ASTs + AbxAdmin
ast_abx <- ast %>%
   mutate(JOIN_START = ORDER_DAY - 30,
          JOIN_END = ORDER_DAY + 30) %>%
   left_join(y = abxDF,
             by = join_by(
                PERSON_ID,
                JOIN_START <= START_DAY,
                JOIN_END >= START_DAY
             )) %>%
   mutate(X = as.integer(START_DAY - ORDER_DAY),
          XT = as.numeric(lubridate::as.duration(START_DATE - ORDER_DATE)) / 86400)# %>% pull(XT) %>% plotBarplot()

# remove cases without enough VAN or DAP before AST: 3,407
ast_abx <- ast_abx %>%
   group_by(PERSON_ID, ORDER_DAY) %>%
   mutate(N7 = sum(unique(X) %in% 4:7),
          N3 = sum(unique(X) %in% 0:3)) %>%
   filter(N3 >= 2L) %>% # 3,407 --> 2,883
   ungroup()

# remove them
ast_abx <- ast_abx %>% group_by(PERSON_ID, ORDER_DAY) %>% filter(!any(X %in% -3:-14)) %>% ungroup() # 2,774


# Define cohorts based on what happened days 0, 1, 2, and 3 in terms of treatment
# just vancomycin in that period? van group
# just daptomycin in that period? dap group
# both? early switch? concurrent? d3dap group (or exclude)
# when is the first day of abx?
ast_abx <- ast_abx %>%
   arrange(PERSON_ID, ORDER_DAY, START_DATE, desc(ABX)) %>%
   group_by(PERSON_ID, ORDER_DAY) %>%
   mutate(ALL_DAP_TIMES = list(XT[X %in% -1:3 & ABX == 'DAPTOMYCIN']),
          ALL_VAN_TIMES = list(XT[X %in% -1:3 & ABX == 'VANCOMYCIN']),,
          FIRST_DAP_DAY = list(X[X %in% -1:3 & ABX == 'DAPTOMYCIN']),
          FIRST_VAN_DAY = list(X[X %in% -1:3 & ABX == 'VANCOMYCIN']),
          NUM_VAN_DAYS = sum(ABX[X %in% -1:3] == 'VANCOMYCIN'),
          NUM_DAP_DAYS = sum(ABX[X %in% -1:3] == 'DAPTOMYCIN'),
          NUM_ABX_DAYS = sum(unique(X) %in% -1:3)) %>%
   ungroup() %>%
   mutate(NUM_VAN_DOSES = lengths(ALL_VAN_TIMES),
          NUM_DAP_DOSES = lengths(ALL_DAP_TIMES)) %>%
   summarise(FIRST_ABX_DAY = min(X[X %in% -1:3]),
             FIRST_ABX_TIME = min(XT[X %in% -1:3]),
             TRT = paste(unique(ABX[X %in% -1:3]), collapse=','),
             TRT7 = paste(unique(ABX[X %in% -1:7]), collapse=','),
             .by = c(PERSON_ID, ORDER_DATE, RESULT_DATE, ORDER_DAY, RESULT_DAY, N3, N7, FIRST_DAP_DAY, FIRST_VAN_DAY, NUM_VAN_DOSES, NUM_DAP_DOSES, 
                     ALL_DAP_TIMES, ALL_VAN_TIMES, NUM_VAN_DAYS, NUM_DAP_DAYS, NUM_ABX_DAYS, GENDER, time, AGE, status, time_censored,
                     FACILITY, DAYS_SINCE_PRV, HospAcq, EARLY_CULTURE, NURSING_HOME, EMERGENCY_DEPT, RESP_CULT_DAY, TIME_BW_ADMIT_ORDER))
ast_abx$FIRST_DAP_DAY  <- sapply(ast_abx$FIRST_DAP_DAY, function(l) ifelse(length(l) == 0L, NA, min(l)))
ast_abx$FIRST_VAN_DAY  <- sapply(ast_abx$FIRST_VAN_DAY, function(l) ifelse(length(l) == 0L, NA, min(l)))
ast_abx$FIRST_DAP_TIME <- sapply(ast_abx$ALL_DAP_TIMES, function(l) ifelse(length(l) == 0L, NA, min(l)))
ast_abx$FIRST_VAN_TIME <- sapply(ast_abx$ALL_VAN_TIMES, function(l) ifelse(length(l) == 0L, NA, min(l)))

# distribution of treatment strategies
ast_abx %>% count(TRT)
ast_abx %>% count(TRT7)
ast_abx %>% count(TRT, TRT7)

ast_abx <- ast_abx %>% filter(TRT != 'DAPTOMYCIN,VANCOMYCIN') # 2,723

# put restrictions on this
ast_abx %>% filter(TRT == 'VANCOMYCIN,DAPTOMYCIN') %>% count(FIRST_VAN_DAY, FIRST_DAP_DAY, sort=TRUE)
# how much time before switch
ast_abx %>% filter(TRT == 'VANCOMYCIN,DAPTOMYCIN') %>% mutate(X = FIRST_DAP_TIME - FIRST_VAN_TIME) %>% pull(X) %>% hist(breaks=24*4)


# how much overlap in time spent on VAN and DAP concurrently? and is the final VAN day LATER than the DAP day?
df <- ast_abx %>% filter(TRT == 'VANCOMYCIN,DAPTOMYCIN') %>% select(ALL_VAN_TIMES, ALL_DAP_TIMES)
keep <- logical(nrow(df))
dap_switch_time_r <- dap_switch_time_a <- rep(NA, nrow(df))
for (i in seq_len(nrow(df))) {
   van <- df$ALL_VAN_TIMES[[i]]
   dap <- df$ALL_DAP_TIMES[[i]]
   overlap <- min(c(max(van), max(dap))) - max(c(min(van), min(dap)))
   if (max(dap) >= max(van) & overlap <= 1) {
      keep[i] <- TRUE
      dap_switch_time_r[i] <- min(dap) - min(van)
      dap_switch_time_a[i] <- min(dap)
   }
}
w <- which(ast_abx$TRT == 'VANCOMYCIN,DAPTOMYCIN')
ast_abx$DAP_SWITCH_TIME_A <- ast_abx$DAP_SWITCH_TIME_R <- NA
ast_abx$DAP_SWITCH_TIME_A[w] <- dap_switch_time_a
ast_abx$DAP_SWITCH_TIME_R[w] <- dap_switch_time_r
ast_abx <- ast_abx[-w[!keep], ]
rm(df, keep, i, van, dap, overlap, w, dap_switch_time_a, dap_switch_time_r)

par(mfrow=c(2,2))
ast_abx %>% filter(TRT == 'VANCOMYCIN,DAPTOMYCIN') %>% pull(FIRST_VAN_TIME) %>% hist(xlim=c(-1.5, 4), main='First VAN time')
ast_abx %>% filter(TRT == 'VANCOMYCIN,DAPTOMYCIN') %>% pull(FIRST_DAP_TIME) %>% hist(xlim=c(-1.5, 4), main='First DAP time')
ast_abx %>% filter(TRT == 'VANCOMYCIN,DAPTOMYCIN') %>% pull(DAP_SWITCH_TIME_A) %>% hist(xlim=c(-1.5, 4), main='Relative to culture')
ast_abx %>% filter(TRT == 'VANCOMYCIN,DAPTOMYCIN') %>% pull(DAP_SWITCH_TIME_R) %>% hist(xlim=c(-1.5, 4), main='Relative to first VAN')

ast_abx %>%
   filter(TRT == 'VANCOMYCIN,DAPTOMYCIN') %>%
   count(DAP_SWITCH_TIME_R <= 1.5, DAP_SWITCH_TIME_A <= 1.5)
ast_abx <- ast_abx %>%
   rename(DAP_SWITCH_TIME = DAP_SWITCH_TIME_R) %>%
   select(-DAP_SWITCH_TIME_A)


ast_abx %>% summarise(n = n(), 
                      ndth = sum(time_censored < 30), 
                      resp = sum(!is.na(RESP_CULT_DAY)) / n,
                      EARLY_CULTURE = sum(EARLY_CULTURE, na.rm=T) / n, 
                      ED = sum(EMERGENCY_DEPT, na.rm=T) / n,
                      HOSP_ACQ = sum(HospAcq, na.rm=T) / n, 
                      NURSE_HOME = sum(NURSING_HOME) / n,
                      FEMALE = sum(GENDER == 'FEMALE') / n,
                      FIRST_ABX_DAY_2_orLater = sum(FIRST_ABX_DAY >= 2) / n,
                      fdth = ndth / n,
                      mean(FIRST_ABX_TIME),
                      .by=TRT)

ast_abx %>% 
   mutate(RESP = !is.na(RESP_CULT_DAY)) %>%
   summarise(n = n(),
             d30 = sum(time_censored < 30) / n, 
             .by=c(RESP, TRT)) %>%
   arrange(TRT, RESP)

ast_abx <- ast_abx %>%
   select(PERSON_ID, ORDER_DAY, ORDER_DATE, FIRST_ABX_DAY, DAP_SWITCH_TIME, FIRST_VAN_TIME, FIRST_DAP_TIME, TRT, GENDER, AGE, status, time, time_censored, FACILITY,
          DAYS_SINCE_PRV, HospAcq, EARLY_CULTURE, NURSING_HOME, EMERGENCY_DEPT, RESP_CULT_DAY, TIME_BW_ADMIT_ORDER)
##### END #####


##### JOIN DX_CODES #####
df <- ast_abx
if (FALSE) {
   source('~/Desktop/EHR/EHR work/config_file.R')
   ids <- unique(df$PERSON_ID)
   chunks <- mapply(':', seq(1, floor(length(ids) / 1000) * 1000 + 1, 1000), c(seq(1000, floor(length(ids) / 1000) * 1000, 1000), length(ids)))
   dx <- tibble()
   for (chunk in chunks) {
      print(chunk[1])
      dx <- rbind(
         dx,
         tbl(conn, in_schema('AMB_ETL', 'LAB_SENS_DX_VW')) %>% 
            filter(PERSON_ID %in% local(ids[chunk]),
                   lubridate::year(DX_FROM_DATE) %in% 2017:2023) %>%
            collect()
      )
   }
   dx <- dx %>%
      mutate(DX_DATE = lubridate::date(DX_FROM_DATE),
             CODE_DESCRIPTION = tolower(CODE_DESCRIPTION)) %>% 
      select(PERSON_ID, DX_DATE, DX_CODE, CODE_DESCRIPTION) %>% 
      distinct()
   rm(result, chunks, chunk, ids)
   save(dx, file='~/Desktop/EHR/EHR work/RdataFiles/causal_prep/MRSA_bacteremia/Diagnoses.Rdata')
}

## Load diagnoses and join with data
load(file='~/Desktop/EHR/EHR work/RdataFiles/causal_prep/MRSA_bacteremia/Diagnoses.Rdata')
# how many actually had diagnosis codes during their stay?
# use this as an indication for whether they didn't have diagnosis vs. missing data...
length(unique(df$PERSON_ID))
length(unique(dx$PERSON_ID))
length(intersect(dx$PERSON_ID, df$PERSON_ID)) # have codes for everyone AT SOME POINT


# join and re-format 1 row per infection, flag for each comorbidity
#### WITHIN WEEK BEFORE BLOOD CULTURE
dfx <- df %>%
   mutate(JOIN_START = ORDER_DAY + FIRST_ABX_DAY - 7,
          JOIN_END = ORDER_DAY + FIRST_ABX_DAY) %>%
   left_join(y = dx,
             by = join_by(
                PERSON_ID,
                JOIN_START <= DX_DATE,
                JOIN_END >= DX_DATE
             )) %>%
   select(-JOIN_START, -JOIN_END) %>%
   mutate(X = as.integer(DX_DATE - ORDER_DAY),
          CODE_DESCRIPTION = ifelse(is.na(CODE_DESCRIPTION), 'none', CODE_DESCRIPTION)) %>%
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
   select(-DX_DATE, -DX_CODE, -CODE_DESCRIPTION, -X) %>%
   distinct()
dfx %>% summarise(n(), across(SepticShock_1w:Respiratory_1w, ~ sum(.) / n()), .by=TRT)


#### WITHIN MONTH OF BLOOD CULTURE
dfx <- dfx %>%
   mutate(JOIN_START = ORDER_DAY - 30,
          JOIN_END = (ORDER_DAY + FIRST_ABX_DAY)) %>%
   left_join(y = dx,
             by = join_by(
                PERSON_ID,
                JOIN_START <= DX_DATE,
                JOIN_END >= DX_DATE
             )) %>%
   select(-JOIN_START, -JOIN_END) %>%
   mutate(X = as.integer(DX_DATE - ORDER_DAY),
          CODE_DESCRIPTION = ifelse(is.na(CODE_DESCRIPTION), 'none', CODE_DESCRIPTION)) %>%
   group_by(PERSON_ID, ORDER_DAY) %>%
   mutate(
      PulmCircDis_1m = any(grepl('^I26|^I27|^I28\\.[089]', DX_CODE)),
      CPD_1m = any(grepl('^I27\\.28|^I27\\.9|J4[0-7]|^J6[0-7]|^J68\\.4|^J70\\.[13]', DX_CODE)),
      MyocInfarc_1m = any(grepl('^I21|^I22|^I25', DX_CODE)),
      MetastSolidTumor_1m = any(grepl('^C7[7-9]|^C80', DX_CODE)),
      Malignancy_1m = any(grepl('^C[01][0-9]|^C2[0-6]|^C3[01234789]|^C4[013]|^C4[5-9]|^C5[0-8]|^C6[0-9]|^C7[0-6]|^C8[123458]|^C9[0-7]', DX_CODE)),
   ) %>%
   ungroup() %>%
   select(-DX_DATE, -DX_CODE, -CODE_DESCRIPTION, -X) %>%
   distinct()
dfx %>% summarise(across(PulmCircDis_1m:Malignancy_1m, ~ sum(.) / n()), .by=TRT)


#### WITHIN 2 YEARS OF BLOOD CULTURE
dfx <- dfx %>%
   mutate(JOIN_START = ORDER_DAY - 730,
          JOIN_END = ORDER_DAY + FIRST_ABX_DAY) %>%
   left_join(y = dx,
             by = join_by(
                PERSON_ID,
                JOIN_START <= DX_DATE,
                JOIN_END >= DX_DATE
             )) %>%
   select(-JOIN_START, -JOIN_END) %>%
   mutate(X = as.integer(DX_DATE - ORDER_DAY),
          CODE_DESCRIPTION = ifelse(is.na(CODE_DESCRIPTION), 'none', CODE_DESCRIPTION)) %>%
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
   select(-DX_DATE, -DX_CODE, -CODE_DESCRIPTION, -X) %>%
   distinct()
dfx %>% summarise(across(OsteoChronic:Smoking, ~ sum(.) / n()), .by=TRT)
##### END #####


##### JOIN OTHER ISOLATES #####
load(file='~/Desktop/EHR/EHR work/RdataFiles/ALL_clean_ASTs.Rdata')
otherBugs <- astDF %>%
   filter(lubridate::year(ORDER_DAY) %in% 2017:2023) %>%
   mutate(MRSA_BSI = BUG == 'Staphylococcus aureus' & OXACILLIN == 1L & BLOOD) %>%
   filter(!MRSA_BSI) %>%
   select(PERSON_ID, ORDER_DAY, BUG, ESBL, VANCOMYCIN) %>%
   distinct() %>%
   rename(oORDER_DAY = ORDER_DAY)

dfx <- dfx %>%
   mutate(JOIN_START = ORDER_DAY - 7,
          JOIN_END = ORDER_DAY + FIRST_ABX_DAY) %>%
   left_join(y = otherBugs,
             by = join_by(
                PERSON_ID,
                JOIN_START <= oORDER_DAY,
                JOIN_END >= oORDER_DAY
             )) %>%
   group_by(PERSON_ID, ORDER_DAY) %>%
   mutate(NumOtherIsolates = lengths(list(BUG)),
          VRE = any(grepl('Enterococcus', BUG) & VANCOMYCIN == 1L),
          MSSA = any(BUG == 'Staphylococcus aureus'),
          ESBLbug = any(ESBL == 1L)) %>%
   ungroup() %>%
   mutate(across(c(VRE, MSSA, ESBLbug), ~ ifelse(is.na(.), FALSE, .))) %>%
   select(-ESBL, -VANCOMYCIN, -BUG, -JOIN_START, -JOIN_END, -oORDER_DAY) %>%
   rename(ESBL = ESBLbug) %>%
   distinct()
rm(astDF, otherBugs); gc()
##### END #####


##### JOIN VAN MICs #####
if (FALSE) {
   load(file = '~/Desktop/EHR/EHR work/RdataFiles/lab_micro_sens_all_vw.Rdata')
   vanDF <- astrDF %>%
      filter(grepl('OXACILLIN|VANCOMYCIN', ANITBIOTIC_NAME)) %>%
      filter(PERSON_ID %in% unique(dfx$PERSON_ID)) %>%
      mutate(BUG = tolower(ORGANISM_NAME)) %>%
      mutate(BUG = gsub(' ', '', BUG)) %>%
      mutate(STATUS = case_when(
         grepl('SENSITIVE|NEGATIVE', SUSCEPTIBILITY_NAME) | grepl('Sensitive', SENSITIVITY_VALUE) ~ 0,
         grepl('INTERMEDIATE', SUSCEPTIBILITY_NAME) | grepl('Intermediate', SENSITIVITY_VALUE) ~ 1,
         grepl('RESISTANT|POSITIVE', SUSCEPTIBILITY_NAME) | grepl('Resistant', SENSITIVITY_VALUE) ~ 2
      )) %>%
      filter(grepl('MRSA|staph(ylococcus)aureus', BUG)) %>%
      filter(ANITBIOTIC_NAME == 'VANCOMYCIN' | (ANITBIOTIC_NAME == 'OXACILLIN' & STATUS == 2L)) %>%
      mutate(RESULT_DATE = as.Date(RESULT_DATE, format='%m/%d/%Y')) %>%
      filter(ANITBIOTIC_NAME == 'VANCOMYCIN') %>%
      select(PERSON_ID, RESULT_DATE, SENSITIVITY_VALUE, BUG, ORGANISM_NAME) %>%
      distinct()
   save(vanDF, file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/MRSA_bacteremia/vancomycin_MIC_data.Rdata')
}
load(file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/MRSA_bacteremia/vancomycin_MIC_data.Rdata')
micDF <- vanDF %>%
   right_join(
      x = .,
      y = dfx %>% 
         select(PERSON_ID, ORDER_DAY) %>%
         mutate(JOIN_START = ORDER_DAY + 1, JOIN_END = ORDER_DAY + 10),
      by = join_by(PERSON_ID, between(RESULT_DATE, JOIN_START, JOIN_END))
   ) %>%
   select(-JOIN_START, -JOIN_END) %>%
   mutate(MOD_SENS_VALUE = gsub(' +Sensitive|S |\\(|\\)|<=', '', SENSITIVITY_VALUE)) %>% #count(SENSITIVITY_VALUE) %>% print(n=25)
   mutate(VAN_MIC = case_when(
      MOD_SENS_VALUE == 'S' ~ NA,
      .default = as.numeric(MOD_SENS_VALUE)
   )) %>%
   select(PERSON_ID, ORDER_DAY, VAN_MIC) %>%
   arrange(PERSON_ID, ORDER_DAY, VAN_MIC) %>%
   distinct() %>%
   group_by(PERSON_ID, ORDER_DAY) %>%
   slice_max(VAN_MIC, with_ties = FALSE) %>%
   ungroup() %>% 
   mutate(VAN_MIC_2 = VAN_MIC >= 2)

dfx <- dfx %>% left_join(micDF, by=join_by(PERSON_ID, ORDER_DAY))

rm(astrDF, micDF, vanDF); gc()
##### END #####


# death rates with and without respiratory cultures?
dfx %>%
   mutate(RESP = !is.na(RESP_CULT_DAY)) %>%
   summarise(n = n(),
             d30 = sum(time_censored < 30) / n(),
             .by = c(TRT, RESP)) %>%
   tidyr::pivot_wider(id_cols=TRT, values_from=d30, names_from=RESP)
dfx <- dfx %>% filter(is.na(RESP_CULT_DAY)) %>% select(-RESP_CULT_DAY)


save(dfx, file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/MRSA_bacteremia/MRSA_BSI_DAP_VAN.Rdata')
# write.table(x = dfx,
#             file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/MRSA_bacteremia/MRSAbactAllAbx.csv',
#             quote = FALSE,
#             row.names = FALSE,
#             sep = '\t')
# load(file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/MRSA_bacteremia/MRSA_BSI_DAP_VAN.Rdata')









####################################################################
###################### PLOT TESTING GROUNDS ########################
####################################################################
library(dplyr)
library(survival)
load(file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/MRSA_bacteremia/MRSA_BSI_DAP_VAN.Rdata')
source('~/Desktop/EHR/EHR-mining/Scripts/AnalysisScripts/KaplanMeierCurveFxns.R')
source('~/Desktop/WrightLab/UsefulRfuncs/newBarplotFxn.R')
source('~/Desktop/EHR/EHR-mining/Scripts/AnalysisScripts/DrawTableOneFxn.R')
source('~/Desktop/EHR/EHR-mining/Scripts/AnalysisScripts/ProcessTableOneFxn.R')

cex <- 1.5
par(cex.main=cex, cex.axis=cex, cex.lab=cex, cex=cex)

site_names <- c("CHP" = 'Childrens', "UPMCALT" = 'Altoona',  "UPMCBED" = 'Bedford', 
                "UPMCCHA" = 'Chatauqua', "UPMCEAS" = 'East', "UPMCHAM" = 'Hamot', 
                "UPMCHZN" = 'Horizon', "UPMCJAM" = 'Jameson', "UPMCMCK" = 'McKeesport', 
                "UPMCMER" = 'Mercy', "UPMCMUN" = '', "UPMCMWH" = 'Magee-W', 
                "UPMCNOR" = 'Northwest', "UPMCPAS" = 'Passavant', "UPMCPUH" = 'Presbyterian', 
                "UPMCSHY" = 'Shadyside', "UPMCSMH" = 'St. Margaret', 
                "UPMCSOL" = 'SOL', "UPMCSUN" = 'SUN',  'UPMCLOC' = 'LOC',
                "UPMCWIL" = 'Williamsport')
site_groups <- list(
   academic = c('Mercy', 'Presbyterian', 'Shadyside', 'Magee-W'),
   regional = c('Hamot', 'Williamsport', 'Jameson', 'Altoona'),
   community = c('Passavant', 'East', 'McKeesport', 'St. Margaret'),
   rural = c('Bedford', 'Northwest', 'Horizon', 'Chatauqua')
)

col_vec <- c('VAN' = "#0000FF", 'iDAP' = "#ef5675", 'sDAP' = "#ffa600", 
             'allDAP' = '#FF0000', 'e1DAP' = '#ffa688', 'e2DAP' = '#ffa688', 'i2DAP' = '#ef5655',
             'eVAN' = '#1111aa', 'lVAN' = '#5555FF')
col_vec2 <- c('VAN' = "#0000FF22", 'iDAP' = "#ef567572", 'sDAP' = "#ffa60066", 
              'allDAP' = '#FF0000', 'e1DAP' = '#ffa68866', 'e2DAP' = '#ffa68866', 'i2DAP' = '#ef565588',
              'eVAN' = '#1111aa', 'lVAN' = '#5555FF')

site_cols <- setNames(c('#000000', '#0000FF', '#FF0000', '#00FF00'), c('academic', 'community', 'regional', 'rural'))

dfx <- dfx %>%
   mutate(
      TRT = case_when(
         TRT == 'VANCOMYCIN' ~ 'VAN',
         TRT == 'DAPTOMYCIN' ~ 'iDAP',
         TRT == 'VANCOMYCIN,DAPTOMYCIN' ~ 'sDAP'
      )
   ) %>%
   mutate(TRTs = ifelse(TRT == 'VAN', 'VAN', 'DAP'),
          FEMALE = GENDER == 'FEMALE',
          year = factor(lubridate::year(ORDER_DAY)),
          FACILITY = site_names[FACILITY]) %>%
   select(-GENDER)


dfx$CATEGORY <- NA
group_indices <- sapply(site_groups, function(x) which(dfx$FACILITY %in% x))
for (i in seq_along(group_indices)) {
   dfx$CATEGORY[group_indices[[i]]] <- names(site_groups)[i]
}

dfx %>% count(CATEGORY)
dfx %>% count(CATEGORY, FACILITY)
dfx <- dfx %>% filter(!is.na(CATEGORY))
dfx <- dfx %>% filter(is.na(time) | time > 3L)

dfx <- dfx %>%
   mutate(RECENT_MRSA = !is.na(DAYS_SINCE_PRV),
          RECENT_MRSA_1y = !is.na(DAYS_SINCE_PRV) & DAYS_SINCE_PRV <= 365L)# %>%
# mutate(DAYS_SINCE_PRV = case_when(
#    is.na(DAYS_SINCE_PRV) ~ -1 * 365 * length(2017:2023),
#    .default = -1 * DAYS_SINCE_PRV
# ))

dfx %>%
   #select(VAN_MIC, TRT) %>% #table() %>% barplot(beside=TRUE)
   summarise(n = n(),
             mic2 = sum(VAN_MIC_2, na.rm=T) / n,
             missing = sum(is.na(VAN_MIC)) / n,
             mean(VAN_MIC, na.rm=T),
             median(VAN_MIC, na.rm=T),
             .by = TRT)

dfx <- dfx %>%
   mutate(VAN_MIC = case_when(
      VAN_MIC == 0.75 ~ 0.5,
      VAN_MIC == 1.5 ~ 1,
      is.na(VAN_MIC) ~ 1,
      .default = VAN_MIC
   )) %>%
   mutate(VAN_MIC_2 = VAN_MIC >= 2)
dfx %>%
   summarise(n = n(),
             mic2 = sum(VAN_MIC_2, na.rm=T) / n,
             mean(VAN_MIC, na.rm=T),
             median(VAN_MIC, na.rm=T),
             .by = TRT)

t <- dfx %>%
   select(VAN_MIC, TRT) %>%
   table()
t <- apply(t, 2, function(x) x / sum(x))
barplot(t(t), beside=TRUE, legend=TRUE)
rm(t)



# cleanup variables
rm(group_indices, site_groups, i, site_names)




dfx %>%
   summarise(n = n(),
             DAP = sum(TRT %in% c('iDAP', 'sDAP')) / n,
             d30 = sum(time_censored < 30) / n,
             .by = c(CATEGORY, FACILITY)) %>%
   arrange(CATEGORY, desc(DAP))
t <- dfx %>%
   summarise(n = n(),
             VAN = sum(TRT == 'VAN'),
             iDAP = sum(TRT == 'iDAP'),
             sDAP = sum(TRT == 'sDAP'),
             DAP = sum(TRT %in% c('iDAP', 'sDAP')),
             .by = c(CATEGORY, FACILITY)) %>%
   arrange(CATEGORY, FACILITY)
write.csv(t, 
          file='~/Desktop/EHR/EHR-mining/Scripts/AnalysisScripts/TreatmentCohorts/MRSA_bacteremia/facility_category_counts.csv',
          quote = FALSE,
          row.names = FALSE)
rm(t)

dfx %>% summarise(n = n(),
                  n30 = sum(time_censored < 30),
                  r30 = n30 / n,
                  age = mean(AGE),
                  .by=CATEGORY)

dfx %>% summarise(n = n(),
                  n30 = sum(time_censored < 30),
                  r30 = n30 / n,
                  age = mean(AGE),
                  .by=c(CATEGORY, TRT)) %>%
   arrange(CATEGORY, TRT)

dfx %>% select(CATEGORY, TRT) %>% table()

pcx <- plr <- numeric(180L)
for (i in seq_along(pcx)) {
   print(i)
   censor_time <- i
   df <- dfx %>%
      filter(TRT != 'sDAP') %>%
      select(TRT, time, time_censored, status, AGE) %>%
      mutate(status = ifelse(is.na(time) | time > censor_time, 0, 1)) %>%
      mutate(time_censored = ifelse(is.na(time) | time > censor_time, censor_time, time))
   cox <- df %>%
      coxph(Surv(time_censored, status) ~ TRT + AGE, data=.) %>%
      summary()
   pcx[i] <- cox$coefficients['TRTVAN', 'Pr(>|z|)']
   lr <- df %>%
      glm(formula=time_censored ~ TRT + AGE, data=.) %>%
      summary()
   plr[i] <- lr$coefficients['TRTVAN', 'Pr(>|t|)']
}
plot(NA, type='l', log='y', ylim=c(0.01, 1), xlim=c(4,180), main='P-values on trt variable under coxph and linear regression')
abline(h = 0.05, lty=2)
abline(v = c(14, 30, 45), lty=2)
lines(pcx)
lines(plr, col='blue')
legend('bottomright', legend=c('coxph', 'linear'), col=c('black', 'blue'))
rm(pcx, plr, i, df, cox, lr)




# when is the first abx administered relative to blood culture??
s <- seq(-2, 4.5, 0.25)
setupPlot <- function(h, main='', leg='') {
   plot(NA, xlim=range(s), ylim=c(0, h), main=main, xlab='Days (relative to culture)', ylab='Density',
        cex.lab=cex, cex.axis=cex, cex.main=cex)
   abline(v=0, lty=3, lwd=1)
   legend('topright', legend=leg, col=col_vec[leg], lwd=cex*1.75, cex=cex)
}
drawDensityLines <- function(v, col) {
   lines(density(v), col=col, lwd=cex*1.75)
   abline(v = median(v), col=col, lwd=cex*1.75, lty=3)
}
par(mfrow=c(3,2), mar=c(5,4,2,2), mgp = c(2, 0.75, 0), tck=-0.015)
plot.new();plot.new();plot.new()
setupPlot(h=2.2, main='First ABX admin', leg=c('VAN', 'sDAP', 'iDAP'))
drawDensityLines(dfx$FIRST_VAN_TIME[dfx$TRT == 'VAN'], col_vec['VAN'])
drawDensityLines(dfx$FIRST_VAN_TIME[dfx$TRT == 'sDAP'], col_vec['sDAP'])
drawDensityLines(dfx$FIRST_DAP_TIME[dfx$TRT == 'iDAP'], col_vec['iDAP'])
setupPlot(h=2.2, main='First VAN admin', leg=c('VAN', 'sDAP'))
drawDensityLines(dfx$FIRST_VAN_TIME[dfx$TRT == 'VAN'], col_vec['VAN'])
drawDensityLines(dfx$FIRST_VAN_TIME[dfx$TRT == 'sDAP'], col_vec['sDAP'])
setupPlot(h=0.75, main='First DAP admin', leg=c('sDAP', 'iDAP'))
drawDensityLines(dfx$FIRST_DAP_TIME[dfx$TRT == 'sDAP'], col_vec['sDAP'])
drawDensityLines(dfx$FIRST_DAP_TIME[dfx$TRT == 'iDAP'], col_vec['iDAP'])
rm(s, setupPlot, drawDensityLines)




# drug distribution stratified by facility category
par(mfrow=c(2,1), mar=c(4,8,2,5), tck=-0.015)
t <- dfx %>%
   summarise(n = n(),
             nD30 = sum(time_censored < 30),
             nVAN = sum(TRT == 'VAN'), 
             niDAP = sum(TRT == 'iDAP'), 
             nsDAP = sum(TRT == 'sDAP'),
             fVAN = nVAN / n,
             fiDAP = niDAP / n,
             fsDAP = nsDAP / n,
             fD30 = nD30 / n,
             .by=CATEGORY) %>% arrange(desc(fVAN))
h <- t(as.matrix(t %>% select(fVAN, fiDAP, fsDAP)))
colnames(h) <- t$CATEGORY
cs <- cumsum(c(0, h[,ncol(h)]))
b <- barplot(h, horiz=TRUE, col=paste0(col_vec, 'aa'), names.arg=rep('', ncol(h)), xaxt='n', cex=cex, cex.axis=cex, cex.lab=cex)
text(x=-0.025, y=b, labels=colnames(h), adj=1, xpd=NA, cex=cex)
text(x=cs[-length(cs)] + diff(cs) / 2, y=b[length(b)]+0.75, labels=gsub('f', '', rownames(h)), xpd=NA, cex=cex)
text(x=1.025, y=b, adj=0, xpd=NA, labels=t$n, cex=cex)
rm(t, h, cs, b)



# drug use trend over time
t <- dfx %>% select(year, TRT) %>% table()
t <- apply(t, 1, function(x) x / sum(x) * 100)
b <- barplot(t, ylab='% of cases', col=paste0(col_vec[rownames(t)], 'aa'), cex.axis=cex, cex.names=cex, cex.lab=cex)
cs <- cumsum(c(0, t[,ncol(t)]))
text(x=b[length(b)] + 1, y=cs[1:3] + diff(cs) / 2, labels=rownames(t), xpd=NA, cex=cex)
rm(t, b, cs)

dfx %>%
   filter(CATEGORY == 'academic', TRT != 'iDAP') %>%
   mutate(year = case_when(
      year %in% 2017:2018 ~ '2017-2019',
      year %in% 2019:2020 ~ '2020-2021',
      year %in% 2021:2023 ~ '2022-2023'
   )) %>%
   summarise(n=n(),
             sep = sum(Sepsis_1w) / n,
             end = sum(Endocarditis_1w) / n,
             ost = sum(Osteomyelitis_1w) / n,
             age = mean(AGE),
             myo = sum(MyocInfarc) / n,
             hyp = sum(Hyperlipid) / n,
             dia = sum(OnDialysis) / n,
             nur = sum(NURSING_HOME) / n,
             dem = sum(Dementia) / n,
             dib = sum(Diabetes) / n,
             aki = sum(AKI_1w) / n,
             obe = sum(Obesity) / n,
             .by=c(TRT, year)) %>%
   select(-n) %>%
   arrange(year, TRT) %>%
   select(-TRT) %>%
   group_by(year) %>%
   mutate(across(sep:obe, ~ round(last(.) - first(.), 1))) %>%
   ungroup() %>%
   distinct()

# 30-day mortality by year?
t <- dfx %>%
   mutate(d30 = time_censored < 30L) %>%
   mutate(year = case_when(
      year %in% 2017:2018 ~ '2017-2019',
      year %in% 2019:2020 ~ '2020-2021',
      year %in% 2021:2023 ~ '2022-2023'
   )) %>%
   summarise(
      n = n(),
      d30 = sum(d30) / n * 100,
      .by = c(year, TRT)
   ) %>%
   tidyr::pivot_wider(
      id_cols = year,
      names_from = TRT,
      values_from = d30
   ) %>%
   arrange(year) %>%
   select(year, VAN, sDAP, iDAP)
par(mar=c(4,8,2,7))
barplot(t %>% select(-year) %>% as.matrix() %>% t, beside=TRUE, col=paste0(col_vec[c('VAN', 'sDAP', 'iDAP')], 'aa'),
        names.arg=t$year, main='30-day mortality', ylim=c(0, 20), ylab='%', legend=TRUE, 
        args.legend=list(x='topright', inset=c(-0.1,0), cex=cex))
rm(t)




# 30-day mortality by site by treatment
n <- dfx %>%
   summarise(n = n(),
             iDAP = sum(TRT == 'iDAP' & time_censored < 30), 
             sDAP = sum(TRT == 'sDAP' & time_censored < 30),
             VAN = sum(TRT == 'VAN' & time_censored < 30),
             .by=CATEGORY) %>% arrange(desc(CATEGORY)) %>% select(-n)
N <- dfx %>%
   summarise(n = n(),
             iDAP = sum(TRT == 'iDAP'),
             sDAP = sum(TRT == 'sDAP'),
             VAN = sum(TRT == 'VAN'),
             .by=CATEGORY) %>% arrange(desc(CATEGORY)) %>% select(-n)
f <- n$CATEGORY
n <- as.matrix(n %>% select(-CATEGORY))
N <- as.matrix(N %>% select(-CATEGORY))
n <- rbind(n, colSums(n))
N <- rbind(N, colSums(N))
rownames(n) <- rownames(N) <- c(f, 'OVERALL')

p <- n / N
se <- qnorm(0.975) * sqrt((1/n) * p * (1-p))

par(mfrow=c(1,1), mar=c(5, 11, 2, 2), mgp=c(2.5, 0.75, 0), tck=-0.015)
b <- barplot(t(p), horiz=TRUE, beside=TRUE, names.arg=rep('', nrow(p)), xlim=c(0, 0.3), args.legend = list(cex=cex),
             col=paste0(col_vec[colnames(p)], 'aa'), legend=TRUE, xlab='30-day mortality', cex.axis=cex, cex.names=cex, cex.lab=cex)
text(y=apply(b, 2, mean), x=-0.035, adj=1, xpd=NA, labels=rownames(p), cex=cex)
text(y=b, x=0.01, adj=0, col='white', labels=t(n), cex=cex)
text(y=b, x=-0.01, adj=1, labels=t(N), xpd=NA, cex=cex)
text(y=mean(b[,ncol(b)])+1.75, x=0.01, labels='30-day deaths', xpd=NA, adj=c(0, 0), cex=cex)
text(y=mean(b[,ncol(b)])+1.75, x=-0.01, labels='Total', xpd=NA, adj=c(1, 0), cex=cex)
rm(n, N, f, p, se, b)




# site distribution stratified by drug
par(mfrow=c(2,1), mar=c(3, 4, 3, 2), mgp=c(2, 0.3, 0), tck=-0.015)
m <- dfx %>%
   filter() %>%
   select(CATEGORY, TRT) %>% 
   table()
m <- apply(m, 2, function(x) x / sum(x) * 100)
cs <- apply(m, 2, function(x) cumsum(c(0, x)))
cs <- apply(cs, 2, function(x) x[-length(x)] + diff(x) / 2)
b <- barplot(m, col=paste0(site_cols, 'aa'), ylab='% of group', cex.names=cex, cex.axis=cex, cex.lab=cex)
text(x=rep(b, each=nrow(m)), y=cs, labels=rownames(m), xpd=NA, col=rep(c('white', 'black'), each=2), cex=cex)
rm(m, cs, b)


# by category and year
t <- dfx %>% select(CATEGORY, year, TRT) %>% table()
t <- aperm(apply(t, c(1, 2), function(x) x / sum(x) * 100), c(2,3,1))
labels <- apply(expand.grid(dimnames(t[,'2023',-3])), 1, paste, collapse=', ')
yvals <- as.vector(t[,'2023',-3])
o <- order(yvals, decreasing = TRUE)
yvals <- yvals[o]
labels <- labels[o]

par(mar=c(3,4,1,15))
plot(NA, xlim=c(2017, 2023), ylim=c(0,35), xlab='Year', ylab='% cases', cex=cex, cex.axis=cex, cex.lab=cex)
for (s in 1:2) { # loop over iDAP, sDAP
   m <- t[,,s]
   for (r in seq_len(nrow(m))) # loop over facility category
      lines(x=2017:2023, y=m[r,], col=site_cols[rownames(m)[r]], lwd=3, lty=s)
}
legend('right', inset=c(-0.56, 0), bty='n', xpd=NA, legend=labels, lwd=3, cex=cex,
       col=site_cols[gsub('(.+), .+', '\\1', labels)], lty=ifelse(grepl('iDAP', labels), 1, 2))
rm(t, labels, yvals, o, s, m, r)



## KAPLAN-MEIER CURVES
par(mfrow=c(1,1), mar=c(16, 14, 16, 14), mgp=c(2.5, 0.7, 0), tck=-0.015)
plotKP(df = dfx %>%
          select(-time) %>% 
          rename(time = time_censored),
       cex = cex,
       cohort='MRSA bacteremia', 
       col_vec=col_vec[1:3])
## END ##




cohorts <- c('VAN vs. iDAP', 'VAN vs. i2DAP', 'VAN vs. sDAP', 'VAN vs. allDAP', 'VAN vs. e1DAP', 'VAN vs. e2DAP', 'VAN vs. i2DAP')#, 'iDAP vs. sDAP', 'eVAN vs. lVAN')

for (cohort in cohorts) {
   print(cohort)
   
   trt <- unlist(strsplit(cohort, ' vs. ', fixed=TRUE))
   if (trt[2] == 'allDAP') {
      df <- dfx %>% mutate(TRT = ifelse(TRT == trt[1], TRT, 'allDAP'))
   } else if (trt[1] == 'eVAN' & trt[2] == 'lVAN') {
      trt <- c('eVAN', 'lVAN')
      df <- dfx %>%
         filter(TRT == 'VAN') %>%
         mutate(TRT = case_when(
            year %in% 2017:2021 ~ 'eVAN',
            year %in% 2022:2023 ~ 'lVAN'
         )) %>%
         filter(!is.na(TRT))
   } else if (trt[2] == 'e1DAP') {
      df <- dfx %>%
         filter(TRT == 'VAN' | (TRT == 'sDAP' & DAP_SWITCH_TIME < 1)) %>%
         mutate(TRT = ifelse(TRT == 'VAN', 'VAN', 'e1DAP'))
   } else if (trt[2] == 'e2DAP') {
      df <- dfx %>%
         filter(TRT == 'VAN' | (TRT == 'sDAP' & DAP_SWITCH_TIME < 2)) %>%
         mutate(TRT = ifelse(TRT == 'VAN', 'VAN', 'e2DAP'))
   } else if (trt[2] == 'i2DAP') {
      df <- dfx %>%
         filter(TRT %in% c('VAN', 'iDAP') | (TRT == 'sDAP' & DAP_SWITCH_TIME < 1)) %>%
         mutate(TRT = ifelse(TRT == 'VAN', 'VAN', 'i2DAP'))
   } else {
      df <- dfx %>% filter(TRT %in% trt)
   }
   df <- df %>%
      select(!c(FIRST_ABX_DAY, time)) %>% 
      rename(time = time_censored) %>%
      mutate(TRT = as.factor(TRT))
   
   
   # prep data for plotting
   d <- df %>%
      summarise(n = n(),
                n14 = sum(time < 14),
                f14 = n14 / n * 100,
                n30 = sum(time < 30),
                f30 = n30 / n * 100,
                .by=TRT)
   
   h <- as.matrix(d %>% select(f14, f30) %>% slice(2:1))
   b <- as.vector(barplot(h, plot=FALSE, beside=TRUE))
   
   
   # PLOT
   {
      layout(mat=matrix(c(1,2,4,3,3,3), nrow=3), widths=c(1, 1.5))
      par(mgp=c(2, 0.6, 0), mar=c(3, 4, 3.5, 3.5), tck=-0.015, oma=c(3,0,2,0), xpd=FALSE)
      
      # 30-day mortality
      b <- barplot(h, beside=TRUE, names.arg=rep('', 2), cex.names=cex, cex.lab=cex, cex.main=cex, ylim=c(0, 25), yaxt='n')
      title(main='All-cause mortality rates', line=1, cex.main=cex)
      text(x=b, y=-1.25, cex=cex, labels=rev(trt), xpd=NA)
      text(x=c(mean(b[1:2]), mean(b[3:4])), y=-3.6, cex=cex, labels=c('14-day', '30-day'), xpd=NA)
      axis(side=2, las=1, cex.axis=cex, tck=-0.015, at=seq(0, 25, 5), labels=paste0(seq(0, 25, 5), '%'))
      text(x=b, y=h+3.5, cex=cex, labels=paste0(round(h, 1), '%'))
      text(x=b, y=h+1.5, cex=cex, labels=paste0('n=', unlist(d %>% select(n14, n30) %>% slice(2:1))))
      
      # survival curve
      plotKP(df=df, cohort='30-day Kaplan-meier curves', trt=trt, cex=cex, col_vec=col_vec[trt])
      
      
      # visualize covariate imbalance
      vars <- names(df)[!names(df) %in% c('PERSON_ID', 'ORDER_DAY', 'ORDER_DATE', 'DAP_SWITCH_TIME', 'FIRST_VAN_TIME', 'FIRST_DAP_TIME',
                                          'status', 'FACILITY', 'TRT', 'TRTs', 'time', 'DAYS_SINCE_PRV')]
      tOg <- tableone::CreateTableOne(vars=vars, strata='TRT', data=df, smd=TRUE)
      table1 <- processTableOne(tOg)
      drawTableOne(table1, trt)
      #print(table1[table1$pvals < 0.1,])
      mtext(text=paste0(trt[2], ' (n = ', sum(df$TRT == trt[2]), '), ', trt[1], ' (n = ', sum(df$TRT == trt[1]), ')'), outer=TRUE, at=0.5, cex=cex-0.3, font=2)
      
      ## FACILITY CATEGORY SPECIFIC PRESCRIBING PATTERNS
      tables1 <- list()
      tables1$overall <- table1
      
      for (cat in unique(df$CATEGORY)) {
         tOg <- tableone::CreateTableOne(vars=vars, strata='TRT', data=df %>% filter(CATEGORY == cat), smd=TRUE)
         table1_cat <- processTableOne(tOg)
         table1_cat <- table1_cat[rownames(table1), ]
         cat(cat, '\n')
         
         tables1[[cat]] <- table1_cat
         
         table1_cat$sicker_group <- ifelse(table1_cat$diffs > 0, trt[1], trt[2])
         #print(table1_cat[table1_cat$pvals < 0.1,])
      }
      
      conf_vars <- unique(unlist(sapply(tables1, function(d) rownames(d[d$pvals < 0.1,]))))
      o <- order(tables1$overall[conf_vars, 'diffs'])
      tables1 <- lapply(tables1, function(d) d[conf_vars,][o,])
      
      
      # Propensity adjusted odds ratio estimates
      vars <- unique(c('AGE', 'CATEGORY', 'year', 'FEMALE', 'Obesity', 'RECENT_MRSA', 'RECENT_MRSA_1y',  # very general
                       'Endocarditis_1w', 'Respiratory_1w', 'AKI_1w', 'Cellulitis_1w', 'Osteomyelitis_1w', 'Peritonitis_1w', 'Sepsis_1w', 'SepticShock_1w', # acute clinical
                       rownames(table1[abs(table1$diffs) > 0.1, ]),
                       rownames(table1[table1$pvals < 0.1, ]),
                       conf_vars))
      write.table(x = vars, quote = FALSE, row.names = FALSE, col.names = FALSE,
                  file = paste0('~/Desktop/EHR/EHR-mining/Scripts/AnalysisScripts/TreatmentCohorts/MRSA_bacteremia/conf_vars_', trt[2], '.txt'))
      
      df <- df %>% mutate(trt1 = as.integer(TRT == trt[1]))
      formula <- as.formula(paste0('Surv(time, status) ~ trt1 + ', paste(vars, collapse=' + ')))
      # setup plot area
      ORplot <- function(trt) {
         par(mar=c(5, 4, 2.75, 3.5))
         plot(NA, xlim=c(0.5,2.5), ylim=c(1/5, 5), log='y', xaxt='n', xlab='', ylab='', yaxt='n', main='30-day cox proportional hazards', cex.main=cex)
         title(ylab='Hazard ratio', line=2.5, cex.lab=cex)
         title(xlab='Model', line=3.5, cex.lab=cex)
         axis(side=2, at=c(1/4, 1/2, 1, 2, 3, 4), labels=as.character(c(1/4, 1/2, 1, 2, 3, 4)), las=1)
         text(x=1:2, y=1/6, labels=c('Unadjusted', 'Inverse propensity\nscore weighted'), adj=c(0.5, 1), xpd=NA, cex=cex)
         abline(h = 1, lty=2)
         arrows(x0=0.5, y0=1/4.1, y1=4.1, length=0.1, code=3)
         text(x=0.45, y=c(1/4.8, 4.8), adj=0, labels=paste0('favors ', c(trt[1], trt[2])), cex=cex)
      }
      ORplot(trt)
      
      s <- summary(coxph(formula = formula, data=df))
      est <- s$coefficients[1, 'coef']
      se <- qnorm(0.975) * s$coefficients[1, 'se(coef)']
      points(x=1, y=exp(est), pch=16, cex=cex)
      arrows(x0=1, y0=exp(est-se), y1=exp(est+se), code=3, angle=90, length=0.05, lwd=cex)
      
      prop_model <- df %>% # predicts trt1
         glm(formula = paste0('trt1 ~ ', paste(vars, collapse=' + ')),
             family = binomial(),
             data = .)
      
      df$prop_score <- predict.glm(prop_model, newdata=df, type='response')
      df <- df %>%
         mutate(prop_weights = case_when(
            trt1 == 1L ~ 1 / prop_score,
            trt1 == 0L ~ 1 / (1 - prop_score)
         ))
      
      s <- summary(coxph(formula = formula, data=df, weights=prop_weights))
      est <- s$coefficients['trt1', 'coef']
      se <- qnorm(0.975) * s$coefficients[1, 'se(coef)']
      points(x=2, y=exp(est), pch=16, cex=cex)
      arrows(x0=2, y0=exp(est-se), y1=exp(est+se), code=3, angle=90, length=0.05, lwd=cex)
      
      #coefs_raw <- glm(formula=paste0('time ~ TRT + ', paste(vars, collapse=' + ')), data=df) %>% summary()
      #coefs_adj <- glm(formula=paste0('time ~ TRT + ', paste(vars, collapse=' + ')), data=df, weights=prop_weights) %>% summary()
      df <- df %>% select(-prop_score, -prop_weights, -trt1)
      
      # facility specifics
      for (t in seq_along(tables1)) {
         tables1[[t]] <- tables1[[t]][!rownames(tables1[[t]]) %in% c('year', 'CATEGORY'), ]
      }
      mat <- matrix(NA, nrow=nrow(tables1[[1]]), ncol=length(tables1), dimnames=list(rownames(tables1[[1]]), names(tables1)))
      for (col in seq_along(tables1)) {
         d <- tables1[[col]]
         cell_cols <- ifelse(d$diffs > 0, trt[1], trt[2])
         cell_cols <- ifelse(d$pvals < 0.1, col_vec[cell_cols], col_vec2[cell_cols])
         if (any(d$diffs == 0)) cell_cols[d$diffs == 0] <- 'white'
         mat[,col] <- cell_cols
      }
      xpos <- rep(1:ncol(mat), each=nrow(mat))
      ypos <- rep(1:nrow(mat), times=ncol(mat))
      
      par(mfrow=c(1,1), mar=c(8, 13, 7, 2))
      plot(NA, xlim=c(0.5,5.5), ylim=c(0.5,nrow(mat)+0.5), xaxs='i', yaxs='i', axes=F, ann=F)
      rect(xleft = xpos - 0.5,
           xright = xpos + 0.5,
           ybottom = ypos - 0.5,
           ytop = ypos + 0.5,
           col = mat,
           xpd = NA)
      text(x=0.4, y=1:nrow(mat), adj=1, xpd=NA, labels = rownames(tables1[[1]]))
      text(x=1:5, y=nrow(mat)+1, xpd=NA, labels=stringr::str_to_sentence(names(tables1)))
      rect(xleft=rep(c(2,3), times=2), xright=rep(c(3,4), times=2), 
           ybottom=rep(c(-1,-2), each=2), ytop=rep(c(0, -1), each=2), 
           col=c(col_vec[trt], col_vec2[trt]), xpd=NA)
      text(x=1.9, y=c(-0.5, -1.5), adj=1, xpd=NA, labels=c('p < 0.1', 'not sig.'))
      text(x=c(2.5, 3.5), y=-2.5, xpd=NA, labels=trt)
      text(x=3, y=-3.2, xpd=NA, labels=c('Greater rate of variable in group:'))
   }
}
##### END #####




















