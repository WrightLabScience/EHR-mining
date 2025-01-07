# This script takes the currently separate antibiotic med_admin txt files
# and combines then into a large tibble and cleans them
library(dplyr)


# first SENS_MED_ADMIN_VW
path <- '~/Desktop/EHR/EHR work/RdataFiles/MED_ADMIN_data/SENS_MED_ADMIN_VW/'
ls <- list.files(path)
abxDF <- tibble()
for (i in seq_along(ls)) {
   cat(i, ls[i], '\n')
   abxDF <- rbind(abxDF, 
                  tibble(read.table(file = paste0(path, ls[i]), header=TRUE, sep='\t')))
}
rm(path, ls, i)


# then SENS_MED_ADMIN_IV_VW
path <- '~/Desktop/EHR/EHR work/RdataFiles/MED_ADMIN_data/SENS_MED_ADMIN_IV_VW/'
ls <- list.files(path)
for (i in seq_along(ls)) {
   cat(i, ls[i], '\n')
   abxDF <- rbind(abxDF, 
                  tibble(read.table(file = paste0(path, ls[i]), header=TRUE, sep='\t')))
}
rm(path, ls, i)


# arrange and remove duplicates, then save
abxDF <- abxDF %>%
   arrange(PERSON_ID, ADMIN_START_DATE) %>%
   distinct()

save(abxDF, file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_2017_2023_AbxAdmin.Rdata')


###################################
######## CLEAN ABX TIBBLE #########
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_2017_2023_AbxAdmin.Rdata')
length(unique(abxDF$PERSON_ID)) # 322,355

# get counts of each medication
meds <- abxDF %>% count(MEDICATION, sort=TRUE)
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/antibiotic_names/abx_names_for_abxDF.Rdata')

# assign each individual antibiotic search term to the respective row of meds
meds$ABX <- character(nrow(meds))
for (i in seq_along(abx_all)) {
   print(i)
   abx <- abx_all[i]
   w <- grep(abx, meds$MEDICATION)
   meds$ABX[w] <- paste(meds$ABX[w], abx, sep=',')
}
rm(i, w, abx_all, abx)

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
      
      .default = ABX
   ))
meds <- setNames(meds$ABX, meds$MEDICATION)

# create antibiotic column from the named vector of antibiotics
start <- Sys.time()
abxDF <- abxDF %>% mutate(ABX = meds[MEDICATION])
print(Sys.time() - start) # 0.3 seconds for ~8.65 million rows
rm(start, meds)
names(abxDF$ABX) <- NULL

# handle AMPHOTERICIN,GENTAMICIN,VANCOMYCIN
w <- which(abxDF$ABX == 'AMPHOTERICIN,GENTAMICIN,VANCOMYCIN')
df1 <- df2 <- df3 <- abxDF[w,]
df1$ABX <- 'AMPHOTERICIN'
df2$ABX <- 'GENTAMICIN'
df3$ABX <- 'VANCOMYCIN'
df <- rbind(df1, df2, df3) %>% arrange(PERSON_ID, ADMIN_START_DATE)
abxDF <- rbind(abxDF[-w,], df)
rm(df1, df2, df3, df, w)


# Remove unnecessary columns (for now) and take distinct()
abxDF <- abxDF %>%
   select(!c(ORDER_ID, ENCOUNTER_ID, MEDICATION_CODE, MEDICATION, ADMIN_ROUTE, ADMIN_DOSAGE, DOSAGE_UNIT)) %>%
   mutate(across(c(ADMIN_START_DATE, ADMIN_END_DATE), ~ strptime(., format='%Y-%m-%d %T')),
          PERSON_ID = as.character(PERSON_ID)) %>%
   rename(START_DATE = ADMIN_START_DATE,
          END_DATE = ADMIN_END_DATE) %>%
   mutate(START_DAY = lubridate::as_date(START_DATE)) %>%
   filter(lubridate::year(START_DATE) %in% 2016:2024) %>%
   distinct() %>% # 8,649,274 --> 8,587,158
   arrange(PERSON_ID, ABX, START_DATE)

abxDF <- abxDF %>% filter(PERSON_ID != '1.001e+09')

# most common drugs
abxDF %>% count(ABX, sort=TRUE)
abxDF %>% count(lubridate::year(START_DAY))
abxDF %>% count(is.na(END_DATE)) # 34,257
abxDF %>% count(START_DATE == END_DATE) # 3,326,883 (~40%)


abxDF %>%
   mutate(ADMIN_DURATION = as.numeric(lubridate::as.duration(END_DATE - START_DATE)) / 86400) %>%
   summarise(max = max(ADMIN_DURATION, na.rm=T),
             median = median(ADMIN_DURATION, na.rm=T),
             .by = ABX) %>%
   arrange(desc(median))

save(abxDF, file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_CLEANED_2017_2023_AbxAdmin.Rdata')






