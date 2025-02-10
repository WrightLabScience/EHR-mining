# This script takes the currently separate antibiotic med_admin txt files
# and combines then into a large tibble and cleans them
library(dplyr)
setwd('~/Desktop/EHR/EHR work/RdataFiles/MED_ADMIN_DATA/')


# first SENS_MED_ADMIN_VW
path <- 'SENS_MED_ADMIN_ABX_DATA/'
ls <- list.files(path)
abxDF <- tibble()
for (i in seq_along(ls)) {
   cat(i, ls[i], '\n')
   abxDF <- rbind(abxDF, 
                  tibble(read.table(file = paste0(path, ls[i]), header=TRUE, sep='\t')))
}
rm(path, ls, i)


# then SENS_MED_ADMIN_IV_VW
path <- 'SENS_MED_ADMIN_IV_ABX_DATA/'
ls <- list.files(path)
for (i in seq_along(ls)) {
   cat(i, ls[i], '\n')
   abxDF <- rbind(abxDF, 
                  tibble(read.table(file = paste0(path, ls[i]), header=TRUE, sep='\t')))
}
rm(path, ls, i)


# arrange and remove duplicates, then save
abxDF <- abxDF %>%
   arrange(PERSON_ID, ADMIN_START_DATE, MEDICATION) %>%
   distinct()

save(abxDF, file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_2015_2024_AbxAdmin.Rdata')




###################################
######## CLEAN ABX TIBBLE #########
library(dplyr)
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_2015_2024_AbxAdmin.Rdata')
length(unique(abxDF$PERSON_ID)) # 398,266

# get counts of each medication
meds <- abxDF %>% count(MEDICATION, sort=TRUE) # 1,114
load(file = '~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/abx_names_for_abxDF.Rdata')

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
      ABX == 'CLAVULANATE,TICARCILLIN' ~ 'TICARCILLIN/CLAVULANATE',
      ABX == 'DURLOBACTAM,SULBACTAM' ~ 'SULBACTAM/DURLOBACTAM',
      .default = ABX
   ))
meds <- setNames(meds$ABX, meds$MEDICATION)

# create antibiotic column from the named vector of antibiotics - shockingly fast
abxDF <- abxDF %>% mutate(ABX = unname(meds[MEDICATION]))
rm(meds)

# handle multi antibiotic administrations that should be considered separately
meds <- abxDF %>% count(ABX, sort=TRUE)
meds %>% filter(grepl(',', ABX))
separateAndExtend <- function(abx_list) {
   w <- which(abxDF$ABX == abx_list)
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
   rbind(abxDF[-w,], df)
}
abxDF <- separateAndExtend('BACITRACIN,NEOMYCIN,POLYMYXIN B')
abxDF <- separateAndExtend('BACITRACIN,POLYMYXIN B')
abxDF <- separateAndExtend('POLYMYXIN B,TRIMETHOPRIM')
abxDF <- separateAndExtend('NEOMYCIN,POLYMYXIN B')
abxDF <- separateAndExtend('COLISTIN,GENTAMICIN')
abxDF <- separateAndExtend('AMPHOTERICIN,GENTAMICIN,VANCOMYCIN')
abxDF <- separateAndExtend('COLISTIN,NEOMYCIN')
rm(meds, separateAndExtend)




# Remove unnecessary columns (for now) and take distinct()
abxDF <- abxDF %>%
   select(!c(ORDER_ID, ENCOUNTER_ID, MEDICATION_CODE, MEDICATION, ADMIN_ROUTE, ADMIN_DOSAGE, DOSAGE_UNIT)) %>%
   mutate(across(c(ADMIN_START_DATE, ADMIN_END_DATE), ~ strptime(., format='%Y-%m-%d %T')),
          PERSON_ID = as.character(PERSON_ID)) %>%
   rename(START_DATE = ADMIN_START_DATE,
          END_DATE = ADMIN_END_DATE) %>%
   mutate(START_DAY = lubridate::as_date(START_DATE)) %>%
   distinct() %>%
   arrange(PERSON_ID, START_DATE, ABX)

abxDF %>% count(lubridate::year(START_DATE)) %>% print(n=15)

abxDF <- abxDF %>% 
   filter(PERSON_ID != '1.001e+09') %>%
   filter(!is.na(START_DATE)) %>%
   filter(lubridate::year(START_DATE) >= 2014L)

# most common drugs
abxDF %>% count(ABX, sort=TRUE)


abxDF %>%
   mutate(ADMIN_DURATION = as.numeric(lubridate::as.duration(END_DATE - START_DATE)) / 86400) %>%
   summarise(max = max(ADMIN_DURATION, na.rm=T),
             median = median(ADMIN_DURATION, na.rm=T),
             .by = ABX) %>%
   arrange(desc(median))

save(abxDF, file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_CLEANED_2015_2024_AbxAdmin.Rdata')






