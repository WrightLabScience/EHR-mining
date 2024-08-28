library(dplyr)

load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_2017_2023_AbxAdmin.Rdata') # NEW - 2017-2023
length(unique(abxDF$PERSON_ID)) # 317,502

# get counts of each medication
meds <- abxDF %>% count(MEDICATION, sort=TRUE)

load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/antibiotic_names/ast_antibiotics.Rdata')
abx_ast <- abx; rm(abx)
abx_ast <- abx_ast[abx_ast != 'PAS']
abx_ast <- c(abx_ast, 'PASER', 'SULFAMETHOXAZOLE/TRIMETHOPRIM', 'QUINUPRISTIN', 'DALFOPRISTIN', 'RELEBACTAM')
abx_abbr <- tibble(read.table(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/antibiotic_names/ABX_ABBR.txt', header = TRUE)) %>%
   select(Antibiotic_Name) %>%
   mutate(Antibiotic_Name = gsub('-', '/', Antibiotic_Name),
          Antibiotic_Name = gsub('_', ' ', Antibiotic_Name),
          Antibiotic_Name = toupper(Antibiotic_Name)) %>% unlist()
names(abx_abbr) <- NULL
abx_gen <- read.table(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/antibiotic_names/AntibioticNames.txt', sep='\t', header = TRUE)
abx_gen <- toupper(abx_gen$Generic)
abx_all <- sort(unique(c(abx_ast, abx_abbr, abx_gen)))
rm(abx_abbr, abx_ast, abx_gen)
abx_all <- unlist(strsplit(abx_all, '/'))
abx_all <- abx_all[!abx_all %in% c('CILASTATIN', 'PENICILLIN V', 'PENICILLIN G', 'STREPTOMYCIN_HIGH_LEVEL', 'STREPTOMYCIN_SYNERGY', 'GENTAMICIN_HIGH_LEVEL', 'GENTAMICIN_SYNERGY', 'CARBAPENEM_INACTIVATION_TEST', 'BETA_LACTAMASE', 'CLINDAMYCIN_INDUCIBLE')]
abx_all <- sort(unique(abx_all))
abx_all <- sapply(abx_all, function(x) grep(x, meds$MEDICATION))
abx_all <- abx_all[lengths(abx_all) > 0L]

# assign each individual antibiotic search term to the respective row of meds
meds$ABX <- character(nrow(meds))
for (i in seq_along(abx_all)) {
   print(i)
   rows <- abx_all[[i]]
   meds$ABX[rows] <- paste0(meds$ABX[rows], ',', names(abx_all)[i])
}
rm(i, rows, abx_all)

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
      .default = ABX
   ))
meds <- setNames(meds$ABX, meds$MEDICATION)

# create antibiotic column from the named vector of antibiotics
start <- Sys.time()
abxDF <- abxDF %>% mutate(ABX = meds[MEDICATION]) # this is large, so takes a minute
print(Sys.time() - start) # 0.3 seconds for ~11 million rows
rm(start, meds)
names(abxDF$ABX) <- NULL

# Remove unnecessary columns (for now) and take distinct()
abxDF <- abxDF %>%
   select(!c(ORDER_ID, ENCOUNTER_ID, MEDICATION_CODE, MEDICATION, ADMIN_ROUTE, ADMIN_DOSAGE, DOSAGE_UNIT)) %>%
   distinct() # 11,004,716 --> 9,681,748

# most common drugs
abxDF %>% count(ABX, sort=TRUE)
abxDF %>% count(substr(ADMIN_START_DATE, 1, 4)) %>% print(n=25)
abxDF %>% count(is.na(ADMIN_END_DATE)) # 2,638
abxDF %>% count(ADMIN_START_DATE == ADMIN_END_DATE) # 4,350,249 (~1/2)

abxDF %>%
   mutate(ADMIN_DURATION = as.numeric(lubridate::as.duration(ADMIN_END_DATE - ADMIN_START_DATE)) / 86400) %>%
   summarise(max = max(ADMIN_DURATION, na.rm=T),
             median = median(ADMIN_DURATION, na.rm=T),
             .by = ABX) %>%
   arrange(desc(median))

# re-arrange by ANTIBIOTIC instead of date
abxDF <- abxDF %>% arrange(PERSON_ID, ABX)

save(abxDF, file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_CLEANED_2017_2023_AbxAdmin.Rdata')






