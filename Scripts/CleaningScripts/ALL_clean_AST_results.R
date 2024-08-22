####################################################
############# CLEAN UP AST RESULT TABLE ############
####################################################
library(dplyr)
library(tidyr)

load(file = '~/Desktop/EHR/EHR work/RdataFiles/lab_micro_sens_all_vw.Rdata')

# simple prep
astrDF <- astrDF %>%
   rename(PATH_NAME = ORGANISM_NAME,
          ANTIBIOTIC = ANITBIOTIC_NAME) %>%
   mutate(PATH_NAME = gsub('  ', ' ', PATH_NAME)) %>%
   mutate(PATH_NAME = stringr::str_to_lower(PATH_NAME)) %>%
   mutate(STATUS = case_when(
      grepl('SENSITIVE|NEGATIVE', SUSCEPTIBILITY_NAME)              | grepl('Sensitive',              SENSITIVITY_VALUE) ~ 0,
      grepl('RESISTANT|POSITIVE|INTERMEDIATE', SUSCEPTIBILITY_NAME) | grepl('Resistant|Intermediate', SENSITIVITY_VALUE) ~ 1
   )) %>%
   select(-ORGANISM_ID, -ANTIBIOTIC_LOINC_CODE, -SUSCEPTIBILITY_NAME, -SENSITIVITY_VALUE, -SENSITIVITY_COMMENT)


# get just the unique pathogen names for cleaning
source('~/Desktop/EHR/EHR work/data/CleanPathogenNames/CleanPathogenNames.R')
path_names <- astrDF %>%
   count(PATH_NAME, sort=TRUE) %>%
   mutate(BUG = NA_character_)
path_names <- cleanPathogenNames(path_names)
path_names <- path_names %>% filter(!is.na(BUG))
path_names <- setNames(object = path_names$BUG, nm = path_names$PATH_NAME)
astrDF <- astrDF %>%
   mutate(BUG = path_names[PATH_NAME]) %>%
   mutate(BUG = ifelse(is.na(BUG), 'Did not match', BUG))
names(astrDF$BUG) <- NULL
rm(path_names)


# Clean up the antibiotic names
astrDF <- astrDF %>%
   mutate(ANTIBIOTIC = case_when(
      grepl('.+ [0-9]+(\\.[0-9])? ?(UG|MCG)/ML', ANTIBIOTIC) ~ gsub('^(.+) [0-9]+(\\.[0-9])? ?(UG|MCG)/ML', '\\1', ANTIBIOTIC),
      grepl('.+ \\([0-9]+\\.[0-9]+\\)', ANTIBIOTIC) ~ gsub('(.+) \\([0-9]+\\.[0-9]+\\)', '\\1', ANTIBIOTIC),
      grepl('PENICILLIN \\((MENINGITIS|NON MENINGITIS|ORAL|PNEUMONIA|OTHER|E-TEST)\\)', ANTIBIOTIC) ~ 'PENICILLIN',
      grepl('CEFOTAXIME \\((MENINGITIS|NON MENINGITIS)\\)', ANTIBIOTIC) ~ 'CEFOTAXIME',
      grepl('CEFTRIAXONE \\((MENINGITIS|NON MENINGITIS)\\)', ANTIBIOTIC) ~ 'CEFTRIAXONE',
      grepl('CEFUROXIME', ANTIBIOTIC) ~ 'CEFUROXIME',                                                                   # no other cefuroxime other than axetil, oral, sodium, parenteral
      ANTIBIOTIC %in% c('ESBL', 'ESBL DRUG SCREEN', 'BETA LACTAMASE', 'CONFIRMATORY ESBL') ~ 'BETA_LACTAMASE',
      ANTIBIOTIC == 'TRIMETHOPRIM/SULFA' ~ 'TRIMETHOPRIM/SULFAMETHOXAZOLE',
      ANTIBIOTIC == 'SULFA/TRIMETHOPRIM' ~ 'TRIMETHOPRIM/SULFAMETHOXAZOLE',
      ANTIBIOTIC == 'SULFA' ~ 'SULFAMETHOXAZOLE',
      ANTIBIOTIC == 'IMIPENEM SUSCEPTIBILITY' ~ 'IMIPENEM',
      ANTIBIOTIC == 'CARBAPENEM INACTIVATION TEST' ~ 'CARBAPENEM_INACTIVATION_TEST',
      ANTIBIOTIC == 'HIGH CONCENTRATION OXACILLIN' ~ 'OXACILLIN',
      grepl('PIPER', ANTIBIOTIC) & grepl('TAZO', ANTIBIOTIC) ~ 'PIPERACILLIN/TAZOBACTAM',
      grepl('(QUIN(U|I)PRISTIN/DALFO)|(SYNERCID)', ANTIBIOTIC) ~ 'SYNERCID',
      grepl('(AMOXICILLIN/CLAV)|(AUGMENTIN)', ANTIBIOTIC) ~ 'AMOXICILLIN/CLAVULANATE',
      grepl('GENTAMICIN SYN', ANTIBIOTIC) ~ 'GENTAMICIN_SYNERGY',
      ANTIBIOTIC %in% c('GENTAMICIN 500', 'HIGH LEVEL GENTAMICIN', 'HIGH GENTAMICIN CONCENTRATION') ~ 'GENTAMICIN_HIGH_LEVEL',
      grepl('STREPTOMYCIN SYN', ANTIBIOTIC) ~ 'STREPTOMYCIN_SYNERGY',
      grepl('STREPTOMYCIN \\d|\\(', ANTIBIOTIC) ~ 'STREPTOMYCIN',
      ANTIBIOTIC == 'HIGH LEVEL STREPTOMYCIN' ~ 'STREPTOMYCIN_HIGH_LEVEL',
      ANTIBIOTIC == 'CEFOXITIN SCREEN' ~ 'CEFOXITIN',
      ANTIBIOTIC == '5-FLUOROCYTOSINE 24 HRS:' ~ 'FLUCYTOSINE',
      ANTIBIOTIC ==  '5-FLUCYSTOSINE' ~ 'FLUCYTOSINE',
      grepl('INDUCIBLE CLINDAMYCIN', ANTIBIOTIC) ~ 'CLINDAMYCIN_INDUCIBLE',
      ANTIBIOTIC == 'IMIPENEM RELEBACTAM' ~ 'IMIPENEM/RELEBACTAM',
      ANTIBIOTIC == 'VANCOMYCIN SCREEN' ~ 'VANCOMYCIN',
      ANTIBIOTIC == 'ACYCLOVIR IC 50' ~ 'ACYCLOVIR',
      ANTIBIOTIC == 'TICARCILLIN/CLAV' ~ 'TICARCILLIN/CLAVULANATE',
      ANTIBIOTIC == 'PYRAZINAMIDE 25.0-50.0UG/ML' ~ 'PYRAZINAMIDE',
      ANTIBIOTIC == 'DAPTIOMYCIN' ~ 'DAPTOMYCIN',
      ANTIBIOTIC == 'AMP/SUBLACT.' ~ 'AMPICILLIN/SULBACTAM',
      grepl('AMOXICILLIN (&|AND) CLAVULANATE', ANTIBIOTIC) ~ 'AMOXICILLIN/CLAVULANATE',
      ANTIBIOTIC == 'CEFTOLOZANE-TAZOBACTAM' ~ 'CEFTOLOZANE/TAZOBACTAM',
      grepl('AMP', ANTIBIOTIC) & grepl('SULBAC', ANTIBIOTIC) ~ 'AMPICILLIN/SULBACTAM',
      grepl('TICARCILLIN', ANTIBIOTIC) & grepl('CLAVULANATE', ANTIBIOTIC) ~ 'TICARCILLIN/CLAVULANATE',
      .default = ANTIBIOTIC
   )) %>%
   filter(!ANTIBIOTIC %in% c('COMMENT:', 'ORGANISM', 'NA', 'ACYCLOVIR'))


# arrange the data
astrDF <- astrDF %>%
   arrange(PERSON_ID, ORDER_PROC_ID, PATH_NAME, ANTIBIOTIC)

# to collapse identical results (ignoring LINE_NUM and RESULT_DATE)
astrDF <- astrDF %>%
   select(-LINE_NUM, -RESULT_DATE) %>%
   distinct() # 21,544,642 --> 21,506,453


# 21,506,453 --> 21,501,852
astrDF1 <- astrDF %>% filter(n() == 1L, .by = c(ORDER_PROC_ID, PATH_NAME, ANTIBIOTIC)) # 21,497,251
astrDF2 <- astrDF %>% filter(n() > 1L, .by = c(ORDER_PROC_ID, PATH_NAME, ANTIBIOTIC))  # 9,202
astrDF2 <- astrDF2 %>% slice_max(STATUS, by = c(ORDER_PROC_ID, PATH_NAME, ANTIBIOTIC)) # 4,601
astrDF <- rbind(astrDF1, astrDF2) %>%
   arrange(PERSON_ID, ORDER_PROC_ID, PATH_NAME, ANTIBIOTIC) # 21,501,852
rm(astrDF1, astrDF2)

sapply(astrDF, function(x) sum(is.na(x)))
astrDF %>% count(is.na(PATH_NAME), is.na(ANTIBIOTIC), is.na(STATUS))

astrDF <- astrDF %>% filter(!is.na(PATH_NAME), !is.na(ANTIBIOTIC)) # 21,501,852 --> 21,478,659
astrDF %>% count(is.na(STATUS)) # 1,650,200 (out of 21,478,659)

# get all antibiotic names
abx <- unique(astrDF$ANTIBIOTIC)
save(abx, file = '~/Desktop/EHR/EHR work/data/ast_antibiotics.Rdata')

astrDF <- astrDF %>% # 21,478,659 --> 1,860,372
   pivot_wider(id_cols = c(PERSON_ID, ORDER_PROC_ID, PATH_NAME, BUG),
               names_from = ANTIBIOTIC,
               values_from = STATUS) %>%
   select(-PATH_NAME)

astrDF <- astrDF %>% distinct() # 1,860,372 --> 1,846,711


save(astrDF, file = '~/Desktop/EHR/EHR work/RdataFiles/AST_results_clean.Rdata')

