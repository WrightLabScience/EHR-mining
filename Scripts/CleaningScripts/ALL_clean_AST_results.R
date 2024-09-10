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
   mutate(RESULT_DATE = as.Date(substr(RESULT_DATE,1,10), format='%m/%d/%Y')) %>%
   mutate(STATUS = case_when(
      grepl('SENSITIVE|NEGATIVE', SUSCEPTIBILITY_NAME)              | grepl('Sensitive',              SENSITIVITY_VALUE) ~ 0,
      grepl('RESISTANT|POSITIVE|INTERMEDIATE', SUSCEPTIBILITY_NAME) | grepl('Resistant|Intermediate', SENSITIVITY_VALUE) ~ 1
   )) %>%
   select(-ORGANISM_ID, -ANTIBIOTIC_LOINC_CODE, -SUSCEPTIBILITY_NAME, -SENSITIVITY_VALUE, -SENSITIVITY_COMMENT)


# get just the unique pathogen names for cleaning
source(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/CleanPathogenNames/CleanPathogenNames.R')
path_names <- astrDF %>% # 76,948
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
      stringi::stri_detect_regex(pattern='.+ [0-9]+(\\.[0-9])? ?(UG|MCG)/ML',                                        str=ANTIBIOTIC) ~ gsub('^(.+) [0-9]+(\\.[0-9])? ?(UG|MCG)/ML', '\\1', ANTIBIOTIC),
      stringi::stri_detect_regex(pattern='.+ \\([0-9]+\\.[0-9]+\\)',                                                 str=ANTIBIOTIC) ~ gsub('(.+) \\([0-9]+\\.[0-9]+\\)', '\\1', ANTIBIOTIC),
      stringi::stri_detect_regex(pattern='PENICILLIN \\((MENINGITIS|NON MENINGITIS|ORAL|PNEUMONIA|OTHER|E-TEST)\\)', str=ANTIBIOTIC) ~ 'PENICILLIN',
      stringi::stri_detect_regex(pattern='CEFOTAXIME \\((MENINGITIS|NON MENINGITIS)\\)',                             str=ANTIBIOTIC) ~ 'CEFOTAXIME',
      stringi::stri_detect_regex(pattern='CEFTRIAXONE \\((MENINGITIS|NON MENINGITIS)\\)',                            str=ANTIBIOTIC) ~ 'CEFTRIAXONE',
      stringi::stri_detect_regex(pattern='CEFUROXIME',                                                               str=ANTIBIOTIC) ~ 'CEFUROXIME',    # no other cefuroxime other than axetil, oral, sodium, parenteral
      stringi::stri_detect_regex(pattern='(QUIN(U|I)PRISTIN/DALFO)|(SYNERCID)',                                      str=ANTIBIOTIC) ~ 'SYNERCID',
      stringi::stri_detect_regex(pattern='(AMOXICILLIN/CLAV)|(AUGMENTIN)',                                           str=ANTIBIOTIC) ~ 'AMOXICILLIN/CLAVULANATE',
      stringi::stri_detect_regex(pattern='GENTAMICIN SYN',                                                           str=ANTIBIOTIC) ~ 'GENTAMICIN_SYNERGY',
      ANTIBIOTIC %in% c('ESBL', 'ESBL DRUG SCREEN', 'BETA LACTAMASE', 'CONFIRMATORY ESBL') ~ 'BETA_LACTAMASE',
      ANTIBIOTIC == 'TRIMETHOPRIM/SULFA' ~ 'TRIMETHOPRIM/SULFAMETHOXAZOLE',
      ANTIBIOTIC == 'SULFA/TRIMETHOPRIM' ~ 'TRIMETHOPRIM/SULFAMETHOXAZOLE',
      ANTIBIOTIC == 'SULFA' ~ 'SULFAMETHOXAZOLE',
      ANTIBIOTIC == 'IMIPENEM SUSCEPTIBILITY' ~ 'IMIPENEM',
      ANTIBIOTIC == 'CARBAPENEM INACTIVATION TEST' ~ 'CARBAPENEM_INACTIVATION_TEST',
      ANTIBIOTIC == 'HIGH CONCENTRATION OXACILLIN' ~ 'OXACILLIN',
      stringi::stri_detect_regex(pattern='PIPER', str=ANTIBIOTIC) & stringi::stri_detect_regex(pattern='TAZO', str=ANTIBIOTIC) ~ 'PIPERACILLIN/TAZOBACTAM',
      ANTIBIOTIC %in% c('GENTAMICIN 500', 'HIGH LEVEL GENTAMICIN', 'HIGH GENTAMICIN CONCENTRATION') ~ 'GENTAMICIN_HIGH_LEVEL',
      stringi::stri_detect_regex(pattern='STREPTOMYCIN SYN', str=ANTIBIOTIC) ~ 'STREPTOMYCIN_SYNERGY',
      stringi::stri_detect_regex(pattern='STREPTOMYCIN \\d|\\(', str=ANTIBIOTIC) ~ 'STREPTOMYCIN',
      ANTIBIOTIC == 'HIGH LEVEL STREPTOMYCIN' ~ 'STREPTOMYCIN_HIGH_LEVEL',
      ANTIBIOTIC == 'CEFOXITIN SCREEN' ~ 'CEFOXITIN',
      ANTIBIOTIC == '5-FLUOROCYTOSINE 24 HRS:' ~ 'FLUCYTOSINE',
      ANTIBIOTIC ==  '5-FLUCYSTOSINE' ~ 'FLUCYTOSINE',
      stringi::stri_detect_regex(pattern='INDUCIBLE CLINDAMYCIN', str=ANTIBIOTIC) ~ 'CLINDAMYCIN_INDUCIBLE',
      ANTIBIOTIC == 'IMIPENEM RELEBACTAM' ~ 'IMIPENEM/RELEBACTAM',
      ANTIBIOTIC == 'VANCOMYCIN SCREEN' ~ 'VANCOMYCIN',
      ANTIBIOTIC == 'ACYCLOVIR IC 50' ~ 'ACYCLOVIR',
      ANTIBIOTIC == 'TICARCILLIN/CLAV' ~ 'TICARCILLIN/CLAVULANATE',
      ANTIBIOTIC == 'PYRAZINAMIDE 25.0-50.0UG/ML' ~ 'PYRAZINAMIDE',
      ANTIBIOTIC == 'DAPTIOMYCIN' ~ 'DAPTOMYCIN',
      ANTIBIOTIC == 'AMP/SUBLACT.' ~ 'AMPICILLIN/SULBACTAM',
      stringi::stri_detect_regex(pattern='AMOXICILLIN (&|AND) CLAVULANATE', str=ANTIBIOTIC) ~ 'AMOXICILLIN/CLAVULANATE',
      ANTIBIOTIC == 'CEFTOLOZANE-TAZOBACTAM' ~ 'CEFTOLOZANE/TAZOBACTAM',
      stringi::stri_detect_regex(pattern='AMP', str=ANTIBIOTIC) & stringi::stri_detect_regex(pattern='SULBAC', str=ANTIBIOTIC) ~ 'AMPICILLIN/SULBACTAM',
      stringi::stri_detect_regex(pattern='TICARCILLIN', str=ANTIBIOTIC) & stringi::stri_detect_regex(pattern='CLAVULANATE', str=ANTIBIOTIC) ~ 'TICARCILLIN/CLAVULANATE',
      .default = ANTIBIOTIC
   )) %>%
   filter(!ANTIBIOTIC %in% c('COMMENT:', 'ORGANISM', 'NA', 'ACYCLOVIR'))


# arrange the data
astrDF <- astrDF %>%
   arrange(PERSON_ID, ORDER_PROC_ID, PATH_NAME, ANTIBIOTIC)

# to collapse identical results (ignoring LINE_NUM and RESULT_DATE)
astrDF <- astrDF %>%
   select(-LINE_NUM) %>%
   # select(-RESULT_DATE) %>%
   distinct() # 21,544,642 --> 21,510,091


# dup_ids <- astrDF$ORDER_PROC_ID[duplicated(astrDF[-3])]
# length(dup_ids) # 3,638



# # 21,506,453 --> 21,501,852
# astrDF1 <- astrDF %>% filter(n() == 1L, .by = c(ORDER_PROC_ID, PATH_NAME, ANTIBIOTIC)) # 21,497,251
# astrDF2 <- astrDF %>% filter(n() > 1L, .by = c(ORDER_PROC_ID, PATH_NAME, ANTIBIOTIC))  # 9,202
# astrDF2 <- astrDF2 %>% slice_max(STATUS, by = c(ORDER_PROC_ID, PATH_NAME, ANTIBIOTIC)) # 4,601
# astrDF <- rbind(astrDF1, astrDF2) %>%
#    arrange(PERSON_ID, ORDER_PROC_ID, PATH_NAME, ANTIBIOTIC) # 21,501,852
# rm(astrDF1, astrDF2)

sapply(astrDF, function(x) sum(is.na(x)))
astrDF %>% count(is.na(PATH_NAME), is.na(ANTIBIOTIC), is.na(STATUS))

astrDF <- astrDF %>% filter(!is.na(PATH_NAME), !is.na(ANTIBIOTIC)) # 21,510,091 --> 21,486,898
astrDF %>% count(is.na(STATUS)) # 1,650,554

# get all antibiotic names
abx <- unique(astrDF$ANTIBIOTIC)
save(abx, file = '~/Desktop/EHR/EHR work/data/ast_antibiotics.Rdata')

astrDF <- astrDF %>% # 21,486,898 --> 1,860,785
   pivot_wider(names_from = ANTIBIOTIC,
               values_from = STATUS,
               values_fn = max)

# astrDF <- astrDF %>% distinct() # 1,860,372 --> 1,846,711

astrDF <- astrDF %>%
   mutate(OXACILLIN = case_when(
      is.na(OXACILLIN) & BUG == 'Staphylococcus aureus' & grepl('cillin resis|mrsa', PATH_NAME) & !grepl('previously reported as mrsa', PATH_NAME) ~ 1,
      is.na(OXACILLIN) & BUG == 'Staphylococcus aureus' & grepl('methicillin sensit', PATH_NAME) ~ 0,
      .default = OXACILLIN
   ))



save(astrDF, file = '~/Desktop/EHR/EHR work/RdataFiles/AST_results_clean.Rdata')

   
   
   
   
   
   