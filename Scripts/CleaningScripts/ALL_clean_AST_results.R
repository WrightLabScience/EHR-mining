####################################################
############# CLEAN UP AST RESULT TABLE ############
####################################################
start <- Sys.time()
library(dplyr)
library(tidyr)

load(file = '~/Desktop/EHR/EHR work/RdataFiles/lab_micro_all_vw.Rdata')


# simple prep
astrDF <- astrDF2 %>%
   rename(PATH_NAME = ORGANISM_NAME,
          ANTIBIOTIC = ANITBIOTIC_NAME) %>%
   mutate(PATH_NAME = gsub('  ', ' ', PATH_NAME)) %>%
   mutate(
      PATH_NAME = stringr::str_to_lower(PATH_NAME),
      RESULT_DATE = as.Date(substr(RESULT_DATE,1,10), format='%m/%d/%Y'),
      STATUS = case_when(
         grepl('SENSITIVE|NEGATIVE', SUSCEPTIBILITY_NAME)              | grepl('Sensitive',              SENSITIVITY_VALUE) ~ 0,
         grepl('RESISTANT|POSITIVE|INTERMEDIATE', SUSCEPTIBILITY_NAME) | grepl('Resistant|Intermediate', SENSITIVITY_VALUE) ~ 1
      )
   ) %>%
   select(-LOINC_CODE, -SUSCEPTIBILITY_NAME, -SENSITIVITY_VALUE)
rm(astrDF2); gc()


# get just the unique pathogen names for cleaning
source(file = '~/Desktop/EHR-mining/UsefulDataForCleaning/CleanPathogenNames/CleanPathogenNames.R')
astrDF <- astrDF %>% filter(!(is.na(PATH_NAME) & is.na(ANTIBIOTIC)))
path_names <- astrDF %>% # 78,777
   count(PATH_NAME, sort=TRUE) %>%
   mutate(BUG = NA_character_)
path_names <- cleanPathogenNames(path_names)
path_names <- path_names %>% filter(!is.na(BUG))
path_names <- setNames(object = path_names$BUG, nm = path_names$PATH_NAME)
astrDF <- astrDF %>%
   mutate(BUG = unname(path_names[PATH_NAME])) %>%
   mutate(BUG = ifelse(is.na(BUG), 'Did not match', BUG))
astrDF$BUG[which(is.na(astrDF$PATH_NAME))] <- NA
rm(path_names)


# Clean up the antibiotic names
# astrDF <- astrDF %>%
abx_names <- astrDF %>% count(ANTIBIOTIC)
abx_names <- abx_names %>%
   mutate(ABX = case_when(
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
   ))
w <- which(abx_names$ABX %in% c('COMMENT:', 'ORGANISM', 'NA', 'ACYCLOVIR') | is.na(abx_names$ANTIBIOTIC)) # 5
abx_names <- abx_names[-w,] %>% select(-n) %>% distinct()
abx_names <- setNames(abx_names$ABX,
                      abx_names$ANTIBIOTIC)
rm(w)
w <- which(astrDF$ANTIBIOTIC %in% c('COMMENT:', 'ORGANISM', 'NA', 'ACYCLOVIR') | is.na(astrDF$ANTIBIOTIC)) # 34
astrDF$ANTIBIOTIC[w] <- 'Unknown'
astrDF$ANTIBIOTIC <- unname(abx_names[astrDF$ANTIBIOTIC])


astrDF <- astrDF %>%
   select(-LINE_NUM) %>%
   distinct()


# if same antibiotic tested against same pathogen on same order, take maximum
astrDF <- astrDF %>% group_by(ORDER_PROC_ID, PATH_NAME, ANTIBIOTIC)
astrDFm <- astrDF %>% filter(n() > 1L)  # 16,461 (8,222 groups)
astrDFm <- astrDFm %>% slice_max(STATUS) # 11,852 (8,222 groups) - still some with much later result_dates
astrDFm <- astrDFm %>% slice_min(RESULT_DATE) # 8,222
astrDFm <- ungroup(astrDFm)
astrDF <- rbind(
   astrDF %>% filter(n() == 1L) %>% ungroup(), 
   astrDFm
) %>%
   arrange(PERSON_ID, ORDER_PROC_ID, PATH_NAME, ANTIBIOTIC) # 21,501,852
rm(astrDFm); gc()


# check for missingness
sapply(astrDF, function(x) sum(is.na(x)))
astrDF %>% count(is.na(PATH_NAME), is.na(ANTIBIOTIC), is.na(STATUS))
astrDF$ANTIBIOTIC[which(is.na(astrDF$ANTIBIOTIC))] <- 'Unknown'
astrDF$ANTIBIOTIC[which(astrDF$ANTIBIOTIC == 'PAS')] <- 'Unknown'


# get all antibiotic names
abx <- sort(unique(astrDF$ANTIBIOTIC))
abx <- abx[abx != 'Unknown']
save(abx, file = '~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/ast_antibiotics.Rdata')

######################################################################################
save(astrDF, file = '~/Desktop/EHR/EHR work/RdataFiles/lab_micro_sens_all_cleaned_long.Rdata')
print(Sys.time() - start) # ~14.5 minutes
######################################################################################



library(dplyr)
library(tidyr)

load(file = '~/Desktop/EHR/EHR work/RdataFiles/lab_micro_sens_all_cleaned_long.Rdata')


# pivot to wide
start <- Sys.time()
astrDF <- astrDF %>% # 21,478,682 --> 1,789,478
   pivot_wider(id_cols = c(PERSON_ID, ORDER_PROC_ID, RESULT_DATE, BUG),
               names_from = ANTIBIOTIC,
               values_from = STATUS,
               values_fn = max,
               unused_fn = ~ paste(unique(.), collapse=', ')) %>%
   relocate(PATH_NAME, .after=BUG)
print(Sys.time() - start) # 1.4 minutes

# Detect MRSA from pathogen_name when oxacillin status is missing
astrDF <- astrDF %>%
   mutate(OXACILLIN = case_when(
      is.na(OXACILLIN) & BUG == 'Staphylococcus aureus' & grepl('cillin resis|mrsa', PATH_NAME) & !grepl('previously reported as mrsa', PATH_NAME) ~ 1,
      is.na(OXACILLIN) & BUG == 'Staphylococcus aureus' & grepl('methicillin sensit', PATH_NAME) ~ 0,
      .default = OXACILLIN
   ))


astrDF <- astrDF %>% 
   select(-Unknown) %>%
   mutate(CEPHALEXIN = NA)

astrDF <- astrDF %>%
   mutate(across(.cols = CEFEPIME:CEPHALEXIN,
                 .fns = as.integer))



# i think some of the enterobacterales I called ESBL-producing may not actually be
# probably the ones that were missing CRO and were tested against any carbapenem
ent <- c("Escherichia coli", 
         "Enterobacter aerogenes", "Enterobacter cloacae",
         "Citrobacter freundii", 
         "Klebsiella aerogenes", "Klebsiella oxytoca", "Klebsiella pneumoniae", "Klebsiella variicola", 
         "Morganella morganii",
         "Proteus mirabilis",
         "Serratia marcescens")
enDF <- astrDF %>% filter(BUG %in% ent)
w <- which(grepl('esbl|extended spectrum beta lactamase', astrDF$PATH_NAME) 
           & !grepl('not an ?esbl|low dilution esbl|possible esbl', astrDF$PATH_NAME)) # 41,203

astrDF$ESBL <- 0L
astrDF$ESBL[w] <- 1L

w <- which(astrDF$CEFTRIAXONE == 1L & astrDF$BUG %in% ent) # 54,193

astrDF$ESBL[w] <- 1L




# astrDF_og <- astrDF
# # w <- which(grepl('esbl|extended spectrum beta lactamase', astrDF$PATH_NAME) & !grepl('not an ?esbl|low dilution esbl|possible esbl', astrDF$PATH_NAME))
# 
# 
# astrDF <- astrDF %>%
#    mutate(ENT = BUG %in% c("Escherichia coli", "Klebsiella pneumoniae", 'Pseudomonas aeruginosa', 'Proteus mirabilis', 'Proteus penneri',
#                            "Enterobacter cloacae", "Klebsiella oxytoca", "Serratia marcescens", "Enterobacter aerogenes", "Klebsiella aerogenes",
#                            "Proteus vulgaris"),
#           ESBL_FLAG = grepl('esbl|extended spectrum beta lactamase', PATH_NAME) & !grepl('not an ?esbl|low dilution esbl|possible esbl', PATH_NAME)) %>%
#    mutate(ESBL = ENT & (ESBL_FLAG | CEFTRIAXONE == 1)) %>%
#    mutate(ESBL = case_when(
#       is.na(ESBL) ~ 0,
#       ESBL ~ 1,
#       !ESBL ~ 0
#    )) %>%
#    select(-ENT, -ESBL_FLAG)
# 
# bugs <- astrDF %>% count(BUG, sort=TRUE)
# 
# astrDF <- astrDF %>%
#    mutate(EKR = BUG %in% c("Klebsiella pneumoniae", "Klebsiella oxytoca", "Klebsiella variicola", "Klebsiella species", "Klebsiella ozaenae", "Klebsiella ornithinolytica",
#                            'Escherichia coli', "Raoultella planticola", "Raoultella ornithinolytica", "Raoultella species")) %>% # E coli, Klebsiella, Raoultella
#    mutate(ESBL2 = EKR & (CEFEPIME == 1L | CEFOTAXIME == 1L | CEFTRIAXONE == 1L | CEFTAZIDIME == 1L)) # resistant to 3rd or 4th gen cephalosporins
# 
# w <- which(astrDF$EKR & astrDF$ESBL2 & astrDF$ESBL == 0L)
# astrDF$ESBL[w] <- 1
# astrDF <- astrDF %>% select(-EKR, -ESBL2)


# 
# en <- astrDF %>%
#    filter(BUG %in% c("Escherichia coli", "Klebsiella pneumoniae", 'Proteus mirabilis', 'Proteus penneri',
#                      "Enterobacter cloacae", "Klebsiella oxytoca", "Serratia marcescens", "Enterobacter aerogenes", "Klebsiella aerogenes",
#                      "Proteus vulgaris"))
# en %>%
#    summarise(esbl = sum(ESBL),
#              not = sum(!ESBL),
#              total = n(),
#              rate = esbl / total * 100,
#              .by=BUG) %>%
#    arrange(desc(rate))
# 
# 
# en <- en %>%
#    mutate(testCarba = !is.na(MEROPENEM) | !is.na(DORIPENEM) | !is.na(ERTAPENEM) | !is.na(IMIPENEM),
#           missCRO = is.na(CEFTRIAXONE),
#           ESBL2 = )
# 
# sum(en$testCarba) / nrow(en) # 87%
# sum(en$missCRO) / nrow(en)   # 26%
# sum(en$ESBL) / nrow(en)      # 4.2%
# 
# en %>% select(testCarba, missCRO, ESBL) %>% table()
# (97571)
# 
# 
# f <- function(x) {
#    total <- nrow(x)
#    x <- sapply(x %>% select(CEFEPIME:DELAFLOXACIN), function(x) c(sum(x, na.rm=T), sum(!is.na(x))))
#    rownames(x) <- c('R', 'Tot')
#    x <- x[, x['Tot',] > 0]
#    x <- data.frame(t(x))
#    x$rate <- round(x$R / x$Tot * 100, 1)
#    x$fracT <- round(x$Tot / total * 100, 1)
#    x <- x[x$fracT > 50, c('rate', 'fracT')]
#    x <- x[order(x$fracT), ]
#    return(x)
# }
# 
# en <- en %>% filter(as.integer(substr(RESULT_DATE,1,4)) > 2020)
# 
# p0 <- f(en %>% filter(ESBL))
# p1 <- f(en %>% filter(!ESBL, missCRO, testCarba))
# p2 <- f(en %>% filter(ESBL, missCRO, testCarba))
# p3 <- f(en %>% filter(!ESBL, !missCRO, testCarba))
# p4 <- f(en %>% filter(ESBL, !missCRO, testCarba))
# p5 <- f(en %>% filter(!ESBL, missCRO, !testCarba))
# p6 <- f(en %>% filter(ESBL, missCRO, !testCarba))
# p7 <- f(en %>% filter(!ESBL, !missCRO, !testCarba))
# p8 <- f(en %>% filter(ESBL, !missCRO, !testCarba))
# 
# # LEVOFLOXACIN                  72.9  60.4
# # MEROPENEM                      3.2  62.2
# # TOBRAMYCIN                    46.2  69.8
# # NITROFURANTOIN                26.8  70.2
# # ERTAPENEM                      5.2  70.6
# # GENTAMICIN                    31.3  90.8
# # CIPROFLOXACIN                 78.2  91.8
# # TRIMETHOPRIM/SULFAMETHOXAZOLE 66.0  93.2
# 
# 
# 
# en %>% select(testCarba, missCRO) %>% table() # 86% of no carbapenem are missing CRO vs. 25% of yes carbapenem are missing CRO, testing against both is most likely, testing neither is 2nd likely
# en %>% filter(testCarba, !missCRO) %>% count(CEFTRIAXONE)  #     tested against carbapenem, has CRO, most are CRO-susceptible
# en %>% filter(!testCarba, !missCRO) %>% count(CEFTRIAXONE) # NOT tested against carbapenem, has CRO, most are CRO-susceptible
# 
# en %>% select(testCarba, ESBL) %>% table() # tested against carbapenem, slightly more likely to be flagged as ESBL
# en %>% filter(testCarba, ESBL) %>% count(CEFTRIAXONE)  #     tested against carbapenem, has CRO, most are CRO-susceptible
# en %>% filter(!testCarba, !) %>% count(CEFTRIAXONE)






save(astrDF, file = '~/Desktop/EHR/EHR work/RdataFiles/AST_results_clean.Rdata')



astrDF <- astrDF %>% mutate(year = as.integer(substr(RESULT_DATE, 1,4)))

x <- astrDF %>% filter(BUG == 'Escherichia coli') # ESBL prevalence has double in last 12 years 4.3% --> 8.5%
x <- astrDF %>% filter(BUG == 'Pseudomonas aeruginosa') # declined??
x <- astrDF %>% filter(BUG == 'Klebsiella pneumoniae') # 7% --> 13% over last 12 years
x %>%
   summarise(esbl = sum(ESBL),
             n = n(),
             rate = esbl / n * 100,
             .by = year) %>%
   arrange(desc(year)) %>%
   print(n=21) # all bugs: doubled from 2.9 --> 5.1% over last 12 years


