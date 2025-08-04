library(dplyr)
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ASTs_AbxAdmin_blood_2017_2023_imputed_encs_flagged_WIDE_survival.Rdata')
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/antibiotic_names/abbr_named.Rdata')
empDF$TARG_ABBR <- sapply(empDF$ABX_TAR, function(x) paste(abbr[x], collapse='+'))

empDF %>% 
   filter(BUG == 'Enterococcus faecium', VANCOMYCIN == 1L) %>%
   count(grepl('DAPTOMYCIN', ABX_TAR) | grepl('DAPTOMYCIN', ABX_BTW),
         grepl('LINEZOLID', ABX_TAR) | grepl('LINEZOLID', ABX_BTW))

empDF %>% filter(BUG == 'Enterococcus faecium', VANCOMYCIN == 1L,
                 #substr(ORDER_DAY,1,4) == '2017',
                 lengths(ABX_TAR) == 1L) %>% 
   count(TARG_ABBR, sort=TRUE)
# 37 DAP
# 19 LZD

empDF %>% filter(BUG == 'Staphylococcus aureus', OXACILLIN == 0L,
                 #substr(ORDER_DAY,1,4) == '2017',
                 lengths(ABX_TAR) == 1L) %>% 
   count(TARG_ABBR, sort=TRUE)
# antibiotic -all -2017
# CFZ        1492   123
# OXA         678    94
# NFC         157    40
# VAN         101    19
# TZP          88     8

empDF %>% filter(BUG == 'Escherichia coli', ESBL == 1L,
                 lengths(ABX_TAR) == 1L) %>%
   count(TARG_ABBR, sort=TRUE)
# 245 - ertapenem
# 206 - meropenem
# 18 - pip/tazo
# 13 - cefepime

empDF %>% filter(BUG == 'Escherichia coli', ESBL == 0L,
                 lengths(ABX_TAR) == 1L) %>%
   count(TARG_ABBR, sort=TRUE)
# CRO - 1344
# CFZ - 488
# TZP - 352
# CIP - 180
# FEP - 154
# CXM - 146 (cefuroxime)


empDF %>% filter(BUG == 'Klebsiella pneumoniae', ESBL == 1L,
                 lengths(ABX_TAR) == 1L) %>%
   count(TARG_ABBR, sort=TRUE)
# MEM - 70
# ETP - 42
# CZA - 10 (ceftaz-avibactam)
 










