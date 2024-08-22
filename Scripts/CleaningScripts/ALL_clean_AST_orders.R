###################################################
############# CLEAN UP AST ORDER TABLE ############
###################################################
library(dplyr)
library(tidyr)

load(file = '~/Desktop/EHR/EHR work/RdataFiles/lab_micro_result_all_vw.Rdata')

# astoDF %>% filter(grepl('CATH', ORDER_NAME)) %>% count(ORDER_NAME, sort=TRUE)
# astoDF %>% filter(grepl('CATH', SPECIMEN_TYPE)) %>% count(SPECIMEN_TYPE, sort=TRUE)
# astoDF %>% filter(grepl('BOTTLE', SPECIMEN_TYPE)) %>% count(SPECIMEN_TYPE, sort=TRUE)
# astoDF %>% filter(grepl('BOTTLE', ORDER_NAME)) %>% count(ORDER_NAME, sort=TRUE)
# astoDF %>% filter(grepl('ANAEROBIC', ORDER_NAME)) %>% count(ORDER_NAME, sort=TRUE)

#### SITE OF CULTURE ####
# on %>%
#    mutate(SITE = case_when(
#       grepl('URIN', ORDER_NAME) ~ 'urine',
#       grepl('BLOOD', ORDER_NAME) ~ 'blood',
#       grepl('WOUND', ORDER_NAME) ~ 'wound',
#       grepl('SPUTUM|TRACH|RESP|THROAT|BRONCH|BAL', ORDER_NAME) ~ 'respiratory',
#       grepl('SURGICAL', ORDER_NAME) ~ 'surgical tissue',
#       grepl('OPHTHAL|EYE|EAR|NOSE|NASAL', ORDER_NAME) ~ 'eye, ear, nose',
#       grepl('STOOL', ORDER_NAME) ~ 'stool',
#       grepl('GENITAL', ORDER_NAME) ~ 'genital',
#       grepl('BONE', ORDER_NAME) ~ 'bone',
#       grepl('CSF', ORDER_NAME) ~ 'csf',
#       .default = ORDER_NAME
#    )) %>%
#    summarise(n = sum(n), .by=SITE) %>%
#    arrange(desc(n)) %>%
#    print(n=100)


orders <- astoDF %>%
   count(ORDER_NAME, sort=TRUE) %>%
   mutate(BLOOD = grepl('BLOOD', ORDER_NAME) & !grepl('NON BLOOD|PERITONEAL DIALYSATE', ORDER_NAME)) %>%
   filter(BLOOD)
orders <- orders$ORDER_NAME

specs <- astoDF %>%
   count(SPECIMEN_TYPE, sort=TRUE) %>%
   mutate(BLOOD = grepl('BLOOD', SPECIMEN_TYPE) & !grepl('ASCITES FLUID IN BLOOD BOTTLES', SPECIMEN_TYPE)) %>%
   filter(BLOOD)
specs <- specs$SPECIMEN_TYPE

values <- astoDF %>%
   count(RESULT_VALUE, sort=TRUE) %>%
   mutate(BLOOD = (stringi::stri_detect_regex(pattern='BLOOD|Peripheral Blood', str=RESULT_VALUE)
                   | stringi::stri_detect_regex(pattern='blood culture', str=RESULT_VALUE, case_insensitive=TRUE))
                  & stringi::stri_detect_regex(pattern='CELL', str=RESULT_VALUE, negate=TRUE)) %>%
   filter(BLOOD)
values <- values$RESULT_VALUE


# create column indicating whether I think it refers to a bloodstream infection
astoDF <- astoDF %>% mutate(BLOOD = ORDER_NAME %in% orders | SPECIMEN_TYPE %in% specs | RESULT_VALUE %in% values)
rm(values, orders, specs)

# make date-time column
# make "pathogen name" column (which occasional is mentioned in the result_value column)
astoDF <- astoDF %>%
   mutate(across(contains('DATE'), ~ strptime(., format='%m/%d/%Y %T')),
          PATH_NAME = stringr::str_to_lower(gsub('  ', ' ', RESULT_VALUE)))

astoDF %>% count(is.na(ORDER_DATE)) # 8,656 (out of 7.7 million)
astoDF <- astoDF %>% 
   filter(!is.na(ORDER_DATE)) %>%
   select(!c(REFERENCE_UNIT, REFERENCE_LOW, REFERENCE_HIGH, COMPONENT_NAME, RESULT_LAB_NAME, SPECIMEN_COLLECTED_DATE, 
             SPECIMEN_RECEIVED_DATE, LOINC_CODE, ORDER_NAME, SPECIMEN_TYPE, RESULT_VALUE)) %>%
   distinct() %>% # 7,702,107 --> 6.6 million
   arrange(PERSON_ID, ORDER_DATE, ORDER_PROC_ID)

# some order_proc_ids are blood, but blood wasn't assigned TRUE for all rows of that ID
blood_order_ids <- unique(astoDF$ORDER_PROC_ID[astoDF$BLOOD])
length(blood_order_ids) # 132,060
length(unique(astoDF$ORDER_PROC_ID)) # 1,550,181
table(astoDF$BLOOD[astoDF$ORDER_PROC_ID %in% blood_order_ids]) # ~21K are FALSE!
astoDF$BLOOD[astoDF$ORDER_PROC_ID %in% blood_order_ids] <- TRUE # no longer
rm(blood_order_ids)

gc()


# get Pathogen Names
source(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/CleanPathogenNames/CleanPathogenNames.R')

# AST ORDERS
path_names <- astoDF %>%
   count(PATH_NAME, sort=TRUE) %>%
   filter(!is.na(PATH_NAME)) %>%
   mutate(BUG = NA_character_) # ~482K !!
path_names <- cleanPathogenNames(path_names)
path_names <- path_names %>% filter(!is.na(BUG))
path_names <- setNames(object = path_names$BUG, nm = path_names$PATH_NAME)
astoDF <- astoDF %>%
   mutate(BUG = path_names[PATH_NAME]) %>%
   select(-PATH_NAME) %>%
   distinct()

rm(path_names, cleanPathogenNames)


# remove instances where ORDER_DATE comes after RESULT_DATE - clear error
astoDF <- astoDF %>%
   filter(RESULT_DATE > ORDER_DATE)

# sometimes everything is identical, but one row has a MUCH later RESULT_DATE, take the minimum
astoDF <- astoDF %>% group_by(PERSON_ID, ORDER_PROC_ID, ORDER_DATE, BLOOD, BUG)
astoDF1 <- astoDF %>% filter(n() == 1L)
astoDF2 <- astoDF %>% filter(n() > 1L)
astoDF2 <- astoDF2 %>% slice_min(RESULT_DATE)
astoDF <- rbind(astoDF1, astoDF2) %>%
   arrange(PERSON_ID, ORDER_PROC_ID, ORDER_DATE, RESULT_DATE)
rm(astoDF1, astoDF2)

# sometimes everything is identical, but one row has a MUCH early ORDER_DATE, take the maximum
astoDF <- astoDF %>% group_by(PERSON_ID, ORDER_PROC_ID, RESULT_DATE, BLOOD, BUG)
astoDF1 <- astoDF %>% filter(n() == 1L) %>% ungroup()
astoDF2 <- astoDF %>% filter(n() > 1L)
astoDF2 <- astoDF2 %>% slice_max(ORDER_DATE) %>% ungroup()
astoDF <- rbind(astoDF1, astoDF2) %>%
   arrange(PERSON_ID, ORDER_PROC_ID, ORDER_DATE, RESULT_DATE)
rm(astoDF1, astoDF2)


save(astoDF, file = '~/Desktop/EHR/EHR work/RdataFiles/AST_orders_clean.Rdata')




