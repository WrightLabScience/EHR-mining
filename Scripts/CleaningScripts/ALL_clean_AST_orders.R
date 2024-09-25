###################################################
############# CLEAN UP AST ORDER TABLE ############ Total run-time = ~ 10 minutes probable
###################################################
library(dplyr)
library(tidyr)

load(file = '~/Desktop/EHR/EHR work/RdataFiles/lab_micro_result_all_vw.Rdata')



# site of infection -------------------------------------------------------

# stn <- astoDF %>% count(SPECIMEN_TYPE, sort=TRUE)
# 
# stn %>%
#    filter(grepl('CATH', SPECIMEN_TYPE))
# stn %>%
#    filter(grepl('CYST', SPECIMEN_TYPE)) %>%
#    print(n=32)
# 
# 
# stn %>%
#    mutate(SITE = case_when(
#       grepl('URIN|UREOSTOMY|URETEROSTOMY|MIDSTREAM|PASSPORT', SPECIMEN_TYPE) ~ 'urine',
#       SPECIMEN_TYPE == 'CATHETER URINE' ~ 'urine',
#       grepl('BILIARY DRAINAGE|G(-| )TUBE|ABDOMINAL|T-TUBE|ASCITES|PERITONEAL DIALYSATE', SPECIMEN_TYPE) ~ 'body fluid',
#       grepl('ASCITES|FLUID|ABDOMEN|BILE|^PUS$|GALLBLADDER', SPECIMEN_TYPE) ~ 'body fluid',
#       grepl('CSF|CNS', SPECIMEN_TYPE) ~ 'csf',
#       grepl('SPUTUM|TRACH|RESP|THROAT|BRONCH|BAL|BRUSH|CHEST TUBE|ASPIRATE|PIGTAIL|LUNG', SPECIMEN_TYPE) ~ 'respiratory',
#       grepl('INCISION|LESION|SUTURE LINE|SEROMA|SURGICAL|GRAFT|JACKSON PRATT|JP SITE', SPECIMEN_TYPE) ~ 'surgical',
#       grepl('GENITAL|VAGINAL|CERVIX|VULVA|LABIA|PENIS|PLACENTA|GROIN|URETHRA|BARTHOLIN|PERINEUM|UMBILICAL|SEMEN', SPECIMEN_TYPE) ~ 'genital',
#       grepl('CATHETER', SPECIMEN_TYPE) ~ 'catheter',
#       SPECIMEN_TYPE == 'PERMA-CATH' ~ 'catheter',
#       grepl('BLOOD', SPECIMEN_TYPE) ~ 'blood',
#       grepl('JP.+DRAINAGE|WOUND|ABSCESS', SPECIMEN_TYPE) ~ 'wound',
#       SPECIMEN_TYPE %in% c('DRAINAGE', 'TISSUE', 'ABCESS') ~ 'wound',
#       grepl('SKIN|BUTTOCKS|TOE|LEG|FOOT|FINGER|KNEE|BREAST|AXILLA|SCALP|ANKLE|ELBOW|HIP|ARM|HEEL|THIGH|HAND|BACK|COCCYX|SHOULDER|NECK|WRIST|CHEST|NAIL|CELLULITIS|PELVIS', SPECIMEN_TYPE) ~ 'skin',
#       grepl('NOSE|NASAL|NASOPHARYNGEAL|SINUS|MEATUS|NARES', SPECIMEN_TYPE) ~ 'nose',
#       grepl('EAR', SPECIMEN_TYPE) ~ 'ear',
#       grepl('OPHTHAL|EYE|CORNEA|CONJUNCTIVA|VITREOUS|LENS|AQUEOUS|LACRIMAL', SPECIMEN_TYPE) ~ 'eye',
#       grepl('STOOL|RECTAL|JEJUNAL|COLOSTOMY|ILEOSTOMY|ULCER|PERI-ANAL|UREOSTOMY|ANUS|RECTUM|APPENDIX', SPECIMEN_TYPE) ~ 'intestinal',
#       grepl('BONE|FEMORAL|SACRAL|TIBIA|SACRUM', SPECIMEN_TYPE) ~ 'bone',
#       grepl('ORAL|MOUTH|DENTAL|TONGUE', SPECIMEN_TYPE) ~ 'oral',
#       SPECIMEN_TYPE %in% c('EXIT SITE', 'PORT SITE') ~ 'catheter',
#       SPECIMEN_TYPE == 'STRAIGHT CATH' ~ 'urine',
#       grepl('PROSTHETIC', SPECIMEN_TYPE) ~ 'prosthetic',
#       SPECIMEN_TYPE == 'CARDIAC PACEMAKER OR LEADS' ~ 'pacemaker',
#       grepl('HEMATOMA', SPECIMEN_TYPE) ~ 'hematoma',
#       grepl('LYMPH', SPECIMEN_TYPE) ~ 'lymph',
#       grepl('SYNOVIAL', SPECIMEN_TYPE) ~ 'joint',
#       SPECIMEN_TYPE %in% c('KIDNEY', 'KIDNEY BIOPSY') ~ 'kidney',
#       SPECIMEN_TYPE == 'CYST' ~ 'cyst',
#       is.na(SPECIMEN_TYPE) | SPECIMEN_TYPE %in% c('SWAB FROM UNKNOWN SITE', 'OTHER SOURCE,DEFINE IN COMMENTS', 'BACTERIAL ISOLATE', 'SWAB', 'UNKNOWN SITE', 'MYCOBACTERIAL ISOLATE', 
#                                                   'DRAW SITE NOT SPECIFIED') ~ 'site not specified',
#       SPECIMEN_TYPE %in% c('DONOR', 'SITE', 'BIOPSY', 'STONE') ~ 'unclassified site specified'
#    )) %>%
#    # filter(grepl('CYST', SPECIMEN_TYPE)) 
#    filter(grepl('CATH', SPECIMEN_TYPE)) %>% print(n=42)
# filter(is.na(SITE)) %>%
#    #mutate(SITE = ifelse(is.na(SITE) & !is.na(SPECIMEN_TYPE), , 'site not specified')) %>%
#    #summarise(n = sum(n), .by=SITE)
#    print(n=20)
# 
# 
# # if it comes down to it...
# astoDF %>% 
#    filter(SPECIMEN_TYPE == 'OTHER SOURCE,DEFINE IN COMMENTS') %>%
#    count(RESULT_VALUE, sort=TRUE)
# 
# 'WND CULT', 'SURG W ANAER CULT'
# 
# ons <- astoDF %>% count(ORDER_NAME, sort=TRUE)
# 
# ons %>%
#    mutate(SITE = case_when(
#       grepl('URIN', ORDER_NAME) ~ 'urine',
#       grepl('BLOOD', ORDER_NAME) ~ 'blood',
#       grepl('WOUND|DRAINAGE|ABSCESS|TISSUE', ORDER_NAME) ~ 'wound',
#       grepl('SKIN', ORDER_NAME) ~ 'skin',
#       grepl('SPUTUM|TRACH|RESP|THROAT|BRONCH|BAL|BRUSH', ORDER_NAME) ~ 'respiratory',
#       ORDER_NAME == 'FLUID/ASPIRATE CULTURE' ~ 'respiratory',
#       grepl('SURGICAL', ORDER_NAME) ~ 'surgical tissue',
#       grepl('OPHTHAL|EYE|EAR|NOSE|NASAL|NASOPHARYNGEAL', ORDER_NAME) ~ 'eye, ear, nose',
#       grepl('STOOL|RECTAL', ORDER_NAME) ~ 'stool',
#       grepl('GENITAL', ORDER_NAME) ~ 'genital',
#       grepl('BONE', ORDER_NAME) ~ 'bone',
#       grepl('CSF', ORDER_NAME) ~ 'csf',
#       grepl('ORAL|MOUTH', ORDER_NAME) ~ 'oral',
#       grepl('CATHETER', ORDER_NAME) ~ 'catheter',
#       grepl('PROSTHESIS', ORDER_NAME) ~ 'prosthetic',
#       ORDER_NAME %in% c('BODY FLUID CULTURE & GRAM STAIN', 'BODY FLUID CULTURE - AEROBIC', 'AMNIOTIC FLUID CULT(MAGEE ONLY)', 'SEMEN CULTURE', 'BREAST MILK CULTURE(MWH)') ~ 'body fluid',
#       ORDER_NAME %in% c('CULTURE, AEROBIC BACTERIA', 'CULTURE,AEROBIC (W/GRAM STAIN)', 'AEROBIC CULTURE (OO)', 'CULTURE, AEROBIC BACTERIA WITH GRAM STAIN', 'AUTOPSY AEROBIC CULTURE (OO)') ~ 'aerobic',
#       ORDER_NAME %in% c('ANAEROBIC CULTURE', 'ANAEROBIC CULTURE (OO)', 'CULTURE, ANAEROBIC BACTERIA WITH GRAM STAIN', 'AUTOPSY ANAEROBIC CULTURE (OO)') ~ 'anaerobic',
#       .default = 'site not specified'
#    )) %>%
#    summarise(n = sum(n), .by=SITE) %>%
#    arrange(desc(n)) %>%
#    print(n=108)



# astoDF %>% filter(grepl('CATH', ORDER_NAME)) %>% count(ORDER_NAME, sort=TRUE)
# astoDF %>% filter(grepl('CATH', SPECIMEN_TYPE)) %>% count(SPECIMEN_TYPE, sort=TRUE)
# astoDF %>% filter(grepl('BOTTLE', SPECIMEN_TYPE)) %>% count(SPECIMEN_TYPE, sort=TRUE)
# astoDF %>% filter(grepl('BOTTLE', ORDER_NAME)) %>% count(ORDER_NAME, sort=TRUE)
# astoDF %>% filter(grepl('ANAEROBIC', ORDER_NAME)) %>% count(ORDER_NAME, sort=TRUE)

# SITE OF CULTURE
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


###############


# make date-time column
# make "pathogen name" column (which occasional is mentioned in the result_value column)
astoDF <- astoDF %>%
   mutate(across(contains('DATE'), ~ strptime(., format='%m/%d/%Y %T')),
          PATH_NAME = stringr::str_to_lower(gsub('  ', ' ', RESULT_VALUE)))


orders <- astoDF %>%
   count(ORDER_NAME, sort=TRUE) %>%
   filter(grepl('BLOOD', ORDER_NAME)) %>%
   filter(!ORDER_NAME %in% c('PERITONEAL DIALYSATE IN BLOOD (OO)', 'FUNGUS CULTURE, NON BLOOD', 'AUTOPSY BLOOD CULTURE (OO)', 'BODY FLUID IN BLOOD CULTURE BOTTLES'))
orders <- orders$ORDER_NAME

# mark the above ORDER_NAME as blood cultures
astoDF <- astoDF %>% mutate(BLOOD = ORDER_NAME %in% orders)

# for the non-blood cultures, look elsewhere to determine if blood culture or not - start with SPECIMEN_TYPE field
specs <- astoDF %>%
   filter(is.na(ORDER_NAME)) %>%
   count(SPECIMEN_TYPE, sort=TRUE)

# if AST ORDER_NAME is missing, but SPECIMEN_TYPE says "PERIPHERAL BLOOD" then we will take it as a blood culture
# it appears that in ~mid 2018 there was a switch from putting "BLOOD CULTURE" in the ORDER_NAME field --> "PERIPHERAL BLOOD" in the SPECIMEN_TYPE field
w <- which(astoDF$SPECIMEN_TYPE %in% c('PERIPHERAL BLOOD', 'BLOOD CULTURE FROM LINE DRAW', 'BLOOD', 'BLOOD FROM BROVIAC', 'BLOOD FROM TRIPLE LUMEN', 'CATHETER BLOOD',
                                       'BLOOD FROM ARTERIAL LINE', 'WHOLE BLOOD', 'HICKMAN BLOOD', 'BLOOD COLLECTED BY UNIT PERSONNEL')
           & (is.na(astoDF$ORDER_NAME) | astoDF$ORDER_NAME == 'CULTURE, AEROBIC BACTERIA'))
astoDF$BLOOD[w] <- TRUE


# Finally, if there is a missing ORDER_NAME and SPECIMEN_TYPE did not indicate blood
# astoDF %>%
#    filter(is.na(ORDER_NAME) & is.na(SPECIMEN_TYPE)) %>%
#    mutate(RESULT_VALUE = gsub('  ', ' ', tolower(RESULT_VALUE))) %>%
#    count(RESULT_VALUE, sort=TRUE) %>%
#    filter(stringi::stri_detect_regex(pattern='blood|peripheral blood|blood culture', str=RESULT_VALUE, case_insensitive=TRUE),
#           stringi::stri_detect_regex(pattern='cell|bloody|coloniz|contamination|abscess|drainage|ascites', str=RESULT_VALUE, case_insensitive=TRUE, negate=TRUE)) %>%
#    filter(n > 10L)
w <- which(is.na(astoDF$ORDER_NAME) & is.na(astoDF$SPECIMEN_TYPE)) # 522,393
w2 <- intersect(which(stringi::stri_detect_regex(pattern='blood|peripheral blood|blood culture',                       str=astoDF$RESULT_VALUE[w], case_insensitive=TRUE)),
                which(stringi::stri_detect_regex(pattern='cell|bloody|coloniz|contamination|abscess|drainage|ascites', str=astoDF$RESULT_VALUE[w], case_insensitive=TRUE, negate=TRUE))) # 4,965
w2 <- w[w2]
astoDF$BLOOD[w2] <- TRUE


rm(orders, specs, w, w2)



astoDF %>% count(is.na(ORDER_DATE)) # 8,656 (out of 7.7 million) - has RESULT_DATE, most SPECIMEN_COLLECTED_DATE, SPECIMEN_RECEIVED_DATE
astoDF <- astoDF %>% 
   filter(!is.na(ORDER_DATE)) %>%
   select(!c(REFERENCE_UNIT, REFERENCE_LOW, REFERENCE_HIGH, COMPONENT_NAME, RESULT_LAB_NAME, SPECIMEN_COLLECTED_DATE, 
             SPECIMEN_RECEIVED_DATE, LOINC_CODE, ORDER_NAME, SPECIMEN_TYPE, RESULT_VALUE)) %>%
   distinct() %>% # 7,702,107 --> 6.6 million
   arrange(PERSON_ID, ORDER_DATE, ORDER_PROC_ID)

# some order_proc_ids are blood, but blood wasn't assigned TRUE for all rows of that ID
blood_order_ids <- unique(astoDF$ORDER_PROC_ID[astoDF$BLOOD])
length(blood_order_ids) # 131,748
length(unique(astoDF$ORDER_PROC_ID)) # 1,550,181
table(astoDF$BLOOD[astoDF$ORDER_PROC_ID %in% blood_order_ids]) # 20,254 are FALSE!
astoDF$BLOOD[astoDF$ORDER_PROC_ID %in% blood_order_ids] <- TRUE # no longer
rm(blood_order_ids)

gc()


# PATHOGEN NAMES
source(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/CleanPathogenNames/CleanPathogenNames.R')
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
   distinct() # 6,560,904 --> 3,023,195
rm(path_names, cleanPathogenNames)

# CLEANING
# remove instances where ORDER_DATE comes after RESULT_DATE - clear error
astoDF <- astoDF %>% filter(RESULT_DATE > ORDER_DATE) # 3,023,195 --> 3,019,615

# sometimes everything is identical, but one row has a MUCH later RESULT_DATE, take the minimum
start <- Sys.time()
astoDF <- astoDF %>% group_by(PERSON_ID, ORDER_PROC_ID, ORDER_DATE, BLOOD, BUG)
astoDF1 <- astoDF %>% filter(n() == 1L)
astoDF2 <- astoDF %>% filter(n() > 1L)
astoDF2 <- astoDF2 %>% slice_min(RESULT_DATE)
astoDF <- rbind(astoDF1, astoDF2) %>% arrange(PERSON_ID, ORDER_PROC_ID, ORDER_DATE, RESULT_DATE)
rm(astoDF1, astoDF2); gc()
print(Sys.time() - start) # ~2 minutes

# sometimes everything is identical, but one row has a MUCH early ORDER_DATE, take the maximum
start <- Sys.time()
astoDF <- astoDF %>% group_by(PERSON_ID, ORDER_PROC_ID, RESULT_DATE, BLOOD, BUG)
astoDF1 <- astoDF %>% filter(n() == 1L) %>% ungroup()
astoDF2 <- astoDF %>% filter(n() > 1L)
astoDF2 <- astoDF2 %>% slice_max(ORDER_DATE) %>% ungroup()
astoDF <- rbind(astoDF1, astoDF2)
rm(astoDF1, astoDF2); gc()
print(Sys.time() - start) # ~1 minute

# ADD BACK A SINGLE ROW PER ORDER WITHOUT PATH_NAME IF IT IS MISSING A PATH_NAME
astoDF <- rbind(astoDF,
                astoDF %>% 
                   filter(all(!is.na(BUG)), 
                          .by = c(PERSON_ID, ORDER_DATE, ORDER_PROC_ID, RESULT_DATE)) %>% 
                   mutate(BUG = NA_character_)) # 3,074,669
astoDF <- astoDF %>% arrange(PERSON_ID, ORDER_PROC_ID, ORDER_DATE, RESULT_DATE)


save(astoDF, file = '~/Desktop/EHR/EHR work/RdataFiles/AST_orders_clean.Rdata')




