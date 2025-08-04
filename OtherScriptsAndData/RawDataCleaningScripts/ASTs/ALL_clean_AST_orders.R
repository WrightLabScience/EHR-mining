###################################################
############# CLEAN UP AST ORDER TABLE ############ Total run-time = ~ 10 minutes probable
###################################################
source('~/Desktop/EHR/EHR work/config_file.R')
astoDF <- tbl(conn, in_schema('AMB_ETL', 'LAB_MICRO_RESULT_ALL_VW')) %>% collect()


# site of infection -------------------------------------------------------

stn <- astoDF %>% count(SPECIMEN_TYPE, sort=TRUE)
stn <- stn %>%
   mutate(SITE = case_when(
      grepl('URIN|UREOSTOMY|URETEROSTOMY|MIDSTREAM|PASSPORT', SPECIMEN_TYPE) ~ 'urine',
      SPECIMEN_TYPE == 'CATHETER URINE' ~ 'urine',
      grepl('BILIARY DRAINAGE|G(-| )TUBE|ABDOMINAL|T-TUBE|ASCITES|PERITONEAL DIALYSATE|PERITONEAL LAVAGE|PERITONEAL TISSUE', SPECIMEN_TYPE) ~ 'body_fluid',
      grepl('ASCITES|FLUID|ABDOMEN|BILE|^PUS$|GALLBLADDER', SPECIMEN_TYPE) ~ 'body_fluid',
      grepl('CSF|CNS', SPECIMEN_TYPE) ~ 'csf',
      grepl('SPUTUM|TRACH|RESP|THROAT|BRONCH|BAL|BRUSH|CHEST TUBE|ASPIRATE|LUNG', SPECIMEN_TYPE) ~ 'respiratory',
      grepl('INCISION|LESION|SUTURE LINE|SEROMA|SURGICAL|GRAFT|JACKSON PRATT|JP SITE', SPECIMEN_TYPE) ~ 'surgical',
      grepl('GENITAL|VAGINAL|CERVIX|VULVA|LABIA|PENIS|PLACENTA|GROIN|URETHRA|BARTHOLIN|PERINEUM|UMBILICAL|SEMEN', SPECIMEN_TYPE) ~ 'genital',
      grepl('CATHETER', SPECIMEN_TYPE) ~ 'catheter',
      SPECIMEN_TYPE == 'PERMA-CATH' ~ 'catheter',
      grepl('BLOOD', SPECIMEN_TYPE) ~ 'blood',
      grepl('JP.+DRAINAGE|WOUND|ABSCESS', SPECIMEN_TYPE) ~ 'wound',
      SPECIMEN_TYPE %in% c('DRAINAGE', 'TISSUE', 'ABCESS') ~ 'wound',
      grepl('SKIN|BUTTOCKS|TOE|LEG|FOOT|FINGER|KNEE|BREAST|AXILLA|SCALP|ANKLE|ELBOW|HIP|ARM|HEEL|THIGH|HAND|BACK|COCCYX|SHOULDER|NECK|WRIST|CHEST|NAIL|CELLULITIS|PELVIS', SPECIMEN_TYPE) ~ 'skin',
      grepl('NOSE|NASAL|NASOPHARYNGEAL|SINUS|MEATUS|NARES', SPECIMEN_TYPE) ~ 'nose',
      grepl('EAR', SPECIMEN_TYPE) ~ 'ear',
      grepl('OPHTHAL|EYE|CORNEA|CONJUNCTIVA|VITREOUS|LENS|AQUEOUS|LACRIMAL', SPECIMEN_TYPE) ~ 'eye',
      grepl('STOOL|RECTAL|JEJUNAL|COLOSTOMY|ILEOSTOMY|ULCER|PERI-ANAL|UREOSTOMY|ANUS|RECTUM|APPENDIX', SPECIMEN_TYPE) ~ 'intestinal',
      grepl('BONE|FEMORAL|SACRAL|TIBIA|SACRUM', SPECIMEN_TYPE) ~ 'bone',
      grepl('ORAL|MOUTH|DENTAL|TONGUE|TONSIL', SPECIMEN_TYPE) ~ 'oral',
      SPECIMEN_TYPE %in% c('EXIT SITE', 'PORT SITE') ~ 'catheter',
      SPECIMEN_TYPE == 'STRAIGHT CATH' ~ 'urine',
      grepl('PROSTHETIC', SPECIMEN_TYPE) ~ 'prosthetic',
      SPECIMEN_TYPE == 'CARDIAC PACEMAKER OR LEADS' ~ 'pacemaker',
      grepl('HEMATOMA', SPECIMEN_TYPE) ~ 'hematoma',
      grepl('LYMPH', SPECIMEN_TYPE) ~ 'lymph',
      grepl('SYNOVIAL', SPECIMEN_TYPE) ~ 'joint',
      SPECIMEN_TYPE %in% c('KIDNEY', 'KIDNEY BIOPSY') ~ 'kidney',
      SPECIMEN_TYPE == 'CYST' ~ 'cyst',
      is.na(SPECIMEN_TYPE) | SPECIMEN_TYPE %in% c('SWAB FROM UNKNOWN SITE', 'OTHER SOURCE,DEFINE IN COMMENTS', 'BACTERIAL ISOLATE', 'SWAB', 'UNKNOWN SITE', 'MYCOBACTERIAL ISOLATE',
                                                  'DRAW SITE NOT SPECIFIED') ~ 'site unknown',
      SPECIMEN_TYPE %in% c('DONOR', 'SITE', 'BIOPSY', 'STONE') ~ 'site unknown'
   )) %>%
   mutate(SITE = case_when(
      SITE %in% c('eye', 'ear', 'nose') ~ 'eye, ear, nose',
      SITE == 'respiratory' & SPECIMEN_TYPE == 'GASTRIC ASPIRATE' ~ 'intestinal',
      is.na(SITE) ~ 'site unkown',
      .default = SITE
   ))
# stn %>% summarise(n=sum(n), .by=SITE) %>% arrange(desc(n))


ons <- astoDF %>% count(ORDER_NAME, sort=TRUE)
ons <- ons %>%
   mutate(SITE = case_when(
      grepl('URIN', ORDER_NAME) ~ 'urine',
      grepl('BLOOD', ORDER_NAME) & !ORDER_NAME %in% c('PERITONEAL DIALYSATE IN BLOOD (OO)', 'FUNGUS CULTURE, NON BLOOD', 'AUTOPSY BLOOD CULTURE (OO)', 'BODY FLUID IN BLOOD CULTURE BOTTLES') ~ 'blood',
      grepl('WOUND|DRAINAGE|ABSCESS', ORDER_NAME) ~ 'wound',
      grepl('SKIN', ORDER_NAME) ~ 'skin',
      grepl('SPUTUM|TRACH|RESP|THROAT|BRONCH|BAL|BRUSH|CHEST TUBE|ASPIRATE|LUNG', ORDER_NAME) ~ 'respiratory',
      ORDER_NAME == 'FLUID/ASPIRATE CULTURE' ~ 'respiratory',
      grepl('SURGICAL', ORDER_NAME) ~ 'surgical tissue',
      grepl('OPHTHAL|EYE|EAR|NOSE|NASAL|NASOPHARYNGEAL', ORDER_NAME) ~ 'eye, ear, nose',
      grepl('STOOL|RECTAL', ORDER_NAME) ~ 'stool',
      grepl('GENITAL', ORDER_NAME) ~ 'genital',
      grepl('BONE', ORDER_NAME) ~ 'bone',
      grepl('CSF', ORDER_NAME) ~ 'csf',
      grepl('ORAL|MOUTH', ORDER_NAME) ~ 'oral',
      grepl('CATHETER', ORDER_NAME) ~ 'catheter',
      grepl('PROSTHESIS', ORDER_NAME) ~ 'prosthetic',
      ORDER_NAME %in% c('BODY FLUID CULTURE & GRAM STAIN', 'BODY FLUID CULTURE - AEROBIC', 'AMNIOTIC FLUID CULT(MAGEE ONLY)', 'SEMEN CULTURE', 'BREAST MILK CULTURE(MWH)') ~ 'body_fluid',
      #ORDER_NAME %in% c('CULTURE, AEROBIC BACTERIA', 'CULTURE,AEROBIC (W/GRAM STAIN)', 'AEROBIC CULTURE (OO)', 'CULTURE, AEROBIC BACTERIA WITH GRAM STAIN', 'AUTOPSY AEROBIC CULTURE (OO)') ~ 'aerobic',
      #ORDER_NAME %in% c('ANAEROBIC CULTURE', 'ANAEROBIC CULTURE (OO)', 'CULTURE, ANAEROBIC BACTERIA WITH GRAM STAIN', 'AUTOPSY ANAEROBIC CULTURE (OO)') ~ 'anaerobic',
      .default = 'site not specified'
   ))
# ons %>% summarise(n=sum(n), .by=SITE) %>% arrange(desc(n))


stn <- setNames(stn$SITE, stn$SPECIMEN_TYPE)
ons <- setNames(ons$SITE, ons$ORDER_NAME)
astoDF <- astoDF %>%
   mutate(SITE_SPECIMEN = unname(stn[SPECIMEN_TYPE]),
          SITE_ORDER = unname(ons[ORDER_NAME]))
# rm(stn, ons)

###############


# make date-time column
# make "pathogen name" column (which occasional is mentioned in the result_value column)
astoDF <- astoDF %>%
   mutate(across(contains('DATE'), ~ strptime(., format='%m/%d/%Y %T')),
          PATH_NAME = stringr::str_to_lower(gsub('  ', ' ', RESULT_VALUE)))

# TODO: provide multiple columns - one for each common isolation source
keep_sites <- c('urine', 'wound', 'blood', 'respiratory', 'body_fluid', 'skin')
astoDF$SITE <- NA_character_

# Start with keep sites as indicated by SITE_ORDER
astoDF$SITE[astoDF$SITE_ORDER %in% keep_sites] <- astoDF$SITE_ORDER[astoDF$SITE_ORDER %in% keep_sites]

# Fill in gaps where SITE_SPECIMEN indicated a keep site
w <- which(is.na(astoDF$SITE_ORDER) & astoDF$SITE_SPECIMEN %in% keep_sites)
astoDF$SITE[w] <- astoDF$SITE_SPECIMEN[w]
rm(w)

# Fill in NAs in SITE with either 'other'
astoDF$SITE[is.na(astoDF$SITE)] <- 'other'


# astoDF <- astoDF %>%
#    mutate(BLOOD = (!is.na(SITE_SPECIMEN) & SITE_SPECIMEN == 'blood') | (!is.na(SITE_ORDER) & SITE_ORDER == 'blood'),
#           RESPIRATORY = (!is.na(SITE_SPECIMEN) & SITE_SPECIMEN == 'respiratory') | (!is.na(SITE_ORDER) & SITE_ORDER == 'respiratory'))
# 
# # Finally, if there is a missing ORDER_NAME and SPECIMEN_TYPE did not indicate blood
# w <- which((astoDF$SITE_ORDER == 'site not specified' | is.na(astoDF$SITE_ORDER)) & (astoDF$SITE_SPECIMEN == 'site unknown' | is.na(astoDF$SITE_SPECIMEN))) # 651,242
# w2 <- intersect(which(stringi::stri_detect_regex(pattern='blood|peripheral blood|blood culture',                       str=astoDF$RESULT_VALUE[w], case_insensitive=TRUE)),
#                 which(stringi::stri_detect_regex(pattern='cell|bloody|coloniz|contamination|abscess|drainage|ascites', str=astoDF$RESULT_VALUE[w], case_insensitive=TRUE, negate=TRUE))) # 5,185
# w2 <- w[w2]
# astoDF$BLOOD[w2] <- TRUE
# rm(w, w2)
# 
# 
# w <- which((astoDF$SITE_ORDER == 'site not specified' | is.na(astoDF$SITE_ORDER)) & (astoDF$SITE_SPECIMEN == 'site unknown' | is.na(astoDF$SITE_SPECIMEN))) # 651,242
# w2 <- intersect(which(stringi::stri_detect_regex(pattern='SPUTUM|TRACH|RESP|THROAT|BRONCH|BAL|BRUSH|CHEST TUBE|ASPIRATE|LUNG', str=astoDF$RESULT_VALUE[w], case_insensitive=TRUE)),
#                 which(stringi::stri_detect_regex(pattern='representative of the lower respiratory tract', str=astoDF$RESULT_VALUE[w], case_insensitive=TRUE, negate=TRUE))) # 3,380
# w2 <- w[w2]
# astoDF$RESPIRATORY[w2] <- TRUE
# rm(w, w2)



astoDF %>% count(is.na(ORDER_DATE)) # 10,312 (out of 8.1 million) - has RESULT_DATE, most SPECIMEN_COLLECTED_DATE, SPECIMEN_RECEIVED_DATE
astoDF <- astoDF %>% 
   filter(!is.na(ORDER_DATE)) %>%
   #select(PERSON_ID, ORDER_DATE, RESULT_DATE, ORDER_PROC_ID, PATH_NAME, SITE) %>%
   #distinct() %>% # 8,096,195 --> 6,895,485
   arrange(PERSON_ID, ORDER_DATE, ORDER_PROC_ID)

# astoDF <- astoDF %>%
#    mutate(X = 1L) %>%
#    tidyr::pivot_wider(
#       names_from = SITE,
#       values_from = X,
#       values_fill = 0L,
#       names_prefix = 'SITE_'
#    )


# PATHOGEN NAMES
source(file = '~/Desktop/EHR-mining/UsefulDataForCleaning/CleanPathogenNames/CleanPathogenNames.R')
path_names <- astoDF %>%
   count(PATH_NAME, sort=TRUE) %>%
   filter(!is.na(PATH_NAME)) %>%
   mutate(BUG = NA_character_) # 519,749
path_names <- cleanPathogenNames(path_names)
path_names <- path_names %>% filter(!is.na(BUG))
path_names <- setNames(object = path_names$BUG, 
                       nm = path_names$PATH_NAME) # 27,818
astoDF <- astoDF %>%
   mutate(BUG = unname(path_names[PATH_NAME])) %>%
   select(-PATH_NAME) %>%
   distinct() # 6,895,485 --> 3,171,098
rm(path_names, cleanPathogenNames)

# CLEANING
# remove instances where ORDER_DATE comes after RESULT_DATE - clear error
astoDF <- astoDF %>% filter(RESULT_DATE > ORDER_DATE) # 3,167,414

# sometimes everything is identical, but one row has a MUCH later RESULT_DATE, take the minimum
start <- Sys.time()
astoDF <- astoDF %>% group_by(PERSON_ID, ORDER_PROC_ID, ORDER_DATE, SITE, BUG)
astoDF2 <- astoDF %>% filter(n() > 1L)
astoDF2 <- astoDF2 %>% slice_min(RESULT_DATE)
astoDF <- rbind(
   astoDF %>% filter(n() == 1L) %>% ungroup(), 
   astoDF2 %>% ungroup()
) %>% 
   arrange(PERSON_ID, ORDER_PROC_ID, ORDER_DATE, RESULT_DATE, BUG)
rm(astoDF2); gc()
print(Sys.time() - start) # ~2 minutes

# sometimes everything is identical, but one row has a MUCH early ORDER_DATE, take the maximum
start <- Sys.time()
astoDF <- astoDF %>% group_by(PERSON_ID, ORDER_PROC_ID, RESULT_DATE, SITE, BUG)
astoDF2 <- astoDF %>% filter(n() > 1L)
astoDF2 <- astoDF2 %>% slice_max(ORDER_DATE) %>% ungroup()
astoDF <- rbind(
   astoDF %>% filter(n() == 1L) %>% ungroup(), 
   astoDF2
) %>%
   arrange(PERSON_ID, ORDER_PROC_ID, ORDER_DATE, RESULT_DATE, BUG)
rm(astoDF2); gc()
print(Sys.time() - start) # ~1 minute

# ADD BACK A SINGLE ROW PER ORDER WITHOUT PATH_NAME IF IT IS MISSING A PATH_NAME
astoDF <- rbind(astoDF,
                astoDF %>% 
                   filter(all(!is.na(BUG)), 
                          .by = c(PERSON_ID, ORDER_DATE, RESULT_DATE, ORDER_PROC_ID, SITE)) %>% 
                   mutate(BUG = NA_character_)) # 3,223,538
astoDF <- astoDF %>% arrange(PERSON_ID, ORDER_DATE, RESULT_DATE, ORDER_PROC_ID, BUG, SITE)


save(astoDF, file = '~/Desktop/EHR/EHR work/RdataFiles/AST_orders_clean.Rdata')




