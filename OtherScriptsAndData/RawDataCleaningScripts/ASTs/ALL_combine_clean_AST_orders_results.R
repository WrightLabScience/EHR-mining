###############################################################
################# JOIN AST ORDERS AND RESULTS ################# # ~5 minute run-time
###############################################################
start <- Sys.time()
library(dplyr)
library(tidyr)

load(file = '~/Desktop/EHR/EHR work/RdataFiles/AST_orders_clean.Rdata')  # 3,074,667
load(file = '~/Desktop/EHR/EHR work/RdataFiles/AST_results_clean.Rdata') # 1,789,478

astrDF <- astrDF %>% rename(RESULT_DAY = RESULT_DATE)
astoDF <- astoDF %>% 
   mutate(RESULT_DAY = as.Date(substr(RESULT_DATE, 1, 10))) %>%
   select(!c(ORDER_NAME, COMPONENT_NAME, RESULT_VALUE, REFERENCE_UNIT, REFERENCE_LOW, REFERENCE_HIGH, RESULT_LAB_NAME, 
             SPECIMEN_COLLECTED_DATE, SPECIMEN_RECEIVED_DATE, SPECIMEN_TYPE, LOINC_CODE, SITE_SPECIMEN, SITE_ORDER)) %>%
   distinct()



# first join on BUG for those orders where bug is present 
astDF <- left_join(x = astrDF,
                   y = astoDF %>% filter(!is.na(BUG)), # order HAS bug - 1/2
                   by = join_by(PERSON_ID, ORDER_PROC_ID, BUG, RESULT_DAY)) %>%
   relocate(ORDER_DATE, RESULT_DATE, RESULT_DAY, SITE, .after=ORDER_PROC_ID)

# then join without bug on the remaining records in ast that do not have an order_date
# with those orders in asto that do not have a bug
astrDF_unM <- astDF %>%
   filter(is.na(ORDER_DATE)) %>%
   select(!c(ORDER_DATE, RESULT_DATE, SITE))
astoDF_unM <- astoDF %>% # 1,545,343 <-- 1,551,609 (out of 3,074,667)
   filter(is.na(BUG)) %>% 
   select(-BUG) %>%
   distinct()

astDF_unM <- left_join(x = astrDF_unM,
                       y = astoDF_unM,
                       by = join_by(PERSON_ID, ORDER_PROC_ID, RESULT_DAY)) %>%
   relocate(ORDER_DATE, RESULT_DATE, RESULT_DAY, SITE, .after=ORDER_PROC_ID)

astDF_unM <- astDF_unM %>% filter(!is.na(ORDER_DATE))
astDF_M   <- astDF     %>% filter(!is.na(ORDER_DATE))

# combine initially matched and unmatched rows
astDF <- rbind(astDF_M, astDF_unM)
# astDF <- astDF %>% arrange(PERSON_ID, ORDER_DATE, RESULT_DATE)

rm(astDF_M, astDF_unM, astoDF_unM, astrDF_unM)
gc()

rm(astoDF, astrDF)

# remove unnecessary variables that may obscure effective duplicates
astDF <- astDF %>% 
   select(-ORDER_PROC_ID, -RESULT_DAY, -PATH_NAME) %>% 
   distinct()


### If same bug, same AST, same day, take minimum order time and minimum result date-time
# but first, how common is this?
astDF <- astDF %>%
   mutate(ORDER_DAY = lubridate::as_date(ORDER_DATE),
          RESULT_DAY = lubridate::as_date(RESULT_DATE)) %>%
   relocate(ORDER_DAY, RESULT_DAY, .after=RESULT_DATE)

# if same order_day
astDF <- astDF %>% # 1,858,510
   group_by_all() %>%
   ungroup(ORDER_DATE, RESULT_DATE, RESULT_DAY) # same order_day

astDFm <- astDF %>% filter(n() > 1L)
astDFm <- astDFm %>%
   summarise(RESULT_DAY = min(RESULT_DAY),
             ORDER_DATE = min(ORDER_DATE),
             RESULT_DATE = min(RESULT_DATE)) %>%
   ungroup()
astDF <- rbind(
   astDF %>% filter(n() == 1L) %>% ungroup(), 
   astDFm
) # 1,823,110 rows
rm(astDFm); gc()


astDF <- astDF %>%
   group_by_all() %>%
   ungroup(ORDER_DATE, RESULT_DATE, ORDER_DAY) # same result_day

# if same result_day, take minimum order date-time and minimum result date-time
astDFm <- astDF %>% filter(n() > 1L) # 8,937 rows if second (59,634 if we do this one first)
astDFm <- astDFm %>% # 4,444 rows if second (28,321 if we do this one first)
   summarise(ORDER_DAY = min(ORDER_DAY),
             ORDER_DATE = min(ORDER_DATE),
             RESULT_DATE = min(RESULT_DATE)) %>%
   ungroup()
astDF <- rbind(
   astDF %>% filter(n() == 1L) %>% ungroup(),
   astDFm
)
rm(astDFm); gc()


nrow(astDF) # 1,819,539

astDF <- astDF %>% arrange(PERSON_ID, ORDER_DATE, RESULT_DATE)


###################################################################################
save(astDF, file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_clean_ASTs.Rdata')
###################################################################################
print(Sys.time() - start) # ~6-8 minutes








