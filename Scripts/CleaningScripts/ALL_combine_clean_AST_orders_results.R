###############################################################
################# JOIN AST ORDERS AND RESULTS ################# # 3.4 minute run-time
###############################################################
start <- Sys.time()
library(dplyr)
library(tidyr)

load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/data_path_name.Rdata')
load(file = paste0(data_path_name, 'AST_orders_clean.Rdata'))  # 3,074,667
load(file = paste0(data_path_name, 'AST_results_clean.Rdata')) # 1,789,478

astrDF <- astrDF %>% rename(RESULT_DAY = RESULT_DATE)
astoDF <- astoDF %>% mutate(RESULT_DAY = as.Date(substr(RESULT_DATE, 1, 10)))

length(astoDF$ORDER_PROC_ID[which(is.na(astoDF$BUG))])
length(astoDF$ORDER_PROC_ID[which(!is.na(astoDF$BUG))])
length(astrDF$ORDER_PROC_ID[which(is.na(astrDF$BUG))])
length(astrDF$ORDER_PROC_ID[which(!is.na(astrDF$BUG))])

length(unique(astrDF$ORDER_PROC_ID[which(is.na(astrDF$BUG))]))
length(intersect(astrDF$ORDER_PROC_ID[which(is.na(astrDF$BUG))],
                 astoDF$ORDER_PROC_ID[which(!is.na(astoDF$BUG))])) # 8,257 ! for AST results that are missing BUG, a subset of AST orders have the BUG NAME!


# first join on BUG for those orders where bug is present 
astDF <- left_join(x = astrDF,
                   y = astoDF %>% filter(!is.na(BUG)), # order HAS bug - 1,523,058
                   by = join_by(PERSON_ID, ORDER_PROC_ID, BUG, RESULT_DAY)) %>%
   relocate(ORDER_DATE, RESULT_DATE, BLOOD, .after=ORDER_PROC_ID) %>%
   relocate(RESULT_DAY, .after=RESULT_DATE)

astDF %>% count(is.na(ORDER_DATE)) # 366,253 still need an order date

# then join without bug on the remaining records in ast that do not have an order_date
# with those orders in asto that do not have a bug
astrDF_unM <- astDF %>% # 366,253
   filter(is.na(ORDER_DATE)) %>%
   select(!c(ORDER_DATE, RESULT_DATE, BLOOD))
astoDF_unM <- astoDF %>% # 1,545,343 <-- 1,551,609 (out of 3,074,667)
   filter(is.na(BUG)) %>% 
   select(-BUG) %>%
   distinct()

astDF_unM <- left_join(x = astrDF_unM,
                       y = astoDF_unM,
                       by = join_by(PERSON_ID, ORDER_PROC_ID, RESULT_DAY)) %>%
   relocate(ORDER_DATE, RESULT_DATE, BLOOD, .after=ORDER_PROC_ID) %>%
   relocate(RESULT_DAY, .after=RESULT_DATE)

astDF_unM <- astDF_unM %>% filter(!is.na(ORDER_DATE)) # 356,225 (from 366,253)
astDF_M   <- astDF     %>% filter(!is.na(ORDER_DATE)) # 1,443,587

# combine initially matched and unmatched rows = 1,846,711
astDF <- rbind(astDF_M, astDF_unM) # 1,850,758
astDF <- astDF %>% arrange(PERSON_ID, ORDER_DATE, RESULT_DATE)

rm(astDF_M, astDF_unM, astoDF_unM, astrDF_unM)
gc()

rm(astoDF, astrDF)

# remove unnecessary variables that may obscure effective duplicates
astDF <- astDF %>% 
   select(-ORDER_PROC_ID, -RESULT_DAY, -PATH_NAME) %>% 
   distinct() # 1,835,382


### If same bug, same AST, same day, take minimum order time and minimum result date-time
# but first, how common is this?
astDF <- astDF %>% # 1,835,382
   mutate(ORDER_DAY = as.Date(substr(ORDER_DATE,1,10)),
          RESULT_DAY = as.Date(substr(RESULT_DATE,1,10))) %>%
   relocate(ORDER_DAY, RESULT_DAY, .after=RESULT_DATE)


# if same order_day
astDF <- astDF %>% # 1,804,069 (1,789,627 groups) if SECOND (1,794,092 groups if FIRST)
   group_by_all() %>%
   ungroup(ORDER_DATE, RESULT_DATE, RESULT_DAY) # same order_day

astDFm <- astDF %>% filter(n() > 1L) # 77,923 rows if first (28,494 if we do this one second)
astDFm <- astDFm %>% # 36,633 rows if first (14,052 if we do this one second)
   summarise(RESULT_DAY = min(RESULT_DAY),
             ORDER_DATE = min(ORDER_DATE),
             RESULT_DATE = min(RESULT_DATE)) %>%
   ungroup()

astDF <- astDF %>% filter(n() == 1L) %>% ungroup()
astDF <- rbind(astDF, astDFm) %>% arrange(PERSON_ID, ORDER_DATE, RESULT_DATE)
rm(astDFm)


astDF <- astDF %>% # 1,835,382 rows (1,804,069 groups) if first (1,789,599 groups if second)
   group_by_all() %>%
   ungroup(ORDER_DATE, RESULT_DATE, ORDER_DAY) # same result_day

# if same result_day, take minimum order date-time and minimum result date-time
astDFm <- astDF %>% filter(n() > 1L) # 8,937 rows if second (59,634 if we do this one first)
astDFm <- astDFm %>% # 4,444 rows if second (28,321 if we do this one first)
   summarise(ORDER_DAY = min(ORDER_DAY),
             ORDER_DATE = min(ORDER_DATE),
             RESULT_DATE = min(RESULT_DATE)) %>%
   ungroup()

astDF <- astDF %>% filter(n() == 1L) %>% ungroup()
astDF <- rbind(astDF, astDFm) %>% arrange(PERSON_ID, ORDER_DATE, RESULT_DATE)
rm(astDFm)


nrow(astDF) 
# 1,789,627 rows (if I do same result_day first, then same order_day)
# 1,789,599 rows (if I do same order_day first, then same result_day)


###################################################################################
save(astDF, file = paste0(data_path_name, '/ALL_clean_ASTs.Rdata'))
###################################################################################
print(Sys.time() - start) # 3.4 minutes








