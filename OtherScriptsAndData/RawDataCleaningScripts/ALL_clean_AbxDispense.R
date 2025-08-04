###############################################################
############ READ IN AND CLEAN MED DISPENSE TABLE #############
###############################################################
load(file = '~/Desktop/EHR/EHR work/data/drug_abx_ndc.Rdata')

abxDF <- tbl(conn, in_schema('AMB_ETL', 'SENS_MED_DISPENSE_VW')) %>%
   filter(DRUG_NAME %in% local(ndc_abx$DRUG_NAME[1:1000])) %>%
   collect()

for (i in c(1001, 2001, 3001)) {
   print(i)
   end <- ifelse(i == 3001, nrow(ndc_abx), i+999)
   abxDF <- rbind(abxDF,
                  tbl(conn, in_schema('AMB_ETL', 'SENS_MED_DISPENSE_VW')) %>%
                     filter(DRUG_NAME %in% local(ndc_abx$DRUG_NAME[i:end])) %>%
                     collect())
}
rm(i, end)

# prep drug_name --> abx "dictionary"
abx <- setNames(ndc_abx$ABX, ndc_abx$DRUG_NAME)

# simple modifications
abxDF <- abxDF %>%
   mutate(ABX = abx[DRUG_NAME],
          LAST_DISPENSE_DATE = paste0(LAST_DISPENSE_DATE, ' 12:00:00'),
          QUANTITY_DISPENSED = as.integer(QUANTITY_DISPENSED)) %>%
   mutate(LAST_DISPENSE_DATE = strptime(LAST_DISPENSE_DATE, format='%m/%d/%Y %T')) %>%
   arrange(PERSON_ID, LAST_DISPENSE_DATE) %>%
   select(-NDC, -AMOUNT_REMAINING, -SUPPLY_ORDERED)

# handle multiple drugs
lens <- lengths(strsplit(abxDF$ABX, split=', '))
abxDF$LEN <- lens; rm(lens)
w3 <- which(abxDF$LEN == 3)
abxDF3 <- abxDF[w3,] %>% select(-LEN)
abxDF3_1 <- abxDF3_2 <- abxDF3_3 <- abxDF3
abxs <- strsplit(abxDF3$ABX, split=', ')
abxDF3_1$ABX <- sapply(abxs, '[', 1)
abxDF3_2$ABX <- sapply(abxs, '[', 2)
abxDF3_3$ABX <- sapply(abxs, '[', 3)
abxDF3 <- rbind(abxDF3_1, abxDF3_2, abxDF3_3)
rm(abxDF3_1, abxDF3_2, abxDF3_3, abxs)

w2 <- which(abxDF$LEN == 2)
abxDF2 <- abxDF[w2,] %>% select(-LEN)
abxDF2_1 <- abxDF2_2 <- abxDF2
abxs <- strsplit(abxDF2$ABX, split=', ')
abxDF2_1$ABX <- sapply(abxs, '[', 1)
abxDF2_2$ABX <- sapply(abxs, '[', 2)
abxDF2 <- rbind(abxDF2_1, abxDF2_2)
rm(abxDF2_1, abxDF2_2, abxs)

abxDF <- rbind(abxDF[-union(w2, w3), ] %>% select(-LEN),
               abxDF2,
               abxDF3)
rm(w2, w3, abxDF2, abxDF3)


# handle abbreviations
abbr <- ndc_abx %>% select(ABX, ABX_ABBR) %>% distinct()
abbr <- setNames(abbr$ABX_ABBR, abbr$ABX)
abxDF <- abxDF %>%
   mutate(ABX_ABBR = abbr[ABX]) %>%
   mutate(ABX_ABBR = ifelse(is.na(ABX_ABBR), ABX, ABX_ABBR)) %>%
   relocate(ABX, ABX_ABBR, .after=DRUG_NAME) %>%
   arrange(PERSON_ID, LAST_DISPENSE_DATE)
rm(abx, abbr, ndc_abx, all_antibiotics, ASTcounts)




##########################################################
############# INFER ABX THERAPY INTERVALS ################
##########################################################
abxDF %>%
   count(FREQUENCY, sort=TRUE) %>%
   print(n=100)

# impute quantity and frequency when missing or text
abxDF <- abxDF %>%
   mutate(QUANTITY_DISPENSED = case_when(
      grepl('ONCE|1XONLY', FREQUENCY) ~ 1,
      QUANTITY_DISPENSED == 0 ~ 1,
      .default = QUANTITY_DISPENSED
   )) %>%
   mutate(TPD = case_when(
      FREQUENCY == 'Daily' ~ 1,
      grepl('ONCE|1XONLY|QPM|QAM|AtBedtime', FREQUENCY) ~ 1,
      grepl('BID', FREQUENCY) ~ 2,
      grepl('TID', FREQUENCY) ~ 3,
      grepl('QID', FREQUENCY) ~ 4,
      FREQUENCY == 'QMndWedFri' ~ 3/7,
      FREQUENCY == 'QTuThSat' ~ 3/7,
      .default = NA
   ))
w <- grep('^Q[0-9]+H$', abxDF$FREQUENCY)
abxDF$TPD[w] <- 24 / as.numeric(gsub('^Q([0-9]+)H$', '\\1', abxDF$FREQUENCY[w]))
w <- grep('^Q[0-9]+HRS$', abxDF$FREQUENCY)
abxDF$TPD[w] <- 24 / as.numeric(gsub('^Q([0-9]+)HRS$', '\\1', abxDF$FREQUENCY[w]))
w <- grep('[0-9]XD.*', abxDF$FREQUENCY)
abxDF$TPD[w] <- as.numeric(gsub('([0-9])XD.*', '\\1', abxDF$FREQUENCY[w]))
w <- grep('^Q[0-9]+Min$', abxDF$FREQUENCY)
abxDF$TPD[w] <- 1440 / as.numeric(gsub('^Q([0-9]+)Min$', '\\1', abxDF$FREQUENCY[w]))
w <- which(grepl('CPOE', abxDF$FREQUENCY) & abxDF$QUANTITY_DISPENSED <= 4)
abxDF$TPD[w] <- 1
rm(w)

# Impute END_DATE from "START_DATE", FREQUENCY AND QUANTITY
abxDF <- abxDF %>%
   mutate(DISPENSE_DATE_2 = case_when(
      TPD == 1 ~ LAST_DISPENSE_DATE + 86400*((QUANTITY_DISPENSED / TPD) - 1),
      TPD != 1 ~ LAST_DISPENSE_DATE + 86400*(QUANTITY_DISPENSED / TPD)
   )) %>%
   relocate(DISPENSE_DATE_2, .after=LAST_DISPENSE_DATE) %>%
   rename(START_DATE = LAST_DISPENSE_DATE,
          END_DATE = DISPENSE_DATE_2)

# FLAG MISSING DATES
#abxDF %>% count(is.na(START_DATE), is.na(END_DATE))
w <- which(is.na(abxDF$END_DATE))
abxDF$END_DATE[w] <- abxDF$START_DATE[w]
rm(w)

abxDF$LEN_OF_THERAPY <- as.numeric(lubridate::as.duration(abxDF$END_DATE - abxDF$START_DATE)) / 86400
abxDF <- abxDF %>% filter(LEN_OF_THERAPY < 150)

# x <- as.numeric(lubridate::as.duration(abxDF$END_DATE - abxDF$START_DATE)) / 86400
# range(x) # 0, 187,597.7
# median(x) # 0
# mean(x) # 2.85
# hist(x, breaks=diff(range(x)), xlim=c(0, 25))
# rm(x)

abxDF <- abxDF %>%
   select(-DRUG_NAME, -QUANTITY_DISPENSED, -FREQUENCY, -TPD, -LEN_OF_THERAPY) %>%
   distinct() %>%
   arrange(PERSON_ID, ABX, START_DATE, END_DATE)

print(object.size(abxDF), units='Mb') # 538.2 Mb

# abxDF <- abxDF %>%
#    mutate(SAME_ABX = ABX == lag(ABX)) %>%
#    mutate(OVERLAP = case_when(
#       !SAME_ABX ~ FALSE,
#       SAME_ABX ~ START_DATE == lag(END_DATE),
#    ))

gc()

save(abxDF, file = '~/Desktop/EHR/EHR work/ALL/abxDF_all.Rdata')







