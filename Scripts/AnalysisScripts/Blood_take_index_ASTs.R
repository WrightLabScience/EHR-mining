library(dplyr)

load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/data_path_name.Rdata')
load(file = paste0(data_path_name, 'ALL_clean_ASTs.Rdata'))
print(object.size(astDF), units='Mb') # 1,790.4 Mb

# Remove any Coagulase-negative Staph species if they appear in blood cultures
astDF <- astDF %>% filter(!(BLOOD & grepl('Staph', BUG) & !grepl('Staphylococcus aureus', BUG))) # 1,734,028 --> 1,697,680
# Remove candida?
astDF <- astDF %>% filter(!grepl('Cryptococcus|Aspergillus|Candida', BUG)) # 1,691,137


# only take 2015 - 2023
astDF <- astDF %>% filter(substr(ORDER_DATE, 1, 4) %in% as.character(2015:2023)) # 1,164,988
length(unique(astDF$PERSON_ID)) # 440,285


##### REMOVE NON-INDEX CULTURES - FLAG OTHER CULTURES ####
# prep data
astDF <- astDF %>%
   mutate(ORDER_DAY  = as.Date(substr(ORDER_DATE, 1, 10)),
          RESULT_DAY = as.Date(substr(RESULT_DATE,1, 10))) %>%
   relocate(ORDER_DAY, RESULT_DAY, .after=RESULT_DATE)

# mark which rows are a patient's first record
# also which rows are part of a calendar day with multiple isolates on the same day
astDF <- astDF %>%
   mutate(FIRST_RECORD = PERSON_ID != lag(PERSON_ID)) %>%
   mutate(MULT_ISO = (ORDER_DAY == lag(ORDER_DAY) & PERSON_ID == lag(PERSON_ID)) | (ORDER_DAY == lead(ORDER_DAY) & PERSON_ID == lead(PERSON_ID))) %>%
   relocate(FIRST_RECORD, MULT_ISO, .after=RESULT_DAY)
astDF$FIRST_RECORD[1] <- TRUE
astDF$MULT_ISO[1] <- (astDF$ORDER_DAY[1] == lead(astDF$ORDER_DAY)[1] & astDF$PERSON_ID[1] == lead(astDF$PERSON_ID)[1])
astDF$MULT_ISO[nrow(astDF)] <- (astDF$ORDER_DAY[nrow(astDF)] == lag(astDF$ORDER_DAY)[nrow(astDF)]  & astDF$PERSON_ID[nrow(astDF)] == lag(astDF$PERSON_ID)[nrow(astDF)])
sum(astDF$MULT_ISO) / nrow(astDF)     # 32.5%
sum(astDF$FIRST_RECORD) / nrow(astDF) # 36.6%

# calculate days since previous culture
astDF$DAYS_SINCE_PRV_CULTURE <- as.integer(astDF$ORDER_DAY - lag(astDF$ORDER_DAY))
astDF$DAYS_SINCE_PRV_CULTURE[astDF$FIRST_RECORD] <- NA
astDF <- astDF %>% relocate(DAYS_SINCE_PRV_CULTURE, .before=BUG)
# days since previous culture for multiple isolates are mis-marked now (if they're not the first row)
# if a multi-isolate row is 2nd or 3rd in a "series",
# give it the same number of days since previous culture
shift <- 1
w <- which(astDF$PERSON_ID == lead(astDF$PERSON_ID) & astDF$ORDER_DAY == lead(astDF$ORDER_DAY)) # lead() means "shift=1"
while(length(w) > 0) {
   cat(shift, length(w), '\n')
   astDF$DAYS_SINCE_PRV_CULTURE[w + shift] <- astDF$DAYS_SINCE_PRV_CULTURE[w]
   shift <- shift + 1
   w <- w[which(astDF$PERSON_ID[w + shift] == astDF$PERSON_ID[w] & astDF$ORDER_DAY[w + shift] == astDF$ORDER_DAY[w])]
}
# patient's first records for multi-isolates are also screwed up, mark them as TRUE
astDF$FIRST_RECORD[is.na(astDF$DAYS_SINCE_PRV_CULTURE)] <- TRUE


# Indicate whether the current culture has subsequent culture with results BEFORE this one
astDF$DAYS_TO_NEXT_RESULT <- as.numeric(lubridate::as.duration(lead(astDF$RESULT_DATE) - astDF$RESULT_DATE)) / 86400
astDF$DAYS_TO_NEXT_RESULT[lead(astDF$PERSON_ID) != astDF$PERSON_ID] <- NA
astDF <- astDF %>% relocate(DAYS_TO_NEXT_RESULT, .before=BUG)

astDF %>% count(DAYS_TO_NEXT_RESULT < 0)         # ~23K (~2%)
astDF %>% count(BLOOD & DAYS_TO_NEXT_RESULT < 0) # ~11K (~1%)

astDF %>% count(DAYS_TO_NEXT_RESULT == 0, MULT_ISO) # are these behaving correctly?
#   `DAYS_TO_NEXT_RESULT == 0` MULT_ISO      n
# 1 FALSE                      FALSE     ~431K - single isolate with some future culture
# 5 NA                         FALSE     ~380K - single isolate with no future culture
# 3 TRUE                       FALSE        90 - single isolate with identical result_date-time as next culture

# 4 TRUE                       TRUE      ~184K - multiple isolates on same order_day - not last in "series"
# 2 FALSE                      TRUE      ~144K - multiple isolates on same order_day - last in "series" - some future culture
# 6 NA                         TRUE       ~61K - multiple isolates on same order-day - last in "series" - no future culture


# Indicate if positive culture within X days. DAYS_TO_NEXT_CULTURE?
astDF$DAYS_TO_NEXT_CULTURE <- lead(astDF$DAYS_SINCE_PRV_CULTURE)
astDF <- astDF %>% relocate(DAYS_TO_NEXT_CULTURE, .after=DAYS_SINCE_PRV_CULTURE)

# Indicate if multiple blood isolates
astDF <- astDF %>%
   group_by(PERSON_ID, ORDER_DAY, BLOOD) %>%
   mutate(MULT_BLOOD_ISO = BLOOD & n() > 1L) %>%
   ungroup() %>%
   relocate(MULT_BLOOD_ISO, .after=MULT_ISO)



### ignore isolates, collapse to 1 culture per day
df <- astDF %>%
   select(PERSON_ID, FIRST_RECORD, BLOOD, ORDER_DAY, DAYS_SINCE_PRV_CULTURE, MULT_BLOOD_ISO, MULT_ISO) %>% 
   distinct()
df$LAST_CULT_BLOOD <- lag(df$BLOOD)
df$LAST_CULT_BLOOD[df$FIRST_RECORD] <- NA

sum(is.na(df$DAYS_SINCE_PRV_CULTURE[df$BLOOD])) / sum(df$BLOOD) # ~42% of blood culture days are a patients first culture day EVER

med <- median(df$DAYS_SINCE_PRV_CULTURE[df$BLOOD], na.rm=T); med # 40 days is the median number of days since the previous culture
sum(df$DAYS_SINCE_PRV_CULTURE[df$BLOOD] < med, na.rm=T) / sum(df$BLOOD)  # ~29% of blood culture days are < median days removed from a previous culture


# Previously, I removed cultures that were within 6 months of ANY previous culture
# This was my way of ensuring that I take the "first positive (index) blood culture per encounter) - Kadri et al (2021)
# This turns out to be ~43% of blood cultures removed, probably for the "sickest" patients (who have more adjacent blood cultures)
# The AMP team suggest I don't do this, but didn't suggest an alternative
# Perhaps I will only take the first culture in a 14-day window? (this could represent "separate" infections)
# That would be 13% of blood cultures

# Blood cultures tend to be clustered more than other cultures
for (days in c(365, 180, 90, 60, 30, 14, 7, 3)) {
   p1 <- round(sum(df$DAYS_SINCE_PRV_CULTURE[df$BLOOD] < days, na.rm=T) / sum(df$BLOOD) * 100)
   p2 <- round(sum(df$LAST_CULT_BLOOD[df$BLOOD & df$DAYS_SINCE_PRV_CULTURE < days], na.rm=T) / sum(df$BLOOD & df$DAYS_SINCE_PRV_CULTURE < days, na.rm=T) * 100)
   cat(p1, '% of cultures are <', days, ' after a previous culture\n', sep='')
   cat(p2, '% of previous cultures were blood\n\n', sep='')
}



t <- table(df$DAYS_SINCE_PRV_CULTURE, df$BLOOD)
t <- t[1:60, 2:1]
t <- t(t)
barplot(t, main = 'Days since previous culture', legend=TRUE, args.legend=list(title='Blood culture'))

t <- df %>% 
   filter(BLOOD) %>% 
   count(DAYS_SINCE_PRV_CULTURE, LAST_CULT_BLOOD) %>% 
   filter(!is.na(LAST_CULT_BLOOD)) %>% 
   tidyr::pivot_wider(values_from=n, names_from=LAST_CULT_BLOOD) %>%
   rename(lastBlood = `TRUE`, lastNonBlood = `FALSE`) %>%
   mutate(total = lastBlood + lastNonBlood) %>%
   select(-lastNonBlood) %>%
   mutate(pcnt = lastBlood / total) %>%
   filter(DAYS_SINCE_PRV_CULTURE <= 180L)

plot(x = t$DAYS_SINCE_PRV_CULTURE, y = t$pcnt, ylim=c(0, 0.65), ylab='', pch=16, cex=0.5,
     xlab = 'Days since previous culture', main='Proportion previous cultures = blood')
lines(stats::filter(t$pcnt, rep(1/5, 5), sides=2), lwd=2)

df %>% count(BLOOD) # ~61K

1 - sum(df$FIRST_RECORD | df$DAYS_SINCE_PRV_CULTURE > 14L) / nrow(df)                   # ~10%
1 - sum((df$FIRST_RECORD | df$DAYS_SINCE_PRV_CULTURE > 14L) & df$BLOOD) / sum(df$BLOOD) # ~21%

1 - sum(df$FIRST_RECORD | df$DAYS_SINCE_PRV_CULTURE > 7L) / nrow(df)                   # ~7%
1 - sum((df$FIRST_RECORD | df$DAYS_SINCE_PRV_CULTURE > 7L) & df$BLOOD) / sum(df$BLOOD) # ~19%

# THIS ALL MEANS THAT FOR BLOODSTREAM
astDF %>%
   filter(BLOOD) %>%              # 70,312 isolates
   group_by(PERSON_ID, ORDER_DAY) # 59,503 infections
astDF %>%
   filter(BLOOD, DAYS_SINCE_PRV_CULTURE > 14L | FIRST_RECORD) %>% # 55,499 isolates
   group_by(PERSON_ID, ORDER_DAY)                                 # 46,967 infections
astDF %>%
   filter(substr(ORDER_DAY,1,4) %in% as.character(2017:2023)) %>%
   filter(BLOOD, DAYS_SINCE_PRV_CULTURE > 14L | FIRST_RECORD) %>% # 44,421 isolates
   group_by(PERSON_ID, ORDER_DAY)                                 # 38,027 infections



df <- astDF %>% filter(BLOOD)

df %>%
   summarise(n = n(),
             pcnt3 = 1 - sum(DAYS_SINCE_PRV_CULTURE > 3L | FIRST_RECORD) / n(),
             pcnt7 = 1 - sum(DAYS_SINCE_PRV_CULTURE > 7L | FIRST_RECORD) / n(),
             pcnt14 = 1 - sum(DAYS_SINCE_PRV_CULTURE > 14L | FIRST_RECORD) / n(),
             pcnt30 = 1 - sum(DAYS_SINCE_PRV_CULTURE > 30L | FIRST_RECORD) / n(),
             .by = BUG) %>%
   arrange(desc(n)) %>%
   filter(n > 200L) %>%
   arrange(desc(pcnt3))

# keep only blood cultures
astDF <- astDF %>% filter(BLOOD)


# what bugs?
df <- astDF %>% filter(DAYS_SINCE_PRV_CULTURE > 14L | FIRST_RECORD)
df %>% count(BUG, sort=TRUE)

df %>% filter(BUG == 'Staphylococcus aureus') %>% count(OXACILLIN, sort=TRUE)

df %>% filter(BUG == 'Staphylococcus aureus') %>% filter(is.na(OXACILLIN))
# 1000005613, 2017-02-11, its there!


source('~/Desktop/EHR/EHR work/config_file.R')
x <- tbl(conn, in_schema('AMB_ETL', 'LAB_MICRO_SENS_ALL_VW')) %>%
   filter(PERSON_ID == '1000005613') %>%
   collect()
x %>%
   select(RESULT_DATE, ORGANISM_NAME, ANITBIOTIC_NAME, SUSCEPTIBILITY_NAME, SENSITIVITY_VALUE, SENSITIVITY_COMMENT) %>% 
   distinct() %>%
   filter(grepl('AUREUS', ORGANISM_NAME), ANITBIOTIC_NAME == 'OXACILLIN') %>%
   mutate(RESULT_DATE = as.Date(RESULT_DATE, format='%m/%d/%Y')) %>%
   arrange(RESULT_DATE) %>% filter(substr(RESULT_DATE,1,4) == '2017')

load(file = paste0(data_path_name, 'ALL_clean_ASTs.Rdata'))

#

save(astDF, file = '~/Desktop/EHR/EHR work/RdataFiles/ASTs_blood_2015_2023.Rdata')









