start <- Sys.time()
library(dplyr)
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/data_path_name.Rdata')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_clean_ASTs.Rdata') # 1,734,071
print(object.size(astDF), units='Mb') # 1,803.7 Mb

# Remove any Coagulase-negative Staph species if they appear in blood cultures
astDF <- astDF %>% filter(!(BLOOD & grepl('Staph', BUG) & !grepl('Staphylococcus aureus|Staphylococcus lugdunensis', BUG))) # 1,698,329
astDF <- astDF %>% filter(!grepl('Cryptococcus|Aspergillus|Candida', BUG)) # 1,691,137


# only take 2015 - 2023
astDF <- astDF %>% filter(substr(ORDER_DATE, 1, 4) %in% as.character(2015:2023)) # 1,165,383
length(unique(astDF$PERSON_ID)) # 440,378
gc()


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
sum(astDF$MULT_ISO) / nrow(astDF)     # 28.9%
sum(astDF$FIRST_RECORD) / nrow(astDF) # 37.8%

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
rm(shift, w)
# patient's first records for multi-isolates are also screwed up, mark them as TRUE
astDF$FIRST_RECORD[is.na(astDF$DAYS_SINCE_PRV_CULTURE)] <- TRUE


# Indicate how many days since previous blood culture
astDF$DAYS_SINCE_LAST_BLOOD <- NA_integer_
astDF <- astDF %>% relocate(DAYS_SINCE_LAST_BLOOD, .after=BLOOD)
ids <- unique(unname(astDF %>% filter(sum(BLOOD) > 1L, .by=PERSON_ID) %>% select(PERSON_ID) %>% unlist())) # 12,200 patients
w <- which(astDF$PERSON_ID %in% ids & astDF$BLOOD) # 36,724
astDF$DAYS_SINCE_LAST_BLOOD[w] <- ifelse(test = astDF$PERSON_ID[w] == lag(astDF$PERSON_ID[w]),
                                         yes = astDF$ORDER_DAY[w] - lag(astDF$ORDER_DAY[w]),
                                         no = NA_integer_)
shift <- 1
w <- w[which(astDF$PERSON_ID[w] == lead(astDF$PERSON_ID[w]) & astDF$ORDER_DAY[w] == lead(astDF$ORDER_DAY[w]))]
w <- w[astDF$BLOOD[w + shift]]
while(length(w) > 0) {
   cat(shift, length(w), '\n')
   astDF$DAYS_SINCE_LAST_BLOOD[w + shift] <- astDF$DAYS_SINCE_LAST_BLOOD[w]
   shift <- shift + 1
   w <- w[which(astDF$PERSON_ID[w + shift] == astDF$PERSON_ID[w] & astDF$ORDER_DAY[w + shift] == astDF$ORDER_DAY[w])]
   w <- w[astDF$BLOOD[w + shift]]
}
rm(shift, w)




# Indicate whether the current culture has subsequent culture with results BEFORE this one
astDF$DAYS_TO_NEXT_RESULT <- as.numeric(lubridate::as.duration(lead(astDF$RESULT_DATE) - astDF$RESULT_DATE)) / 86400
astDF$DAYS_TO_NEXT_RESULT[lead(astDF$PERSON_ID) != astDF$PERSON_ID] <- NA
astDF <- astDF %>% relocate(DAYS_TO_NEXT_RESULT, .before=BUG)

astDF %>% count(DAYS_TO_NEXT_RESULT < 0)         # ~23K (~2%)
astDF %>% count(BLOOD & DAYS_TO_NEXT_RESULT < 0) # ~11K (~1%)

astDF %>% count(DAYS_TO_NEXT_RESULT == 0, MULT_ISO) # are these behaving correctly?
#   `DAYS_TO_NEXT_RESULT == 0` MULT_ISO      n
# 1 FALSE                      FALSE     ~443K - single isolate with some future culture
# 5 NA                         FALSE     ~385K - single isolate with no future culture
# 3 TRUE                       FALSE        90 - single isolate with identical result_date-time as next culture

# 4 TRUE                       TRUE      ~148K - multiple isolates on same order_day - not last in "series"
# 2 FALSE                      TRUE      ~133K - multiple isolates on same order_day - last in "series" - some future culture
# 6 NA                         TRUE       ~55K - multiple isolates on same order-day - last in "series" - no future culture


# Indicate if positive culture within X days. DAYS_TO_NEXT_CULTURE?
astDF$DAYS_TO_NEXT_CULTURE <- lead(astDF$DAYS_SINCE_PRV_CULTURE)
astDF <- astDF %>% relocate(DAYS_TO_NEXT_CULTURE, .after=DAYS_SINCE_PRV_CULTURE)

# # Indicate if multiple blood isolates
# astDF <- astDF %>%
#    group_by(PERSON_ID, ORDER_DAY, BLOOD) %>%
#    mutate(MULT_BLOOD_ISO = BLOOD & n() > 1L) %>%
#    ungroup() %>%
#    relocate(MULT_BLOOD_ISO, .after=MULT_ISO)
# 
# astDF %>% count(MULT_ISO, MULT_BLOOD_ISO)
# astDF %>% count(BLOOD)



### ignore isolates, collapse to 1 culture per day
df <- astDF %>%
   select(PERSON_ID, FIRST_RECORD, BLOOD, ORDER_DAY, DAYS_SINCE_PRV_CULTURE, MULT_ISO, DAYS_SINCE_LAST_BLOOD) %>% 
   distinct()
df$LAST_CULT_BLOOD <- lag(df$BLOOD)
df$LAST_CULT_BLOOD[df$FIRST_RECORD] <- NA

sum(is.na(df$DAYS_SINCE_PRV_CULTURE[df$BLOOD])) / sum(df$BLOOD) # ~42% of blood culture days are a patients first culture day EVER

med <- median(df$DAYS_SINCE_PRV_CULTURE[df$BLOOD], na.rm=T); med # 40 days is the median number of days since previous culture
sum(df$DAYS_SINCE_PRV_CULTURE[df$BLOOD] <= 30, na.rm=T) / sum(df$BLOOD) # ~27% of blood culture days are < 30 days removed from a previous culture

med <- median(df$DAYS_SINCE_LAST_BLOOD[df$BLOOD], na.rm=T); med # 18 days is the median number of days since previous BLOOD culture
sum(df$DAYS_SINCE_LAST_BLOOD[df$BLOOD] <= 30, na.rm=T) / sum(df$BLOOD)  # ~17% of blood culture days are < 30 days removed from a previous BLOOD culture

# Previously, I removed cultures that were within 6 months of ANY previous culture
# This was my way of ensuring that I take the "first positive (index) blood culture per encounter) - Kadri et al (2021)
# == ensuring independence between infections
# This turns out to be ~43% of blood cultures removed, probably for the "sickest" patients (who have more adjacent blood cultures)
# The AMP team suggest I don't do this, but didn't suggest an alternative
# Perhaps I will only take the first culture in a 30-day window? (this could represent "separate" infections)
# That would still be removing 27% of blood cultures

# Blood cultures tend to be clustered more than other cultures
for (days in c(365, 180, 90, 60, 30, 14, 7, 3)) {
   x <- round(sum(df$DAYS_SINCE_PRV_CULTURE[df$BLOOD] < days, na.rm=T) / sum(df$BLOOD) * 100)
   cat(x, '% of blood cultures are <', days, ' after any previous culture\n', sep='')
   x <- round(sum(df$DAYS_SINCE_LAST_BLOOD[df$BLOOD] < days, na.rm=T) / sum(df$BLOOD) * 100)
   cat(x, '% of blood cultures are <', days, ' after previous BLOOD culture\n\n', sep='')
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
t2 <- df %>% 
   filter(BLOOD) %>% 
   count(DAYS_SINCE_LAST_BLOOD, LAST_CULT_BLOOD) %>% 
   filter(!is.na(LAST_CULT_BLOOD)) %>% 
   tidyr::pivot_wider(values_from=n, names_from=LAST_CULT_BLOOD) %>%
   rename(lastBlood = `TRUE`, lastNonBlood = `FALSE`) %>%
   mutate(total = lastBlood + lastNonBlood) %>%
   select(-lastNonBlood) %>%
   mutate(pcnt = lastBlood / total) %>%
   filter(DAYS_SINCE_LAST_BLOOD <= 180L)

plot(x = t$DAYS_SINCE_PRV_CULTURE, y = t$pcnt, ylim=c(0, 1), ylab='', pch=16, cex=0.5,
     xlab = 'Days since previous culture', main='Proportion previous cultures = blood')
lines(stats::filter(t$pcnt, rep(1/6, 6), sides=2), lwd=2)
points(x=t2$DAYS_SINCE_LAST_BLOOD, y=t2$pcnt, col='red', pch=16, cex=0.5)
lines(stats::filter(t2$pcnt, rep(1/6, 6), sides=2), lwd=2, col='red')


df %>% count(BLOOD) # ~60K

1 - sum(df$FIRST_RECORD | df$DAYS_SINCE_PRV_CULTURE > 14L) / nrow(df)                   # ~10%
1 - sum((df$FIRST_RECORD | df$DAYS_SINCE_PRV_CULTURE > 14L) & df$BLOOD) / sum(df$BLOOD) # ~21%

1 - sum(df$FIRST_RECORD | df$DAYS_SINCE_PRV_CULTURE > 7L) / nrow(df)                   # ~7%
1 - sum((df$FIRST_RECORD | df$DAYS_SINCE_PRV_CULTURE > 7L) & df$BLOOD) / sum(df$BLOOD) # ~19%

# THIS ALL MEANS THAT FOR BLOODSTREAM
astDF %>%
   filter(substr(ORDER_DAY,1,4) %in% as.character(2017:2023)) %>%
   filter(BLOOD) %>%              # ~52K isolates
   group_by(PERSON_ID, ORDER_DAY) # ~48K infections

# ONLY KEEP BLOOD - 2017-2023 - at least 30 days removed from previous blood culture
astDF <- astDF %>%                                                      # 1.1-1.2 mil
   filter(substr(ORDER_DAY,1,4) %in% as.character(2017:2023)) %>%       #   ~900K
   filter(BLOOD) %>%                                                    #    ~53K
   filter(is.na(DAYS_SINCE_LAST_BLOOD) | DAYS_SINCE_LAST_BLOOD >= 30L)  #    43,990







astDF %>% filter(BUG == 'Staphylococcus aureus') %>% count(OXACILLIN, sort=TRUE)
astDF %>% filter(BUG == 'Escherichia coli') %>% count(ESBL)
t <- astDF %>% filter(BUG == 'Enterococcus faecium') %>% select(VANCOMYCIN, AMPICILLIN) %>% table()
t; fisher.test(t)


save(astDF, file = '~/Desktop/EHR/EHR work/RdataFiles/ASTs_blood_2017_2023.Rdata')
print(Sys.time() - start) # 38 seconds








