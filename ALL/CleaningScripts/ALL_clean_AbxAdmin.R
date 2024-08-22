library(dplyr)
library(tidyr)

load(file = '~/Desktop/EHR/EHR work/RdataFiles/sens_med_admin_and_iv_vw.Rdata')
length(unique(abxDF$PERSON_ID)) # 150,225


# make a few changes to the medication column to coerce particular antibiotic names
abxDF <- abxDF %>%
   mutate(MEDICATION = gsub('-', '/', MEDICATION)) %>%
   mutate(MEDICATION = case_when(
      MEDICATION == 'SULFAMETHOXAZOLE/TRIMETHOPRIM' ~ 'TRIMETHOPRIM/SULFAMETHOXAZOLE',
      MEDICATION == 'DALFOPRISTIN/QUINUPRISTIN' ~ 'SYNERCID',
      .default = MEDICATION
   ))

# get counts of each medication
meds <- abxDF %>% count(MEDICATION, sort=TRUE)



# AST antibiotics
load(file = '~/Desktop/EHR/EHR work/data/ast_antibiotics.Rdata')
abx_ast <- abx; rm(abx)
abx_ast <- abx_ast[abx_ast != 'PAS']
abx_ast <- c(abx_ast, 'PASER')

# Antibiotic names, abbreviations, and classes
abx_abbr <- tibble(read.table(file = '~/Desktop/EHR/EHR work/data/ABX_ABBR.txt', header = TRUE)) %>%
   select(Antibiotic_Name) %>%
   mutate(Antibiotic_Name = gsub('-', '/', Antibiotic_Name),
          Antibiotic_Name = gsub('_', ' ', Antibiotic_Name),
          Antibiotic_Name = toupper(Antibiotic_Name)) %>%
   unlist()
names(abx_abbr) <- NULL

# Antibiotic
abx_gen <- read.table(file = '~/Desktop/EHR/EHR work/data/AntibioticNames.txt', sep='\t', header = TRUE)
abx_gen <- abx_gen$Generic

# COMBINE ALL
abx_all <- sort(unique(toupper(c(abx_ast, abx_abbr, abx_gen))))
length(abx_all) # 216
rm(abx_gen, abx_ast, abx_abbr)

# look through actaul medication names for hits from this list
s <- sapply(abx_all, function(a) grep(a, meds$MEDICATION, value=TRUE, ignore.case=TRUE))
s[lengths(s) > 0L]
rm(s)


# get rows of full data.frame
med_rows <- sort(unique(unlist(sapply(abx_all, function(a) grep(a, meds$MEDICATION, ignore.case=TRUE)))))
meds <- meds[med_rows, ]
meds$ABX <- NA_character_
for (i in seq_len(nrow(meds))) {
   w <- which(meds$MEDICATION[i] == abx_all)
   if (length(w) == 1) {
      meds$ABX[i] <- meds$MEDICATION[i]
   } else {
      meds$ABX[i] <- paste(abx_all[which(stringi::stri_detect_regex(str = meds$MEDICATION[i], pattern = abx_all))],
                           collapse = ', ')
   }
}
rm(med_rows, i, w, abx_all)


# clean up the >1 match instances
meds %>% filter(grepl(',', ABX)) # 20
meds <- meds %>%
   mutate(ABX = case_when(
      ABX == 'PENICILLIN G' ~ 'BENZYLPENICILLIN',
      ABX == 'PENICILLIN, PENICILLIN G' ~ 'BENZYLPENICILLIN',
      ABX == 'CIPROFLOXACIN, OFLOXACIN' ~ 'CIPROFLOXACIN',
      ABX == 'PENICILLIN, PENICILLIN V' ~ 'PENICILLIN V',
      .default = ABX
   ))
meds %>% filter(grepl(',', ABX)) # 12 real ones


# create antibiotic column from the named vector of antibiotics
meds <- setNames(meds$ABX, nm = meds$MEDICATION)
start <- Sys.time()
abxDF <- abxDF %>% mutate(ABX = meds[MEDICATION]) # this is large, so takes a minute
print(Sys.time() - start) # 3-4.4 seconds for 24 million rows
rm(start, meds)

# how many of the medications were antibiotics
abxDF %>% count(is.na(ABX)) # 1.3 million (out of 23.7 million) are antibiotics
medsDF <- abxDF %>% filter(is.na(ABX))
abxDF <- abxDF %>% filter(!is.na(ABX))


# Check for remaining brand name antibiotics
abx_b <- read.table(file = '~/Desktop/EHR/EHR work/data/AntibioticNames.txt', sep='\t', header = TRUE)
abx_b <- abx_b$Brand
abx_b <- unlist(strsplit(abx_b, split=', '))
abx_b <- abx_b[abx_b != ' ']
abx_b <- toupper(abx_b)
meds2 <- unique(medsDF$MEDICATION)
s <- sapply(abx_b, function(a) grep(a, meds2))
any(lengths(s) > 0L) # FALSE - GOOD!!
rm(medsDF, abx_b, s, meds2)
###########




abxDF <- abxDF %>%
   select(-ORDER_ID, -ENCOUNTER_ID, -MEDICATION_CODE, -MEDICATION)


abxDF %>%
   count(MEDICATION, sort=TRUE) %>%
   filter(grepl('TOP(ICAL)?$', MEDICATION)) %>%
   print(n=50)
c('TOPICAL', 'OPHTH', 'OTIC', 'TOP', 'OPHTHALMIC')


abxDF %>% count(ABX, sort=TRUE)
length(unique(abxDF$PERSON_ID)) # 85,713 (only ~1/2 of individuals with AST results in 2017 received antibiotics)


abxDF %>% count(is.na(ADMIN_END_DATE)) # ~25,748

abxDF %>%
   mutate(ADMIN_DURATION = as.numeric(lubridate::as.duration(ADMIN_END_DATE - ADMIN_START_DATE)) / 86400) %>%
   summarise(median = median(ADMIN_DURATION, na.rm=T),
             .by = ABX) %>%
   arrange(desc(median)) # SYNERCID is the only antibiotic for which the duration > 0 (1 hour in this case)

abxDF <- abxDF %>%
   select(-ADMIN_ROUTE, -ADMIN_DOSAGE, -DOSAGE_UNIT)


# how long are therapy courses?
abxDF <- abxDF %>% arrange(PERSON_ID, ABX)

names(abxDF$ABX) <- NULL


save(abxDF, file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_2017_AbxAdmin.Rdata')






