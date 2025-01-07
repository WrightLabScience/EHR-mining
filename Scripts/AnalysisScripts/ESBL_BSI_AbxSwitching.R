library(dplyr)
source('~/Desktop/WrightLab/UsefulRfuncs/newBarplotFxn.R')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ASTs_AbxAdmin_blood_2017_2023_imputed_encs_flagged_LONG.Rdata')

# Get bacteremia ESBL ASTs and antibiotics
bugs <- unique(empDF$BUG[empDF$ESBL == 1L])
cohort <- empDF %>%
   #filter(ESBL == 1L) %>%
   filter(BUG %in% bugs) %>%
   select(PERSON_ID, BUG, ORDER_DAY, RESULT_DAY, ABX, START_DAY, FLAG, ESBL) %>%
   distinct()
cohort <- cohort %>%
   group_by_all() %>%
   ungroup(RESULT_DAY) %>%
   slice_min(RESULT_DAY) %>%
   ungroup()
cohort$PATHOGEN <- paste0(ifelse(cohort$ESBL == 1L, 'ESBL - ', ''), gsub('^([A-Z])[a-z]+ ([a-z]+)$', '\\1. \\2', cohort$BUG))

asts <- cohort %>% 
   select(PERSON_ID, ORDER_DAY, PATHOGEN) %>% distinct() %>% 
   summarise(BUG = paste(sort(PATHOGEN), collapse=', '),
             .by=c(PERSON_ID, ORDER_DAY))
asts %>% count(grepl('ESBL', BUG)) # 2,327 ESBL + 17,647 non = 19,974 cases total (11.7% ESBL)
abx <- cohort %>%
   select(PERSON_ID, ORDER_DAY, ABX, START_DAY, FLAG) %>%
   distinct() %>%
   mutate(PROX = as.integer(START_DAY - ORDER_DAY)) %>%
   filter(PROX %in% -10:30)
df <- left_join(x = asts, 
                y = abx,
                by = join_by(PERSON_ID, ORDER_DAY)) %>%
   mutate(ESBL = grepl('ESBL', BUG))
rm(asts, abx)


# how many have 0 antibiotics? ~10% which bugs? ~a good mix - ESBL missing 10.1%, non-ESBL missing 8.24%
# remove them...
df %>% mutate(hasAbx = !is.na(ABX)) %>% select(PERSON_ID, ORDER_DAY, BUG, ESBL, hasAbx) %>% distinct() %>% count(hasAbx) # 236 / 2,327 (~10%)
df %>% mutate(hasAbx = !is.na(ABX)) %>% select(PERSON_ID, ORDER_DAY, BUG, ESBL, hasAbx) %>% distinct() %>% summarise(n=n(), missingAbx = sum(!hasAbx) / n * 100, .by=ESBL)
df %>% mutate(hasAbx = !is.na(ABX)) %>% select(PERSON_ID, ORDER_DAY, BUG, ESBL, hasAbx) %>% distinct() %>% summarise(n=n(), missingAbx = sum(!hasAbx) / n * 100, .by=BUG) %>% filter(n > 200) %>% arrange(desc(missingAbx))
df <- df %>% filter(!is.na(ABX))
df %>% group_by(PERSON_ID, ORDER_DAY) # 18,284 cases

# how many antibiotics prescribed early? before culture? 18-22%
# flag them...
df %>% mutate(X = as.integer(START_DAY - ORDER_DAY)) %>% pull(X) %>% plotBarplot()
df %>% group_by(PERSON_ID, ORDER_DAY) %>% filter(any(PROX < -2)) # 1,748 / 18,284 (18%) are on antibiotics before blood culture
df <- df %>% group_by(PERSON_ID, ORDER_DAY) %>% mutate(earlyAbx = any(PROX < -2)) %>% ungroup()
df %>% summarise(sum(earlyAbx) / n(), .by=ESBL) # ESBLs - 35.9%! vs. 20.4%!


# do patients who receive antibiotics early have higher rates of discordant empiric therapy?
df %>% 
   filter(PROX == 0) %>% # count(FLAG) # 56.7% discordance
   summarise(FLAG = paste(sort(unique(FLAG)), collapse=', '), 
             .by=c(PERSON_ID, ORDER_DAY, BUG, ESBL, START_DAY, PROX, earlyAbx)) %>%
   mutate(FLAG = case_when(
      grepl('CONCORDANT', FLAG) ~ 'CONCORDANT',
      grepl('DISCORDANT', FLAG) ~ 'DISCORDANT',
      .default = FLAG
   )) %>%
   summarise(n = n(),
             sum(FLAG == 'CONCORDANT') / n,
             sum(FLAG == 'DISCORDANT') / n,
             sum(FLAG == 'NOTTESTED') / n,
             .by=c(earlyAbx, ESBL)) %>%
   arrange(ESBL, earlyAbx)
# early antibiotic use highly associated with ESBLs
# rates of discordance are much higher among ESBLs
# early antibiotic use is not associated with greater rates of empiric discordance among ESBL
# early antibiotic use IS associated with greater rates of empiric discordance among non-ESBL

tn <- df %>%
   filter(PROX == 0, ESBL == 0L, FLAG == 'DISCORDANT') %>%
   filter(!earlyAbx) %>% count(ABX, sort=TRUE) %>%
   rename(noEarlyUse = n)
ty <- df %>%
   filter(PROX == 0, ESBL == 0L, FLAG == 'DISCORDANT') %>%
   filter(earlyAbx) %>% count(ABX, sort=TRUE) %>%
   rename(yesEarlyUse = n)

t <- full_join(ty, tn, by=join_by(ABX)) %>% 
   mutate(across(c(yesEarlyUse, noEarlyUse), ~ . / sum(., na.rm=T) * 100)) %>%
   mutate(across(c(yesEarlyUse, noEarlyUse), ~ ifelse(is.na(.), 0, .)))
t$diff <- t$yesEarlyUse - t$noEarlyUse
t <- t %>% arrange(diff)

load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/antibiotic_names/abbr_named.Rdata')
par(mar=c(3.5, 4.5, 1, 1), mgp=c(2.5, 0.5, 0), tck=-0.015)
b <- barplot(t$diff, xlim=c(-10, 10), xaxt='n', xlab='% difference', main='early - not early', horiz=TRUE)
axis(side=1, las=1)
w <- which(abs(t$diff) > 1)
text(y=b[w], x=-10.5, labels=abbr[t$ABX][w], adj=1, xpd=NA)
text(x=-5, y=mean(b), labels='Used less in patients\ntreated prior to culture')
text(x=5, y=mean(b), labels='Used more in patients\ntreated prior to culture')



# how many different (2-5) and which antibiotics are given? (VAN!!, MER, TZP, ERT, FEP, CRO)
df %>% select(PERSON_ID, ORDER_DAY, BUG, ABX) %>% distinct() %>% group_by(PERSON_ID, ORDER_DAY, BUG) %>% tally() %>% pull(n) %>% plotBarplot()
df %>% select(PERSON_ID, ORDER_DAY, BUG, ABX) %>% distinct() %>% count(ABX, sort=T)
df %>% filter(lubridate::year(ORDER_DAY) %in% 2021:2023, PROX == 1) %>% select(PERSON_ID, ORDER_DAY, BUG, ABX) %>% distinct() %>% count(ABX, sort=T)
df %>% filter(lubridate::year(ORDER_DAY) < 2021, PROX == 1) %>% select(PERSON_ID, ORDER_DAY, BUG, ABX) %>% distinct() %>% count(ABX, sort=T)



# join with encounters
ids <- unique(df$PERSON_ID)
source('~/Desktop/EHR/EHR work/config_file.R')
chunks <- mapply(FUN = ':',
                 seq(1, floor(length(ids) / 1000) * 1000 + 1, 1000),
                 c(seq(1000, floor(length(ids) / 1000) * 1000, 1000), length(ids)))
encs <- tibble()
for (chunk in chunks) {
   print(chunk[1])
   encs <- rbind(encs,
                 tbl(conn, in_schema('AMB_ETL', 'SENS_ENCOUNTER_VW')) %>%
                    filter(PERSON_ID %in% local(ids[chunk])) %>%
                    collect())
}
# fix date times
encs <- encs %>%
   mutate(across(c(ADMIT_DATE, DISCHARGE_DATE), ~ strptime(., format='%m/%d/%Y %T'))) %>%
   mutate(ADMIT_DAY = lubridate::as_date(ADMIT_DATE),
          DISCHARGE_DAY = lubridate::as_date(DISCHARGE_DATE)) %>%
   arrange(PERSON_ID, ADMIT_DATE, DISCHARGE_DATE)

# join with df
x <- df %>%
   select(PERSON_ID, ORDER_DAY, BUG, ESBL, earlyAbx) %>% distinct() %>%
   mutate(JOIN_START = ORDER_DAY - 14,
          JOIN_END = ORDER_DAY + 7) %>%
   left_join(y = encs,
             by = join_by(
                PERSON_ID,
                JOIN_START <= ADMIT_DAY,
                JOIN_END >= ADMIT_DAY
             )) %>%
   mutate(PROX = as.integer(ORDER_DAY - ADMIT_DAY)) %>% #pull(X) %>% plotBarplot()
   select(-JOIN_START, -JOIN_END) %>%
   relocate(ADMIT_DAY, DISCHARGE_DAY, .after=ORDER_DAY) %>%
   group_by(PERSON_ID, ORDER_DAY) %>%
   mutate(BEFORE = any(ORDER_DAY > DISCHARGE_DAY)) %>%
   ungroup() %>%
   mutate(WITHIN = ORDER_DAY >= ADMIT_DAY & ORDER_DAY <= DISCHARGE_DAY) %>%
   filter(WITHIN) %>%
   mutate(HOSP_ACQ = PROX >= 2)

fisher.test(x %>% select(ESBL, BEFORE) %>% table())
fisher.test(x %>% select(ESBL, earlyAbx) %>% table())
fisher.test(x %>% select(ESBL, HOSP_ACQ) %>% table())

fisher.test(x %>% select(earlyAbx, HOSP_ACQ) %>% table()) # 29.4
fisher.test(x %>% select(BEFORE, HOSP_ACQ) %>% table())   # 1.3 - you wouldn't expect recent visits to predispose to hospital-acquired infections
fisher.test(x %>% select(BEFORE, earlyAbx) %>% table())   # 10.7

x %>% summarise(n=n(), sum(ESBL) / n(), .by=earlyAbx) # 2X ESBL % for earlyAbx
x %>% summarise(n=n(), sum(ESBL) / n(), .by=HOSP_ACQ) # 1.5X ESBL for hospital-acquired
x %>% summarise(n=n(), sum(ESBL) / n(), .by=BEFORE) # 1.8X ESBL for hospital-acquired
x %>% summarise(n=n(), sum(ESBL) / n(), .by=c(HOSP_ACQ, earlyAbx, BEFORE)) %>% arrange(HOSP_ACQ, earlyAbx, BEFORE)
# among community-acquired, earlyAbx,     2.5X more ESBL if recent visit
# among community-acquired, no earlyAbx,  1.4X more ESBL if recent visit

summary(glm(ESBL ~ earlyAbx, family=binomial(link='logit'), data=x))
summary(glm(ESBL ~ HOSP_ACQ, family=binomial(link='logit'), data=x))
summary(glm(ESBL ~ BEFORE, family=binomial(link='logit'), data=x))
summary(glm(ESBL ~ BEFORE + earlyAbx + HOSP_ACQ, family=binomial(link='logit'), data=x))

# among hospital-acquired, earlyAbx, 2.5X more ESBL if recent visit

# among community-acquired, 1.8X ESBL for earlyAbx
# among hospital-acquired, 2X ESBL for earlyAbx

# add flag for if this culture+visit were preceded recently by another
df %>%
   select(PERSON_ID, ORDER_DAY, BUG, ESBL, earlyAbx) %>% distinct() %>%
   mutate(JOIN_START = ORDER_DAY - 14,
          JOIN_END = ORDER_DAY + 7) %>%
   left_join(y = encs,
             by = join_by(
                PERSON_ID,
                JOIN_START <= ADMIT_DAY,
                JOIN_END >= ADMIT_DAY
             )) %>%
   mutate(PROX = as.integer(ORDER_DAY - ADMIT_DAY)) %>% #pull(X) %>% plotBarplot()
   select(-JOIN_START, -JOIN_END) %>%
   relocate(ADMIT_DAY, DISCHARGE_DAY, .after=ORDER_DAY) %>%
   mutate(WITHIN = ORDER_DAY >= ADMIT_DAY & ORDER_DAY <= DISCHARGE_DAY) %>%
   group_by(PERSON_ID, ORDER_DAY) %>%
   mutate(BEFORE = any(ORDER_DAY > DISCHARGE_DAY)) %>%
   ungroup() %>%
   filter(WITHIN) %>%
   mutate(HOSP_ACQ = PROX >= 2) %>%
   summarise(n=n(), 
             esbl_rate = sum(ESBL) / n(), 
             .by=BEFORE)




