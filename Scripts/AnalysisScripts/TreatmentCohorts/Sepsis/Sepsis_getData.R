source('~/Desktop/EHR/EHR work/config_file.R')
source('~/Desktop/WrightLab/UsefulRfuncs/newBarplotFxn.R')

##### GET SEPSIS DIAGNOSIS DATA #####
dx <- tibble(dbGetQuery(conn, "select * from AMB_ETL.LAB_SENS_DX_VW where CODE_DESCRIPTION like '%sepsis%' or CODE_DESCRIPTION like '%Sepsis%'"))
dx <- dx %>%
   filter(substr(DX_FROM_DATE,1,4) %in% as.character(2017:2023)) %>% 
   mutate(DX_FROM_DATE = lubridate::as_date(DX_FROM_DATE),
          DX_TO_DATE   = lubridate::as_date(DX_TO_DATE)) %>%
   distinct() %>% 
   arrange(PERSON_ID, DX_FROM_DATE, CODE_DESCRIPTION)

dx <- dx %>% 
   mutate(DX_TO_DATE = case_when(is.na(DX_TO_DATE) ~ DX_FROM_DATE, .default=DX_TO_DATE)) %>%
   summarise(DX_TO_DATE = max(DX_TO_DATE),
             across(.cols = c(DX_CODE, CODE_TYPE, CODE_DESCRIPTION, PRIMARY_DX_IND),
                    .fns = ~ paste(sort(unique(.)), collapse='; ')),
             .by=c(PERSON_ID, DX_FROM_DATE)) %>%
   mutate(DURATION = as.integer(DX_TO_DATE - DX_FROM_DATE)) %>%
   group_by(PERSON_ID) %>%
   mutate(DAYS_SINCE_PRV = as.integer(DX_FROM_DATE - lag(DX_FROM_DATE))) %>%
   ungroup()

dx <- dx %>% filter(DURATION >= 0L)

dx %>% pull(DURATION) %>% plotBarplot(xmax=60)

dx <- dx %>% filter(is.na(DAYS_SINCE_PRV) | DAYS_SINCE_PRV > 90L)
dx %>% group_by(PERSON_ID, DX_FROM_DATE)

save(dx, file='~/Desktop/EHR/EHR work/RdataFiles/causal_prep/Sepsis/Sepsis_Dx.Rdata')
##### END #####


##### GET ENCOUNTERS DATA #####
if (FALSE) {
   load(file='~/Desktop/EHR/EHR work/RdataFiles/causal_prep/Sepsis/Sepsis_Dx.Rdata')
   ids <- unique(dx$PERSON_ID)
   chunks <- mapply(FUN = ':', 
                    seq(1, floor(length(ids) / 1000) * 1000 + 1, 1000), 
                    c(seq(1000, floor(length(ids) / 1000) * 1000 + 1, 1000), length(ids)))
   encs <- tibble()
   for (chunk in chunks) {
      print(chunk[1])
      encs <- rbind(encs,
                    tbl(conn, in_schema('AMB_ETL', 'SENS_ENCOUNTER_VW')) %>%
                       filter(PERSON_ID %in% local(ids[chunk])) %>%
                       collect())
   }
   rm(chunk, chunks, ids)
   encs <- encs %>%
      mutate(across(.cols = c(ADMIT_DATE, DISCHARGE_DATE),
                    .fns = ~ strptime(x=., format='%m/%d/%Y %T')))
   source('~/Desktop/EHR/EHR-mining/Scripts/CleaningScripts/CombineOverlappingEncountersFxn.R')
   encs_og <- encs
   w <- c(which(encs_og$PERSON_ID == '1000605689' & encs_og$ADMIT_DATE == '2021-09-09 01:34:00'),
          which(encs_og$PERSON_ID == '1004454851' & encs_og$ADMIT_DATE == '2022-03-24 12:31:00'))
   encs_og <- encs_og[-w, ]; rm(w)
   encs <- combineOverlappingEncounters(encs_og) # 5-10 minutes
   rm(combineOverlappingEncounters)
   
   encs_og <- encs_og %>%
      mutate(ADMIT_DAY = lubridate::as_date(ADMIT_DATE),
             DISCHARGE_DAY = lubridate::as_date(DISCHARGE_DATE))
   
   save(encs, file='~/Desktop/EHR/EHR work/RdataFiles/causal_prep/Sepsis/EncountersCleaned.Rdata')
   save(encs_og, file='~/Desktop/EHR/EHR work/RdataFiles/causal_prep/Sepsis/EncountersRaw.Rdata')
}

library(dplyr)
source('~/Desktop/WrightLab/UsefulRfuncs/newBarplotFxn.R')
load(file='~/Desktop/EHR/EHR work/RdataFiles/causal_prep/Sepsis/Sepsis_Dx.Rdata')
load(file='~/Desktop/EHR/EHR work/RdataFiles/causal_prep/Sepsis/EncountersCleaned.Rdata')
load(file='~/Desktop/EHR/EHR work/RdataFiles/causal_prep/Sepsis/EncountersRaw.Rdata')
encs <- encs %>%
   select(-COMBINED_ADJ, -ENCOUNTER_TYPE) %>%
   group_by(PERSON_ID, ADMIT_DATE, DISCHARGE_DATE, ADMIT_DAY, DISCHARGE_DAY) %>%
   reframe(across(.cols = c(FACILITY, ADMIT_SOURCE, ADMIT_TYPE),
                  .fns = ~ paste(unique(sort(unlist(strsplit(., ', ')))), collapse=', '))) %>% 
   mutate(JOIN_START = ADMIT_DAY - 1)


# we have DIAGNOSES: dx
# we have ENCOUNTERS: encs, encs_og
# JOIN THEM
df <- dx %>% 
   rename(SEPSIS_DX_DATE = DX_FROM_DATE, 
          SEPSIS_DESCRIPTION = CODE_DESCRIPTION) %>%
   select(!c(DX_TO_DATE, DX_CODE, PRIMARY_DX_IND, DURATION)) %>%
   #filter(PERSON_ID != '1001922498') %>%
   left_join(
      x = .,
      y = encs,
      by = join_by(
         PERSON_ID,
         between(SEPSIS_DX_DATE, JOIN_START, DISCHARGE_DAY)
      )
   ) %>%
   select(-JOIN_START) %>%
   mutate(SEPDX_m_ADMIT = as.integer(SEPSIS_DX_DATE - ADMIT_DAY))
df %>% group_by(PERSON_ID, SEPSIS_DX_DATE) %>% filter(n() > 1L) # 0 good
df %>% count(is.na(ADMIT_DATE))
df %>% count(is.na(DAYS_SINCE_PRV) | DAYS_SINCE_PRV >= 60L)
df <- df %>% filter(!is.na(ADMIT_DATE))

df$NURSING_HOME <- grepl('Nursing', df$ADMIT_SOURCE) # much higher mortality rate, nearly double 30-day
df$EMERGENCY_DEPT <- grepl('ED|Emergency|Trauma|Urgent', df$ADMIT_TYPE) # 6% higher 30-day mortality rate
df$TRANSFER <- grepl(', ', df$FACILITY)

df <- df %>% mutate(DAYS_BW_SEP_DISCH = as.integer(DISCHARGE_DAY - SEPSIS_DX_DATE))
df %>% pull(DAYS_BW_SEP_DISCH) %>% plotBarplot(xmax=14)
df %>% count(DAYS_BW_SEP_DISCH > 4L)
df <- df %>% filter(DAYS_BW_SEP_DISCH > 4L)

df %>% count(TRANSFER) # 3,636
df <- df %>%
   select(!c(ADMIT_DATE, ADMIT_DAY, DISCHARGE_DATE, DISCHARGE_DAY)) %>%
   rename(FACILITIES = FACILITY) %>%
   left_join(y = encs_og %>% 
                select(PERSON_ID, ADMIT_DAY, DISCHARGE_DAY, FACILITY) %>% 
                mutate(JOIN_START = ADMIT_DAY - 1), 
             by = join_by(
                PERSON_ID, 
                between(SEPSIS_DX_DATE, JOIN_START, DISCHARGE_DAY)
             )) %>%
   distinct() %>%
   group_by(PERSON_ID, SEPSIS_DX_DATE) %>%
   slice_max(DISCHARGE_DAY) %>%
   filter(n() == 1L) %>%
   ungroup() %>%
   select(-ADMIT_DAY, -DISCHARGE_DAY, -JOIN_START)
any(is.na(df$FACILITY))

# which facilities are represented and do they use the primary sepsis diagnosis code at vastly different rates?
df %>% 
   summarise(n = n(), 
             x = sum(grepl('Sepsis, unspecified organism', SEPSIS_DESCRIPTION)) / n, 
             .by=FACILITY) %>% 
   arrange(desc(n)) %>% 
   filter(n > 500L)

# how many patients have sepsis diagnosis within 1 day of admission?
sum(df$SEPDX_m_ADMIT == 0, na.rm=T) / nrow(df)      # 95.3%
sum(df$SEPDX_m_ADMIT %in% -1:1, na.rm=T) / nrow(df) # 97.0%


# Ryan Shields wanted to grab the people who were diagnosed with sepsis within 48 hours of admission
# and find out how many of them were given each antibiotic
# this remove MOSTLY patients without ANY ENCOUNTERS data,
# and also small number that developed sepsis at some point in their visit
df <- df %>% select(!c(ADMIT_SOURCE, ADMIT_TYPE, CODE_TYPE, FACILITIES)) # ~50K
save(df, file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/Sepsis/Sepsis_Dx_Enc.Rdata')
##### END #####



### load dx + enc data ###
library(dplyr)
load(file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/Sepsis/Sepsis_Dx_Enc.Rdata')
source('~/Desktop/WrightLab/UsefulRfuncs/newBarplotFxn.R')
### END ###



##### GET DEMOGRAPHIC/DEATH DATA #####
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_DEMO.Rdata')
df <- df %>%
   left_join(x = .,
             y = dth,
             by = join_by(PERSON_ID)) %>%
   mutate(time = as.integer(DEATH_DATE - SEPSIS_DX_DATE),
          AGE = as.integer(SEPSIS_DX_DATE - DOB) / 365) %>%
   filter(is.na(time) | time > 4L)

df %>% filter(time < 3L) %>% count(time)
df %>% summarise(sum(time < 30, na.rm=T) / n())
df %>% summarise(n=n(), sum(time < 30, na.rm=T) / n, .by=FACILITY) %>% filter(n > 500L)
rm(dth)
##### END #####



##### GET AST DATA #####
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_clean_ASTs.Rdata')
astDF <- astDF %>%
   filter(PERSON_ID %in% unique(df$PERSON_ID),
          lubridate::year(ORDER_DAY) %in% 2017:2023)

w <- which(astDF$BLOOD & grepl('Staph', astDF$BUG) & !grepl('Staphylococcus aureus|Staphylococcus lugdunensis', astDF$BUG)) # 49,167
astDF <- astDF[-w, ]; rm(w)
astDF <- astDF %>% select(PERSON_ID, ORDER_DAY, RESULT_DAY, BLOOD, BUG, OXACILLIN, VANCOMYCIN, ESBL) %>% distinct()
astDF$FUNGUS <- grepl('Cryptococcus|Aspergillus|Candida', astDF$BUG)

df <- df %>%
   mutate(JOIN_START = SEPSIS_DX_DATE - 1,
          JOIN_END = SEPSIS_DX_DATE + 1) %>%
   left_join(y = astDF,
             multiple = 'all',
             by = join_by(
                PERSON_ID,
                between(y$ORDER_DAY, x$JOIN_START, x$JOIN_END)
             )) %>% #mutate(X = as.integer(ORDER_DAY - SEPSIS_DX_DATE)) %>% count(is.na(X))
   select(-JOIN_START, -JOIN_END) %>%
   mutate(ORDER_DAY = case_when(
      is.na(ORDER_DAY) ~ as.Date('1900-01-01', format='%Y-%m-%d'),
      .default = ORDER_DAY
   )) %>%
   group_by(PERSON_ID, SEPSIS_DX_DATE) %>%
   mutate(
      ISOLATE = any(!is.na(BUG)),
      BUGs = list(unique(sort(BUG))),
      BLOOD = any(BLOOD),
      FUNGI = any(FUNGUS),
      BACTERIA = any(!FUNGUS),
      ESBL = case_when(
         all(is.na(BUG)) ~ NA,
         any(!is.na(BUG) & ESBL == 1L) ~ TRUE,
         .default = FALSE
      ),
      MRSA = case_when(
         all(is.na(BUG)) ~ NA,
         any(BUG == 'Staphylococcus aureus' & !is.na(OXACILLIN) & OXACILLIN == 1L) ~ TRUE,
         .default = FALSE
      ),
      VRE = case_when(
         all(is.na(BUG)) ~ NA,
         any(BUG %in% c('Enterococcus faecium', 'Enterococcus faecalis') & !is.na(VANCOMYCIN) & VANCOMYCIN == 1L) ~ TRUE,
         .default = FALSE
      ),
      ORDER_DAY = min(ORDER_DAY),
      RESULT_DAY = min(RESULT_DAY)
   ) %>%
   select(-OXACILLIN, -VANCOMYCIN, -FUNGUS, -BUG) %>% 
   distinct() %>%
   ungroup() %>%
   mutate(ORDER_DAY = case_when(
      ORDER_DAY == as.Date('1900-01-01', format='%Y-%m-%d') ~ NA,
      .default = ORDER_DAY
   ))
df %>% count(MRSA, VRE, ESBL)
df %>% count(FUNGI, BACTERIA)

x <- sum(!is.na(df$ORDER_DAY))
x / nrow(df) # only 55.2% have AST results
sum(df$BLOOD, na.rm=T) / x  # 53.4% had bacteremia
sum(df$FUNGI, na.rm=T) / x  # 1.1% had fungus
sum(df$VRE, na.rm=T) / x    # 1.6% had VRE 
sum(df$MRSA, na.rm=T) / x   # 10.8% had MRSA
sum(df$ESBL, na.rm=T) / x   #  8.4% had an ESBL
df %>% filter(ISOLATE) %>% count(FUNGI, BLOOD) # most instances of fungus are from blood
rm(x)

w <- grepl('Methicillin resistant Staph', df$SEPSIS_DESCRIPTION, ignore.case=TRUE)
x <- df$MRSA
x[is.na(x)] <- -1
table(x, w)
#     w
# x   FALSE    TRUE
# -1  21,934    341  - no AST collected, code says MRSA
# 0   24,371     97  - AST = MSSA, code sometimes says MRSA
# 1    1,599  1,360  - AST = MRSA, code agrees ~1/2 the time
rm(w, x, astDF)
##### END #####



##### JOIN WITH ANTIBIOTIC ADMINISTRATION #####
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_CLEANED_2017_2023_AbxAdmin.Rdata')
abxDF <- abxDF %>% 
   filter(PERSON_ID %in% unique(df$PERSON_ID)) %>% # 9.7 mil --> 2.2 mil
   select(PERSON_ID, START_DAY, ABX) %>%
   distinct()

df <- df %>%
   mutate(JOIN_START = SEPSIS_DX_DATE - 1,
          JOIN_END = SEPSIS_DX_DATE + 7) %>%
   left_join(y = abxDF,
             by = join_by(
                PERSON_ID,
                between(y$START_DAY, x$JOIN_START, x$JOIN_END)
             )) %>%
   select(-JOIN_START, -JOIN_END) %>%
   mutate(ABX_m_DX = as.integer(START_DAY - SEPSIS_DX_DATE)) #%>% pull(ABX_m_DX) %>% plotBarplot()

df <- df %>% filter(!is.na(ABX))

df %>% group_by(PERSON_ID, SEPSIS_DX_DATE) # ~47K cases
df %>% summarise(n3 = sum(unique(ABX_m_DX) %in% 0:3), .by=c(PERSON_ID, SEPSIS_DX_DATE)) %>% count(n3)
df <- df %>% group_by(PERSON_ID, SEPSIS_DX_DATE) %>% filter(sum(unique(ABX_m_DX) %in% 0:3) >= 3L) %>% ungroup()


# most patients receive 3 or more different antibiotics in the days 0 - 3
df %>% 
   filter(ABX_m_DX %in% 0:3) %>% 
   select(PERSON_ID, SEPSIS_DX_DATE, ABX) %>%
   distinct() %>%
   summarise(n = n(),
             .by = c(PERSON_ID, SEPSIS_DX_DATE)) %>% 
   pull(n) %>% 
   plotBarplot(main='Number of unique abx administered within 2 days of diagnosis')


# which antibiotics are given?
# which antibiotics given on a per infection basis by day
abx <- c('PIPERACILLIN/TAZOBACTAM', 'VANCOMYCIN', 'AZITHROMYCIN', 'METRONIDAZOLE', 
         'CEFEPIME', 'CEFTRIAXONE', 'MEROPENEM', 'CEFAZOLIN', 'AMPICILLIN/SULBACTAM')
d <- df %>% filter(ABX_m_DX %in% 0:4) %>% select(PERSON_ID, SEPSIS_DX_DATE, ABX_m_DX) %>% distinct() %>% select(ABX_m_DX) %>% table()
t <- df %>%
   filter(ABX_m_DX %in% 0:4, ABX %in% abx) %>% 
   select(PERSON_ID, SEPSIS_DX_DATE, ABX_m_DX, ABX) %>% 
   distinct() %>%
   select(ABX_m_DX, ABX) %>%
   table()
for (i in 1:5) {
   t[i, ] <- t[i, ] / d[i] * 100
}

load(file='~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/antibiotic_names/abbr_named.Rdata')
set.seed(5L)
col_vec <- paletteer_d("tidyquant::tq_light")
col_vec <- setNames(col_vec[sample(1:length(col_vec))], abx)
plot(NA, xlim=c(0, 4), ylim=c(0, 55), ylab='% of sepsis cases', xlab='Days after diganosis')
for (i in seq_len(ncol(t))) {
   lines(x=0:4, y=t[,i], lwd=3, col=col_vec[colnames(t)[i]])
}
legend('topright', legend=abbr[abx], col=col_vec, lwd=3)
rm(d, t, col_vec, i)


# what do treatment patterns and trajectories look like??
# which concurrent antibiotic combos are common?
X <- df %>% 
   filter(ABX_m_DX %in% 0:4, ABX %in% abx) %>% 
   select(PERSON_ID, SEPSIS_DX_DATE, ABX_m_DX, ABX) %>% 
   distinct() %>%
   mutate(Rx = 1L) %>% 
   tidyr::pivot_wider(names_from=ABX, values_from=Rx, values_fill=0)

N <- D <- setNames(vector('list', length(0:4)), paste0('d', 0:4))
npairs <- (length(abx)^2 - length(abx)) / 2

for (d in 0:4) {
   print(d)
   x <- X %>% filter(ABX_m_DX == d)
   
   # directional dependency
   t <- matrix(NA, nrow=length(abx), ncol=length(abx), dimnames=list(abbr[abx], abbr[abx]))
   for (a in seq_along(abx)) {
      for (b in seq_along(abx)[-a]) {
         t0 <- table(x[[abx[a]]], x[[abx[b]]])
         t[a,b] <- (t0[2,2] / sum(t0[2,])) / (t0[1,2] / sum(t0[1,])) # how much more/less b is used when a used
      }
   }
   N[[d+1]] <- t
   
   # statistical dependency between pairs
   m <- matrix(NA, nrow=2, ncol=npairs, dimnames=list(c('odds_ratio', 'p_value'), rep('', npairs)))
   i <- 1
   for (a in seq_along(abx)[-length(abx)]) {
      for (b in seq_along(abx)[(a+1):length(abx)]) {
         f <- fisher.test(table(x[[abx[a]]], x[[abx[b]]])) # odds ratio
         m['odds_ratio', i] <- f$estimate
         m['p_value', i] <- f$p.value
         colnames(m)[i] <- paste0(abbr[abx[a]], '_', abbr[abx[b]])
         i <- i + 1
      }
   }
   D[[d+1]] <- data.frame(t(m))
}
rm(X, d, x, t, a, b, t0, f, m, npairs, i)

# x <- numeric(length = (length(abx) * length(abx) - length(abx)) / 2)
# i <- 1
# for (a in seq_along(abx)) {
#    for (b in seq_along(abx)[-a]) {
#       x[i] <- N$d0[a,b] / N$d0[b,a]
#       i <- i + 1
#    }
# }
# rm(x, i, a, b)

par(mfrow=c(2,3))
sapply(N, hist, breaks=seq(0,4.5,0.1), xlim=c(0, 4.5))

par(mfrow=c(2,3))
sapply(D, function(x) hist(x['odds_ratio', ], breaks=seq(0,4.5,0.1), xlim=c(0,4.5)))



# what are the pairs of abx that are co-prescribed 
# more or less frequently than expected by chance?
D$d0[which(D$d0$odds_ratio < 1/3),]
D$d1[which(D$d1$odds_ratio < 1/3),]
D$d2[which(D$d2$odds_ratio < 1/3),]
D$d3[which(D$d3$odds_ratio < 1/3),]
D$d4[which(D$d4$odds_ratio < 1/3),]

D$d0[which(D$d0$odds_ratio > 1.5),]
D$d1[which(D$d1$odds_ratio > 1.5),]
D$d2[which(D$d2$odds_ratio > 1.5),]
D$d3[which(D$d3$odds_ratio > 1.5),]
D$d4[which(D$d4$odds_ratio > 1.5),]


# how does this change over time?
# it weakens slightly over time
# meaning that 
par(mfrow=c(1,1))
plot(NA, xlim=c(0,4), ylim=c(-3, 3), ylab='log odds', xlab='Days after sepsis dx')
or <- sapply(D, function(x) x$odds_ratio)
rownames(or) <- rownames(D$d0)
apply(or, 1, function(y) {
   lines(x=0:4, y=log(y), col=ifelse(log(y) < 0, 'red', 'blue'))
})
wilcox.test(or[,1], or[,5], paired=TRUE)
rm(D, N, or)


# what are patients being switched to/from?
# this is a mess, so let's start with where people end up
# do people end up with fewer antibiotics than they started with?


# number of antibiotics given each day, increases initially, then decreases over time
x <- df %>%
   filter(ABX %in% abx, is.na(time) | time > 6L, ABX_m_DX %in% 0:6) %>%
   select(PERSON_ID, SEPSIS_DX_DATE, ABX, ABX_m_DX) %>%
   distinct() %>%
   group_by(PERSON_ID, SEPSIS_DX_DATE) %>%
   reframe(
      d0 = length(unique(ABX[ABX_m_DX == 0L])),
      d1 = length(unique(ABX[ABX_m_DX == 1L])),
      d2 = length(unique(ABX[ABX_m_DX == 2L])),
      d3 = length(unique(ABX[ABX_m_DX == 3L])),
      d4 = length(unique(ABX[ABX_m_DX == 4L])),
      d5 = length(unique(ABX[ABX_m_DX == 5L])),
      d6 = length(unique(ABX[ABX_m_DX == 6L])),
      n0 = length(unique(ABX[ABX_m_DX == 0L])),
      n1 = length(unique(ABX[ABX_m_DX %in% 0:1])),
      n2 = length(unique(ABX[ABX_m_DX %in% 0:2])),
      n3 = length(unique(ABX[ABX_m_DX %in% 0:3])),
      n4 = length(unique(ABX[ABX_m_DX %in% 0:4])),
      n5 = length(unique(ABX[ABX_m_DX %in% 0:5])),
      n6 = length(unique(ABX[ABX_m_DX %in% 0:6]))
   )


par(mfrow=c(2,1))
vioplot::vioplot(x %>% select(d0:d6), ylab='Number of distinct abx given', xlab='Day after sepsis dx')
lines(x=1:7, y=sapply(x %>% select(d0:d6), mean), lwd=4)
plot(NA, xlim=c(0,6), ylim=c(0,6), ylab='Number of distinct abx given', xlab='Day after sepsis dx')
for (i in 1:6) {
   d <- tibble(d0 = unlist(x %>% select(d0:d6) %>% select(i)), d1 = unlist(x %>% select(d0:d6) %>% select(i+1L)))
   d <- d %>% count(d0, d1) %>% mutate(f = n / max(n) * 10)
   for (j in seq_len(nrow(d))) {
      segments(x0=i-1, x1=i, y0=d$d0[j], y1=d$d1[j], lwd=d$f[j])
   }
}

# cumulative number of distinct antibiotics prescribed over time
vioplot::vioplot(x %>% select(n0:n6), ylab='Number of total abx given up to this point', xlab='Day after sepsis dx')
lines(x=1:7, y=sapply(x %>% select(n0:n6), mean), lwd=4)
plot(NA, xlim=c(0,6), ylim=c(0,6), ylab='Number of total abx given up to this point', xlab='Day after sepsis dx')
for (i in 1:6) {
   d <- tibble(d0 = unlist(x %>% select(n0:n6) %>% select(i)), d1 = unlist(x %>% select(n0:n6) %>% select(i+1L)))
   d <- d %>% count(d0, d1) %>% mutate(f = n / max(n) * 10)
   for (j in seq_len(nrow(d))) {
      segments(x0=i-1, x1=i, y0=d$d0[j], y1=d$d1[j], lwd=d$f[j])
   }
}
### LESSON LEARNED: THE SWITCHES MUST BE HAPPENING IN THE FIRST 1-3 DAYS


# we can tell if there was a switch IF:
# the total increased and the number given stayed the same
x %>% filter(d0 == 2L, d6 == 1L)
df %>% filter(PERSON_ID == '1000000486') %>% select(ABX, ABX_m_DX) %>% distinct %>% arrange(ABX_m_DX, ABX) # cut VAN
df %>% filter(PERSON_ID == '1000000623') %>% select(ABX, ABX_m_DX) %>% distinct %>% arrange(ABX_m_DX, ABX)
df %>% filter(PERSON_ID == '1000001270') %>% select(ABX, ABX_m_DX) %>% distinct %>% arrange(ABX_m_DX, ABX) # cut everything to just CFZ, then total switch to TZP
df %>% filter(PERSON_ID == '1000001856') %>% select(ABX, ABX_m_DX) %>% distinct %>% arrange(ABX_m_DX, ABX) %>% print(n=30)
df %>% filter(PERSON_ID == '1000003826') %>% select(ABX, ABX_m_DX) %>% distinct %>% arrange(ABX_m_DX, ABX)
df %>% filter(PERSON_ID == '1000001942') %>% select(ABX, ABX_m_DX) %>% distinct %>% arrange(ABX_m_DX, ABX)

# need to distinguish between:
# responsible: starting broad to cover, then removing some to narrow
# total switches: 





# for patients who end up on just a single drug, what is it and what did they start on?
df %>%
   filter(ABX %in% abx, ABX_m_DX == 0L) %>%
   select(PERSON_ID, SEPSIS_DX_DATE, ABX, ABX_m_DX) %>%
   distinct() %>%
   filter(n() == 1L, .by=c(PERSON_ID, SEPSIS_DX_DATE)) %>%
   count(ABX, sort=TRUE)
   
x <- df %>%
   filter(ABX %in% c('CEFTRIAXONE', 'PIPERACILLIN/TAZOBACTAM'),
          ABX_m_DX %in% c(0, 5L)) %>%
   select(PERSON_ID, SEPSIS_DX_DATE, ABX, ABX_m_DX) %>%
   distinct() %>%
   arrange(PERSON_ID, SEPSIS_DX_DATE, ABX_m_DX, ABX) %>%
   group_by(PERSON_ID, SEPSIS_DX_DATE) %>%
   filter(any(ABX_m_DX == 0L) & any(ABX_m_DX == 5L)) %>%
   reframe(d0 = paste(sort(unique(ABX[ABX_m_DX == 0L])), collapse=', '),
           d1 = paste(sort(unique(ABX[ABX_m_DX == 5L])), collapse=', ')) %>%
   #filter(!any(ABX[ABX_m_DX == 0L] == 'CEFTRIAXONE')) %>%
   #ungroup() %>% filter(PERSON_ID == '1000060603')
   # summarise(d0 = ABX[],
   #           d1 = last(ABX),
   #           .by=c(PERSON_ID, SEPSIS_DX_DATE))# %>%
   select(-PERSON_ID, -SEPSIS_DX_DATE) %>%
   filter(!grepl(',', d0), !grepl(',', d1))

x %>% filter(d0 == 'PIPERACILLIN/TAZOBACTAM') %>% count(d1, sort=TRUE)
x %>% filter(d0 == 'CEFTRIAXONE') %>% count(d1, sort=TRUE)
x %>% count(d0, d1, sort=TRUE)


df



# given 2 days, what are the switch probabilities for each pair?
# switch meaning the early day abx is no longer present at later and vice versa

   
#


# how many received one of the 4 antibiotics of interest?
df %>% group_by(PERSON_ID, SEPSIS_DX_DATE) # 47,891
df %>%
   filter(ABX_m_DX %in% -1:2) %>%
   filter(ABX %in% c('PIPERACILLIN/TAZOBACTAM', 'MEROPENEM', 'CEFEPIME', 'CEFTRIAXONE')) %>%
   group_by(PERSON_ID, SEPSIS_DX_DATE) # 38,317 received one of the above in the first 3 days after sepsis diagnosis

# collapse to 1 row per isolate
trts <- df %>%
   filter(ABX_m_DX %in% -1:2) %>%
   filter(ABX %in% c('PIPERACILLIN/TAZOBACTAM', 'MEROPENEM', 'CEFEPIME', 'CEFTRIAXONE')) %>%
   arrange(PERSON_ID, SEPSIS_DX_DATE, ABX_m_DX, ABX) %>%
   summarise(ABX_DAYS = list(unique(ABX_m_DX)),
             TZP = any(TZP),
             CRO = any(CRO),
             FEP = any(FEP),
             MEM = any(MEM),
             TRT = paste(sort(unique(ABX)), collapse=','),
             .by = c(PERSON_ID, SEPSIS_DX_DATE, time, AGE, GENDER, ISOLATE, ESBL, VRE, MRSA, SEPDX_m_ADMIT, 
                     FACILITY, NURSING_HOME, EMERGENCY_DEPT, SEPSIS_DESCRIPTION, CULTURE_m_SEPDX)) # ADD MORE VARIABLES TO KEEP FOR FUTURE STEPS

# fix up a few variables...
trts$NUM_ABX_DAYS <- lengths(trts$ABX_DAYS)
getMin <- function(l) {
   w <- which(lengths(l) > 0L)
   new_l <- rep(NA, length(l))
   new_l[w] <- sapply(l[w], min)
   return(new_l)
}
trts$FIRST_ABX_DAY <- getMin(trts$ABX_DAYS)
rm(getMin)
trts <- trts %>% select(-ABX_DAYS)

# early death rates differ by antibiotics
# could use an indication to understand treatment propensity based on really acute deaths
trts %>% filter(TZP) %>% summarise(sum(time <= 2L, na.rm=T) / n()) # 4.88%
trts %>% filter(CRO) %>% summarise(sum(time <= 2L, na.rm=T) / n()) # 2.41%
trts %>% filter(FEP) %>% summarise(sum(time <= 2L, na.rm=T) / n()) # 4.89%
trts %>% filter(MEM) %>% summarise(sum(time <= 2L, na.rm=T) / n()) # 5.28%
trts %>% filter(TRT == 'PIPERACILLIN/TAZOBACTAM') %>% summarise(sum(time <= 2L, na.rm=T) / n()) # 5.46%
trts %>% filter(TRT == 'CEFTRIAXONE')             %>% summarise(sum(time <= 2L, na.rm=T) / n()) # 2.04%
trts %>% filter(TRT == 'CEFEPIME')                %>% summarise(sum(time <= 2L, na.rm=T) / n()) # 5.47%
trts %>% filter(TRT == 'MEROPENEM')               %>% summarise(sum(time <= 2L, na.rm=T) / n()) # 5.28%

trts %>%
   filter(TRT %in% c('PIPERACILLIN/TAZOBACTAM', 'CEFTRIAXONE', 'CEFEPIME', 'MEROPENEM')) %>%
   summarise(n = n(),
             mAGE = mean(AGE),
             fMALE = sum(GENDER == 'MALE') / n,
             fISOLATE = sum(ISOLATE, na.rm=T) / n,
             fESBL = sum(ESBL, na.rm=T) / n,
             fMRSA = sum(MRSA, na.rm=T) / n,
             fVRE = sum(VRE, na.rm=T) / n,
             fNH = sum(NURSING_HOME, na.rm=T) / n,
             fED = sum(EMERGENCY_DEPT, na.rm=T) / n,
             first_abx_day = mean(FIRST_ABX_DAY),
             num_abx_days = mean(NUM_ABX_DAYS),
             d30 = sum(time < 30, na.rm=T) / n,
             d14 = sum(time < 14, na.rm=T) / n, 
             .by = TRT)

trts <- trts %>% # 38,317 --> 30,830
   filter(is.na(time) | time > 2L,
          SEPDX_m_ADMIT %in% 0:2, 
          FIRST_ABX_DAY %in% -1:1, 
          NUM_ABX_DAYS >= 2L)

### UPSET PLOT ###
{
   keep <- which((is.na(trts$time) | trts$time > 2L) & trts$SEPDX_m_ADMIT %in% 0:2) # 33,021
   t <- trts %>% 
      slice(keep) %>%
      count(TZP, CRO, FEP, MEM) %>% 
      arrange(desc(n))
   m <- sort(table(unlist(strsplit(trts$TRT[keep], ','))), decreasing=TRUE)
   names(m)[names(m) == 'PIPERACILLIN/TAZOBACTAM'] <- 'Pip-tazo'
   names(m) <- stringr::str_to_sentence(names(m))
   # setup plot
   pdf(file = '~/Desktop/UpSetSepsisAbx.pdf', height=6.8, width=10)
   par(mfrow=c(1,1), mar=c(9, 14, 1, 1), mgp=c(2, 0.4, 0))
   # intersect set barplot
   xvals <- seq_len(nrow(t))
   plot(NA, xlim=range(1, max(xvals)), ylim=c(0, max(t$n)*1.05), bty='n', axes=F, ann=F)
   axis(side=2, las=1, tck=-0.01)
   bar_width <- 0.45
   rect(xleft = xvals - bar_width,
        xright = xvals + bar_width,
        ybottom = 0, ytop = t$n,
        col = ifelse(test = apply(t %>% select(-n), 1, sum) == 1L, yes = '#555555', no = 'gray'))
   text(x=xvals, y=t$n+250, labels=prettyNum(t$n, big.mark=','), xpd=NA)
   text(x=xvals[apply(t %>% select(-n), 1, sum) == 1L]-0.49, y=t$n[apply(t %>% select(-n), 1, sum) == 1L]+750,
        labels=names(m), xpd=NA, adj=c(0, 0.5))
   # set indicator points
   yvals <- seq(-800, -1 * max(t$n) / 2.5, length.out=length(m))
   points(x=rep(xvals, each=length(yvals)), y=rep(yvals, length(xvals)), xpd=NA, pch=16, cex=3, col='gray')
   text(x=0, y=yvals, labels=names(m), xpd=NA, adj=1)
   for (i in seq_len(nrow(t))) {
      w <- which(unlist(t[i, ] %>% select(-n)))
      points(x=rep(i, length(w)), y=yvals[w], pch=16, cex=3, col='black', xpd=NA)
      if (length(w) > 1L)
         segments(x0=i, y0=min(yvals[w]), y1=max(yvals[w]), xpd=NA, lwd=5)
   }
   # the hardest part...marginal barplot
   bar_bot_x <- xvals[1]-3.5
   bar_top_max_x <- bar_bot_x-2.5
   bar_top_x <- (-1*m / max(m) * 2.5) - 2.5
   bar_width <- 450
   rect(xright = bar_bot_x,
        xleft = bar_top_x,
        ybottom = yvals - bar_width,
        ytop = yvals + bar_width,
        xpd = NA,
        col = '#555555')
   ticks <- (-2.5 * seq(0,15000,5000) / max(m)) - 2.5
   segments(x0 = min(ticks), x1 = max(ticks), y0 = yvals[1]+750, xpd=NA)
   segments(x0 = ticks, y0 = yvals[1]+750, y1 = yvals[1]+850, xpd=NA)
   text(x =ticks, y = yvals[1]+950, 
        labels = seq(0, 15000, 5000), xpd=NA, srt=60, adj=0)
   text(x = mean(ticks), y = yvals[1]+2500, labels='Abx counts', xpd=NA)
   
   dev.off()
   rm(bar_bot_x, bar_top_max_x, bar_top_x, bar_width, i, m, ticks, xvals, yvals, t, w, keep)
}

# get rid of early deaths and later in visit, only keep monotherapies
censor_time <- 30
trts <- trts %>%
   mutate(status = ifelse(is.na(time) | time > censor_time, 0L, 1L)) %>%
   mutate(time_censored = ifelse(status == 0L, censor_time, time))



# relationship between BUGs and 30-day mortality? any BUG and MRSA only
trts %>% mutate(d30 = time < 30) %>% select(MRSA, d30) %>% table() %>% fisher.test() # sig
trts %>% mutate(d30 = time < 30) %>% select(VRE, d30) %>% table() %>% fisher.test() # not sig
trts %>% mutate(d30 = time < 30) %>% select(ESBL, d30) %>% table() %>% fisher.test() # not sig
trts %>% mutate(d30 = time < 30) %>% select(ISOLATE, d30) %>% table() %>% fisher.test() # sig

# out of 4 days [-1, +2] from diagnosis, how many days were antibiotics adminstered?
# ceftriaxone given fewer days, implying those patients cleared faster, required less aggressive to begin with?
par(mfrow=c(4,2))
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/antibiotic_names/abbr_named.Rdata')
for (abx in c('PIPERACILLIN/TAZOBACTAM', 'CEFTRIAXONE', 'CEFEPIME', 'MEROPENEM')) {
   x <- trts %>% filter(is.na(time) | time > 2L, TRT == abx)
   frac3om   <- round(sum(x$NUM_ABX_DAYS  %in% 3:5)  / nrow(x) * 100, 1)
   fracEarly <- round(sum(x$FIRST_ABX_DAY %in% -1:1) / nrow(x) * 100, 1)
   plotBarplot(x$NUM_ABX_DAYS, main=abx)
   plotBarplot(x$FIRST_ABX_DAY, main=abx)
   cat(abbr[abx], frac3om, fracEarly, '\n')
}
rm(abx, frac3om, fracEarly, x, abbr)

w <- which(grepl('with septic shock', trts$SEPSIS_DESCRIPTION) & !grepl('without septic shock', trts$SEPSIS_DESCRIPTION))
trts$SEPTIC_SHOCK <- FALSE
trts$SEPTIC_SHOCK[w] <- TRUE
##### END #####
df <- trts
save(df, file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/Sepsis/Sepsis_Dx_Enc_Ast_Abx.Rdata')



### load dataset up to this point = dx + enc + ast + abx + demo
library(dplyr)
load(file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/Sepsis/Sepsis_Dx_Enc_Ast_Abx.Rdata')
source('~/Desktop/WrightLab/UsefulRfuncs/newBarplotFxn.R')
###


##### GET DX CODES #####
source('~/Desktop/EHR/EHR work/config_file.R')
ids <- unique(df$PERSON_ID)
chunks <- mapply(':',
                 seq(1, floor(length(ids) / 1000) * 1000 + 1, 1000),
                 c(seq(1000, floor(length(ids) / 1000) * 1000, 1000), length(ids)))
dx <- tibble()
for (chunk in chunks) {
   print(chunk[1])
   dx <- rbind(
      dx,
      tbl(conn, in_schema('AMB_ETL', 'LAB_SENS_DX_VW')) %>%
         filter(PERSON_ID %in% local(ids[chunk]),
                lubridate::year(SEPSIS_DX_DATE) %in% 2017:2023) %>%
         collect()
   )
}
dx <- dx %>%
   mutate(DX_FROM_DATE = lubridate::date(DX_FROM_DATE),
          DX_TO_DATE = lubridate::date(DX_TO_DATE),
          CODE_DESCRIPTION = tolower(CODE_DESCRIPTION)) %>% 
   select(PERSON_ID, DX_FROM_DATE, DX_TO_DATE, DX_CODE, CODE_DESCRIPTION) %>% 
   distinct()
rm(result, chunks, chunk, ids)
save(dx, file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/Sepsis/Sepsis_OtherDxCodes.Rdata')
##### END #####


##### JOIN DX CODES #####
# join and re-format 1 row per infection, flag for each comorbidity
# load just other dx codes
load(file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/Sepsis/Sepsis_OtherDxCodes.Rdata')

#### WITHIN WEEK BEFORE BLODO CULTURE
dfx <- df %>%
   mutate(JOIN_START = SEPSIS_DX_DATE + FIRST_ABX_DAY - 7,
          JOIN_END = SEPSIS_DX_DATE + FIRST_ABX_DAY) %>%
   left_join(y = dx,
             multiple = 'all',
             by = join_by(
                PERSON_ID,
                between(y$DX_FROM_DATE, x$JOIN_START, x$JOIN_END)
             )) %>%
   select(-JOIN_START, -JOIN_END) %>%
   mutate(X = as.integer(DX_FROM_DATE - SEPSIS_DX_DATE),#) %>% pull(X) %>% plotBarplot()
          CODE_DESCRIPTION = ifelse(is.na(CODE_DESCRIPTION), 'none', CODE_DESCRIPTION)) %>%
   group_by(PERSON_ID, SEPSIS_DX_DATE) %>%
   mutate(
      AKI_1w = any(grepl('^N17.9', DX_CODE)), 
      Endocarditis_1w = any(grepl('^I33.0', DX_CODE)),
      Osteomyelitis_1w = any(grepl('^M86.9', DX_CODE)),
      Cellulitis_1w = any(grepl('^L03.90', DX_CODE)),
      Peritonitis_1w = any(grepl('^K65.9', DX_CODE))
   ) %>%
   ungroup() %>%
   select(-DX_FROM_DATE, -DX_TO_DATE, -DX_CODE, -CODE_DESCRIPTION, -X) %>%
   distinct()
dfx %>% summarise(across(AKI_1w:Peritonitis_1w, ~ sum(.) / n()))
dfx %>% summarise(across(AKI_1w:Peritonitis_1w, ~ sum(.) / n()), .by=TRT)


#### WITHIN MONTH OF BLOOD CULTURE
dfx <- dfx %>%
   mutate(JOIN_START = SEPSIS_DX_DATE - 30,
          JOIN_END = (SEPSIS_DX_DATE + FIRST_ABX_DAY)) %>%
   left_join(y = dx,
             multiple = 'all',
             by = join_by(
                PERSON_ID,
                between(y$DX_FROM_DATE, x$JOIN_START, x$JOIN_END)
             )) %>%
   select(-JOIN_START, -JOIN_END) %>%
   mutate(X = as.integer(DX_FROM_DATE - SEPSIS_DX_DATE),#) %>% pull(X) %>% plotBarplot()
          CODE_DESCRIPTION = ifelse(is.na(CODE_DESCRIPTION), 'none', CODE_DESCRIPTION)) %>%
   group_by(PERSON_ID, SEPSIS_DX_DATE) %>%
   mutate(
      MyocInfarc_1m = any(grepl('^I21|^I22|^I25', DX_CODE)),
      CPD_Pneum_1m = any(grepl('^I26|^I27|^I28\\.[089]|^J4[0-7]|^J6[0-8]|^J70', DX_CODE)),
      MetastSolidTumor_1m = any(grepl('^C7[7-9]|^C80', DX_CODE)),
      Malignancy_1m = any(grepl('^C[01][0-9]|^C2[0-6]|^C3[01234789]|^C4[013]|^C4[5-9]|^C5[0-8]|^C6[0-9]|^C7[0-6]|^C8[123458]|^C9[0-7]', DX_CODE)),
   ) %>%
   ungroup() %>%
   select(-DX_FROM_DATE, -DX_TO_DATE, -DX_CODE, -CODE_DESCRIPTION, -X) %>%
   distinct()
dfx %>% summarise(across(MyocInfarc_1m:Malignancy_1m, ~ sum(.) / n()))
dfx %>% summarise(across(MyocInfarc_1m:Malignancy_1m, ~ sum(.) / n()), .by=TRT)


#### WITHIN 2 YEARS OF BLOOD CULTURE
dfx <- dfx %>%
   mutate(JOIN_START = SEPSIS_DX_DATE - 730,
          JOIN_END = (SEPSIS_DX_DATE + FIRST_ABX_DAY)) %>%
   left_join(y = dx,
             multiple = 'all',
             by = join_by(
                PERSON_ID,
                between(y$DX_FROM_DATE, x$JOIN_START, x$JOIN_END)
             )) %>%
   select(-JOIN_START, -JOIN_END) %>%
   mutate(X = as.integer(DX_FROM_DATE - SEPSIS_DX_DATE),#) %>% pull(X) %>% plotBarplot()
          CODE_DESCRIPTION = ifelse(is.na(CODE_DESCRIPTION), 'none', CODE_DESCRIPTION)) %>%
   group_by(PERSON_ID, SEPSIS_DX_DATE) %>%
   mutate(
      OnDialysis = any(CODE_DESCRIPTION == 'dependence on renal dialysis'),
      Obesity = any(grepl('^E66', DX_CODE)),
      WeightLoss = any(grepl('^E4[0-6]|^R63\\.4|^R64', DX_CODE)),
      Anemia = any(grepl('^D50\\.[089]|^D5[1-3]', DX_CODE)),
      Hypothyroid = any(grepl('^E0[0-3]|^E89', DX_CODE) & !grepl('E000\\.0|E030|E015\\.2|E001\\.0', DX_CODE)),
      FluidElectroDis = any(grepl('^E22\\.2|^E86|^E87', DX_CODE)),
      Coagulopathy = any(grepl('^D6[5-8]|^D69\\.[16789]', DX_CODE)),
      Alcohol = any(grepl('^F10|^E52|^G62\\.1|^I42\\.6|^K29\\.2|^K70\\.0|^K70\\.3|^K70\\.9|^T51|^Z50\\.2|^Z71\\.4|^Z72\\.1', DX_CODE)),
      Drugs = any(grepl('^F1[12345689]|Z71\\.5|Z72\\.2', DX_CODE)),
      Psychoses = any(grepl('^F20|^F2[234589]|F30\\.2|F31\\.2|F31\\.5', DX_CODE)),
      Depression = any(grepl('^F20\\.4|^F31\\.[3-5]|^F32|^F33|^F34\\.1|^F41\\.2|^F43\\.2', DX_CODE)),
      NeuroDisease = any(grepl('^G1[0-3]|^G2[0-2]|^G25\\.4|^G25\\.5|^G31\\.[289]|^G32|^G3[5-7]|^G4[01]|^G93\\.[14]|^R47\\.0|^R56', DX_CODE)),
      CardiacArrythm = any(grepl('^I44\\.[1-3]|^I45\\.6|^I45\\.9|^I4[7-9]|^R00\\.[018]|^T82\\.1|^Z45\\.0|^Z95\\.0', DX_CODE)),
      MyocInfarc = any(grepl('^I21|^I22|^I25', DX_CODE)),
      CompHypertension = any(grepl('^I1[135]', DX_CODE)),
      UncompHypertension = any(grepl('^I10', DX_CODE)),
      CongHeartFailure = any(grepl('^I42|^I43|^I50', DX_CODE)),
      PeriphVasDis = any(grepl('^I70|^I71|^I73|^I79|^K55|^Z95', DX_CODE)),
      CereVasDis = any(grepl('^G45|^G46|^I6[0-9]', DX_CODE)),
      Dementia = any(grepl('^F0[01235]|^G30', DX_CODE)),
      CPD_Pneum = any(grepl('^I26|^I27|^I28\\.[089]|^J4[0-7]|^J6[0-8]|^J70', DX_CODE)), # this includes additional codes from Elixhauser
      RheumaticDis = any(grepl('^L94\\.[013]|^M05|^M06|^M08|^M12\\.[03]|^M3[0-6]|^M45|M46\\.[189]', DX_CODE)), # this includes additional codes from Elixhauser
      PepticUlcerDis = any(grepl('^K2[5-8]', DX_CODE)),
      MildLiverDis = any(grepl('^B18|^K7[1346]|^Z94', DX_CODE)),
      Diabetes = any(grepl('^E1[1-4]', DX_CODE)),
      HemiParaplegia = any(grepl('^G8[0-4]', DX_CODE)),
      RenalDis = any(grepl('^N03|^N05|^N18|^N19|^Z49', DX_CODE)),
      Malignancy = any(grepl('^C[01][0-9]|^C2[0-6]|^C3[01234789]|^C4[013]|^C4[5-9]|^C5[0-8]|^C6[0-9]|^C7[0-6]|^C8[123458]|^C9[0-7]', DX_CODE)),
      ModSevLivDis = any(grepl('I85|K72', DX_CODE)),
      MetastSolidTumor = any(grepl('^C7[7-9]|^C80', DX_CODE)),
      AIDS_HIV = any(grepl('^B2[0124]', DX_CODE)),
      Hyperlipid = any(grepl('(^| )hyperlipid', CODE_DESCRIPTION)),
      Smoking = any(grepl('nicotine', CODE_DESCRIPTION))
   ) %>%
   ungroup() %>%
   select(-DX_FROM_DATE, -DX_TO_DATE, -DX_CODE, -CODE_DESCRIPTION, -X) %>%
   distinct()
dfx %>% summarise(across(OnDialysis:Smoking, ~ sum(.) / n()))
dfx %>% summarise(across(OnDialysis:Smoking, ~ sum(.) / n()), .by=TRT)
##### END #####



# transfers tend to use TZP more often rather than CRO and MEM
# slightly higher 30-day mortality
dfx %>% 
   mutate(mult = grepl(',', FACILITY)) %>%
   filter(TRT %in% c('PIPERACILLIN/TAZOBACTAM', 'CEFTRIAXONE', 'CEFEPIME', 'MEROPENEM')) %>%
   mutate(TRT = case_when(
      TRT == 'PIPERACILLIN/TAZOBACTAM' ~ 'TZP',
      TRT == 'CEFEPIME' ~ 'FEP',
      TRT == 'CEFTRIAXONE' ~ 'CRO',
      TRT == 'MEROPENEM' ~ 'MEM'
   )) %>%
   summarise(n = n(),
             nD30 = sum(time_censored < 30),
             nTZP = sum(TRT == 'TZP'),
             nCRO = sum(TRT == 'CRO'),
             nFEP = sum(TRT == 'FEP'),
             nMEM = sum(TRT == 'MEM'),
             fTZP = nTZP / n,
             fCRO = nCRO / n,
             fFEP = nFEP / n,
             fMEM = nMEM / n,
             fD30 = nD30 / n,
             .by=mult)




save(dfx, file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/Sepsis/Sepsis_AllData.Rdata')
write.table(x = dfx,
            file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/Sepsis/SepsisAllAbx.csv',
            quote = FALSE,
            row.names = FALSE,
            sep = '\t')






####################################################################
###################### PLOT TESTING GROUNDS ########################
####################################################################

##### END #####















