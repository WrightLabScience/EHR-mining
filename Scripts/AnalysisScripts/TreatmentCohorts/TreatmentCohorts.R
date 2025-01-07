# The purpose of this script is to get aggregate counts for potential cohorts
# of empiric / definitive "trials" comparison groups
# 
# Version 1: choose pathogen - link AST orders to antibiotic administration events
#        --ESBL - TZP vs. Carbapenems (lots treated outpatient presumably...)
#        --MDRO Pseudomonas - ceftolozane-tazobactam vs. ceftazidime-avibactam
#        --E. coli - Fluoroquinolones vs. others
#        --MRSA (BSI) - VAN vs. DAP
# 
# Version 2: choose condition - get encounter - link ASTs and AbxAdmin
#        --Sepsis - FEP vs. MER
#        --VAP - Zosyn, FEP, MER
#        Outpatient UTIs - SXT, NIT, etc.
#        --C. diff - VAN vs. fedaxomycin (Clostridium difficile-associated disease)
writeTable <- function(x, fname) {
   write.table(x, quote=F, row.names=F, sep='\t',
               file=paste0('~/Desktop/EHR/EHR-mining/Scripts/AnalysisScripts/TreatmentCohorts/counts_txts/', fname, '.txt'))
}
#


##### ESBL - TZP vs. carbapenems #####
library(dplyr)
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_clean_ASTs.Rdata')

esbl <- astDF %>% filter(ESBL == 1L)
rm(astDF); gc()

esbl <- esbl %>% 
   filter(as.integer(RESULT_DAY - ORDER_DAY) < 14) %>%
   filter(substr(ORDER_DAY,1,4) %in% as.character(2017:2023)) %>% # 44,820
   select(PERSON_ID, ORDER_DAY, RESULT_DAY) %>%
   mutate(JOIN_START = ORDER_DAY - 2,
          JOIN_END = RESULT_DAY + 7) %>%
   distinct()
length(unique(esbl$PERSON_ID)) # 24K

# get anbtiotic administration, restrict to same cohort
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_CLEANED_2017_2023_AbxAdmin.Rdata')
abxDF <- abxDF %>% filter(PERSON_ID %in% unique(esbl$PERSON_ID)) # 1.73 million

# Join
df <- esbl %>%
   left_join(y = abxDF,
             by = join_by(
                PERSON_ID,
                JOIN_START <= START_DAY,
                JOIN_END >= START_DAY
             )) %>%
   select(-JOIN_START, -JOIN_END, -START_DATE, -END_DATE) %>%
   distinct()

# bin into empiric and definitive
df <- df %>%
   mutate(WINDOW = case_when(
      START_DAY <= ORDER_DAY + 1 ~ 'EMPIRIC',
      START_DAY >= RESULT_DAY ~ 'DEFINITIVE'
   ))
df %>% count(WINDOW) # 124K def, 72K emp, 51K none

# how many infections are missing antibiotics? about 1/2
df %>% mutate(hasAbx = !is.na(ABX)) %>% select(PERSON_ID, ORDER_DAY, RESULT_DAY, hasAbx) %>% distinct() %>% count(hasAbx)

# restrict to treated cohorts of interest - TZP vs. Carbapenem
# at this point, we are excluding ESBLs treated outpatient, because medication administration data only!
df <- df %>%
   filter(!is.na(ABX)) %>%
   mutate(TRT = case_when(
      ABX == 'PIPERACILLIN/TAZOBACTAM' ~ 'TZP',
      ABX %in% c('MEROPENEM', 'ERTAPENEM', 'IMIPENEM', 'DORIPENEM') ~ 'CAR',
      .default = 'otherAbx'
   )) %>%
   select(PERSON_ID, ORDER_DAY, TRT, WINDOW) %>%
   distinct()

# did patients receive just one? how often both TZP and CARB? Did they switch?
df <- df %>%
   group_by(PERSON_ID, ORDER_DAY, WINDOW) %>%
   mutate(TRTC = case_when(
      any(TRT == 'TZP') & !any(TRT == 'CAR') ~ 'Pip-tazo',
      !any(TRT == 'TZP') &  any(TRT == 'CAR') ~ 'Carbapenem',
      any(TRT == 'TZP') &  any(TRT == 'CAR') ~ 'both',
      !any(TRT == 'TZP') & !any(TRT == 'CAR') ~ 'otherAbx'
   )) %>%
   ungroup() %>%
   select(PERSON_ID, ORDER_DAY, WINDOW, TRTC) %>%
   distinct()
df <- df %>% mutate(WINDOW = ifelse(is.na(WINDOW), 'BETWEEN', WINDOW))
df <- df %>%
   tidyr::pivot_wider(
      values_from = TRTC,
      names_from = WINDOW,
      values_fill = 'noAbx'
   ) %>%
   relocate(BETWEEN, .after=EMPIRIC)

e <- df %>% count(EMPIRIC) %>% rename(EmpiricTherapy = EMPIRIC)
d <- df %>% count(DEFINITIVE) %>% rename(DefinitiveTherapy = DEFINITIVE)
t <- df %>% 
   count(EMPIRIC, BETWEEN, DEFINITIVE, sort=TRUE) %>%
   filter(n > 100L) %>%
   rename(EmpiricTherapy = EMPIRIC,
          BetweenTherapy = BETWEEN,
          DefinitiveTherapy = DEFINITIVE)
writeTable(e, 'esbl_TZP_CAR_empiric')
writeTable(d, 'esbl_TZP_CAR_definitive')
writeTable(t, 'esbl_TZP_CAR_transitions')
rm(e, d, t, writeTable)

# crude survival analysis
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_DEMO.Rdata')
df <- df %>% 
   filter(DEFINITIVE %in% c('Carbapenem', 'Pip-tazo')) %>%
   select(-EMPIRIC, -BETWEEN) %>%
   rename(TRT = DEFINITIVE)

df <- df %>%
   left_join(y = dth,
             by = join_by(PERSON_ID))
df %>% summarise(median(as.numeric(lubridate::as.duration(ORDER_DAY - DOB)) / (365 * 24 * 60 * 60)), .by=TRT) # extremely similar median ages
df %>% summarise(sum(GENDER == 'MALE') / n(), .by=TRT) # women slightly more likely to pip-tazo

censor_time <- 90
dfs <- df %>%
   select(-GENDER, -DOB, -PATIENT_STATUS) %>%
   mutate(time = as.integer(DEATH_DATE - ORDER_DAY)) %>%
   mutate(status = case_when(
      is.na(time) | time > censor_time ~ 0,
      .default = 1
   )) %>%
   mutate(time = case_when(
      is.na(time) | time > censor_time ~ censor_time,
      .default = time
   )) %>%
   select(-DEATH_DATE) %>%
   filter(time > 2)

sample_sizes <- table(dfs$TRT)
fit <- survfit(Surv(time, status) ~ TRT, data = dfs)
coxph(Surv(time, status) ~ TRT, data = dfs)
pval <- survdiff(Surv(time, status) ~ TRT, data = dfs)$pvalue
d <- setNames(lapply(unique(summary(fit)$strata), function(x) which(summary(fit)$strata == x)),
              as.character(unique(summary(fit)$strata)))
s <- summary(fit)
xvals <- s$time
yvals <- s$surv
ste <- s$std.err
col_vec <- c('Carbapenems' = '#FF0000', 'Pip-tazo' = '#0000FF')
ypos <- seq(0.575, 0.52, length.out=3)

plot(NA, ylim=c(0.5,1), xlim=c(0, censor_time), xaxt='n', yaxt='n', xaxs='i', yaxs='i', xpd=NA, xlab='Days', ylab='Surviving fraction')
title(main = 'ESBL - Carbapenems vs. Pip-tazo')
plotSurv(xvals, yvals, ste, d, col_vec)
axis(side=1, at=seq(0, censor_time, 15))
axis(side=2, las=1)
abline(v=30, lty=2)
text(x=2, y=ypos[1], adj=0, labels=paste0('p = ', prettyNum(pval, digits=3, scientific=-1)), font=ifelse(pval<0.05, 2, 1))
text(x=2, y=ypos[2:3], adj=0, col=col_vec,
     labels=paste0(names(col_vec), ' (n = ', prettyNum(sample_sizes, big.mark=','), ')'))
##### END #####



##### ESBL Bacteremia - TZP vs. carbapenems #####
library(dplyr)
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ASTs_AbxAdmin_blood_2017_2023_imputed_encs_flagged_WIDE_survival.Rdata')
df <- empDF %>% 
   filter(ESBL == 1L) %>%
   select(PERSON_ID, ORDER_DAY, AGE, ABX_TAR, ABX_TARl, SURV_TIME) %>%
   distinct()

tzp <- 'PIPERACILLIN/TAZOBACTAM'
car <- 'MEROPENEM|ERTAPENEM|DORIPENEM|IMIPENEM'

df <- df %>%
   summarise(ABX = paste(sort(unique(unlist(ABX_TAR), unlist(ABX_TARl))), collapse=','),
             .by = c(PERSON_ID, ORDER_DAY, AGE, SURV_TIME)) %>%
   mutate(TRT = case_when(
      ABX == '' ~ 'noAbx',
       grepl(tzp, ABX) & !grepl(car, ABX) ~ 'TZP',
      !grepl(tzp, ABX) &  grepl(car, ABX) ~ 'CAR',
       grepl(tzp, ABX) &  grepl(car, ABX) ~ 'both',
      ABX != '' & !grepl(tzp, ABX) & !grepl(car, ABX) ~ 'otherAbx',
   ))

library(survival)
censor_time <- 90
df <- df %>%
   rename(time = SURV_TIME) %>%
   mutate(status = case_when(
      is.na(time) | time > censor_time ~ 0,
      .default = 1
   )) %>%
   mutate(time = case_when(
      is.na(time) | time > censor_time ~ censor_time,
      .default = time
   )) %>%
   filter(TRT %in% c('TZP', 'CAR'))

sample_sizes <- table(df$TRT)
fit <- survfit(Surv(time, status) ~ TRT, data = df)
cox <- coxph(Surv(time, status) ~ TRT + AGE, data = df)
pval <- summary(cox)[['coefficients']][,'Pr(>|z|)'][1]
#pval <- survdiff(Surv(time, status) ~ TRT, data = df)$pvalue
d <- setNames(lapply(unique(summary(fit)$strata), function(x) which(summary(fit)$strata == x)),
              as.character(unique(summary(fit)$strata)))
s <- summary(fit)
xvals <- s$time
yvals <- s$surv
ste <- s$std.err
col_vec <- c('CAR' = '#FF0000', 'TZP' = '#0000FF')
ypos <- seq(0.575, 0.52, length.out=3)

plot(NA, ylim=c(0.5,1), xlim=c(0, censor_time), xaxt='n', yaxt='n', xaxs='i', yaxs='i', xpd=NA, xlab='Days', ylab='Surviving fraction')
plotSurv(xvals, yvals, ste, d, col_vec)
axis(side=1, at=seq(0, censor_time, 15))
axis(side=2, las=1)
abline(v=30, lty=2)
text(x=2, y=ypos[1], adj=0, labels=paste0('p = ', prettyNum(pval, digits=3, scientific=-1)), font=ifelse(pval<0.05, 2, 1))
text(x=2, y=ypos[2:3], adj=0, col=col_vec,
     labels=paste0(names(col_vec), ' (n = ', prettyNum(sample_sizes, big.mark=','), ')'))
#####




##### Sepsis - FEP vs. MER #####
library(dplyr)
source('~/Desktop/EHR/EHR work/config_file.R')

# first get all sepsis diagnosis dates
dx <- tibble(dbGetQuery(conn, "select * from AMB_ETL.LAB_SENS_DX_VW where CODE_DESCRIPTION like '%Sepsis%'"))
dx <- dx %>% 
   filter(substr(DX_FROM_DATE,1,4) %in% as.character(2017:2023)) %>% 
   mutate(DX_FROM_DATE = as.Date(substr(DX_FROM_DATE,1,10))) %>%
   select(PERSON_ID, CODE_DESCRIPTION, DX_FROM_DATE) %>% 
   distinct() %>% 
   arrange(PERSON_ID, DX_FROM_DATE, CODE_DESCRIPTION) %>% 
   summarise(CODE_DESCRIPTION = paste(sort(unique(CODE_DESCRIPTION)), collapse=' + '), 
             .by=c(PERSON_ID, DX_FROM_DATE))

# now get ASTs
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_clean_ASTs.Rdata')
astDF <- astDF %>%
   filter(PERSON_ID %in% unique(dx$PERSON_ID),
          substr(ORDER_DAY,1,4) %in% as.character(2017:2023)) # 1.7 mil --> 215K

df <- dx %>%
   mutate(JOIN_START = DX_FROM_DATE - 4,
          JOIN_END = DX_FROM_DATE + 4) %>%
   left_join(y = astDF %>% select(PERSON_ID, ORDER_DAY, RESULT_DAY, BUG) %>% distinct(),
             multiple = 'first',
             by = join_by(
                PERSON_ID,
                between(y$ORDER_DAY, x$JOIN_START, x$JOIN_END)
             )) %>%
   select(-JOIN_START, -JOIN_END)
df %>% count(is.na(ORDER_DAY)) # ~1/2 have bugs
df <- df %>% mutate(X = as.integer(ORDER_DAY - DX_FROM_DATE))
barplot(table(df$X))

# join dx+AST with abxDF
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_CLEANED_2017_2023_AbxAdmin.Rdata')
abxDF <- abxDF %>% 
   filter(PERSON_ID %in% unique(df$PERSON_ID)) %>% # 9.7 mil --> 1.9 mil
   select(PERSON_ID, START_DAY, ABX) %>%
   distinct()

df$JOIN_START <- df$DX_FROM_DATE - 2
w <- which((df$ORDER_DAY - 2) < df$JOIN_START)
df$JOIN_START[w] <- df$ORDER_DAY[w] - 2
df$JOIN_END <- df$DX_FROM_DATE + 7
w <- which((df$RESULT_DAY + 7) > df$JOIN_END)
df$JOIN_END[w] <- df$RESULT_DAY[w] + 7

df2 <- df %>%
   mutate(JOIN_START = DX_FROM_DATE - 2,
          JOIN_END = DX_FROM_DATE + 14) %>%
   left_join(y = abxDF,
             by = join_by(
                PERSON_ID,
                JOIN_START <= START_DAY,
                JOIN_END >= START_DAY
             )) %>%
   select(-JOIN_START, -JOIN_END) %>% #count(hasBug = !is.na(BUG), hasAbx = !is.na(ABX))
   mutate(FEP = ABX == 'CEFEPIME',
          MER = ABX == 'MEROPENEM')

e <- df2 %>%
   filter(!is.na(BUG)) %>%
   filter(START_DAY >= ORDER_DAY - 2 & START_DAY <= ORDER_DAY + 1) %>%
   select(PERSON_ID, DX_FROM_DATE, ORDER_DAY, FEP, MER) %>%
   distinct() %>%
   group_by(PERSON_ID, ORDER_DAY, DX_FROM_DATE) %>%
   mutate(TRT = case_when(
      !any(FEP) & any(MER) ~ 'MER',
      any(FEP) & any(MER) ~ 'both',
      any(FEP) & !any(MER) ~ 'FEP',
      is.na(FEP) & is.na(MER) ~ 'noAbx',
      !any(FEP) & !any(MER) ~ 'otherAbx'
   )) %>%
   ungroup() %>%
   select(PERSON_ID, DX_FROM_DATE, ORDER_DAY, TRT) %>%
   distinct() %>%
   count(TRT) %>%
   rename(EmpiricTherapy = TRT)
writeTable(e, 'sepsis_FEP_MER_empiric_bug')

d <- df2 %>%
   filter(!is.na(BUG)) %>%
   filter(START_DAY >= RESULT_DAY - 1 & START_DAY <= RESULT_DAY + 7) %>%
   select(PERSON_ID, DX_FROM_DATE, ORDER_DAY, FEP, MER) %>% 
   distinct() %>%
   group_by(PERSON_ID, ORDER_DAY, DX_FROM_DATE) %>%
   mutate(TRT = case_when(
      !any(FEP) & any(MER) ~ 'MER',
      any(FEP) & any(MER) ~ 'both',
      any(FEP) & !any(MER) ~ 'FEP',
      is.na(FEP) & is.na(MER) ~ 'noAbx',
      !any(FEP) & !any(MER) ~ 'otherAbx'
   )) %>%
   ungroup() %>%
   select(PERSON_ID, DX_FROM_DATE, ORDER_DAY, TRT) %>%
   distinct() %>%
   count(TRT) %>%
   rename(DefinitiveTherapy = TRT)
writeTable(d, 'sepsis_FEP_MER_definitive_bug')

a <- df2 %>%
   filter(!is.na(BUG)) %>%
   select(PERSON_ID, DX_FROM_DATE, ORDER_DAY, FEP, MER) %>% 
   distinct() %>%
   group_by(PERSON_ID, ORDER_DAY, DX_FROM_DATE) %>%
   mutate(TRT = case_when(
      !any(FEP) & any(MER) ~ 'MER',
      any(FEP) & any(MER) ~ 'both',
      any(FEP) & !any(MER) ~ 'FEP',
      is.na(FEP) & is.na(MER) ~ 'noAbx',
      !any(FEP) & !any(MER) ~ 'otherAbx'
   )) %>%
   ungroup() %>%
   select(PERSON_ID, DX_FROM_DATE, ORDER_DAY, TRT) %>%
   distinct() %>%
   count(TRT) %>%
   rename(PathogenIDedTherapy = TRT)
writeTable(a, 'sepsis_FEP_MER_all_bug')

n <- df2 %>%
   filter(is.na(BUG)) %>%
   select(PERSON_ID, DX_FROM_DATE, FEP, MER) %>% 
   distinct() %>%
   group_by(PERSON_ID, DX_FROM_DATE) %>%
   mutate(TRT = case_when(
      !any(FEP) & any(MER) ~ 'MER',
      any(FEP) & any(MER) ~ 'both',
      any(FEP) & !any(MER) ~ 'FEP',
      is.na(FEP) & is.na(MER) ~ 'noAbx',
      !any(FEP) & !any(MER) ~ 'otherAbx'
   )) %>%
   ungroup() %>%
   select(PERSON_ID, DX_FROM_DATE, TRT) %>%
   distinct() %>%
   count(TRT) %>%
   rename(NoPathogenTherapy = TRT)
writeTable(n, 'sepsis_FEP_MER_nobug')

o <- df2 %>%
   select(PERSON_ID, DX_FROM_DATE, FEP, MER) %>% 
   distinct() %>%
   group_by(PERSON_ID, DX_FROM_DATE) %>%
   mutate(TRT = case_when(
      !any(FEP) & any(MER) ~ 'MER',
      any(FEP) & any(MER) ~ 'both',
      any(FEP) & !any(MER) ~ 'FEP',
      is.na(FEP) & is.na(MER) ~ 'noAbx',
      !any(FEP) & !any(MER) ~ 'otherAbx'
   )) %>%
   ungroup() %>%
   select(PERSON_ID, DX_FROM_DATE, TRT) %>%
   distinct() %>%
   count(TRT) %>%
   rename(OverallTherapy = TRT)
writeTable(o, 'sepsis_FEP_MER_overall')
rm(o, e, d, n, a, writeTable, df, w, abxDF, astDF, dx)


# crude survival analysis
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_DEMO.Rdata')
df <- df2 %>%
   select(PERSON_ID, DX_FROM_DATE, FEP, MER) %>% 
   distinct() %>%
   group_by(PERSON_ID, DX_FROM_DATE) %>%
   mutate(TRT = case_when(
      !any(FEP) & any(MER) ~ 'MER',
      any(FEP) & any(MER) ~ 'both',
      any(FEP) & !any(MER) ~ 'FEP',
      is.na(FEP) & is.na(MER) ~ 'noAbx',
      !any(FEP) & !any(MER) ~ 'otherAbx'
   )) %>%
   ungroup() %>%
   select(PERSON_ID, DX_FROM_DATE, TRT) %>%
   distinct() %>%
   filter(TRT %in% c('MER', 'FEP'))

df <- df %>%
   left_join(y = dth,
             by = join_by(PERSON_ID)) %>% 
   filter(!is.na(DOB))
df %>% summarise(median(as.numeric(lubridate::as.duration(DX_FROM_DATE - DOB)) / (365 * 24 * 60 * 60)), .by=TRT) # similar median ages
df %>% summarise(sum(GENDER == 'MALE') / n(), .by=TRT) # similar

censor_time <- 90
dfs <- df %>%
   select(-GENDER, -DOB, -PATIENT_STATUS) %>%
   mutate(time = as.integer(DEATH_DATE - DX_FROM_DATE)) %>%
   mutate(status = case_when(
      is.na(time) | time > censor_time ~ 0,
      .default = 1
   )) %>%
   mutate(time = case_when(
      is.na(time) | time > censor_time ~ censor_time,
      .default = time
   )) %>%
   select(-DEATH_DATE) %>%
   filter(time > 2)

sample_sizes <- table(dfs$TRT)
fit <- survfit(Surv(time, status) ~ TRT, data = dfs)
coxph(Surv(time, status) ~ TRT, data = dfs)
pval <- survdiff(Surv(time, status) ~ TRT, data = dfs)$pvalue
d <- setNames(lapply(unique(summary(fit)$strata), function(x) which(summary(fit)$strata == x)),
              as.character(unique(summary(fit)$strata)))
s <- summary(fit)
xvals <- s$time
yvals <- s$surv
ste <- s$std.err
col_vec <- c('DAP' = '#FF0000', 'VAN' = '#0000FF')
ypos <- seq(0.575, 0.52, length.out=3)

plot(NA, ylim=c(0.5,1), xlim=c(0, censor_time), xaxt='n', yaxt='n', xaxs='i', yaxs='i', xpd=NA, xlab='Days', ylab='Surviving fraction')
title(main = 'Sepsis - Cefepime vs. Meropenem')
plotSurv(xvals, yvals, ste, d, col_vec)
axis(side=1, at=seq(0, censor_time, 15))
axis(side=2, las=1)
abline(v=30, lty=2)
text(x=2, y=ypos[1], adj=0, labels=paste0('p = ', prettyNum(pval, digits=3, scientific=-1)), font=ifelse(pval<0.05, 2, 1))
text(x=2, y=ypos[2:3], adj=0, col=col_vec,
     labels=paste0(names(col_vec), ' (n = ', prettyNum(sample_sizes, big.mark=','), ')'))
##### END #####



##### E. coli Fluoroquinolones vs. others #####
library(dplyr)
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_clean_ASTs.Rdata')
ecol <- astDF %>% filter(BUG == 'Escherichia coli')
ecol <- ecol %>% filter(LEVOFLOXACIN == 0L | CIPROFLOXACIN == 0L) # 600K --> 444K

# load abx admin
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_CLEANED_2017_2023_AbxAdmin.Rdata')
abxDF <- abxDF %>% filter(PERSON_ID %in% unique(ecol$PERSON_ID)) # 9.7 mil --> 3.7 mil

df <- ecol %>%
   select(PERSON_ID, ORDER_DAY, RESULT_DAY) %>%
   distinct() %>%
   mutate(JOIN_START = ORDER_DAY - 2,
          JOIN_END = RESULT_DAY + 7) %>%
   left_join(y = abxDF,
             by = join_by(PERSON_ID,
                          JOIN_START <= START_DAY,
                          JOIN_END >= START_DAY)) %>%
   mutate(WINDOW = case_when(
      START_DAY >= ORDER_DAY - 2 & START_DAY <= ORDER_DAY + 1 ~ 'EMPIRIC',
      START_DAY >= RESULT_DAY ~ 'DEFINITIVE'
   )) %>%
   filter(!is.na(WINDOW))

df <- df %>% mutate(FLQ = ABX %in% grep('OXACIN', unique(ABX), value=TRUE)) 
df %>% filter(FLQ) %>% count(ABX) # 39K CIP, 15K LVX, 1K MOX, few GAT and OFL
df %>% group_by(PERSON_ID, ORDER_DAY) # ~72K infections

df <- df %>%
   select(PERSON_ID, ORDER_DAY, WINDOW, FLQ) %>%
   distinct() %>%
   group_by(PERSON_ID, ORDER_DAY, WINDOW) %>%
   mutate(TRT = case_when(
      any(FLQ) & !any(!FLQ) ~ 'FLQ',
      any(FLQ) & any(!FLQ) ~ 'both',
      !any(FLQ) & any(!FLQ) ~ 'other'
   )) %>%
   ungroup() %>%
   select(PERSON_ID, ORDER_DAY, WINDOW, TRT) %>%
   distinct()

e <- df %>% filter(WINDOW == 'EMPIRIC') %>% count(TRT) %>% rename(EmpiricTherapy = TRT)
d <- df %>% filter(WINDOW == 'DEFINITIVE') %>% count(TRT) %>% rename(DefinitiveTherapy = TRT)
t <- df %>%
   tidyr::pivot_wider(
      values_from = TRT,
      names_from = WINDOW
   ) %>% 
   relocate(EMPIRIC, .before=DEFINITIVE) %>%
   count(EMPIRIC, DEFINITIVE, sort=TRUE) %>%
   mutate(across(.cols = c(EMPIRIC, DEFINITIVE),
                 .fns = ~ case_when(
                    is.na(.) ~ ' ',
                    .default = .
                 ))) %>%
   rename(EmpiricTherapy = EMPIRIC,
          DefinitiveTherapy = DEFINITIVE)

writeTable(e, 'ecoli_FLQ_empiric')
writeTable(d, 'ecoli_FLQ_definitive')
writeTable(t, 'ecoli_FLQ_transitions')
rm(e, d, t, writeTable, ecol, abxDF, astDF)
gc()


# crude survival analysis
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_DEMO.Rdata')
df <- df %>% 
   filter(WINDOW == 'DEFINITIVE') %>% 
   filter(TRT != 'both')

df <- df %>%
   left_join(y = dth,
             by = join_by(PERSON_ID)) %>% 
   filter(!is.na(DOB))
df %>% summarise(median(as.numeric(lubridate::as.duration(ORDER_DAY - DOB)) / (365 * 24 * 60 * 60)), .by=TRT) # similar median ages
df %>% summarise(sum(GENDER == 'MALE') / n(), .by=TRT) # FL1 = 81% female, other = 71% female

censor_time <- 90
dfs <- df %>%
   mutate(AGE = as.integer(ORDER_DAY - DOB)) %>%
   select(-DOB, -PATIENT_STATUS) %>%
   mutate(time = as.integer(DEATH_DATE - ORDER_DAY)) %>%
   mutate(status = case_when(
      is.na(time) | time > censor_time ~ 0,
      .default = 1
   )) %>%
   mutate(time = case_when(
      is.na(time) | time > censor_time ~ censor_time,
      .default = time
   )) %>%
   select(-DEATH_DATE) %>%
   filter(time > 2)

sample_sizes <- table(dfs$TRT)
fit <- survfit(Surv(time, status) ~ TRT, data = dfs)
coxph(Surv(time, status) ~ TRT, data = dfs)
pval <- survdiff(Surv(time, status) ~ TRT, data = dfs)$pvalue
d <- setNames(lapply(unique(summary(fit)$strata), function(x) which(summary(fit)$strata == x)),
              as.character(unique(summary(fit)$strata)))
s <- summary(fit)
xvals <- s$time
yvals <- s$surv
ste <- s$std.err
col_vec <- c('FLQ' = '#FF0000', 'other' = '#0000FF')
ypos <- seq(0.575, 0.52, length.out=3)

plot(NA, ylim=c(0.5,1), xlim=c(0, censor_time), xaxt='n', yaxt='n', xaxs='i', yaxs='i', xpd=NA, xlab='Days', ylab='Surviving fraction')
title(main = 'E. coli - Fluoroquinolones vs. others')
plotSurv(xvals, yvals, ste, d, col_vec)
axis(side=1, at=seq(0, censor_time, 15))
axis(side=2, las=1)
abline(v=30, lty=2)
text(x=2, y=ypos[1], adj=0, labels=paste0('p = ', prettyNum(pval, digits=3, scientific=-1)), font=ifelse(pval<0.05, 2, 1))
text(x=2, y=ypos[2:3], adj=0, col=col_vec,
     labels=paste0(names(col_vec), ' (n = ', prettyNum(sample_sizes, big.mark=','), ')'))
##### END #####



##### MDRO Pseudomonas - Ceftolozane-tazobactams vs. Ceftazidime-avibactam #####
# nonsusceptibility to ≥1 antibiotic in ≥3 classes: 
#     antipseudomonal penicillins [e.g., piperacillin–tazobactam]
#     cephalosporins [ceftazidime, cefepime]
#     fluoroquinolones
#     aminoglycosides
#     carbapenems [meropenem, imipenem]
library(dplyr)
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_clean_ASTs.Rdata')

pseu <- astDF %>%
   filter(any(BUG == 'Pseudomonas aeruginosa'), .by=c(PERSON_ID, ORDER_DAY)) %>%
   group_by(PERSON_ID, ORDER_DAY) %>%
   mutate(MULT_ISO = n() > 1L) %>%
   filter(BUG == 'Pseudomonas aeruginosa') %>%
   ungroup() %>%
   mutate(APP_R = case_when(`PIPERACILLIN/TAZOBACTAM` == 1L ~ 1L, .default = 0L),
          CEF_R = case_when(CEFTAZIDIME == 1L | CEFEPIME == 1L ~ 1L, .default = 0L),
          FLQ_R = case_when(CIPROFLOXACIN == 1L | LEVOFLOXACIN == 1L ~ 1L, .default = 0L),
          AMI_R = case_when(TOBRAMYCIN == 1L | GENTAMICIN == 1L | AMIKACIN == 1L ~ 1L, .default = 0L),
          CAR_R = case_when(MEROPENEM == 1L | IMIPENEM == 1L ~ 1L, .default = 0L)) %>% 
   mutate(MDR_SCORE = APP_R + CEF_R + FLQ_R + AMI_R + CAR_R)

mdr <- pseu %>%
   filter(MDR_SCORE >= 3) %>%
   select(PERSON_ID, ORDER_DAY, RESULT_DAY, MULT_ISO, APP_R, CEF_R, FLQ_R, AMI_R, CAR_R, MDR_SCORE) %>%
   distinct()

# join with antibiotics
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_CLEANED_2017_2023_AbxAdmin.Rdata')
abxDF <- abxDF %>% filter(PERSON_ID %in% unique(mdr$PERSON_ID)) # 9.7 mil --> 600K


# join
df <- mdr %>% # 123K
   mutate(JOIN_START = ORDER_DAY - 2,
          JOIN_END = RESULT_DAY + 7) %>%
   left_join(y = abxDF,
             by = join_by(
                PERSON_ID,
                JOIN_START <= START_DAY,
                JOIN_END >= START_DAY
             )) %>%
   select(-JOIN_START, -JOIN_END, -START_DATE, -END_DATE) %>%
   mutate(WINDOW = case_when(
      START_DAY >= ORDER_DAY - 2 & START_DAY <= ORDER_DAY + 1 ~ 'EMPIRIC',
      START_DAY >= RESULT_DAY ~ 'DEFINITIVE'
   )) %>%
   filter(!is.na(WINDOW)) %>% # this step removes patients without definitive inpatient treatment
   select(-START_DAY) %>% 
   distinct()

x <- df %>%
   mutate(ABX = ifelse(ABX %in% c('CEFTOLOZANE/TAZOBACTAM', 'CEFTAZIDIME/AVIBACTAM'), ABX, 'other')) %>%
   distinct() %>%
   mutate(CT = ABX == 'CEFTOLOZANE/TAZOBACTAM',
          CZA = ABX == 'CEFTAZIDIME/AVIBACTAM') %>%
   group_by(PERSON_ID, ORDER_DAY, WINDOW) %>%
   mutate(TRT = case_when(
      any(CT) & any(CZA) ~ 'both',
      !any(CT) & any(CZA) ~ 'CZA',
      any(CT) & !any(CZA) ~ 'CT',
      .default = 'neither'
   )) %>%
   ungroup() %>%
   select(PERSON_ID, ORDER_DAY, WINDOW, TRT) %>%
   distinct()

e <- x %>% filter(WINDOW == 'EMPIRIC') %>% count(TRT) %>% rename(EmpiricTherapy = TRT)
d <- x %>% filter(WINDOW == 'DEFINITIVE') %>% count(TRT) %>% rename(DefinitiveTherapy = TRT)
t <- x %>% 
   tidyr::pivot_wider(values_from=TRT, names_from=WINDOW) %>% 
   count(EMPIRIC, DEFINITIVE, sort=TRUE) %>%
   filter(EMPIRIC %in% c('CT', 'CZA') | DEFINITIVE %in% c('CT', 'CZA')) %>%
   rename(EmpiricTherapy = EMPIRIC, DefinitiveTherapy = DEFINITIVE)

writeTable(e, 'MDRPA_CT_CZA_empiric')
writeTable(d, 'MDRPA_CT_CZA_definitive')
writeTable(t, 'MDRPA_CT_CZA_transitions')
rm(e, d, t, writeTable, abxDF, astDF, mdr, pseu)



# crude survival analysis
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_DEMO.Rdata')
df <- x %>% 
   filter(WINDOW == 'DEFINITIVE') %>% 
   filter(TRT %in% c('CT', 'CZA'))
rm(x)

df <- df %>%
   left_join(y = dth,
             by = join_by(PERSON_ID)) %>% 
   filter(!is.na(DOB))
df %>% summarise(median(as.numeric(lubridate::as.duration(ORDER_DAY - DOB)) / (365 * 24 * 60 * 60)), .by=TRT) # similar median ages
df %>% summarise(sum(GENDER == 'MALE') / n(), .by=TRT) # FL1 = 81% female, other = 71% female

censor_time <- 90
dfs <- df %>%
   mutate(AGE = as.integer(ORDER_DAY - DOB)) %>%
   select(-DOB, -PATIENT_STATUS) %>%
   mutate(time = as.integer(DEATH_DATE - ORDER_DAY)) %>%
   mutate(status = case_when(
      is.na(time) | time > censor_time ~ 0,
      .default = 1
   )) %>%
   mutate(time = case_when(
      is.na(time) | time > censor_time ~ censor_time,
      .default = time
   )) %>%
   select(-DEATH_DATE) %>%
   filter(time > 2)

sample_sizes <- table(dfs$TRT)
fit <- survfit(Surv(time, status) ~ TRT, data = dfs)
coxph(Surv(time, status) ~ TRT + GENDER + AGE, data = dfs)
pval <- survdiff(Surv(time, status) ~ TRT, data = dfs)$pvalue
d <- setNames(lapply(unique(summary(fit)$strata), function(x) which(summary(fit)$strata == x)),
              as.character(unique(summary(fit)$strata)))
s <- summary(fit)
xvals <- s$time
yvals <- s$surv
ste <- s$std.err
col_vec <- c('Ceftolozane-tazobactam' = '#FF0000', 'Ceftazidime-avibactam' = '#0000FF')
ypos <- seq(0.575, 0.52, length.out=3)

plot(NA, ylim=c(0.5,1), xlim=c(0, censor_time), xaxt='n', yaxt='n', xaxs='i', yaxs='i', xpd=NA, xlab='Days', ylab='Surviving fraction')
title(main = 'MDRPA\nCeftolozane-tazobactam vs. Ceftazidime-avibactam')
plotSurv(xvals, yvals, ste, d, col_vec)
axis(side=1, at=seq(0, censor_time, 15))
axis(side=2, las=1)
abline(v=30, lty=2)
text(x=2, y=ypos[1], adj=0, labels=paste0('p = ', prettyNum(pval, digits=3, scientific=-1)), font=ifelse(pval<0.05, 2, 1))
text(x=2, y=ypos[2:3], adj=0, col=col_vec,
     labels=paste0(names(col_vec), ' (n = ', prettyNum(sample_sizes, big.mark=','), ')'))
##### END #####



##### C. diff - vancomycin vs. fidaxomicin (DIFICID) #####
# ICD-10: A04.71, A04.72, A04.7
# ICD-9: 008.45
# 506 +toxin +code
# 156 +toxin -code
# 239 -toxin +code
# 
# ~1/2 do not have any antibiotics administered in the days surrounding their Cdiff diagnosis
library(dplyr)
source('~/Desktop/EHR/EHR work/config_file.R')
dx10 <- tibble(dbGetQuery(conn, "select * from AMB_ETL.LAB_SENS_DX_VW where DX_CODE like '%A04.7%'"))
dx9 <- tibble(dbGetQuery(conn, "select * from AMB_ETL.LAB_SENS_DX_VW where DX_CODE like '%008.45%'"))
dx <- rbind(dx9, dx10) %>% 
   filter(substr(DX_FROM_DATE,1,4) %in% as.character(2017:2023)) %>% 
   arrange(PERSON_ID, DX_FROM_DATE) %>%
   mutate(DAY = as.Date(substr(DX_FROM_DATE,1,10)))
rm(dx9, dx10)
dx %>% count(PRIMARY_DX_IND) # 6203 Y, 12346 N

# join with ASTs? do they match? exclude these later to determine treatment for C diff alone??
load('~/Desktop/EHR/EHR work/RdataFiles/ALL_clean_ASTs.Rdata')
astDF <- astDF %>% 
   filter(PERSON_ID %in% unique(dx$PERSON_ID)) %>% 
   select(PERSON_ID, ORDER_DAY) %>% 
   distinct()

dx <- dx %>%
   mutate(JOIN_START = DAY - 2, JOIN_END = DAY + 2) %>%
   left_join(y = astDF,
             multiple = 'first',
             by = join_by(
                PERSON_ID,
                between(y$ORDER_DAY, x$JOIN_START, x$JOIN_END)
             ))# %>%
   #mutate(X = as.integer(ORDER_DAY - DAY))
#barplot(table(dx2$X))
dx %>% count(hasAST = !is.na(ORDER_DAY)) # only ~4K matched to AST, ~15K did not
dx <- dx %>% mutate(hasAST = !is.na(ORDER_DAY)) %>% select(-ORDER_DAY)

# join AbxAdmin
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_CLEANED_2017_2023_AbxAdmin.Rdata')
abxDF <- abxDF %>% 
   filter(PERSON_ID %in% unique(dx$PERSON_ID)) %>% # 9.7 mil --> 1 mil
   select(-START_DATE, -END_DATE) %>%
   distinct()

df <- dx %>%
   mutate(DAY = as.Date(substr(DX_FROM_DATE,1,10))) %>%
   mutate(JOIN_START = DAY, 
          JOIN_END = DAY + 14) %>%
   left_join(y = abxDF,
             by = join_by(
                PERSON_ID,
                JOIN_START <= START_DAY,
                JOIN_END >= START_DAY
             )) %>%
   #mutate(hasAbx = !is.na(ABX)) %>% select(PERSON_ID, DX_FROM_DATE, hasAbx) %>% distinct() %>% count(hasAbx) # 1/2 do not have antibiotics administered
   mutate(X = as.integer(START_DAY - DAY))
barplot(table(df$X))
barplot(table(df$X[df$ABX == 'FIDAXOMICIN']))

x <- df %>%
   select(PERSON_ID, DAY, ABX, hasAST) %>% 
   mutate(ABX = case_when(
      is.na(ABX) ~ NA,
      ABX == 'FIDAXOMICIN' ~ ABX,
      ABX == 'VANCOMYCIN' ~ ABX,
      .default = 'other'
   )) %>%
   distinct() %>%
   group_by(PERSON_ID, DAY, hasAST) %>%
   mutate(TRT = case_when(
      any(ABX == 'VANCOMYCIN') & !any(ABX == 'FIDAXOMICIN') ~ 'VAN',
      !any(ABX == 'VANCOMYCIN') & any(ABX == 'FIDAXOMICIN') ~ 'FID',
      any(ABX == 'VANCOMYCIN') & any(ABX == 'FIDAXOMICIN') ~ 'both',
      all(is.na(ABX)) ~ 'no abx',
      .default = 'other'
   )) %>%
   ungroup() %>%
   select(-ABX) %>%
   distinct()

x %>% filter(!hasAST) %>% count(TRT)
x %>% filter(hasAST) %>% count(TRT)

t <- x %>% count(TRT) %>% rename(Therapy = TRT)
t2 <- x %>% filter(!hasAST) %>% count(TRT) %>% rename(Therapy = TRT)
writeTable(t, 'Cdiff_VAN_FID')
writeTable(t2, 'Cdiff_VAN_FID_noBug')



# crude survival analysis
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_DEMO.Rdata')
df <- x %>% filter(TRT %in% c('VAN', 'FID'))

df <- df %>%
   left_join(y = dth,
             by = join_by(PERSON_ID)) %>% 
   filter(!is.na(DOB))
df %>% summarise(median(as.numeric(lubridate::as.duration(DAY - DOB)) / (365 * 24 * 60 * 60)), .by=TRT) # similar median ages
df %>% summarise(sum(GENDER == 'MALE') / n(), .by=TRT) # VAN = 59% female, FID = 71% female

censor_time <- 90
dfs <- df %>%
   mutate(AGE = as.integer(DAY - DOB)) %>%
   select(-DOB, -PATIENT_STATUS) %>%
   mutate(time = as.integer(DEATH_DATE - DAY)) %>%
   mutate(status = case_when(
      is.na(time) | time > censor_time ~ 0,
      .default = 1
   )) %>%
   mutate(time = case_when(
      is.na(time) | time > censor_time ~ censor_time,
      .default = time
   )) %>%
   select(-DEATH_DATE) %>%
   filter(time > 2)

sample_sizes <- table(dfs$TRT)
fit <- survfit(Surv(time, status) ~ TRT, data = dfs)
coxph(Surv(time, status) ~ TRT + GENDER + AGE, data = dfs)
pval <- survdiff(Surv(time, status) ~ TRT, data = dfs)$pvalue
d <- setNames(lapply(unique(summary(fit)$strata), function(x) which(summary(fit)$strata == x)),
              as.character(unique(summary(fit)$strata)))
s <- summary(fit)
xvals <- s$time
yvals <- s$surv
ste <- s$std.err
col_vec <- c('Fidaxomicin' = '#FF0000', 'Vancomycin' = '#0000FF')
ypos <- seq(0.575, 0.52, length.out=3)

plot(NA, ylim=c(0.5,1), xlim=c(0, censor_time), xaxt='n', yaxt='n', xaxs='i', yaxs='i', xpd=NA, xlab='Days', ylab='Surviving fraction')
title(main = 'C. diff\n ', font.main=4)
title(main = ' \nFidaxomicin vs. Vancomycin')
plotSurv(xvals, yvals, ste, d, col_vec)
axis(side=1, at=seq(0, censor_time, 15))
axis(side=2, las=1)
abline(v=30, lty=2)
text(x=2, y=ypos[1], adj=0, labels=paste0('p = ', prettyNum(pval, digits=3, scientific=-1)), font=ifelse(pval<0.05, 2, 1))
text(x=2, y=ypos[2:3], adj=0, col=col_vec,
     labels=paste0(names(col_vec), ' (n = ', prettyNum(sample_sizes, big.mark=','), ')'))
##### END #####



##### VAP - Pip-tazo vs. Cefepime vs. Meropenem ######
# ICD9: 997.31
# ICD10: J95.851
source('~/Desktop/EHR/EHR work/config_file.R')
dx <- tibble(dbGetQuery(conn, "select * from AMB_ETL.LAB_SENS_DX_VW where DX_CODE = 'J95.851'"))

load('~/Desktop/EHR/EHR work/RdataFiles/ALL_clean_ASTs.Rdata')
astDF <- astDF %>% filter(PERSON_ID %in% unique(dx$PERSON_ID)) # 14K
astDF <- astDF %>% select(PERSON_ID, ORDER_DAY, RESULT_DAY) %>% distinct() # 9K

load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_CLEANED_2017_2023_AbxAdmin.Rdata')
abxDF <- abxDF %>%
   filter(PERSON_ID %in% unique(dx$PERSON_ID)) %>% # 9.7 mil --> 220K
   select(-START_DATE, -END_DATE) %>%
   distinct()

# join with ASTs
df <- dx %>%
   mutate(DAY = as.Date(substr(DX_FROM_DATE, 1, 10))) %>%
   mutate(JOIN_START = DAY - 1,
          JOIN_END = DAY + 14) %>%
   left_join(y = astDF,
             multiple = 'first',
             by = join_by(
                PERSON_ID,
                between(y$ORDER_DAY, x$JOIN_START, x$JOIN_END)
             )) %>%
   mutate(X = as.integer(ORDER_DAY - DAY)) %>%
   select(-JOIN_START, -JOIN_END)
df %>% count(is.na(ORDER_DAY)) # 409 are missing ASTs
barplot(table(df$X))

# join with AbxAdmin
df$JOIN_START <- df$DAY
w <- which((df$ORDER_DAY - 2) < df$JOIN_START)
df$JOIN_START[w] <- df$ORDER_DAY[w] - 2
df$JOIN_END <- df$DAY
w <- which((df$RESULT_DAY + 7) > df$JOIN_END)
df$JOIN_END[w] <- df$RESULT_DAY[w] + 7
rm(w)

df <- df %>%
   left_join(y = abxDF,
             by = join_by(
                PERSON_ID,
                JOIN_START <= START_DAY,
                JOIN_END >= START_DAY
             )) %>%
   mutate(X = as.integer(START_DAY - DAY)) %>%
   mutate(WINDOW = case_when(
      START_DAY >= ORDER_DAY - 2 & START_DAY <= ORDER_DAY + 1 ~ 'EMPIRIC',
      START_DAY >= RESULT_DAY ~ 'DEFINITIVE'
   ))
barplot(table(df$X))

df <- df %>%
   select(PERSON_ID, DAY, ABX, START_DAY, WINDOW) %>%
   distinct()

df %>% count(ABX, sort=TRUE)

df <- df %>%
   mutate(TRT = case_when(
      ABX == 'PIPERACILLIN/TAZOBACTAM' ~ 'TZP',
      ABX == 'CEFEPIME' ~ 'FEP',
      ABX == 'MEROPENEM' ~ 'MER',
      is.na(ABX) ~ 'NoAbx',
      .default = 'otherAbx'
   )) %>%
   select(PERSON_ID, DAY, TRT, WINDOW) %>%
   distinct()

df %>% select(PERSON_ID, DAY, TRT) %>% distinct() %>% count(TRT)

df <- df %>%
   group_by(PERSON_ID, DAY, WINDOW) %>%
   mutate(TRT2 = case_when(
      any(TRT == 'TZP') & !any(TRT == 'FEP') & !any(TRT == 'MER') ~ 'TZPo',
      !any(TRT == 'TZP') &  any(TRT == 'FEP') & !any(TRT == 'MER') ~ 'FEPo',
      !any(TRT == 'TZP') & !any(TRT == 'FEP') &  any(TRT == 'MER') ~ 'MERo',
      
      any(TRT == 'TZP') &  any(TRT == 'FEP') & !any(TRT == 'MER') ~ 'TZP+FEP',
      !any(TRT == 'TZP') &  any(TRT == 'FEP') &  any(TRT == 'MER') ~ 'FEP+MER',
      any(TRT == 'TZP') & !any(TRT == 'FEP') &  any(TRT == 'MER') ~ 'MER+TZP',
      
      any(TRT == 'TZP') &  any(TRT == 'FEP') &  any(TRT == 'MER') ~ 'TZP+FEP+MER',
      
      !any(TRT == 'TZP') & !any(TRT == 'FEP') & !any(TRT == 'MER') & any(TRT == 'otherAbx') ~ 'otherAbx',
      .default = 'noAbx'
   )) %>%
   ungroup() %>%
   select(PERSON_ID, DAY, TRT2, WINDOW) %>%
   distinct()

e <- df %>% filter(WINDOW == 'EMPIRIC') %>% count(TRT2) %>% slice(c(7,2,4, 5,1,3,6, 8)) %>% rename(EmpiricTherapy = TRT2)
d <- df %>% filter(WINDOW == 'DEFINITIVE') %>% count(TRT2) %>% slice(c(7,2,4, 5,1,3,6, 8)) %>% rename(DefinitiveTherapy = TRT2)
t <- df %>% filter(!is.na(WINDOW)) %>% 
   tidyr::pivot_wider(values_from = TRT2,
                      names_from = WINDOW,
                      values_fill = 'noAbx') %>%
   count(EMPIRIC, DEFINITIVE, sort=TRUE) %>%
   slice(1:10) %>%
   rename(EmpiricTherapy = EMPIRIC,
          DefinitiveTherapy = DEFINITIVE)

writeTable(e, 'VAP_TZP_FEP_MER_empiric')
writeTable(d, 'VAP_TZP_FEP_MER_definitive')
writeTable(t, 'VAP_TZP_FEP_MER_transitions')



# crude survival analysis
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_DEMO.Rdata')
df <- df %>% 
   filter(WINDOW == 'DEFINITIVE') %>% 
   filter(TRT2 %in% c('TZPo', 'FEPo', 'MERo')) %>% 
   rename(TRT = TRT2) %>% 
   select(-WINDOW)

df <- df %>%
   left_join(y = dth,
             by = join_by(PERSON_ID)) %>% 
   filter(!is.na(DOB))
df %>% summarise(median(as.numeric(lubridate::as.duration(DAY - DOB)) / (365 * 24 * 60 * 60)), .by=TRT) # similar median ages
df %>% summarise(sum(GENDER == 'MALE') / n(), .by=TRT) # similar

censor_time <- 180
dfs <- df %>%
   mutate(AGE = as.integer(DAY - DOB)) %>%
   select(-DOB, -PATIENT_STATUS) %>%
   mutate(time = as.integer(DEATH_DATE - DAY)) %>%
   mutate(status = case_when(
      is.na(time) | time > censor_time ~ 0,
      .default = 1
   )) %>%
   mutate(time = case_when(
      is.na(time) | time > censor_time ~ censor_time,
      .default = time
   )) %>%
   select(-DEATH_DATE) %>%
   filter(time > 2)

library(survival)
sample_sizes <- table(dfs$TRT)
fit <- survfit(Surv(time, status) ~ TRT, data = dfs)
coxph(Surv(time, status) ~ TRT, data = dfs)
pval <- c(0.0223, 0.407, 0.0765)
d <- setNames(lapply(unique(summary(fit)$strata), function(x) which(summary(fit)$strata == x)),
              as.character(unique(summary(fit)$strata)))
s <- summary(fit)
xvals <- s$time
yvals <- s$surv
ste <- s$std.err
col_vec <- c('Cefepime' = '#FF0000', 'Meropenem' = '#0000FF', 'Pip-tazo' = '#888888')
ypos <- seq(0.24, 0.02, length.out=6)

plot(NA, ylim=c(0,1), xlim=c(0, censor_time), xaxt='n', yaxt='n', xaxs='i', yaxs='i', xpd=NA, xlab='Days', ylab='Surviving fraction')
title(main = 'VAP - Cefepime vs. Meropenem vs. Pip-tazo')
plotSurv(xvals, yvals, ste, d, col_vec)
axis(side=1, at=seq(0, censor_time, 15))
axis(side=2, las=1)
abline(v=30, lty=2)
text(x=2, y=ypos[1:3], adj=0, labels=paste0('p = ', prettyNum(pval, digits=3, scientific=-1), c(' (FEP vs. MER)', ' (FEP vs. TZP)', ' (TZP vs. MER)')), font=ifelse(pval<0.05, 2, 1))
text(x=2, y=ypos[4:6], adj=0, col=col_vec,
     labels=paste0(names(col_vec), ' (n = ', prettyNum(sample_sizes, big.mark=','), ')'))
#### END #####

























