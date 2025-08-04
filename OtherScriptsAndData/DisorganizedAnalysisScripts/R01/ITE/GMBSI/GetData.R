#### LOAD AND PREPROCESS DATA ####
library(dplyr)
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_clean_ASTs.Rdata')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_CLEANED_2015_2024_AbxAdmin.Rdata')
load(file = '~/Desktop/EHR-mining/UsefulDataForCleaning/bug_groups.Rdata')
source('~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/CreateNamedAbxAbbreviations.R')
Rdata_file_path <- '~/Desktop/EHR/EHR work/RdataFiles/R01/'


# Filter ASTs to 
#     blood cultures with gram-negatives
#     index culture by at least 42 days
astDF <- astDF %>%
   filter(BLOOD) %>%
   group_by(PERSON_ID) %>%
   mutate(TIME_SINCE_PRV_BSI = as.integer(ORDER_DAY - lag(ORDER_DAY))) %>%
   filter(is.na(TIME_SINCE_PRV_BSI) | TIME_SINCE_PRV_BSI >= 42L) %>%
   ungroup() %>%
   filter(
      BUG %in% unlist(bug_groups[c('Enterobacterales', 'NonFermGN')]),
      lubridate::year(ORDER_DAY) %in% 2017:2023
   )
rm(bug_groups)
save(astDF, file = paste0(Rdata_file_path, 'astDF_raw_gramnegBSIs.Rdata'))

astDF %>% summarise(sum(ESBL) / n()) # only 10.6% are ESBL, will probably exclude from this analysis


# encounters
source('~/Desktop/EHR-mining/Scripts/CleaningScripts/getRawEncountersData.R')
encsDF <- getRawEncountersData(person_ids=unique(astDF$PERSON_ID), years=2015:2024)
save(encsDF, file = paste0(Rdata_file_path, 'encsDF_raw_gramnegBSIs.Rdata'))

# comorbidities
source('~/Desktop/EHR-mining/Scripts/CleaningScripts/GetRawICDsData.R')
icdsDF <- getRawICDsData(person_ids=unique(astDF$PERSON_ID), years=2015:2024)
save(icdsDF, file = paste0(Rdata_file_path, 'icdsDF_raw_gramnegBSIs.Rdata'))

rm(conn, result, encsDF, icdsDF, getRawICDsData, getRawEncountersData)


# Filter ABXs to those patients
abxDF <- abxDF %>%
   filter(
      PERSON_ID %in% unique(astDF$PERSON_ID),
      lubridate::year(START_DAY) %in% 2017:2023
   ) %>%
   filter(!ABX %in% c('CASPOFUNGIN', 'FLUCONAZOLE', 'METRONIDAZOLE', 'VORICONAZOLE', 'POSACONAZOLE', 'CLOTRIMAZOLE', 
                      'MICAFUNGIN', 'KETOCONAZOLE', 'AMPHOTERICIN', 'TERBINAFINE', 'ITRACONAZOLE', 'FLUCYTOSINE')) %>%
   select(-END_DATE) %>%
   distinct()
w <- which(abxDF$ABX == 'PENICILLIN G')
abxDF$ABX[w] <- 'BENZYLPENICILLIN'
save(abxDF, file = paste0(Rdata_file_path, 'abxDF_raw_gramnegBSIs.Rdata'))
rm(w)
#### END LOAD/PREPROCESS ####



library(dplyr)
Rdata_file_path <- '~/Desktop/EHR/EHR work/RdataFiles/R01/'
load(file = paste0(Rdata_file_path, 'astDF_raw_gramnegBSIs.Rdata'))
load(file = paste0(Rdata_file_path, 'encsDF_raw_gramnegBSIs.Rdata'))
load(file = paste0(Rdata_file_path, 'abxDF_raw_gramnegBSIs.Rdata'))
load(file = paste0(Rdata_file_path, 'icdsDF_raw_gramnegBSIs.Rdata'))
load(file = '~/Desktop/EHR-mining/UsefulDataForCleaning/bug_groups.Rdata')
load(file = '~/Desktop/EHR-mining/UsefulDataForCleaning/UPMC_site_names.Rdata')
source('~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/CreateNamedAbxAbbreviations.R')



#### JOIN ENCOUNTER DATA ####
source('~/Desktop/EHR-mining/Scripts/CleaningScripts/featurizeEncounters.R')
df <- featurizeEncounters(df = astDF %>% select(PERSON_ID, ORDER_DAY, ORDER_DATE, BUG, ESBL) %>% distinct(),
                          encsDF = encsDF,
                          preprocess = FALSE,
                          rename = FALSE)
rm(featurizeEncounters, getFacilityNames, site_names)
df_og <- df


# 9% are missing encounters data !!
df %>% summarise(sum(MISSING_ENCOUNTER_DATA) / n() * 100)
df <- df %>% filter(!MISSING_ENCOUNTER_DATA) %>% select(-MISSING_ENCOUNTER_DATA)
df <- df %>% filter(FACILITY1 != '' & FACILITY2 != '')
rm(encsDF)
#### END ENCOUNTERS ####


#### JOIN DEMOGRAPHCS ####
source('~/Desktop/EHR-mining/Scripts/CleaningScripts/featurizeDemographics.R')
df <- getDemographics(df, preprocess = FALSE)
df <- df %>% 
   filter(is.na(mortality_time) | mortality_time >= 4L) %>%
   relocate(mortality_time, AGE, FEMALE, .before=HOSP_ACQ)
rm(getDemographics)
#### END ####


#### Join ASTs + ABXs + define treatment groups ####
#     only empiric CRO + TZP administrations
# df <- df_og
df <- df %>%
   mutate(
      JOIN_START = ORDER_DATE - 2 * 86400,
      JOIN_END = ORDER_DATE + 4 * 86400
   ) %>%
   left_join(
      x = .,
      y = abxDF,
      by = join_by(
         PERSON_ID,
         between(y$START_DATE, x$JOIN_START, x$JOIN_END)
      )
   ) %>%
   mutate(
      ABX_TIME = as.numeric(lubridate::as.duration(START_DATE - ORDER_DATE)) / 86400,
      ABX_DAY = as.integer(START_DAY - ORDER_DAY)
   ) %>%
   select(!c(JOIN_START, JOIN_END, START_DATE, START_DAY)) %>%
   relocate(ABX, ABX_TIME, ABX_DAY, .after=ESBL)
df_abx <- df
# df <- df_abx



# how many did not receive empiric therapy?
pcntNoAbxInWindow <- function(end_time_window) {
   df %>%
      mutate(empABX = !is.na(ABX) & ABX_TIME < end_time_window) %>%
      group_by(PERSON_ID, ORDER_DAY) %>%
      mutate(empABX = any(empABX)) %>%
      ungroup() %>%
      select(PERSON_ID, ORDER_DAY, empABX) %>%
      distinct() %>%
      summarise(sum(!empABX) / n())
}
pcntNoAbxInWindow(0)   # 87%
pcntNoAbxInWindow(0.5) # 13.5%
pcntNoAbxInWindow(1)   # 6.5%
pcntNoAbxInWindow(2)   # 3.7%
pcntNoAbxInWindow(4)   # 3.2%
rm(pcntNoAbxInWindow)


# keep only patients with emiric AND later ABX
# how many patients are left?
df %>% group_by(PERSON_ID, ORDER_DAY) %>% n_groups()                            # 22,745 - total gram-negative blood cultures
df %>% filter(!is.na(ABX)) %>% group_by(PERSON_ID, ORDER_DAY) %>% n_groups()    # 20,024 - received abx within 4 days
df %>% filter(ABX_TIME < 0.5) %>% group_by(PERSON_ID, ORDER_DAY) %>% n_groups() # 19,683 - received empiric ABX
df %>% group_by(PERSON_ID, ORDER_DAY) %>%                                       # 18,883 - received empiric ABX and ABX ≥ 48h
   filter(any(ABX_TIME < 0.5), any(ABX_TIME > 2)) %>% n_groups()
df %>%                                                                          # 13,185 - received empiric TZP/CRO
   filter(ABX_TIME < 0.5 & ABX %in% c('PIPERACILLIN/TAZOBACTAM', 'CEFTRIAXONE')) %>%
   group_by(PERSON_ID, ORDER_DAY) %>% n_groups()
df %>% group_by(PERSON_ID, ORDER_DAY) %>%                                       # 12,636 - received empiric TZP/CRO + still receiving ABX 48 hours later
   filter(any(ABX_TIME < 0.5 & ABX %in% c('PIPERACILLIN/TAZOBACTAM', 'CEFTRIAXONE')),
          any(ABX_TIME > 2)) %>%
   n_groups()
df %>% group_by(PERSON_ID, ORDER_DAY) %>%                                       # 3,414 - received empiric TZP/CRO + still receiving ABX 48 hours later (but not TZP or CRO)
   filter(any(ABX_TIME < 0.5 & ABX %in% c('PIPERACILLIN/TAZOBACTAM', 'CEFTRIAXONE')),
          any(ABX_TIME > 2)) %>%
   filter(!any(ABX_TIME > 2 & ABX %in% c('PIPERACILLIN/TAZOBACTAM', 'CEFTRIAXONE'))) %>% 
   n_groups()


# filter down to empiric TZP or CRO + still receiving ABX ≥ 48h later
df <- df %>%
   filter(!is.na(ABX)) %>%
   mutate(
      empTC = ABX_TIME < 0.5 & ABX %in% c('PIPERACILLIN/TAZOBACTAM', 'CEFTRIAXONE'),
      defABX = ABX_TIME > 2
   ) %>%
   group_by(PERSON_ID, ORDER_DAY) %>%
   filter(any(empTC), any(defABX)) %>%
   ungroup() %>%
   select(-empTC, -defABX)

df %>%
   filter(ABX_TIME > 2) %>%
   group_by(PERSON_ID, ORDER_DAY) %>% filter(!any(ABX %in% c('PIPERACILLIN/TAZOBACTAM', 'CEFTRIAXONE'))) %>% #n_groups() # 3,414 - discontinued TZP or CRO
   ungroup() %>%
   select(PERSON_ID, ORDER_DAY, ABX) %>%
   distinct() %>%
   count(ABX, sort=TRUE) # FEP, MEM, CFZ, VAN, ETP, CIP, SXT, SAM, CXM, LVX



# how many days / doses of abx after empiric?
# get variable - how many distinct abx each infection was treated with empirically
df <- df %>%
   filter(ABX_TIME < 0.5) %>%
   group_by(PERSON_ID, ORDER_DAY) %>%
   mutate(
      numAbx = length(unique(ABX)),
      firstAbxTime = min(ABX_TIME[ABX %in% c('PIPERACILLIN/TAZOBACTAM', 'CEFTRIAXONE')]),
      firstAnyAbxTime = min(ABX_TIME)
   ) %>%
   ungroup() %>%
   filter(ABX %in% c('PIPERACILLIN/TAZOBACTAM', 'CEFTRIAXONE')) %>%
   mutate(TRT = unname(abbr[ABX])) %>%
   select(-ABX_TIME, -ABX_DAY, -ABX) %>%
   distinct() %>%
   group_by(PERSON_ID, ORDER_DAY) %>%
   mutate(TRT = paste(sort(TRT), collapse=',')) %>%
   ungroup() %>%
   distinct() %>%
   relocate(TRT, numAbx, firstAbxTime, firstAnyAbxTime, .after=BUG) %>%
   mutate(
      x = ifelse(grepl(',', TRT), 2L, 1L),
      numAbx = numAbx - x
   ) %>%
   select(-x)

# how many received both?
df %>% count(TRT)

# is TZP or CRO generally the first abx given? yes.
df %>% mutate(d = round((firstAbxTime - firstAnyAbxTime) * 24)) %>% count(d > 4)


# how many non-empiric abx days for these patients?
df <- left_join(
   x = df,
   y = df_abx %>% 
      filter(ABX_TIME >= 1) %>% 
      summarise(
         n_ABX_DAY = length(unique(ABX_DAY)),
         .by=c(PERSON_ID, ORDER_DAY)
      ),
   by = join_by(PERSON_ID, ORDER_DAY)
) %>%
   relocate(n_ABX_DAY, .after=TRT)

# some don't have any! maybe they died?
df %>% count(n_ABX_DAY)

df <- df %>% filter(n_ABX_DAY > 1L)

# how many received both?
df %>% count(TRT)

# narrow down to those that received either one or the other empirically, but not both
df <- df %>% filter(!grepl(',', TRT))


# add back other empiric abx...
df <- left_join(
   x = df,
   y = df_abx %>% 
      filter(ABX_TIME < 0.5) %>% 
      select(PERSON_ID, ORDER_DAY, ABX) %>% 
      distinct() %>% 
      mutate(ABX = unname(abbr[ABX])) %>% 
      filter(ABX %in% c('VAN', 'FEP', 'CFZ', 'SAM', 'MEM', 'ETP', 'AZM', 'GEN', 'LZD', 'SXT', 'DAP')) %>%
      mutate(ABX = paste0('EMP_', ABX)),
   by = join_by(PERSON_ID, ORDER_DAY)
) %>%
   mutate(X = 1L) %>%
   tidyr::pivot_wider(
      names_from = ABX,
      values_from = X,
      values_fill = 0L
   ) %>%
   select(-`NA`) %>%
   relocate(contains('EMP_'), .before=ESBL)
#### END JOIN ASTs + ABXs ####


#### JOIN DX CODES ####
icdsDF <- icdsDF %>% filter(PERSON_ID %in% unique(df$PERSON_ID))
source('~/Desktop/EHR-mining/Scripts/CleaningScripts/featurizeICDs.R')
df <- featurizeICDs(df=df, icdsDF=icdsDF)
rm(featurizeICDs, icdsDF)
#### END DX CODES ####


#### JOIN LAB VALUES ####
labsDF <- labsDF %>% filter(PERSON_ID %in% unique(df$PERSON_ID))
source('~/Desktop/EHR-mining/Scripts/CleaningScripts/featurizeLabs.R')
df <- featurizeLabs(df=df, labsDF=labsDF)
#### END LAB VALUES ####


df$BUG <- gsub('^([A-Z])[a-z]+ ([a-z]+)$', '\\1_\\2', df$BUG)
df$FACILITY <- gsub('[^A-z]', '', df$FACILITY1)
df <- df %>% select(-FACILITY1, -FACILITY2)

# one-hot encoding
for (var in c('BUG', 'FACILITY')) {
   vals <- unique(df[[var]])
   for (v in vals) {
      df[[v]] <- df[[var]] == v
   }
   df <- df[names(df) != var]
}

# combine empiric carbapenems, ceph
df$EMP_CEF <- as.integer(df$EMP_CFZ | df$EMP_FEP); df <- df %>% select(-EMP_CFZ, -EMP_FEP)
df$EMP_CAR <- as.integer(df$EMP_MEM | df$EMP_ETP); df <- df %>% select(-EMP_MEM, -EMP_ETP)
df$EMP_GRP <- as.integer(df$EMP_LZD | df$EMP_DAP); df <- df %>% select(-EMP_DAP, -EMP_LZD)
df <- df %>% select(-EMP_SAM)
save(df, file=paste0(Rdata_file_path, 'finalDF_gramneg.Rdata'))



#### PLOT STUFF


# plotting variables
Rdata_file_path <- '~/Desktop/EHR/EHR work/RdataFiles/R01/'
plot_path <- '~/Desktop/EHR-mining/Scripts/AnalysisScripts/R01/ITE/GMBSI/Plots/'
col_vec <- c('TZP' = '#3399FF', 'CRO' = '#FFaa66')
pmax <- 80
ptcex <- 1.3
lwd <- 2
leg <- function(pch, lwd=NA, pos='topright', pt.cex=ptcex) {
   legend(x=pos, col=col_vec, pch=pch, lwd=lwd, legend=names(col_vec), xpd=NA, pt.cex=ptcex)
}

library(dplyr)
load(file = paste0(Rdata_file_path, 'finalDF_gramneg.Rdata'))




df %>% summarise(
   n = n(),
   num_abx_days = mean(n_ABX_DAY), 
   esbl_rate = sum(ESBL) / n, 
   abxBeforeCulture = sum(firstAbxTime < -0/24) / n,
   anyAbxBeforeCulture = sum(firstAnyAbxTime < 0) / n,
   firstAbxTime_hours = mean(firstAbxTime) * 24, 
   firstAnyAbxTime_hours = mean(firstAnyAbxTime) * 24,
   .by=TRT
)



#### PLOT 001ab - Empiric prescriptions ####
{
   # TZP group receives: 
   #     VAN way more often than CRO group !!!
   #     greater number of distinct ABX empirically
   t <- df %>%
      summarise(
         across(contains('EMP_'), ~ sum(.) / n() * 100),
         .by = TRT
      ) %>%
      as.data.frame()
   row.names(t) <- t$TRT
   t <- t %>% select(-TRT)
   t <- t[c('CRO','TZP'), ]
   t <- as.matrix(t, rownames.force=TRUE)
   
   
   pdf(file = paste0(plot_path, '001a_other_empiric_abx.pdf'), width=6.5, height=9)
   layout(matrix(1:2, nrow=2), heights=c(1.5,1))
   par(cex=1, mar=c(3,4,1,1), mgp=c(1.5, 0.25, 0), tck=-0.015, cex.main=1)
   
   b <- barplot(t, horiz=TRUE, beside=TRUE, names.arg=rep(NA, ncol(t)), xlab='% of infections treated with X', col=rev(col_vec), xlim=c(0, pmax))
   leg(pch=15)
   text(x=-1, y=apply(b, 2, mean), labels=gsub('EMP_', '', colnames(t)), adj=1, xpd=NA)
   title(main='Other empiric abx', line=-1)
   text(x=t+1, y=b, adj=0, labels=paste0(round(t, 1), '%'))
   
   par(mar=c(2.5,3,2.5,1))
   plot(NA, xlim=c(0,8), ylim=c(0, pmax), xlab='Num other ABX', xpd=NA, ylab='% of infections', yaxt='n')
   title(main='Number of distinct other ABX administered', line=1.2, xpd=NA)
   axis(side=2, las=1)
   leg(pch=16, lwd=lwd)
   plotLinesPoints <- function(x, abx) {
      col <- col_vec[abx]
      y <- x[df$TRT == abx]
      abline(v=mean(y), lty=2, col=col, lwd=lwd)
      if (abx == 'TZP') text(x=mean(y)-0.2, y=86, labels='mean', xpd=NA)
      y <- table(y)
      y <- y / sum(y) * 100
      x <- as.integer(names(y))
      lines(x, y, col=col, lwd=lwd)
      points(x, y, col=col, pch=16, cex=ptcex)
   }
   plotLinesPoints(df$numAbx, 'TZP')
   plotLinesPoints(df$numAbx, 'CRO')
   dev.off()
   
   
   # timing of first abx Rx
   pdf(file = paste0(plot_path, '001b_abx_timing.pdf'), width=5.5)
   # par(mfrow=c(2,1), mar=c(2,3,2,1), mgp=c(1, 0.2, 0))
   par(mfrow=c(2,1), cex=1, mar=c(3,2.5,2,1), mgp=c(1.5, 0.25, 0), tck=-0.015, cex.main=1)
   abxTimeDens <- function(x='firstAbxTime', abx) {
      if (abx == 'TZP') {
         plot(NA, xlim=c(-2, 0.5), ylim=c(0, 7), xlab='Days since blood culture order', ylab='', yaxt='n')
         title(main=ifelse(x == 'firstAbxTime', 'Timing of first TZP or CRO admin', 'Timing of first any ABX admin'),
               line=0.3)
         title(ylab='Density', line=0.8)
      }
      if (abx == 'CRO') {
         leg(pch=NA, lwd=lwd, pos='topleft')
      }
      col <- col_vec[abx]
      x <- df[[x]][df$TRT == abx]
      abline(v=mean(x), lty=2, col=col, lwd=lwd)
      lines(density(x), col=col, lwd=lwd)
   }
   abxTimeDens(abx='TZP')
   abxTimeDens(abx='CRO')
   abxTimeDens(x='firstAnyAbxTime', abx='TZP')
   abxTimeDens(x='firstAnyAbxTime', abx='CRO')
   
   dev.off()
   rm(t, b, plotLinesPoints, abxTimeDens)
}
#### END PLOT 001ab ####


#### PLOT 002 - Encounter, bug, facility, age differences ####
{
   t <- df %>%
      summarise(
         across(.cols = c(ESBL, FEMALE, HOSP_ACQ, NURSING_HOME_HOSPICE, EMERGENCY_DEPT, TRANSFER),
                .fns = ~ sum(.) / n() * 100),
         .by = TRT
      ) %>%
      as.data.frame()
   row.names(t) <- t$TRT
   t <- t %>% select(-TRT)
   t <- t[c('CRO','TZP'), ]
   t <- as.matrix(t, rownames.force=TRUE)
   colnames(t) <- stringr::str_to_sentence(gsub('_', ' ', colnames(t)))
   colnames(t)[colnames(t) == 'Nursing home hospice'] <- 'Nursing\nhome/hospice'
   
   pmax <- 100
   
   pdf(file = paste0(plot_path, '002_encounter_bug_facility_info.pdf'), height=7.5)
   layout(matrix(c(1,2,1,3), nrow=2), heights=c(1,1.3))
   par(cex=1, mar=c(3,7,1.25,3), mgp=c(1.5, 0.25, 0), tck=-0.015, cex.main=1)
   
   
   # Enouncter info
   b <- barplot(t, horiz=TRUE, beside=TRUE, names.arg=rep(NA, ncol(t)), xlab='% of infections with X', col=rev(col_vec), xlim=c(0, pmax))
   title(main='Encounter info', line=0)
   leg(pch=15, pos='bottomright', pt.cex=ptcex+0.2)
   text(x=-26, y=apply(b, 2, mean), labels=colnames(t), adj=0, xpd=NA)
   text(x=t+1, y=b, adj=0, labels=paste0(round(t, 1), '%'), xpd=NA)
   
   
   # bug and facility distributions
   getTable <- function(col) {
      tab <- inner_join(
         x = df %>%
            filter(TRT == 'TZP') %>%
            count(get(col), sort=TRUE) %>%
            rename(tzp = n),
         y = df %>%
            filter(TRT == 'CRO') %>%
            count(get(col), sort=TRUE) %>%
            rename(cro = n),
         by = join_by(`get(col)`)
      ) %>%
         mutate(
            n = tzp + cro,
            tzp = tzp / n * 100,
            cro = cro / n * 100
         ) %>%
         arrange(desc(tzp))
      
      h <- t(as.matrix(tab %>% select(tzp, cro)))
      b <- barplot(height=h, horiz=TRUE, col=col_vec, xlab='% of (bug) infections')
      if (col == 'BUG') {
         tab[[1]] <- gsub('^([A-Z])[a-z]+ ([a-z]+)', '\\1. \\2', tab[[1]])
      }
      text(x=-3, y=b, labels=tab[[1]], adj=1, xpd=NA, font=ifelse(col == 'BUG', 3, 1))
      text(x=103, y=c(b, b[length(b)]+1), labels=c(tab$n, 'N'), adj=0, xpd=NA)
      title(main=stringr::str_to_sentence(gsub('[0-9]', '', col)))
   }
   
   getTable('BUG')
   getTable('FACILITY1')
   # getTable('FACILITY2')
   
   dev.off()
   rm(t,b,getTable)
   
   
   # # AGE 
   # pdf(file = paste0(plot_path, '002b_age.pdf'), height=3.5, width=3.5)
   # par(cex=1, mar=c(3,2.5,1.25,1), mgp=c(1.5, 0.25, 0), tck=-0.015, cex.main=1)
   # ageDens <- function(x='AGE', abx, setplot=FALSE, addlegend=FALSE) {
   #    if (setplot) {
   #       plot(NA, xlim=c(0,100), ylim=c(0, 0.05), xlab='Days since blood culture order', ylab='', yaxt='n')
   #       title(main='Age', line=0.3)
   #       title(ylab='Density', line=0.8)
   #    }
   #    if (addlegend) {
   #       leg(pch=NA, lwd=lwd, pos='topleft', pt.cex=NA)
   #    }
   #    col <- col_vec[abx]
   #    x <- df[[x]][df$TRT == abx]
   #    abline(v=mean(x), lty=2, col=col, lwd=lwd)
   #    lines(density(x), col=col, lwd=lwd)
   # }
   # ageDens(abx='TZP', setplot=TRUE)
   # ageDens(abx='CRO', addlegend=TRUE)
   # dev.off()
   # rm(ageDens)
}
#### END PLOT 002 ####


#### PLOT 003 - Survival ####
{
   library(survival)
   source('~/Desktop/EHR-mining/Scripts/AnalysisScripts/KaplanMeierCurveFxns.R')
   
   pdf(file = paste0(plot_path, '003_survival_curve.pdf'), height=5, width=5)
   par(cex=1, mar=c(4,4,2.5,1), mgp=c(1.5, 0.4, 0), tck=-0.015, cex.main=1)
   plotKP(df = df %>% rename(time = mortality_time),
          strat = 'TRT',
          addstatus = TRUE,
          cex = ptcex,
          cohort = 'Gram-negative bacteremia',
          col_vec = col_vec)
   # plotKP(df = df %>% rename(time = mortality_time),
   #        strat = 'ESBL',
   #        addstatus = TRUE,
   #        cex = ptcex,
   #        cohort = 'Gram-negative bacteremia',
   #        col_vec = c('0' = "#3399FF", '1' = "#FFaa66"))
   dev.off()
}
#### END PLOT 003 ####


#### PLOT 004 - var distribution ####
source('~/Desktop/EHR-mining/Scripts/AnalysisScripts/DrawTableOneFxn.R')
source('~/Desktop/EHR-mining/Scripts/AnalysisScripts/ProcessTableOneFxn.R')

# draw table 1
vars <- df %>% select(!c(PERSON_ID:TRT, mortality_time)) %>% names()
t1 <- tableone::CreateTableOne(vars=vars, strata='TRT', smd=TRUE,
                               data=df %>% 
                                  mutate(across(.cols = c(matches('^(EMP|ICD)_'), FEMALE, HOSP_ACQ, NURSING_HOME_HOSPICE, EMERGENCY_DEPT, MATERNAL, TRANSFER), 
                                                .fns = as.logical)))
t1 <- processTableOne(t1)

pdf(file = paste0(plot_path, '004a_covar_distributions.pdf'), width=9, height=5.5)
par(mfrow=c(1,3))
drawTableOne(t1[grep('^ICD_', rownames(t1)),], trt=names(col_vec), mar_left = 12, addlegend=TRUE, xpos=1.1)
drawTableOne(t1[grep('^(LAB|EMP)_', rownames(t1)),], trt=names(col_vec), mar_left = 12, xpos=1.1)
drawTableOne(t1[-grep('^(ICD|LAB|EMP)_', rownames(t1)),], trt=names(col_vec), mar_left = 12, xpos=1.1)
dev.off()


# propensity score overlap
prop_model <- df %>%
   mutate(TZP = TRT == 'TZP') %>%
   select(TZP, vars) %>%
   glm(TZP ~ ., family=binomial(), data=.)
prop_scores <- predict(prop_model, df[vars], type='response')

pdf(file = paste0(plot_path, '004b_propensity_score_overlap.pdf'), height=5, width=5)
par(mgp=c(1.75, 0.4, 0), tck=-0.01, mar=c(3,3,2,1))
plot(NA, xlim=c(-0.05, 1.05), ylim=c(0, 3.25), xlab='Propensity score', ylab='Density')
lines(density(prop_scores[df$TRT == 'TZP']), col=col_vec['TZP'], lwd=2)
lines(density(prop_scores[df$TRT == 'CRO']), col=col_vec['CRO'], lwd=2)
legend('topright', legend=names(col_vec), col=col_vec, lwd=2)
dev.off()
#### END PLOT 004 ####


#### PLOT 005 - ATE ####
source('~/Desktop/EHR-mining/Scripts/AnalysisScripts/CoxModelFxns.R')

coefs <- modelCoxIPTWcoefs(
   df = df %>%
      select(TRT, mortality_time, BUG, AGE, FEMALE, ESBL, HOSP_ACQ:TRANSFER, FACILITY1, contains('ICD'), firstAnyAbxTime, EMP_VAN, EMP_MEM) %>%
      mutate(TZP = TRT == 'TZP') %>%
      select(-MATERNAL, -TRT),
   time_var = 'mortality_time',
   status_var = NULL,
   treatment_var = 'TZP',
   censor_time = 30L,
   conf_vars = NULL,
   print_conf_vars = FALSE,
   normalize_weights = TRUE
)

pdf(file = paste0(plot_path, '005_cox_HR.pdf'), height=5, width=5)
plotCoxCoefs(coefs)
dev.off()

rm(time_var, status_var, treatment_var, censor_time, conf_vars, print_conf_vars, normalize_weights)
rm(mult, uni, mult_ATE, mult_ATO, mult_SW, ubugs, b, ufac, f, var_order)
#### END PLOT 005 ####












