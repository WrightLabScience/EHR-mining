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
      BUG == 'Staphylococcus aureus',
      OXACILLIN == 1L,
      lubridate::year(ORDER_DAY) %in% 2017:2023
   )
rm(bug_groups)
save(astDF, file = paste0(Rdata_file_path, 'astDF_raw_MRSAbact.Rdata'))


# # encounters
# source('~/Desktop/EHR-mining/Scripts/CleaningScripts/getRawEncountersData.R')
# encsDF <- getRawEncountersData(person_ids=unique(astDF$PERSON_ID), years=2015:2024)
# save(encsDF, file = paste0(Rdata_file_path, 'encsDF_raw_gramnegBSIs.Rdata'))
# 
# # comorbidities
# source('~/Desktop/EHR-mining/Scripts/CleaningScripts/GetRawICDsData.R')
# icdsDF <- getRawICDsData(person_ids=unique(astDF$PERSON_ID), years=2015:2024)
# save(icdsDF, file = paste0(Rdata_file_path, 'icdsDF_raw_gramnegBSIs.Rdata'))
# 
# rm(conn, result, encsDF, icdsDF, getRawICDsData, getRawEncountersData)


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
save(abxDF, file = paste0(Rdata_file_path, 'abxDF_raw_MRSAbact.Rdata'))
rm(w)
#### END LOAD/PREPROCESS ####


library(dplyr)
Rdata_file_path <- '~/Desktop/EHR/EHR work/RdataFiles/R01/'
load(file = paste0(Rdata_file_path, 'OutcomePredictionModel/encsDF_all_raw.Rdata'))
load(file = paste0(Rdata_file_path, 'OutcomePredictionModel/labsDF_all_cleaned.Rdata'))
load(file = paste0(Rdata_file_path, 'OutcomePredictionModel/icdsDF_all_raw.Rdata'))

load(file = paste0(Rdata_file_path, 'astDF_raw_MRSAbact.Rdata'))
load(file = paste0(Rdata_file_path, 'abxDF_raw_MRSAbact.Rdata'))

load(file = '~/Desktop/EHR-mining/UsefulDataForCleaning/bug_groups.Rdata')
load(file = '~/Desktop/EHR-mining/UsefulDataForCleaning/UPMC_site_names.Rdata')
source('~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/CreateNamedAbxAbbreviations.R')



#### JOIN ENCOUNTER DATA ####
source('~/Desktop/EHR-mining/Scripts/CleaningScripts/featurizeEncounters.R')
df <- featurizeEncounters(df = astDF %>% select(PERSON_ID, ORDER_DAY, ORDER_DATE) %>% distinct(),
                          encsDF = encsDF,
                          preprocess = FALSE,
                          rename = FALSE)
rm(featurizeEncounters, getFacilityNames, site_names)
df_og <- df


# 9.7% are missing encounters data !!
df %>% summarise(sum(MISSING_ENCOUNTER_DATA) / n() * 100)
df <- df %>% filter(!MISSING_ENCOUNTER_DATA) %>% select(-MISSING_ENCOUNTER_DATA)
df <- df %>% filter(FACILITY1 != '', FACILITY2 != '')
df <- df %>% select(-MATERNAL)
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
   relocate(ABX, ABX_TIME, ABX_DAY, .after=ORDER_DATE)
df_abx <- df
df <- df %>% filter(is.na(ABX) | ABX == 'VANCOMYCIN' | ABX == 'DAPTOMYCIN')
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
pcntNoAbxInWindow(0)   # 92%
pcntNoAbxInWindow(0.5) # 32%
pcntNoAbxInWindow(1)   # 14%
pcntNoAbxInWindow(2)   # 8.3%
pcntNoAbxInWindow(4)   # 6%
rm(pcntNoAbxInWindow)


# keep only patients with emiric AND later ABX
# how many patients are left?
df %>% group_by(PERSON_ID, ORDER_DAY) %>% n_groups()                            # 22,745 - total gram-negative blood cultures
df %>% filter(!is.na(ABX)) %>% group_by(PERSON_ID, ORDER_DAY) %>% n_groups()    # 20,024 - received abx within 4 days
df %>% filter(ABX_TIME < 0.5) %>% group_by(PERSON_ID, ORDER_DAY) %>% n_groups() # 19,683 - received empiric ABX
df %>% group_by(PERSON_ID, ORDER_DAY) %>%                                       # 18,883 - received empiric ABX and ABX ≥ 48h
   filter(any(ABX_TIME < 0.5), any(ABX_TIME > 2)) %>% n_groups()
df %>%                                                                          # 13,185 - received empiric TZP/CRO
   filter(ABX_TIME < 0.5 & ABX %in% c('VANCOMYCIN', 'DAPTOMYCIN')) %>%
   group_by(PERSON_ID, ORDER_DAY) %>% n_groups()
df %>% group_by(PERSON_ID, ORDER_DAY) %>%                                       # 12,636 - received empiric TZP/CRO + still receiving ABX 48 hours later
   filter(any(ABX_TIME < 0.5 & ABX %in% c('VANCOMYCIN', 'DAPTOMYCIN')),
          any(ABX_TIME > 2)) %>%
   n_groups()



# how many days / doses of abx after empiric?
# get variable - how many distinct abx each infection was treated with empirically
df <- df %>%
   group_by(PERSON_ID, ORDER_DAY) %>%
   mutate(
      numDaysAbx = length(unique(ABX_DAY)),
      firstAbxTime = min(ABX_TIME),
   ) %>%
   ungroup()


# filter down to at least 2 days of abx trt in 5-day window
df <- df %>%
   filter(!is.na(ABX), ABX_TIME <= 3) %>%
   group_by(PERSON_ID, ORDER_DAY) %>%
   filter(length(unique(ABX_DAY)) >= 2L) %>%
   ungroup()

# remove DAP --> VAN
df <- df %>% 
   group_by(PERSON_ID, ORDER_DAY) %>%
   mutate(TRT = paste(unique(abbr[ABX]), collapse=',')) %>%
   filter(TRT != 'DAP,VAN') %>%
   ungroup()

# remove if >24h overlap
w <- which(df$TRT == 'VAN,DAP')
x <- df[w,] %>%
   group_by(PERSON_ID, ORDER_DAY) %>%
   mutate(
      max_van = max(ABX_TIME[ABX == 'VANCOMYCIN']),
      min_van = min(ABX_TIME[ABX == 'VANCOMYCIN']),
      max_dap = max(ABX_TIME[ABX == 'DAPTOMYCIN']),
      min_dap = min(ABX_TIME[ABX == 'DAPTOMYCIN']),
      min1 = min(max_van, max_dap),
      max1 = max(min_van, min_dap),
      overlap = min1 - max1
   ) %>%
   ungroup() %>%
   filter(max_dap >= max_van, overlap <= 1) %>%
   select(!c(max_van, min_van, max_dap, min_dap, min1, max1, overlap, ABX, ABX_TIME, ABX_DAY)) %>%
   distinct()
df <- rbind(
   df[-w,] %>% select(-ABX, -ABX_TIME, -ABX_DAY) %>% distinct(),
   x
) %>%
   relocate(TRT, .after=ORDER_DATE) %>%
   arrange(PERSON_ID, ORDER_DATE)
rm(w, x)
df %>% count(TRT)


# how many non-empiric abx days for these patients?
df <- left_join(
   x = df,
   y = df_abx %>%
      filter(!ABX %in% c('VANCOMYCIN','DAPTOMYCIN')) %>%
      filter(ABX_TIME >= -2, ABX_TIME < 3) %>% 
      summarise(
         n_ABX_DAY = length(unique(ABX_DAY)),
         .by=c(PERSON_ID, ORDER_DAY)
      ),
   by = join_by(PERSON_ID, ORDER_DAY)
) %>%
   relocate(numDaysAbx, n_ABX_DAY, firstAbxTime, .after=TRT)

# some don't have any! maybe they died?
df %>% count(n_ABX_DAY)
df$n_ABX_DAY[is.na(df$n_ABX_DAY)] <- 0L

# how many received both?
df %>% count(TRT)

# add back other empiric abx...
df <- left_join(
   x = df,
   y = df_abx %>%
      filter(!ABX %in% c('VANCOMYCIN','DAPTOMYCIN')) %>% 
      filter(ABX_TIME < 0.5) %>% 
      select(PERSON_ID, ORDER_DAY, ABX) %>% 
      distinct() %>% 
      mutate(ABX = unname(abbr[ABX])) %>% 
      filter(ABX %in% c('TZP', 'FEP', 'CRO')) %>%#, 'AZM', 'SAM', 'CFZ', 'LVX', 'MEM', 'CLI', 'CIP', 'ATM', 'CPT', 'LZD')) %>% count(ABX, sort=TRUE)
      mutate(EMP_BLI = 1L) %>%
      select(-ABX) %>%
      distinct(),
   by = join_by(PERSON_ID, ORDER_DAY)
) %>%
   relocate(EMP_BLI, .after=firstAbxTime)
df %>% count(EMP_BLI)
df$EMP_BLI[is.na(df$EMP_BLI)] <- 0L
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
rm(labsDF, featurizeLabs)
#### END LAB VALUES ####


df$FACILITY <- gsub('[^A-z]', '', df$FACILITY1)
df <- df %>% select(-FACILITY1, -FACILITY2)

# one-hot encoding facilities as groups
site_groups <- list(
   academic = c("Mercy", "Presbyterian", "Shadyside", "MageeW"), 
   regional = c("Hamot", "Williamsport", "Jameson", "Altoona"), 
   community = c("Passavant", "East", "McKeesport", "StMargaret", 'SOL', 'LOC'), 
   rural = c("Bedford", "Northwest", "Horizon", "Chatauqua")
)
for (g in seq_along(site_groups)) {
   df[[paste0('FAC_', names(site_groups)[g])]] <- as.integer(df$FACILITY %in% site_groups[[g]])
}
df <- df %>% select(-FACILITY)


save(df, file=paste0(Rdata_file_path, 'finalDF_MRSAbact.Rdata'))



#### PLOT STUFF


# plotting variables
Rdata_file_path <- '~/Desktop/EHR/EHR work/RdataFiles/R01/'
plot_path <- '~/Desktop/EHR-mining/Scripts/AnalysisScripts/R01/ITE/MRSA/Plots/'
col_vec <- c('VAN' = '#3399FF', 'DAP' = '#FFaa66')
pmax <- 80
ptcex <- 1.3
lwd <- 2
leg <- function(pch, lwd=NA, pos='topright', pt.cex=ptcex) {
   legend(x=pos, col=col_vec, pch=pch, lwd=lwd, legend=names(col_vec), xpd=NA, pt.cex=ptcex)
}

library(dplyr)
load(file = paste0(Rdata_file_path, 'finalDF_MRSAbact.Rdata'))




df %>% summarise(
   n = n(),
   num_abx_days = mean(n_ABX_DAY), 
   abxBeforeCulture = sum(firstAbxTime < -0/24) / n,
   firstAbxTime_hours = mean(firstAbxTime) * 24, 
   .by=TRT
)

df <- df %>%
   mutate(TRT = ifelse(TRT == 'VAN,DAP', 'DAP', TRT))



#### PLOT 002 - Encounter, bug, facility, age differences ####
{
   t <- df %>%
      summarise(
         across(.cols = c(FEMALE, HOSP_ACQ, NURSING_HOME_HOSPICE, EMERGENCY_DEPT, TRANSFER),
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
          cohort = 'MRSA bacteremia',
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
                                  mutate(across(.cols = c(matches('^(EMP|ICD|FAC)_'), FEMALE, HOSP_ACQ, NURSING_HOME_HOSPICE, EMERGENCY_DEPT, TRANSFER), 
                                                .fns = as.logical)))
t1 <- processTableOne(t1)

pdf(file = paste0(plot_path, '004a_covar_distributions.pdf'), width=9, height=8)
par(mfrow=c(1,2))
xpos <- 0.52
drawTableOne(t1[grep('^ICD_', rownames(t1)),], trt=names(col_vec), mar_left = 11.5, addlegend=TRUE, xpos=xpos)
drawTableOne(t1[-grep('^ICD_', rownames(t1)),], trt=names(col_vec), mar_left = 11.5, xpos=xpos)
# drawTableOne(t1[-grep('^(ICD|LAB|EMP)_', rownames(t1)),], trt=names(col_vec), mar_left = 12, xpos=1.1)
dev.off()



source('~/Desktop/EHR-mining/Scripts/AnalysisScripts/CoxModelFxns.R')

coefs <- modelCoxIPTWcoefs(
   df = df %>% 
      mutate(VAN = TRT == 'VAN') %>%
      select(!c(PERSON_ID:TRT)),
   time_var = 'mortality_time',
   status_var = NULL,
   treatment_var = 'VAN',
   censor_time = 30L,
   conf_vars = NULL,
   print_conf_vars = FALSE,
   normalize_weights = TRUE
)

# propensity score overlap
prop_model <- df %>%
   mutate(DAP = TRT == 'DAP') %>%
   select(DAP, vars) %>%
   glm(DAP ~ ., family=binomial(), data=.)
prop_scores <- predict(prop_model, df[vars], type='response')

{
   pdf(file = paste0(plot_path, '004b_propensity_score_overlap.pdf'), height=4.5, width=9.5)
   layout(matrix(c(1,2,3,3), nrow=2))
   par(mgp=c(1.75, 0.4, 0), tck=-0.01, mar=c(3,3,2,1))
   plotDens <- function(with_weights=FALSE) {
      wV <- rep(1 / sum(df$TRT == 'VAN'), sum(df$TRT == 'VAN'))
      wD <- rep(1 / sum(df$TRT == 'DAP'), sum(df$TRT == 'DAP'))
      if (with_weights) {
         wV <- (prop_scores[df$TRT == 'VAN']) / sum(prop_scores[df$TRT == 'VAN'])
         wD <- (1 - prop_scores[df$TRT == 'DAP']) / sum(1 - prop_scores[df$TRT == 'DAP'])
      }
      plot(NA, xlim=c(-0.05, 1.05), ylim=c(0, 6.1), xlab='Propensity score', ylab='Density')
      title(main=paste0(ifelse(with_weights, 'After', 'Before'), ' propensity score weighting (overlap)'), cex.main=1, line=0.5)
      lines(density(prop_scores[df$TRT == 'VAN'], weights=wV), col=col_vec['VAN'], lwd=2)
      lines(density(prop_scores[df$TRT == 'DAP'], weights=wD), col=col_vec['DAP'], lwd=2)
      legend('topright', legend=names(col_vec), col=col_vec, lwd=2)
   }
   plotDens()
   plotDens(TRUE)
   
   # ATE coefs
   par(mar=c(4,3,2,1))
   plotCoxCoefs(coefs, treatment_var = 'VAN', par_given=TRUE)
   dev.off()
}
#### END PLOT 004 ####




# write to csv for econML
df <- df %>% select(-ICD_OnDialysis_2y, -ICD_Peritonitis_1w, -ICD_AIDS_HIV_2y)
write.csv(x = df %>%
             mutate(
                treatment = as.integer(TRT == 'DAP'),
                outcome = ifelse(!is.na(mortality_time) & mortality_time < 30L, 1L, 0L),
                time = lubridate::decimal_date(ORDER_DATE) - 2017
             ) %>%
             select(!c(TRT, mortality_time, PERSON_ID, ORDER_DAY, ORDER_DATE)),
          file = '~/Desktop/EHR/EHR work/RdataFiles/R01/finalDF_MRSAbact.csv',
          quote = FALSE, 
          row.names = FALSE)














