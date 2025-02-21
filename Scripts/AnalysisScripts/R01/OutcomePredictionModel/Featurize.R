##### START FEATURE ENGINEERING #####
# load required data and functions
start <- Sys.time()
library(dplyr)
Rdata_file_path <- '~/Desktop/EHR/EHR work/RdataFiles/R01/OutcomePredictionModel/'
load(file = paste0(Rdata_file_path, 'cohortDF_all.Rdata'))
load(file = paste0(Rdata_file_path, 'encsDF_all_raw.Rdata'))
load(file = paste0(Rdata_file_path, 'labsDF_all_cleaned.Rdata'))
load(file = paste0(Rdata_file_path, 'icdsDF_all_raw.Rdata'))
load(file = paste0(Rdata_file_path, 'abxtDF_all.Rdata'))
source('~/Desktop/EHR-mining/Scripts/CleaningScripts/featurizeLabs.R')
source('~/Desktop/EHR-mining/Scripts/CleaningScripts/featurizeEncounters.R')
source('~/Desktop/EHR-mining/Scripts/CleaningScripts/featurizeICDs.R')
source('~/Desktop/EHR-mining/Scripts/CleaningScripts/featurizeAcutePathogen.R')
source('~/Desktop/EHR-mining/Scripts/CleaningScripts/featurizeDemographics.R')
source('~/Desktop/EHR-mining/Scripts/CleaningScripts/featurizeXdayOutcome.R')
source('~/Desktop/EHR-mining/Scripts/CleaningScripts/featurizeReadmissionLOS.R')
source('~/Desktop/EHR-mining/Scripts/CleaningScripts/featurizeBSIrecurrence.R')
w <- which(labsDF$LAB == 'PCT_NEUTROPHILS')
labsDF <- labsDF[-w[which.max(labsDF$RESULT_VALUE[w])],]
rm(w)

# featurization steps
pathDF <- featurizeAcutePathRes(df) # acute pathogen/resistance
labsDF <- featurizeLabs(df, labsDF) # lab results
admtDF <- featurizeEncounters(df, encsDF) # encounters data
icdsDF <- featurizeICDs(df, icdsDF) # ICD codes
demoDF <- getDemographics(df) # demographics

# encode outcomes
mortDF <- XdayOutcome(df=demoDF, col_name='mortality_time') %>% select(-mortality_time)
readDF <- featurizeReadmissionLOS(df, encsDF)
recrDF <- getBSIrecurrence(df)


# combine into one large dataset
data <- inner_join(
   x = inner_join(
      x = inner_join(
         x = inner_join(
            x = pathDF,
            y = abxtDF,
            by = join_by(PERSON_ID, ORDER_DAY)
         ),
         y = inner_join(
            x = labsDF,
            y = admtDF,
            by = join_by(PERSON_ID, ORDER_DAY)
         ),
         by = join_by(PERSON_ID, ORDER_DAY)
      ),
      y = inner_join(
         x = inner_join(
            x = icdsDF,
            y = demoDF,
            by = join_by(PERSON_ID, ORDER_DAY)
         ),
         y = inner_join(
            x = mortDF,
            y = readDF,
            by = join_by(PERSON_ID, ORDER_DAY)
         ),
         by = join_by(PERSON_ID, ORDER_DAY)
      ),
      by = join_by(PERSON_ID, ORDER_DAY)
   ),
   y = recrDF,
   by = join_by(PERSON_ID, ORDER_DAY)
)

# create composite event
minEventTime <- function(x, y, z) {
   if (all(is.na(c(x, y, z))))
      return(NA_integer_)
   return(min(c(x, y, z), na.rm=T))  
}

data <- data %>%
   rowwise() %>%
   mutate(event_time = minEventTime(mortality_time, readmit_time, BSIrecur_time)) %>%
   relocate(event_time, .after=BSIrecur_time) %>%
   ungroup()

data <- XdayOutcome(data, col_name='event_time', includes_pid = FALSE)

save(data, file = paste0(Rdata_file_path, 'data_feat_outcomes_all.Rdata'))
print(Sys.time() - start)





#### PREPROCESS DATA ####
library(dplyr)
Rdata_file_path <- '~/Desktop/EHR/EHR work/RdataFiles/R01/OutcomePredictionModel/'
load(file = paste0(Rdata_file_path, 'data_feat_outcomes_all.Rdata'))

# remove rows with high missingness - MISSING_ENCOUNTER_DATA column
w <- which(data$ENC_MISSING_ENCOUNTER_DATA) # 2781
if (length(w) > 0L)
   data <- data[-w,]
data <- data %>% select(-ENC_MISSING_ENCOUNTER_DATA)

w <- which(is.na(data$EmpDisc)) # 1353
if (length(w) > 0L)
   data <- data[-w,]

w <- which(data$ENC_FACILITY1 == '' | data$ENC_FACILITY2 == '') # 73
if (length(w) > 0L)
   data <- data[-w,]


# remove patients that died within 7 days
w <- which(data$mortality_time <= 4) # 1810
if (length(w) > 0L)
   data <- data[-w,]
rm(w)


# remove variables with less than 80% completion rate
missDF <- data.frame(
   pcnt_complete = sort(sapply(data %>% select(-PERSON_ID, -ORDER_DAY),
                               FUN = function(x) sum(!is.na(x)) / length(x)))
)
remove_vars <- rownames(missDF)[missDF$pcnt_complete < 0.8] # 20 variables
remove_vars <- remove_vars[!remove_vars %in% c('BSIrecur_time', 'mortality_time', 'readmit_time', 'event_time')]
if (length(remove_vars) > 0L) {
   cat('Removing', length(remove_vars), 'variables from data for high missingness.')
   data <- data %>% select(!c(!!remove_vars))
} else {
   cat('No variables with > 20% missingness!')
}
rm(missDF, remove_vars)

rows_comp_rate <- apply(data %>% select(-PERSON_ID, -ORDER_DAY), 
                        MARGIN = 1, 
                        FUN = function(x) sum(!is.na(x)) / length(x))
hist(rows_comp_rate, xlim=c(0, 1), ylim=c(0,400))
w <- which(rows_comp_rate < 0.8)
if (length(w) > 0L) {
   cat(length(w), 'patients with fewer than 80% variables complete.\n')
   data <- data[-w,]
} else {
   cat('No patient with > 20% missingness!')
}
rm(w, rows_comp_rate)


# variables with low/no variance ??
no_var_vars <- names(which(sapply(data %>% select(-PERSON_ID, -ORDER_DAY),
                                  FUN = function(x) length(unique(x)) == 1L)))
if (length(no_var_vars) > 0L) {
   cat(length(no_var_vars), 'variables removed for being constant')
   data <- data %>% select(!c(!!no_var_vars))
} else {
   cat('No variables low/no variance!')
}


# ensure target binary variables and facility strings are interpreted as factors
data <- data %>%
   mutate(
      across(.cols = c(ENC_FACILITY1, ENC_FACILITY2, matches('^d[0-9]+')),
             .fns = as.factor)
   )


save(data, file = paste0(Rdata_file_path, 'data_feat_outcomes_processed_all.Rdata'))




