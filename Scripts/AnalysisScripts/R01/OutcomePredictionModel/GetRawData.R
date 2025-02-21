library(dplyr)
Rdata_file_path <- '~/Desktop/EHR/EHR work/RdataFiles/R01/OutcomePredictionModel/'

##### START get cohort + abx treatment encoding #####
source('~/Desktop/EHR-mining/Scripts/CleaningScripts/getIndexASTsAbxImputeFlag.R')
years <- 2017:2023
df <- getIndexASTsABXsImputeFlag(year = years, 
                                 index_culture_time_filter = 42L)

abxtDF <- df %>% 
   select(-ORDER_DATE) %>% 
   mutate(across(c(EmpDisc, DefDisc), 
                 as.integer))
df <- df %>% 
   select(PERSON_ID, ORDER_DAY, ORDER_DATE)
rm(rs, bug_groups, abx_class, abbr, getAbxClass, getAbxClassResistanceFlags, getBugClass, getDiscordanceFlag, getIndexASTsABXsImputeFlag, imputeASTs); gc()
save(df, file = paste0(Rdata_file_path, 'cohortDF_all.Rdata'))
save(abxtDF, file = paste0(Rdata_file_path, 'abxtDF_all.Rdata'))


##### GET RAW DATA #####

# lab results
source('~/Desktop/EHR-mining/Scripts/CleaningScripts/getRawLabsData.R')
source('~/Desktop/EHR-mining/Scripts/CleaningScripts/CleanLabResults.R')
labsDF <- getRawLabsData(person_ids=unique(df$PERSON_ID), years=years)
save(labsDF, file = paste0(Rdata_file_path, 'labsDF_all_raw.Rdata'))
labsDF <- cleanLabResults(labsDF)
rm(getRawLabsData, cleanLabResults)
save(labsDF, file = paste0(Rdata_file_path, 'labsDF_all_cleaned.Rdata'))

# encounters
source('~/Desktop/EHR-mining/Scripts/CleaningScripts/getRawEncountersData.R')
encsDF <- getRawEncountersData(person_ids=unique(df$PERSON_ID), years=years)
rm(getRawEncountersData)
save(encsDF, file = paste0(Rdata_file_path, 'encsDF_all_raw.Rdata'))

# comorbidities
source('~/Desktop/EHR-mining/Scripts/CleaningScripts/GetRawICDsData.R')
icdsDF <- getRawICDsData(person_ids=unique(df$PERSON_ID), years=years)
rm(getRawICDsData)
save(icdsDF, file = paste0(Rdata_file_path, 'icdsDF_all_raw.Rdata'))







##### SEPSIS ######



library(dplyr)
Rdata_file_path <- '~/Desktop/EHR/EHR work/RdataFiles/R01/OutcomePredictionModel/'

##### START get cohort + abx treatment encoding #####
source('~/Desktop/EHR-mining/Scripts/CleaningScripts/getIndexASTsAbxImputeFlag.R')
years <- 2017:2023
df <- getIndexASTsABXsImputeFlag(year = years, 
                                 index_culture_time_filter = 42L)

abxtDF <- df %>% 
   select(-ORDER_DATE) %>% 
   mutate(across(c(EmpDisc, DefDisc), 
                 as.integer))
df <- df %>% 
   select(PERSON_ID, ORDER_DAY, ORDER_DATE)
rm(rs, bug_groups, abx_class, abbr, getAbxClass, getAbxClassResistanceFlags, getBugClass, getDiscordanceFlag, getIndexASTsABXsImputeFlag, imputeASTs); gc()
save(df, file = paste0(Rdata_file_path, 'cohortDF_all.Rdata'))
save(abxtDF, file = paste0(Rdata_file_path, 'abxtDF_all.Rdata'))


##### GET RAW DATA #####

# lab results
source('~/Desktop/EHR-mining/Scripts/CleaningScripts/getRawLabsData.R')
source('~/Desktop/EHR-mining/Scripts/CleaningScripts/CleanLabResults.R')
labsDF <- getRawLabsData(person_ids=unique(df$PERSON_ID), years=years)
save(labsDF, file = paste0(Rdata_file_path, 'labsDF_all_raw.Rdata'))
labsDF <- cleanLabResults(labsDF)
rm(getRawLabsData, cleanLabResults)
save(labsDF, file = paste0(Rdata_file_path, 'labsDF_all_cleaned.Rdata'))

# encounters
source('~/Desktop/EHR-mining/Scripts/CleaningScripts/getRawEncountersData.R')
encsDF <- getRawEncountersData(person_ids=unique(df$PERSON_ID), years=years)
rm(getRawEncountersData)
save(encsDF, file = paste0(Rdata_file_path, 'encsDF_all_raw.Rdata'))

# comorbidities
source('~/Desktop/EHR-mining/Scripts/CleaningScripts/GetRawICDsData.R')
icdsDF <- getRawICDsData(person_ids=unique(df$PERSON_ID), years=years)
rm(getRawICDsData)
save(icdsDF, file = paste0(Rdata_file_path, 'icdsDF_all_raw.Rdata'))




