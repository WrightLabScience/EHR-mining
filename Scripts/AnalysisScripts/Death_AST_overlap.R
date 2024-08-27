library(dplyr)

load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_clean_ASTs.Rdata')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/SENS_DeathAgeGender.Rdata')

length(unique(astDF$PERSON_ID)) # 613,854
length(unique(dth$PERSON_ID))   # 93,953
length(intersect(astDF$PERSON_ID, dth$PERSON_ID)) # 93,747

length(unique(astDF$PERSON_ID[astDF$BLOOD])) # 77,977
length(intersect(astDF$PERSON_ID[astDF$BLOOD], dth$PERSON_ID)) # 26,921