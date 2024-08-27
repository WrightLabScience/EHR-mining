source('~/Desktop/EHR/EHR work/config_file.R')

# AST results view
astrDF <- tbl(conn, in_schema('AMB_ETL', 'LAB_MICRO_SENS_ALL_VW')) %>% collect()
print(object.size(astrDF), units='Mb') # 1957.4 Mb - August, 2024
save(astrDF, file = '~/Desktop/EHR/EHR work/RdataFiles/lab_micro_sens_all_vw.Rdata')

# AST orders view
astoDF <- tbl(conn, in_schema('AMB_ETL', 'LAB_MICRO_RESULT_ALL_VW')) %>% collect()
print(object.size(astoDF), units='Mb') # 1396.5 Mb - August, 2024
save(astoDF, file = '~/Desktop/EHR/EHR work/RdataFiles/lab_micro_result_all_vw.Rdata')


# DEATH DATES, AGE, GENDER
getDeathView <- function(view_name) {
   tbl(conn, in_schema('AMB_ETL', view_name)) %>% 
      collect() %>% 
      mutate(TABLE = as.integer(gsub('SENS_DEATH_20([0-9]{2})_VW', '\\1', view_name)))
}
dth <- lapply(grep('DEATH', result$VIEW_NAME, value=TRUE), 
              function(view_name) getDeathView(view_name))
dth <- do.call(rbind, dth)
dth <- dth %>% select(-RESULT_DATE, -TABLE) %>% distinct()
dth %>% filter(n() > 1L, .by=PERSON_ID)
dth <- dth %>% mutate(across(.cols = c(DOB, DEATH_DATE), .fns = ~ as.Date(., format='%m/%d/%Y')))
dth <- dth %>% mutate(AGE_AT_DEATH = as.numeric(lubridate::as.duration(DEATH_DATE - DOB), 'years'))
barplot(table(as.integer(dth$AGE_AT_DEATH)))
barplot(table(as.integer(substr(dth$DOB, 1, 4))))
save(dth, file = '~/Desktop/EHR/EHR work/RdataFiles/SENS_DeathAgeGender.Rdata')




# Med admin (IV too) views
meds    <- tbl(conn, in_schema('AMB_ETL', 'SENS_MED_ADMIN_VW'))    %>% count(MEDICATION, sort=TRUE) %>% collect()
meds_iv <- tbl(conn, in_schema('AMB_ETL', 'SENS_MED_ADMIN_IV_VW')) %>% count(MEDICATION, sort=TRUE) %>% collect()
meds$IV_FLAG <- FALSE
meds_iv$IV_FLAG <- TRUE
meds <- rbind(meds, meds_iv)
meds %>% count(IV_FLAG)


# AST antibiotics
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/antibiotic_names/ast_antibiotics.Rdata')
abx_ast <- abx; rm(abx)
abx_ast <- abx_ast[abx_ast != 'PAS']
abx_ast <- c(abx_ast, 'PASER', 'SULFAMETHOXAZOLE/TRIMETHOPRIM', 'QUINUPRISTIN', 'DALFOPRISTIN', 'RELEBACTAM')

# Antibiotic names, abbreviations, and classes
abx_abbr <- tibble(read.table(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/antibiotic_names/ABX_ABBR.txt', header = TRUE)) %>%
   select(Antibiotic_Name) %>%
   mutate(Antibiotic_Name = gsub('-', '/', Antibiotic_Name),
          Antibiotic_Name = gsub('_', ' ', Antibiotic_Name),
          Antibiotic_Name = toupper(Antibiotic_Name)) %>%
   unlist()
names(abx_abbr) <- NULL

# Antibiotic
abx_gen <- read.table(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/antibiotic_names/AntibioticNames.txt', sep='\t', header = TRUE)
abx_gen <- toupper(abx_gen$Generic)

# COMBINE ALL
abx_all <- sort(unique(c(abx_ast, abx_abbr, abx_gen)))
abx_all <- abx_all[!abx_all %in% c('STREPTOMYCIN_HIGH_LEVEL', 'STREPTOMYCIN_SYNERGY', 'GENTAMICIN_HIGH_LEVEL', 'GENTAMICIN_SYNERGY', 
                                   'CARBAPENEM_INACTIVATION_TEST', 'BETA_LACTAMASE', 'CLINDAMYCIN_INDUCIBLE')]
abx_all <- unlist(strsplit(abx_all, '/'))
abx_all <- sort(unique(abx_all))
write.table(x = paste0('"%', paste(abx_all, collapse='%", "%'), '%"'),
            file='~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/antibiotic_names/Single_abx_names.txt',
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE)
length(abx_all) # 233
rm(abx_gen, abx_ast, abx_abbr)


# search for each antibiotic search pattern in the list of drug_names
# assign the generic antibiotic name to it
meds$ABX <- NA
for (i in seq_along(abx_all)) {
   print(i)
   w <- which(stringi::stri_detect_regex(str = meds$MEDICATION, pattern = abx_all[i]))
   if (grepl('/', abx_all[i])) {
      w <- union(w, which(stringi::stri_detect_regex(str = meds$MEDICATION, pattern = gsub('/', '-', abx_all[i])))) 
   }
   meds$ABX[w] <- paste0(meds$ABX[w], ', ', abx_all[i])
}
meds %>% count(is.na(ABX)) # 1096 matches

# narrow down the list to only the antibiotics and clean up the names
meds <- meds %>%
   filter(!is.na(ABX)) %>%
   mutate(ABX = gsub('NA, ', '', ABX)) %>%
   mutate(ABX = case_when(
      ABX == 'SULFAMETHOXAZOLE, SULFAMETHOXAZOLE/TRIMETHOPRIM, TRIMETHOPRIM, TRIMETHOPRIM/SULFAMETHOXAZOLE' ~ 'TRIMETHOPRIM/SULFAMETHOXAZOLE',
      ABX == 'DALFOPRISTIN, QUINUPRISTIN' ~ 'SYNERCID',
      ABX == 'IMIPENEM, IMIPENEM/CILASTATIN' ~ 'IMIPENEM',
      ABX == 'IMIPENEM, IMIPENEM/CILASTATIN, RELEBACTAM' ~ 'IMIPENEM/RELEBACTAM',
      ABX == 'DOXYCYCLIN, DOXYCYCLINE' ~ 'DOXYCYCLINE',
      ABX == 'SULFAMETHOXAZOLE, TRIMETHOPRIM' ~ 'TRIMETHOPRIM/SULFAMETHOXAZOLE',
      ABX == 'SULFAMETHOXAZOLE, SULFAMETHOXAZOLE/TRIMETHOPRIM, TRIMETHOPRIM' ~ 'TRIMETHOPRIM/SULFAMETHOXAZOLE',
      ABX == 'CIPROFLOXACIN, OFLOXACIN' ~ 'CIPROFLOXACIN',
      ABX == 'PIPERACILLIN, PIPERACILLIN/TAZOBACTAM' ~ 'PIPERACILLIN/TAZOBACTAM',
      ABX == 'AMPICILLIN, AMPICILLIN/SULBACTAM' ~ 'AMPICILLIN/SULBACTAM',
      ABX == 'TRIMETHOPRIM, TRIMETHOPRIM/SULFAMETHOXAZOLE' ~ 'TRIMETHOPRIM/SULFAMETHOXAZOLE',
      ABX == 'LEVOFLOXACIN, OFLOXACIN' ~ 'LEVOFLOXACIN',
      ABX == 'SULFAMETHOXAZOLE, TRIMETHOPRIM/SULFAMETHOXAZOLE' ~ 'TRIMETHOPRIM/SULFAMETHOXAZOLE',
      ABX == 'POLYMYXIN B, POLYMYXIN' ~ 'POLYMYXIN B',
      ABX == 'BACITRACIN, POLYMYXIN B, POLYMYXIN' ~ 'BACITRACIN, POLYMYXIN B',
      ABX == 'CEFTAZIDIME, CEFTAZIDIME/AVIBACTAM' ~ 'CEFTAZIDIME/AVIBACTAM',
      ABX == 'POLYMYXIN B, POLYMYXIN, TRIMETHOPRIM' ~ 'POLYMYXIN B, TRIMETHOPRIM',
      ABX == 'TICARCILLIN, TICARCILLIN/CLAVULANATE' ~ 'TICARCILLIN/CLAVULANATE',
      ABX == 'SULFAMETHOXAZOLE, TRIMETHOPRIM, TRIMETHOPRIM/SULFAMETHOXAZOLE' ~ 'TRIMETHOPRIM/SULFAMETHOXAZOLE',
      ABX == 'BACITRACIN, NEOMYCIN, POLYMYXIN B, POLYMYXIN' ~ 'BACITRACIN, NEOMYCIN, POLYMYXIN B',
      ABX == 'IMIPENEM, IMIPENEM/RELEBACTAM' ~ 'IMIPENEM/RELEBACTAM',
      ABX == 'MEROPENEM, MEROPENEM/VABORBACTAM' ~ 'MEROPENEM/VABORBACTAM',
      ABX == 'CLOXACILLIN, DICLOXACILLIN, OXACILLIN' ~ 'DICLOXACILLIN',
      ABX == 'AMOXICILLIN, AMOXICILLIN/CLAVULANATE' ~ 'AMOXICILLIN/CLAVULANATE',
      ABX == 'PENICILLIN, PENICILLIN G' ~ 'PENICILLIN',
      ABX == 'PENICILLIN, PENICILLIN V' ~ 'PENICILLIN',
      .default = ABX
   ))

save(meds, file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/antibiotic_names/med_admin_abx_names.Rdata')




load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/antibiotic_names/med_admin_abx_names.Rdata')
write.table(x = paste0("'", paste(sort(unique(unlist(strsplit(unique(meds$ABX[meds$IV_FLAG]), split=', |/')))), collapse='|'), "'"),
            file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/antibiotic_names/Med_IV_admin_abx_names.txt',
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
write.table(x = paste0("regexp_like(MEDICATION, '", sort(unique(unlist(strsplit(unique(meds$ABX[!meds$IV_FLAG]), split=', |/')))), "') or"),
            file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/antibiotic_names/Med_admin_abx_names.txt',
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)






