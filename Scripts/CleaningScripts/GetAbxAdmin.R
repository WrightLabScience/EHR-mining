source(file = '~/Desktop/EHR/EHR work/config_file.R')

# GET ANTIBIOTIC SEARCH TERMS
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/antibiotic_names/ast_antibiotics.Rdata')
abx_ast <- abx; rm(abx)
abx_ast <- abx_ast[abx_ast != 'PAS']
abx_ast <- c(abx_ast, 'PASER', 'SULFAMETHOXAZOLE/TRIMETHOPRIM', 'QUINUPRISTIN', 'DALFOPRISTIN', 'RELEBACTAM')
abx_abbr <- tibble(read.table(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/antibiotic_names/ABX_ABBR.txt', header = TRUE)) %>%
   select(Antibiotic_Name) %>%
   mutate(Antibiotic_Name = gsub('-', '/', Antibiotic_Name),
          Antibiotic_Name = gsub('_', ' ', Antibiotic_Name),
          Antibiotic_Name = toupper(Antibiotic_Name)) %>% unlist()
names(abx_abbr) <- NULL
abx_gen <- read.table(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/antibiotic_names/AntibioticNames.txt', sep='\t', header = TRUE)
abx_gen <- toupper(abx_gen$Generic)
abx_all <- sort(unique(c(abx_ast, abx_abbr, abx_gen)))
abx_all <- abx_all[!abx_all %in% c('STREPTOMYCIN_HIGH_LEVEL', 'STREPTOMYCIN_SYNERGY', 'GENTAMICIN_HIGH_LEVEL', 'GENTAMICIN_SYNERGY', 'CARBAPENEM_INACTIVATION_TEST', 'BETA_LACTAMASE', 'CLINDAMYCIN_INDUCIBLE')]
abx_all <- unlist(strsplit(abx_all, '/'))
abx_all <- sort(unique(abx_all))

# Narrow down to which terms are actually found in sens_med_admin_vw
meds <- tbl(conn, in_schema('AMB_ETL', 'SENS_MED_ADMIN_VW')) %>% count(MEDICATION) %>% collect()
abx_meds <- abx_all[sapply(abx_all, function(x) any(grepl(x, meds$MEDICATION)))]

# loop over each search term to query the view
for (i in seq_along(abx_meds)) {
   abx <- abx_meds[i]
   query <- paste0("SELECT * FROM AMB_ETL.SENS_MED_ADMIN_VW WHERE MEDICATION LIKE '%", abx, "%'")
   start <- Sys.time()
   x <- tibble(dbGetQuery(conn, query))
   write.table(x, file = paste0('~/Desktop/EHR/EHR work/RdataFiles/Sens_Med_Admin_txtfiles/', abx, '.txt'), sep='\t', quote=FALSE, row.names=FALSE)
   end <- Sys.time()
   cat(i, abx, end-start, '\n')
}
print('DONE WITH SENS_MED_ADMIN_VW!')



meds_iv <- tbl(conn, in_schema('AMB_ETL', 'SENS_MED_ADMIN_IV_VW')) %>% count(MEDICATION) %>% collect()
abx_meds_iv <- abx_all[sapply(abx_all, function(x) any(grepl(x, meds_iv$MEDICATION)))]

# loop over each search term to query the view
print('Moving onto SENS_MED_ADMIN_IV_VW')
for (i in seq_along(abx_meds_iv)) {
   abx <- abx_meds_iv[i]
   query <- paste0("SELECT * FROM AMB_ETL.SENS_MED_ADMIN_IV_VW WHERE MEDICATION LIKE '%", abx, "%'")
   start <- Sys.time()
   x <- tibble(dbGetQuery(conn, query))
   write.table(x, file = paste0('~/Desktop/EHR/EHR work/RdataFiles/Sens_Med_Admin_txtfiles/', abx, '_IV.txt'), sep='\t', quote=FALSE, row.names=FALSE)
   end <- Sys.time()
   cat(i, abx, end-start, '\n')
}
print('DONE WITH ALL!')


path <- '~/Desktop/EHR/EHR work/RdataFiles/Sens_Med_Admin_txtfiles/'
ls <- list.files(path)

abxDF <- tibble(read.table(file = paste0(path, ls[1]), header=TRUE, sep='\t'))
for (i in 2:length(ls)) {
   cat(i, ls[i], '\n')
   r <- tibble(read.table(file = paste0(path, ls[i]), header=TRUE, sep='\t'))
   abxDF <- rbind(abxDF, r)
}

abxDF <- abxDF %>%
   distinct() %>%
   arrange(PERSON_ID, ADMIN_START_DATE)

save(abxDF, file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_2017_2023_AbxAdmin.Rdata')



