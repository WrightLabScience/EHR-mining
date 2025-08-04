source(file = '~/Desktop/EHR/EHR work/config_file.R')

setwd('~/Desktop/EHR/EHR work/RdataFiles/MED_ADMIN_DATA/')

# SENS_MED_ADMIN_VW
load(file='~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/abx_names_in_med_admin_vw.Rdata')
for (i in seq_along(abx_meds)) {
   abx <- abx_meds[i]
   query <- paste0("SELECT * FROM AMB_ETL.SENS_MED_ADMIN_VW WHERE MEDICATION LIKE '%", abx, "%'")
   start <- Sys.time()
   x <- tibble(dbGetQuery(conn, query))
   write.table(x, file = paste0('SENS_MED_ADMIN_ABX_DATA/', abx, '.txt'), sep='\t', quote=FALSE, row.names=FALSE)
   end <- Sys.time()
   cat(i, abx, end-start, '\n')
}

print('DONE WITH SENS_MED_ADMIN_VW!')
print('Moving onto SENS_MED_ADMIN_IV_VW')

# SENS_MED_ADMIN_IV_VW
load(file='~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/abx_names_in_med_admin_iv_vw.Rdata')
for (i in seq_along(abx_meds_iv)) {
   abx <- abx_meds_iv[i]
   query <- paste0("SELECT * FROM AMB_ETL.SENS_MED_ADMIN_IV_VW WHERE MEDICATION LIKE '%", abx, "%'")
   start <- Sys.time()
   x <- tibble(dbGetQuery(conn, query))
   write.table(x, file = paste0('SENS_MED_ADMIN_IV_ABX_DATA/', abx, '_IV.txt'), sep='\t', quote=FALSE, row.names=FALSE)
   end <- Sys.time()
   cat(i, abx, end-start, '\n')
}

print('DONE WITH ALL!')
