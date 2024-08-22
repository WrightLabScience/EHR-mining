source('~/Desktop/EHR/EHR work/config_file.R')

# AST results view
astrDF <- tbl(conn, in_schema('AMB_ETL', 'LAB_MICRO_SENS_ALL_VW')) %>% collect()
print(object.size(astrDF), units='Mb') # 1957.4 Mb - August, 2024
save(astrDF, file = '~/Desktop/EHR/EHR work/RdataFiles/lab_micro_sens_all_vw.Rdata')

# AST orders view
astoDF <- tbl(conn, in_schema('AMB_ETL', 'LAB_MICRO_RESULT_ALL_VW')) %>% collect()
print(object.size(astoDF), units='Mb') # 1396.5 Mb - August, 2024
save(astoDF, file = '~/Desktop/EHR/EHR work/RdataFiles/lab_micro_result_all_vw.Rdata')

# Med admin (IV too) views
abxDF <- tbl(conn, in_schema('AMB_ETL', 'SENS_MED_ADMIN_VW')) %>% collect()
abxDF_iv <- tbl(conn, in_schema('AMB_ETL', 'SENS_MED_ADMIN_IV_VW')) %>% collect()
print(object.size(abxDF), units='Mb') # 2111.1 Mb
print(object.size(abxDF_iv), units='Mb') # 616.5 Mb
abxDF <- rbind(abxDF, abxDF_iv); rm(abxDF_iv)
save(abxDF, file = '~/Desktop/EHR/EHR work/RdataFiles/sens_med_admin_and_iv_vw.Rdata')