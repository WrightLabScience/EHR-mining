source('~/Desktop/EHR/EHR work/config_file.R')

# AST results view
astrDF <- tbl(conn, in_schema('AMB_ETL', 'LAB_MICRO_SENS_ALL_VW')) %>% collect()
astrDF2 <- tbl(conn, in_schema('AMB_ETL', 'LAB_MICRO_ALL_VW')) %>% collect()
print(object.size(astrDF), units='Mb')  # 1527.3 Mb - Feb, 2025
print(object.size(astrDF2), units='Mb') # 1697.2 Mb - Feb, 2025
length(unique(astrDF$PERSON_ID)) # 493,983
length(unique(astrDF2$PERSON_ID)) # 637,906
length(intersect(astrDF$PERSON_ID, astrDF2$PERSON_ID)) # 493,983
save(astrDF, file = '~/Desktop/EHR/EHR work/RdataFiles/lab_micro_sens_all_vw.Rdata')
save(astrDF2, file = '~/Desktop/EHR/EHR work/RdataFiles/lab_micro_all_vw.Rdata')

# AST orders view
astoDF <- tbl(conn, in_schema('AMB_ETL', 'LAB_MICRO_RESULT_ALL_VW')) %>% collect()
print(object.size(astoDF), units='Mb') # 1464.4 Mb - Feb, 2025
save(astoDF, file = '~/Desktop/EHR/EHR work/RdataFiles/lab_micro_result_all_vw.Rdata')

