source(file = '~/Desktop/EHR/EHR work/config_file.R')

# GET ANTIBIOTIC SEARCH TERMS
load(file = '~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/ast_antibiotics.Rdata')
abx_ast <- abx; rm(abx)
abx_ast <- abx_ast[abx_ast != 'PAS']
abx_ast <- abx_ast[!is.na(abx_ast)]
# add missing antibiotics
abx_ast <- c(abx_ast, 'PASER', 'SULFAMETHOXAZOLE/TRIMETHOPRIM', 'QUINUPRISTIN', 'DALFOPRISTIN', 'RELEBACTAM', 'FIDAXOMICIN', 'AMPHOTERICIN')
abx_abbr <- tibble(read.table(file = '~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/ABX_ABBR.txt', header = TRUE)) %>%
   select(Antibiotic_Name) %>%
   mutate(Antibiotic_Name = gsub('-', '/', Antibiotic_Name),
          Antibiotic_Name = gsub('_', ' ', Antibiotic_Name),
          Antibiotic_Name = toupper(Antibiotic_Name)) %>% unlist()
names(abx_abbr) <- NULL
abx_gen <- read.table(file = '~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/AntibioticNames.txt', sep='\t', header = TRUE)
abx_gen <- toupper(abx_gen$Generic)
abx_all <- sort(unique(c(abx_ast, abx_abbr, abx_gen)))
abx_all <- abx_all[!abx_all %in% c('STREPTOMYCIN_HIGH_LEVEL', 'STREPTOMYCIN_SYNERGY', 'GENTAMICIN_HIGH_LEVEL', 'GENTAMICIN_SYNERGY', 'CARBAPENEM_INACTIVATION_TEST', 'BETA_LACTAMASE', 'CLINDAMYCIN_INDUCIBLE')]
abx_all <- unlist(strsplit(abx_all, '/'))
abx_all <- sort(unique(abx_all))
save(abx_all, file='~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/abx_names_for_abxDF.Rdata')
rm(abx_gen, abx_ast, abx_abbr)

# Narrow down to which terms are actually found in med admin views:

# SENS_MED_ADMIN_VW
meds <- tbl(conn, in_schema('AMB_ETL', 'SENS_MED_ADMIN_VW')) %>% count(MEDICATION) %>% collect()
abx_meds <- abx_all[sapply(abx_all, function(x) any(grepl(x, meds$MEDICATION)))]
length(abx_meds) # 108

save(abx_meds, file='~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/abx_names_in_med_admin_vw.Rdata')


# SENS_MED_ADMIN_IV_VW
meds_iv <- tbl(conn, in_schema('AMB_ETL', 'SENS_MED_ADMIN_IV_VW')) %>% count(MEDICATION) %>% collect()
abx_meds_iv <- abx_all[sapply(abx_all, function(x) any(grepl(x, meds_iv$MEDICATION)))]
length(abx_meds_iv) # 12

save(abx_meds_iv, file='~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/abx_names_in_med_admin_iv_vw.Rdata')






