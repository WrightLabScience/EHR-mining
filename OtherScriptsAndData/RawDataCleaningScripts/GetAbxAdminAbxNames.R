source(file = '~/Desktop/EHR/EHR work/config_file.R')

# GET ANTIBIOTIC SEARCH TERMS
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_clean_ASTs.Rdata')
ast_abx <- astDF %>% select(CEFEPIME:CEPHALEXIN) %>% names
ast_abx <- ast_abx[!ast_abx %in% c('STREPTOMYCIN_HIGH_LEVEL', 'STREPTOMYCIN_SYNERGY', 'GENTAMICIN_HIGH_LEVEL', 'GENTAMICIN_SYNERGY', 'CARBAPENEM_INACTIVATION_TEST', 'BETA_LACTAMASE', 'CLINDAMYCIN_INDUCIBLE')]

# add missing antibiotics
ast_abx <- c(ast_abx, 'PASER', 'SULFAMETHOXAZOLE/TRIMETHOPRIM', 'QUINUPRISTIN', 'DALFOPRISTIN', 'RELEBACTAM', 'FIDAXOMICIN', 'AMPHOTERICIN')

# From the abbreviations file that I took from Andrew's Lancet paper and extended to include additional abx.
abx_abbr <- tibble(read.table(file = '~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/ABX_ABBR.txt', header = TRUE)) %>%
   select(Antibiotic_Name) %>%
   mutate(Antibiotic_Name = gsub('-', '/', Antibiotic_Name),
          Antibiotic_Name = gsub('_', ' ', Antibiotic_Name),
          Antibiotic_Name = toupper(Antibiotic_Name)) %>% 
   unlist() %>%
   unname()

# From the text file of generic and brand antibiotic names I made.
abx_gen <- read.table(file = '~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/AntibioticNames.txt', sep='\t', header = TRUE)
abx_gen <- toupper(abx_gen$Generic)

# Combine, split on /, unique, sort.
abx_all <- sort(unique(c(ast_abx, abx_abbr, abx_gen)))
abx_all <- unlist(strsplit(abx_all, '/'))
abx_all <- sort(unique(abx_all))

# Remove antifungals, anti-parasitics
abx_all <- abx_all[!abx_all %in% c('ANIDULOFUNGIN', 'CASPOFUNGIN', 'CLOTRIMAZOLE', 'FLUCONAZOLE', 'ISAVUCONAZOLE', 'ITRACONAZOLE', 
                                   'KETOCONAZOLE', 'MICAFUNGIN', 'POSACONAZOLE', 'TINIDAZOLE', 'VORICONAZOLE', 'TERBINAFINE', 
                                   'GRISEOFULVIN', 'FLUCYTOSINE', 'MILTEFOSINE', 'AMPHOTERICIN', 'AMPHOTEROCIN B',
                                   'ETHAMBUTOL', 'BEDAQUILINE', 'ISONIAZID', 'RIFABUTIN', 'CLOFAZIMIN', 'PYRAZINAMIDE', 'RIFAPENTINE')]


save(abx_all, file = '~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/abx_names_for_abxDF.Rdata')
rm(abx_gen, ast_abx, abx_abbr)


#### Narrow down to which terms are actually found in med admin views:

# Get all unique medication names from medication administration views.
meds <- tbl(conn, in_schema('AMB_ETL', 'SENS_MED_ADMIN_VW')) %>% count(MEDICATION) %>% collect()
meds_iv <- tbl(conn, in_schema('AMB_ETL', 'SENS_MED_ADMIN_IV_VW')) %>% count(MEDICATION) %>% collect()

# Unique antibiotic names from medication admin views.
abx_admin <- unique(unlist(sapply(abx_all, function(x) grepv(x, meds$MEDICATION))))
abx_admin_iv <- unique(unlist(sapply(abx_all, function(x) grepv(x, meds_iv$MEDICATION))))

save(abx_admin, file = '~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/abx_names_in_med_admin_vw.Rdata')
save(abx_admin_iv, file = '~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/abx_names_in_med_admin_iv_vw.Rdata')

# meds$ABX <- ''
# for (abx in abx_all) {
#    w <- grep(abx, meds$MEDICATION)
#    if (length(w) == 0L) next
#    meds$ABX[w] <- paste(meds$ABX[w], abx, sep=', ')
# }
# meds %>%
#    filter(ABX != '') %>%
#    mutate(ABX = gsub('^, ', '', ABX)) %>%
#    summarise(n=sum(n), .by=ABX) %>%
#    arrange(desc(n)) %>%
#    View()



# Get vasopressor med names from med_admin and med_admin_IV views.
vas_all <- c('PHENYLEPHRINE', 'EPINEPHRINE', 'NOREPINEPHRINE', 'VASOPRESSIN', 'DOPAMINE')
vas_admin <- unique(unlist(sapply(vas_all, function(x) grepv(x, meds$MEDICATION))))
vas_admin_iv <- unique(unlist(sapply(vas_all, function(x) grepv(x, meds_iv$MEDICATION))))

vas_admin <- grepv('OPHTH|CAINE|RACE|TOPIC|SUPP|NASAL|MINERAL|IRRIG', vas_admin, invert=TRUE)
vas_admin_iv <- grepv('OPHTH|CAINE|RACE|TOPIC|SUPP|NASAL|MINERAL|IRRIG', vas_admin_iv, invert=TRUE)

save(vas_admin, file = '~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/vasopressor_names_in_med_admin_vw.Rdata')
save(vas_admin_iv, file = '~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/vasopressor_names_in_med_admin_iv_vw.Rdata')


# Save and find appropriate vasopressors in medications
meds_list <- list(meds, meds_iv)
save(meds_list, file='~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/all_med_names.Rdata')





