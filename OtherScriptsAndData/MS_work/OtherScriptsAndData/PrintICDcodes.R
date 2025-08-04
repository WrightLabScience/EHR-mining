# Get all ICD-10 diagnosis codes used for each comorbidity

library(dplyr)
load(file = '~/Desktop/EHR/EHR work/RdataFiles/BuildCohortDatasets/RawData/icdsDF_raw_2015_2024_bacteremia.Rdata')
x <- icdsDF %>% count(DX_CODE, CODE_TYPE, CODE_DESCRIPTION)
x$CODE_DESCRIPTION <- tolower(x$CODE_DESCRIPTION)

a <- c(
   "SepticShock_1w = any(grepl('with septic shock', CODE_DESCRIPTION))",
   "Sepsis_1w = any(grepl('sepsis', CODE_DESCRIPTION))", 
   "AKI_1w = any(grepl('^N17.9', DX_CODE))",
   "Endocarditis_1w = any(grepl('^I33.0', DX_CODE))",
   "Cellulitis_1w = any(grepl('^L03.90', DX_CODE))",
   "Peritonitis_1w = any(grepl('^K65.9', DX_CODE))",
   "PulmCircDis_1m = any(grepl('^I26|^I27|^I28\\.[089]', DX_CODE))",
   "CPD_1m = any(grepl('^I27\\.28|^I27\\.9|J4[0-7]|^J6[0-7]|^J68\\.4|^J70\\.[13]', DX_CODE))",
   "MyocInfarc_1m = any(grepl('^I21|^I22|^I25', DX_CODE))",
   "MetastSolidTumor_1m = any(grepl('^C7[7-9]|^C80', DX_CODE))",
   "Malignancy_1m = any(grepl('^C[01][0-9]|^C2[0-6]|^C3[01234789]|^C4[013]|^C4[5-9]|^C5[0-8]|^C6[0-9]|^C7[0-6]|^C8[123458]|^C9[0-7]', DX_CODE))",
   "Obesity_2y = any(grepl('^E66', DX_CODE))",
   "WeightLoss_2y = any(grepl('^E4[0-6]|^R63\\.4|^R64', DX_CODE))",
   "Anemia_2y = any(grepl('^D50\\.[089]|^D5[1-3]', DX_CODE))",
   "FluidElectroDis_2y = any(grepl('^E22\\.2|^E86|^E87', DX_CODE))",
   "Coagulopathy_2y = any(grepl('^D6[5-8]|^D69\\.[16789]', DX_CODE))",
   "Alcohol_2y = any(grepl('^F10|^E52|^G62\\.1|^I42\\.6|^K29\\.2|^K70\\.0|^K70\\.3|^K70\\.9|^T51|^Z50\\.2|^Z71\\.4|^Z72\\.1', DX_CODE))",
   "Drugs_2y = any(grepl('^F1[12345689]|Z71\\.5|Z72\\.2', DX_CODE))",
   "Psychoses_2y = any(grepl('^F20|^F2[234589]|F30\\.2|F31\\.2|F31\\.5', DX_CODE))",
   "Depression_2y = any(grepl('^F20\\.4|^F31\\.[3-5]|^F32|^F33|^F34\\.1|^F41\\.2|^F43\\.2', DX_CODE))",
   "NeuroDisease_2y = any(grepl('^G1[0-3]|^G2[0-2]|^G25\\.4|^G25\\.5|^G31\\.[289]|^G32|^G3[5-7]|^G4[01]|^G93\\.[14]|^R47\\.0|^R56', DX_CODE))",
   "CardiacArrythm_2y = any(grepl('^I44\\.[1-3]|^I45\\.6|^I45\\.9|^I4[7-9]|^R00\\.[018]|^T82\\.1|^Z45\\.0|^Z95\\.0', DX_CODE))",
   "MyocInfarc_2y = any(grepl('^I21|^I22|^I25', DX_CODE))",
   "CompHypertension_2y = any(grepl('^I1[135]', DX_CODE))",
   "UncompHypertension_2y = any(grepl('^I10', DX_CODE))",
   "CongHeartFailure_2y = any(grepl('^I42|^I43|^I50', DX_CODE))",
   "PeriphVasDis_2y = any(grepl('^I70|^I71|^I73|^I79|^K55|^Z95', DX_CODE))",
   "CereVasDis_2y = any(grepl('^G45|^G46|^I6[0-9]', DX_CODE))",
   "Dementia_2y = any(grepl('^F0[01235]|^G30', DX_CODE))",
   "CPD_Pneum_2y = any(grepl('^I26|^I27|^I28\\.[089]|^J4[0-7]|^J6[0-8]|^J70', DX_CODE))", # this includes additional codes from Elixhauser
   "RheumaticDis_2y = any(grepl('^L94\\.[013]|^M05|^M06|^M08|^M12\\.[03]|^M3[0-6]|^M45|M46\\.[189]', DX_CODE))", # this includes additional codes from Elixhauser
   "PepticUlcerDis_2y = any(grepl('^K2[5-8]', DX_CODE))",
   "MildLiverDis_2y = any(grepl('^B18|^K7[1346]|^Z94', DX_CODE))",
   "Diabetes_2y = any(grepl('^E1[1-4]', DX_CODE))",
   "HemiParaplegia_2y = any(grepl('^G8[0-4]', DX_CODE))",
   "RenalDis_2y = any(grepl('^N03|^N05|^N18|^N19|^Z49', DX_CODE))",
   "Malignancy_2y = any(grepl('^C[01][0-9]|^C2[0-6]|^C3[01234789]|^C4[013]|^C4[5-9]|^C5[0-8]|^C6[0-9]|^C7[0-6]|^C8[123458]|^C9[0-7]', DX_CODE))",
   "ModSevLivDis_2y = any(grepl('I85|K72', DX_CODE))",
   "MetastSolidTumor_2y = any(grepl('^C7[7-9]|^C80', DX_CODE))",
   "AIDS_HIV_2y = any(grepl('^B2[0124]', DX_CODE))",
   "Hyperlipid_2y = any(grepl('(^| )hyperlipid', CODE_DESCRIPTION))",
   "Smoking_2y = any(grepl('nicotine', CODE_DESCRIPTION))"
)
a <- sapply(1:3, function(x) gsub("^(.+) = any\\(grepl\\('(.+)', (.+)\\)\\)", paste0('\\', x), a))

mat <- matrix(NA, nrow=nrow(a), ncol=2, dimnames=list(NULL, c('Covariate', 'ICD-10 code(s)')))

for (i in seq_len(nrow(a))) {
   w <- grepl(a[i,2], x[[a[i,3]]])
   mat[i, 'Covariate'] <- a[i,1]
   mat[i, 'ICD-10 code(s)'] <- paste(x$DX_CODE[w], collapse=', ')
}



b <- c(
   "Respiratory_1w = any(grepl('^J', DX_CODE) & grepl('infect', CODE_DESCRIPTION))",
   "Osteomyelitis_1w = any(grepl('^M86', DX_CODE) & !grepl('chronic', CODE_DESCRIPTION))",
   "OsteoChronic_2y = any(grepl('M86', DX_CODE) & grepl('chronic', CODE_DESCRIPTION))",
   "Hypothyroid_2y = any(grepl('^E0[0-3]|^E89', DX_CODE) & !grepl('E000\\.0|E030|E015\\.2|E001\\.0', DX_CODE))"
)
b <- sapply(1:6, function(x) gsub("^(.+) = any\\(grepl\\('(.+)', (.+)\\) & (!)?grepl\\('(.+)', (.+)\\)\\)", paste0('\\', x), b))


mat_b <- matrix(NA, nrow=nrow(b), ncol=2, dimnames=dimnames(mat))
for (i in seq_len(nrow(b))) {
   if (b[i,4] == '')  w <- grepl(b[i,2], x[[b[i,3]]]) & grepl(b[i,5], x[[b[i,6]]])
   if (b[i,4] == '!') w <- grepl(b[i,2], x[[b[i,3]]]) & !grepl(b[i,5], x[[b[i,6]]])
   mat_b[i, 'Covariate'] <- b[i,1]
   mat_b[i, 'ICD-10 code(s)'] <- paste(x$DX_CODE[w], collapse=', ')
}

mat <- rbind(mat, mat_b)


mat <- rbind(mat,
             matrix(c('OnDialysis_2y', x$DX_CODE[x$CODE_DESCRIPTION == 'dependence on renal dialysis']),
                    nrow=1,
                    dimnames=dimnames(mat)))

mat <- mat[c(grep('1w', mat[,'Covariate']), grep('1m', mat[,'Covariate']), grep('2y', mat[,'Covariate'])), ]

colnames(mat) <- c('Covariates', 'Codes')

save(mat, file='~/Desktop/EHR-mining/UsefulDataForCleaning/ICDcodes_covariates_matrix.Rdata')








Rdata_file_path <- '~/Desktop/EHR/EHR work/RdataFiles/BuildCohortDatasets/CleanedData/'
source(file='~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/CreateNamedAbxAbbreviations.R')
load(file = paste0('~/Desktop/EHR-mining/Scripts/AnalysisScripts/BuildCohortDatasets/info/treatment/cohort_info.Rdata'))
treatment_var <- 'treatment'
library(twang)

table1s <- all_vars <- vector('list', length(cohort_info))

for (cohort in seq_along(cohort_info)) {
   
   # Cohort information
   cohort_bug <- cohort_info[[cohort]]$bug
   cohort_abx <- cohort_info[[cohort]]$abx
   cohort_name <- paste0(paste(cohort_bug, collapse='_'),
                         ifelse('switch' %in% names(cohort_info[[cohort]]), '_switch', ''), '_',
                         paste(unname(abbr[unlist(strsplit(cohort_abx, split='+', fixed=TRUE))]), collapse='_'), '_',
                         treatment_var)
   cat(cohort_name, '\n')
   names(table1s)[cohort] <- names(all_vars)[cohort] <- gsub('_treatment', '', cohort_name)
   
   # Load clean data.
   load(file = paste0('~/Desktop/EHR/EHR work/RdataFiles/BuildCohortDatasets/CleanedData/treatment/', cohort_name, '.Rdata'))
   all_vars[[cohort]] <- names(df)[!grepl('^OUT_', names(df))]
   
   # Load propensity score data.
   load(file = paste0('~/Desktop/EHR-mining/Scripts/AnalysisScripts/BuildCohortDatasets/info/treatment/propensity_score_info/', cohort_name, '.Rdata'))
   ps <- ps[[cohort_info[[cohort]]$estimand]]
   tab <- bal.table(ps)
   
   cnames <- c('tx.mn', 'ct.mn', 'std.eff.sz', 'pval')
   tab <- lapply(tab, function(t) t[,cnames])
   tab <- merge(x=tab[[1]], y=tab[[2]], by=0)
   
   cnames_mat <- matrix(c("Row.names", "Covariates", "tx.mn.x", "Treatment mean (Un)", 
                          "ct.mn.x", "Control mean (Un)", "std.eff.sz.x", "SMD (Un)", "pval.x", 
                          "P-value (Un)", "tx.mn.y", "Treatment mean (Adj)", "ct.mn.y", 
                          "Control mean (Adj)", "std.eff.sz.y", "SMD (Adj)", "pval.y", 
                          "P-value (Adj)"), ncol=2, byrow=T)
   
   for (i in seq_len(nrow(cnames_mat))) {
      colnames(tab)[colnames(tab) == cnames_mat[i,1]] <- cnames_mat[i,2]
   }
   
   w <- which(abs(tab$`SMD (Adj)`) > 0.1)
   if (length(w) > 0L) {
      tab$Covariates[w] <- paste0('**', tab$Covariates[w], '**')
   }
   
   table1s[[cohort]] <- tab
}



vars <- unique(unlist(all_vars))
var_mat <- matrix('', nrow=length(all_vars), ncol=length(vars), dimnames=list(names(all_vars), vars))
for (i in seq_along(all_vars)) {
   var_mat[i, all_vars[[i]]] <- 'X'
}
var_mat <- t(var_mat)

names(table1s)[1:5] <- colnames(var_mat)[1:5] <- gsub('([A-Z])_([a-z]+)_([A-Z]+)_([A-Z]+)', '*\\1. \\2* - \\3 vs. \\4', colnames(var_mat)[1:5])
names(table1s)[6] <- colnames(var_mat)[6] <- '*S. aureus* (MRSA) - VAN vs. switch to DAP'
names(table1s)[7] <- colnames(var_mat)[7] <- '*E. faecalis* - VAN vs. switch to AMP'
names(table1s)[8] <- colnames(var_mat)[8] <- '*E. faecium* (VRE) - DAP vs. LZD'

var_mat <- var_mat[c(
   grep('^DEM_', rownames(var_mat)),
   grep('^ABX_', rownames(var_mat)),
   grep('^BUG_', rownames(var_mat)),
   grep('^ENC_', rownames(var_mat)),
   grep('^LAB_', rownames(var_mat)),
   grep('^ICD_', rownames(var_mat))
), ]

save(var_mat, file = '~/Desktop/EHR-mining/Scripts/AnalysisScripts/BuildCohortDatasets/OtherScripts/covariates_cohort_matrix.Rdata')


sapply(X = table1s,
       FUN = function(tab) {
          w <- which(abs(tab$`SMD (Adj)`) > 0.1)
          if (length(w) > 0L) {
             length(tab$Covariates[w])
          }
       })

save(table1s, file = '~/Desktop/EHR-mining/Scripts/AnalysisScripts/BuildCohortDatasets/OtherScripts/table1s.Rdata')



abx_mat <- gsub('^.+ - (.+)$', '\\1', names(table1s))
abx_mat <- gsub('switch to ', '', abx_mat)
abx_mat <- unique(unlist(strsplit(abx_mat, split=' vs. ', fixed=TRUE)))
abx_mat <- setNames(names(abbr), abbr)[abx_mat]
abx_mat <- sapply(strsplit(abx_mat, split='/'), function(x) paste(stringr::str_to_sentence(x), collapse='-'))
abx_mat <- sort(abx_mat)
abx_mat <- matrix(c(abx_mat, names(abx_mat)), ncol=2, dimnames=list(NULL, c('Name', 'Abbreviation')))
save(abx_mat, file='~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/abx_abbr_matrix.Rdata')




source('~/Desktop/EHR-mining/UsefulDataForCleaning/handle_UPMC_site_names.R')


x <- split(names(fac2cat), fac2cat)

fac_mat <- matrix(c(names(x), sapply(x, paste, collapse=', ')), nrow=length(x), ncol=2L, dimnames = list(NULL, c('Category', 'UPMC facilities')))

save(fac_mat, file = '~/Desktop/EHR-mining/UsefulDataForCleaning/facility_category_matrix.Rdata')








load(file='~/Desktop/EHR-mining/UsefulDataForCleaning/bug_groups.Rdata')
bug_mat <- matrix(c(paste0('*', names(bug_groups), '*'), sapply(bug_groups, function(x) paste(paste0('*', x, '*'), collapse=', '))), nrow=length(bug_groups), dimnames=list(NULL, c('Group', 'Pathogen genus and species')))


# matrix(unlist(bug_groups, paste, collapse=', ')), nrow=length(bug_groups), dimnames=list(NULL, c('Group', 'Pathogen genus and species'))))

save(bug_mat, file='~/Desktop/EHR-mining/UsefulDataForCleaning/bug_group_matrix.Rdata')






source(file='~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/CreateNamedAbxAbbreviations.R')
abbr_ <- setNames(names(abbr), abbr)

abx <- list(
   'CAR' = c('MEM', 'ETP', 'DOR', 'IPM', 'MVB'),
   'OTH' = c('BAC', 'PB', 'NEO', 'NIT', 'FOF', 'RIF', 'TMP', 'EMB', 'CST', 'LEX'),
   'CFZ' = 'CFZ', 
   'ATM' = 'ATM', 
   'SXT' = 'SXT', 
   'CPT' = 'CPT', 
   'CLI' = 'CLI',
   'FLQ' = c('LVX', 'CIP', 'MXF', 'OFX', 'GAT'),
   'AMG' = c('GEN', 'TOB', 'AMK'),
   'MAC' = c('AZM', 'ERY', 'CAM'),
   'TET' = c('DOX', 'TGC', 'MIN', 'TET'),
   'BLI' = c('TZP', 'SAM', 'AMC'),
   'BLA' = c('AMP', 'OXA', 'NFC', 'AMX'),
   'GRP' = c('LZD', 'DAP'),
   'CF2' = c('CXM', 'FOX', 'CTT'),
   'CF3' = c('CRO', 'CAZ', 'POD', 'CTX'),
   'CF4' = c('FEP', 'CZA', 'CFT'),
   'GLP' = c('VAN', 'TLV')
)

abx <- sapply(abx, function(x) unname(abbr_[x]))
abx <- sapply(abx, function(x) {
   y <- strsplit(x, split='/')
   sapply(y, function(y) paste(stringr::str_to_sentence(y), collapse='-'))
})

abx_mat <- matrix(c(names(abx), unname(sapply(abx, paste, collapse=', '))), nrow=length(abx), dimnames = list(NULL, c('Group', 'Antibiotic names')))
save(abx_mat, file='~/Desktop/EHR-mining/Scripts/AnalysisScripts/BuildCohortDatasets/info/abx_groupings_matrix.Rdata')



sapply(cohort_info, '[[', 'bug')

cohort_mat <- c(
   'E. coli', 'CRO vs. TZP', 
   'K. pneumoniae', 'CRO vs. TZP', 
   'P. mirabilis', 'CRO vs. TZP',
   'P. aeruginosa', 'FEP vs. TZP',
   'S. aureus', 'OXA vs. CFZ',
   'S. aureus (MRSA)', 'VAN vs. switch to DAP',
   'E. faecalis', 'VAN vs. switch to AMP',
   'E. faecium (VRE)', 'DAP vs. LZD'
)
#cohort_mat <- 
matrix(c(cohort_mat[seq(1,length(cohort_mat),2)],
         #paste0('*', cohort_mat[seq(1,length(cohort_mat),2)], '*'),
         cohort_mat[seq(2,length(cohort_mat),2)]),
       nrow=length(cohort_mat)/2,
       dimnames = list(NULL, c('Pathogen', 'Antibiotic Treatments')))


path <- '~/Desktop/EHR-mining/Scripts/AnalysisScripts/BuildCohortDatasets/info/treatment/cat_output_df_building/'
lf <- c(
   "E_coli_CRO_TZP_treatment.txt", 
   "K_pneumoniae_CRO_TZP_treatment.txt", 
   "P_mirabilis_CRO_TZP_treatment.txt", 
   "P_aeruginosa_FEP_TZP_treatment.txt", 
   "S_aureus_OXA_CFZ_treatment.txt", 
   "S_aureus_MRSA_switch_VAN_DAP_treatment.txt", 
   "E_faecalis_switch_VAN_AMP_treatment.txt", 
   "E_faecium_VRE_DAP_LZD_treatment.txt"
)

for (i in seq_along(lf)) {
   cat(gsub('([A-Z]_[a-z]+)_.+', '\\1', lf[i]), '- ')
   cat(readLines(con=paste0(path, lf[i]))[12:13], '\n')
}






#### TABLE OF ATE / CATE ODDS RATIOS 












