bug_groups <- list(
   Staphylococci = c("Staphylococcus aureus", "Staphylococcus lugdunensis"),
   
   Enterococci = c("Enterococcus faecalis", "Enterococcus faecium"),
   
   Streptococci = c("Streptococcus agalactiae", "Streptococcus pneumoniae", "Streptococcus pyogenes"),
   
   Enterobacterales = c("Escherichia coli", 
                        "Enterobacter aerogenes", "Enterobacter cloacae",
                        "Citrobacter freundii", 
                        "Klebsiella aerogenes", "Klebsiella oxytoca", "Klebsiella pneumoniae", "Klebsiella variicola", 
                        "Morganella morganii",
                        "Proteus mirabilis",
                        "Serratia marcescens"),
   
   NonFermGN = c("Pseudomonas aeruginosa", "Stenotrophomonas maltophilia", "Acinetobacter baumannii")
)
save(bug_groups, file = '~/Desktop/EHR-mining/UsefulDataForCleaning/bug_groups.Rdata')

# function to get antibiotic class name and pathogen group names
getBugClass <- function(df) {
   bug_groups_named <- setNames(
      gsub('[0-9]+$', '', names(unlist(bug_groups))),
      unname(unlist(bug_groups))
   )
   
   df <- df %>% 
      mutate(BUGC = unname(bug_groups_named[BUG])) %>%
      relocate(BUGC, .after=BUG)
   
   return(df)
}
abx_class <- list(
   penicillin = c('AMPICILLIN', 'OXACILLIN', 'NAFCILLIN', 'AMOXICILLIN', 'PENICILLIN V', 'BENZYLPENICILLIN', 'DICLOXACILLIN'),
   penicillin_blinhibitor = c('AMOXICILLIN/CLAVULANATE', 'AMPICILLIN/SULBACTAM', 'PIPERACILLIN/TAZOBACTAM'),
   monobactam = c('AZTREONAM'),
   cephalosporin_1g = c('CEFAZOLIN', 'CEPHALEXIN', 'CEFADROXIL'),
   cephalosporin_2g = c('CEFUROXIME', 'CEFOXITIN', 'CEFOTETAN', 'CEFPROZIL'),
   cephalosporin_3g = c('CEFTRIAXONE', 'CEFTAZIDIME', 'CEFOTAXIME', 'CEFPODOXIME', 'CEFDINIR', 'CEFIDEROCOL', 'CEFIXIME'),
   cephalosporin_4g = c('CEFEPIME'),
   cephalosporin_5g = c('CEFTAROLINE'),
   cephalosporin_blinhibitor = c('CEFTOLOZANE/TAZOBACTAM', 'CEFTAZIDIME/AVIBACTAM'),
   carbapenem = c('ERTAPENEM',  'MEROPENEM', 'IMIPENEM', 'DORIPENEM'),
   carbapenems_blinhibitor = c('MEROPENEM/VABORBACTAM', 'IMIPENEM/RELEBACTAM'),
   sulfonamide_dihydrofolate = c('TRIMETHOPRIM/SULFAMETHOXAZOLE'), 
   dihydrofolate = c('TRIMETHOPRIM'),
   macrolide = c('AZITHROMYCIN', 'ERYTHROMYCIN', 'CLARITHROMYCIN'),
   lincosamide = c('CLINDAMYCIN'),
   aminoglycoside = c('GENTAMICIN', 'TOBRAMYCIN', 'NEOMYCIN', 'AMIKACIN'),
   tetracycline = c('DOXYCYCLINE', 'TIGECYCLINE', 'TETRACYCLINE', 'MINOCYCLINE', 'ERAVACYCLINE', 'DEMECLOCYCLINE'),
   fluoroquinolone = c('CIPROFLOXACIN', 'LEVOFLOXACIN', 'MOXIFLOXACIN', 'GATIFLOXACIN', 'OFLOXACIN', 'DELAFLOXACIN'),
   glycopeptide = c('VANCOMYCIN', 'TELAVANCIN', 'DALBAVANCIN', 'ORITAVANCIN'), 
   oxazolidinone = c('LINEZOLID'), 
   lipopeptide = c('DAPTOMYCIN'), 
   polypeptide = c('BACITRACIN'), 
   nitrofuran = c('NITROFURANTOIN'),
   polymyxin = c('POLYMYXIN B', 'COLISTIN'),
   phosphonic = c('FOSFOMYCIN'),
   ansamycin = c('RIFAMPIN'),
   streptogramin = c('SYNERCID', 'PRISTINAMYCIN', 'VIRGINIAMYCIN'),
   tiacumicin = c('FIDAXOMICIN'),
   antiTB = c('ETHAMBUTOL', 'ISONIAZID', 'RIFAMPIN', 'RIFABUTIN', 'PYRAZINAMIDE'),
   antiFungal = c('CASPOFUNGIN', 'FLUCONAZOLE', 'METRONIDAZOLE', 'VORICONAZOLE', 'POSACONAZOLE', 'CLOTRIMAZOLE', 'MICAFUNGIN', 
                   'KETOCONAZOLE', 'AMPHOTERICIN', 'TERBINAFINE', 'ITRACONAZOLE', 'FLUCYTOSINE', 'GRISEOFULVIN', 'NATAMYCIN')
)
save(abx_class, file = '~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/abx_classes.Rdata')

getAbxClass <- function(df) {
   abx_class_named <- setNames(gsub(pattern = '[0-9]+$', 
                                    replacement = '', 
                                    x = names(unlist(abx_class))),
                               unname(unlist(abx_class)))
   df <- df %>% 
      mutate(ABXC = unname(abx_class_named[ABX])) %>%
      relocate(ABX, ABXC, .after=ORDER_DAY)
   return(df)
}

getAbxClassResistanceFlags <- function(df) {
   # check if resistance to any antibiotic in a class
   for (i in seq_along(abx_class)) {
      class <- names(abx_class)[i]
      cat(class, i, '\n')
      abxs <- abx_class[[i]]
      abxs <- abxs[abxs %in% names(df)]
      m <- as.matrix(df[abxs])
      m[is.na(m)] <- 0
      m <- apply(m, 1, function(r) as.integer(any(r == 1L)))
      df[paste0('ResC_', class)] <- m
   }
   return(df)
}











