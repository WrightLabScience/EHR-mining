# function to get antibiotic class name and pathogen group names
getBugClass <- function(df) {
   bug_class <- list(
      Enterobacterales = c('Escherichia coli', 'Klebsiella pneumoniae', 'Proteus mirabilis', 
                           'Enterobacter cloacae', 'Klebsiella oxytoca', 'Serratia marcescens', 
                           'Enterobacter aerogenes',
                           'Klebsiella aerogenes', 'Klebsiella variicola', 'Citrobacter freundii', 'Morganella morganii'),
      Staphylococci = grep('Staph', unique(df$BUG), value=TRUE),
      Streptococci = grep('Strep', unique(df$BUG), value=TRUE),
      Enterococci = grep('Enterococcus', unique(df$BUG), value=TRUE),
      NonFermGN = c('Acinetobacter baumannii', 'Pseudomonas aeruginosa', 'Stenotrophomonas maltophilia')
   )
   bug_class_named <- setNames(
      gsub('[0-9]+$', '', names(unlist(bug_class))),
      unname(unlist(bug_class))
   )
   
   df <- df %>% 
      mutate(BUGC = unname(bug_class_named[BUG])) %>%
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
   fluoroquinolone = c('CIPROFLOXACIN', 'LEVOFLOXACIN', 'MOXIFLOXACIN', 'GATIFLOXACIN', 'OFLOXACIN'),
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











