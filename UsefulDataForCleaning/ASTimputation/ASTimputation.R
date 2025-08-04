# AST imputation function
imputeASTs <- function(df) {
   load(file = '~/Desktop/EHR-mining/UsefulDataForCleaning/bug_groups.Rdata')
   
   df <- df %>% filter(BUG %in% unname(unlist(bug_groups)))
   
   source('~/Desktop/EHR-mining/UsefulDataForCleaning/getAbxBugClassFxn.R')
   df <- getBugClass(df)
   
   # check for throw away rows
   w <- which(df$BUG == 'Staphylococcus aureus' & is.na(df$OXACILLIN))
   if (length(w) > 0L) {
      cat('Removed', length(w), 'missing SA rows.\n')
      df <- df[-w,]
   }
   
   enterobacterales <- bug_groups$Enterobacterales
   source('~/Desktop/EHR-mining/UsefulDataForCleaning/ASTimputation/EUCAST_expected_phenotypes.R')
   rs <- rs %>% 
      filter(BUG %in% unique(df$BUG),
             ABX %in% names(df)) %>%
      distinct()
   # let's do this differently, not just for what is missing an empiric flag...
   bug_rules <- split(rs$BUG, paste0(rs$ABX, rs$value))
   abxs <- gsub('^([A-Z/ ]+)(0|1)$', '\\1', names(bug_rules))
   labels <- as.integer(gsub('^([A-Z/ ]+)(0|1)$', '\\2', names(bug_rules)))
   for (i in seq_along(bug_rules)) {
      w <- which(df$BUG %in% bug_rules[[i]] & is.na(df[[abxs[i]]]))
      # cat(i, length(w), '\n')
      df[[abxs[i]]][w] <- labels[i]
   }
   cat('Finished EUCAST expected phenotypes.\n')
   
   
   ### EXPERT RULES ###
   
   ### Staphylococcus aureus###
   # BETA LACTAM ANTIBIOTICS
   # if MRSA, resistant to all beta-lactams (except few: CEFTAROLINE in this case)
   mrsa_R_abx <- c('AMOXICILLIN/CLAVULANATE', 'AMPICILLIN', 'AMPICILLIN/SULBACTAM', 
                   'CEFAZOLIN', 'CEFEPIME', 'CEFOXITIN', 'CEFTRIAXONE', 'CEFUROXIME', 'CEPHALEXIN',
                   'ERTAPENEM', 'MEROPENEM', 
                   'OXACILLIN', 'PIPERACILLIN/TAZOBACTAM')
   mrsa_S_abx <- c('CEFTAROLINE')
   for (abx in mrsa_R_abx) {
      # get index of MRSA that are missing SR status of antibiotic of interest
      w <- which(is.na(df[[abx]]) & df$BUG == 'Staphylococcus aureus' & df$OXACILLIN == 1L)
      df[[abx]][w] <- 1L
   }
   for (abx in mrsa_S_abx) {
      # get index of MRSA that are missing SR status of antibiotic of interest
      w <- which(is.na(df[[abx]]) & df$BUG == 'Staphylococcus aureus' & df$OXACILLIN == 1L)
      df[[abx]][w] <- 0L
   }
   
   # if MSSA, susceptible to all beta-lactams
   # if MSSA but resistance to benzylpenicillin/beta_lactamase producer, resistant to beta-lactams except isoxazolyl-penicillins and in combo with bl inhibitor
   mssa_S_abx <- c('AMOXICILLIN/CLAVULANATE', 'AMPICILLIN', 'AMPICILLIN/SULBACTAM', 'CEFAZOLIN', 'CEFEPIME', 'CEFOXITIN', 'CEFTRIAXONE', 'CEFUROXIME', 'CEFTAROLINE', 'CEPHALEXIN',
                   'ERTAPENEM', 'MEROPENEM', 'OXACILLIN', 'PIPERACILLIN/TAZOBACTAM')
   for (abx in mssa_S_abx) {
      if (abx == 'AMPICILLIN') {
         # look for beta_lactamase producers, make resistant to ampicillin if so
         w <- which(is.na(df[[abx]]) & df$BUG == 'Staphylococcus aureus' & df$OXACILLIN == 0L & (df$BETA_LACTAMASE == 1L | df$BENZYLPENICILLIN == 1L))
         df[[abx]][w] <- 1L
      }
      # get index of MRSA that are missing SR status of antibiotic of interest
      w <- which(is.na(df[[abx]]) & df$BUG == 'Staphylococcus aureus' & df$OXACILLIN == 0L)
      df[[abx]][w] <- 0L
   }
   rm(mrsa_R_abx, mrsa_S_abx, abx, w, mssa_S_abx)
   
   
   ### Staphylococcus - EUCAST rules: ###
   {
      # CIPRO, LEVO (fluoroquinolone)
      # AZITHRO (macrolide)
      # CLINDA, ERYTH (lincosamide)
      # DOXYCYCLINE (tetracycline)
      # GENTAMICIN, TOBRAMYCIN (aminoglycoside)
      
      # we don't have the norfloxacin screening test...so trouble imputing missing here
      # next, for staph spp. missing levofloxacin or cipro, make R if R to LEVO or MOXI (EUCAST told me so)
      w <- which(grepl('Staphylococcus', df$BUG) & (df$LEVOFLOXACIN == 1L | df$MOXIFLOXACIN == 1L)) # 1,499
      df$LEVOFLOXACIN[w] <- 1L
      df$CIPROFLOXACIN[w] <- 1L
      rm(w)
      
      # IF susceptible to erythromycin and clindamycin,
      # THEN report as susceptible to all macrolides and lincosamides
      # affects (AZITH)
      w <- which(grepl('Staphylococcus', df$BUG) & df$CLINDAMYCIN == 0L & df$ERYTHROMYCIN == 0L) # 4,152
      df$AZITHROMYCIN[w] <- 0L # all were missing (except 4 were already S)
      
      # DOXY 
      # if S-TET, then S-DOX, S-MINO, S-TIGE
      # if R-TET, then either R-DOX, R-MINO    OR    report DOX and MINO individually
      w <- which(grepl('Staphylococcus', df$BUG) & is.na(df$DOXYCYCLINE) & df$TETRACYCLINE == 0L) # 10,813 missing (53 need it (received DOXY empirically))
      df$DOXYCYCLINE[w] <- 0L
      # df %>% filter(grepl('Staphylococcus', BUG)) %>% count(TETRACYCLINE, DOXYCYCLINE)
      w <- which(grepl('Staphylococcus', df$BUG) & df$TETRACYCLINE == 1L) # 873 (table(df$DOXYCYCLINE[w])) - some DOX-S
      df$DOXYCYCLINE[w] <- 1L
      
      # GENTAMICIN, TOBRAMYCIN
      # don't have any other aminoglycoside susceptibility data...
      # df %>% filter(grepl('Staphylococcus', BUG), ABX == 'GENTAMICIN') %>% count(GENTAMICIN) # 9 missing
      # df %>% filter(grepl('Staphylococcus', BUG), ABX == 'GENTAMICIN', is.na(GENTAMICIN)) %>% count(STREPTOMYCIN_SYNERGY, STREPTOMYCIN_HIGH_LEVEL, STREPTOMYCIN, NEOMYCIN, KANAMYCIN, AMIKACIN, TOBRAMYCIN, GENTAMICIN, GENTAMICIN_HIGH_LEVEL, GENTAMICIN_SYNERGY)
      # df %>% filter(grepl('Staphylococcus', BUG), ABX == 'TOBRAMYCIN') %>% count(TOBRAMYCIN, GENTAMICIN, STREPTOMYCIN_SYNERGY, STREPTOMYCIN_HIGH_LEVEL, STREPTOMYCIN, NEOMYCIN, KANAMYCIN, AMIKACIN, TOBRAMYCIN, GENTAMICIN, GENTAMICIN_HIGH_LEVEL, GENTAMICIN_SYNERGY) # all 14 missing
      
      # NITROFURANTION
      # df %>% filter(grepl('Staphylococcus', BUG), ABX == 'NITROFURANTOIN') %>% count(NITROFURANTOIN) # all 8 missing
      
      # BACITRACIN
      # df %>% filter(grepl('Staphylococcus', BUG), ABX == 'BACITRACIN') %>% count(BACITRACIN) # all 47 missing
      
      
      # STAPHYLOCOCCUS LUGDUNENSIS
      # Beta-lactams + AZITHRO + CIPRO + DOXY + ERTAPENEM + ERYTHRO
      # df[unlist(wmiss_list),] %>% filter(BUG == 'Staphylococcus lugdunensis') %>% count(ABX, sort=TRUE)
      # df %>% filter(BUG == 'Staphylococcus lugdunensis') %>% count(OXACILLIN) # 35 (out of 250 have mecA?)
      
      mrsl_R_abx <- c('AMOXICILLIN/CLAVULANATE', 'AMPICILLIN', 'AMPICILLIN/SULBACTAM', 'CEFOXITIN',
                      'CEFAZOLIN', 'CEFEPIME', 'CEFOXITIN', 'CEFTRIAXONE', 'CEFUROXIME', 'CEPHALEXIN',
                      'ERTAPENEM', 'MEROPENEM', 'OXACILLIN', 'PIPERACILLIN/TAZOBACTAM')
      for (abx in mrsl_R_abx) {
         # get index of MRSL that are missing SR status of antibiotic of interest
         w <- which(is.na(df[[abx]]) & df$BUG == 'Staphylococcus lugdunensis' & df$OXACILLIN == 1L)
         df[[abx]][w] <- 1L
      }
      for (abx in mrsl_R_abx) {
         # get index of MSSL that are missing SR status of antibiotic of interest
         w <- which(is.na(df[[abx]]) & df$BUG == 'Staphylococcus lugdunensis' & df$OXACILLIN == 0L)# & df$ABX == abx)
         df[[abx]][w] <- 0L
      }
      rm(abx, w, mrsl_R_abx)
   }
   
   cat('Finished EUCAST expert rules for:\n')
   cat('\tStaphylococci\n')
   
   
   ##### STREPTOCOCCI #####
   # pneumococcus - beta-lactams
   # if PEN-S, report as susceptible to beta-lactams?
   # if PEN-R, report as R to all beta-lactams?
   # if missing PEN, report as suscepible?
   w <- which(df$BUG %in% c('Streptococcus pneumoniae', 'Streptococcus agalactiae', 'Streptococcus pyogenes') & is.na(df$PENICILLIN)) # 492
   df$PENICILLIN[w] <- 0L
   strep_penS_expS <- c("AMPICILLIN", "AMPICILLIN/SULBACTAM", "PIPERACILLIN/TAZOBACTAM", 
                        "AZTREONAM", "CEFAZOLIN", "CEPHALEXIN", "CEFUROXIME", "CEFOXITIN", 
                        "CEFTRIAXONE", "CEFEPIME", "CEFTAROLINE", 'ERTAPENEM', 'MEROPENEM')
   for (abx in strep_penS_expS) {
      w <- which(df$BUG %in% c('Streptococcus pneumoniae', 'Streptococcus agalactiae', 'Streptococcus pyogenes') & df$PENICILLIN == 0L & is.na(df[[abx]]))
      df[[abx]][w] <- 0L
   }
   rm(w, abx, strep_penS_expS)
   
   # missing CIPRO, would be R if LEVO-R or MOXI-R, but none
   # df[unlist(wmiss_list),] %>% filter(BUG == 'Streptococcus pneumoniae', ABX == 'CIPROFLOXACIN') %>% count(CIPROFLOXACIN, LEVOFLOXACIN, MOXIFLOXACIN) # 15 rows
   
   # DOX-S if TET-S
   w <- which(df$BUG == 'Streptococcus pneumoniae' & df$TETRACYCLINE == 0L) # 1093
   df$DOXYCYCLINE[w] <- 0L
   w <- which(df$BUG == 'Streptococcus pneumoniae' & df$TETRACYCLINE == 1L) # 217
   df$DOXYCYCLINE[w] <- 1L
   
   # STILL MISSING: AZTREONAM, AMP, AMP/SUL, AMOX/CLAV, PIP/TAZ, CEFEPIME, CEFTRIAXONE, CEFAZOLIN, CEPHALEXIN, AZITH, CIPRO, DOXY, GENT
   # df[unlist(wmiss_list),] %>% 
   #    filter(BUG %in% c('Streptococcus pneumoniae', 'Streptococcus agalactiae', 'Streptococcus pyogenes'), 
   #           PENICILLIN == 1L) %>% # 192
   #    count(ABX, sort=TRUE)
   
   cat('\tStreptococci\n')
   
   
   
   ##### ENTEROBACTERIACEAE #####
   # miss_bugs[enterobacterales, colSums(miss_bugs[enterobacterales,]) != 0]
   
   df$ESBL <- NA
   df$ESBL[df$BUG %in% enterobacterales] <- 0L
   cro_res <- !is.na(df$CEFTRIAXONE) & df$CEFTRIAXONE == 1L
   carb_tested <- !is.na(df$ERTAPENEM) | !is.na(df$MEROPENEM) | !is.na(df$IMIPENEM) | !is.na(df$DORIPENEM)
   bl_notrep <- is.na(df$CEFTRIAXONE) & is.na(df$`PIPERACILLIN/TAZOBACTAM`) & is.na(df$CEFEPIME)
   # bug is enterobacterales AND ( CRO-R | (CARB tested and BL not reported))
   w <- which(df$BUG %in% enterobacterales & (cro_res | (carb_tested & bl_notrep)))
   df$ESBL[w] <- 1L
   rm(cro_res, carb_tested, bl_notrep)
   
   
   # ESBL producing?
   # resistant to β-lactam antibiotics including penicillins, third-generation cephalosporins, and monobactams (aztreonam)
   # Ceftazidime, Ceftriaxone, Cefotaxime, Cefdinir, Cefpodoxime, Cefixime, Ceftizoxime, Cefoperazone
   esbl_R_set <- c('AMPICILLIN', 'AMPICILLIN/SULBACTAM', 'AMOXICILLIN/CLAVULANATE', 'PIPERACILLIN/TAZOBACTAM', 
                   'AZTREONAM', 
                   'CEFAZOLIN', 'CEPHALEXIN', 'CEFUROXIME', 'CEFOTAXIME', 'CEFTAZIDIME', 'CEFTRIAXONE',
                   'CEFOXITIN', 'CEFTAROLINE', 'CEFEPIME')
   
   for (abx in esbl_R_set) {
      # find ESBLs, impute these as resistant if missing
      w <- which(df$BUG %in% enterobacterales & df$ESBL == 1L & is.na(df[[abx]]))
      df[[abx]][w] <- 1L
      
      # find non-ESBLs, impute these as susceptible if missing
      w <- which(df$BUG %in% enterobacterales & df$ESBL == 0L & is.na(df[[abx]]))
      df[[abx]][w] <- 0L
   }
   
   
   ### CARBAPENEMS
   # surely if AST for carbapenems are missing for enterobacterales
   # they must be susceptible, RIGHT??
   w <- which(df$BUG %in% enterobacterales & is.na(df$MEROPENEM) & (is.na(df$ERTAPENEM) | df$ERTAPENEM == 0L) & (is.na(df$IMIPENEM) | df$IMIPENEM == 0L) & (is.na(df$DORIPENEM) | df$DORIPENEM == 0L))
   df$MEROPENEM[w] <- 0L # 12,477
   w <- which(df$BUG %in% enterobacterales & is.na(df$ERTAPENEM) & (is.na(df$MEROPENEM) | df$MEROPENEM == 0L) & (is.na(df$IMIPENEM) | df$IMIPENEM == 0L) & (is.na(df$DORIPENEM) | df$DORIPENEM == 0L))
   df$ERTAPENEM[w] <- 0L  # 6,629
   
   w <- which(df$BUG %in% enterobacterales & is.na(df$MEROPENEM) & ((!is.na(df$ERTAPENEM) & df$ERTAPENEM == 1L) | (!is.na(df$IMIPENEM) & df$IMIPENEM == 1L) | (!is.na(df$DORIPENEM) & df$DORIPENEM == 1L)))
   df$MEROPENEM[w] <- 1L # 180ish
   w <- which(df$BUG %in% enterobacterales & is.na(df$ERTAPENEM) & ((!is.na(df$MEROPENEM) & df$MEROPENEM == 1L) | (!is.na(df$IMIPENEM) & df$IMIPENEM == 1L) | (!is.na(df$DORIPENEM) & df$DORIPENEM == 1L)))
   df$ERTAPENEM[w] <- 1L  # 76
   
   
   cat('\tEnterobacterales\n')
   
   
   
   ##### ENTEROCOCCUS #####
   enterococci <- c('Enterococcus faecium', 'Enterococcus faecalis')
   # miss_bugs[enterococci,] # missing amp/sul, pip/tazo, carbapenems, doxy, fluoros, bacitracin
   
   # Beta-lactams - EUCAST expert rules are confusing
   # if missing ampicillin, assume resistant if E. faecium, susceptible if E. faecalis
   w <- which(df$BUG == 'Enterococcus faecium' & is.na(df$AMPICILLIN)) # 33
   df$AMPICILLIN[w] <- 1L
   w <- which(df$BUG == 'Enterococcus faecalis' & is.na(df$AMPICILLIN)) # 143
   df$AMPICILLIN[w] <- 0L
   
   # can be inferred from ampicillin - according to EUCAST breakpoints document
   abx_set <- c('AMOXICILLIN/CLAVULANATE', 'AMPICILLIN/SULBACTAM', 'PIPERACILLIN/TAZOBACTAM')
   for (abx in abx_set) {
      w <- which(df$BUG %in% enterococci & is.na(df[[abx]]))
      df[[abx]][w] <- df$AMPICILLIN[w]
   }
   
   
   # E faecalis/faecium - if AMP-R, then ureidopen-R, imipenem-R
   ampR_expR <- c('AZLOCILLIN', 'PIPERACILLIN', 'MEZLOCILLIN', 'IMIPENEM')#, 'PIPERACILLIN/TAZOBACTAM')
   for (abx in ampR_expR) {
      w <- which(df$BUG %in% enterococci & is.na(df[[abx]]) & df$AMPICILLIN == 1L)
      df[[abx]][w] <- 1L
   }
   rm(ampR_expR, abx, w)
   
   # E faecalis ONLY!!! - if AMP-S, mostly TZP-S, PIP-S, AMOX-S
   ampS_expS <- c('AMOXICILLIN', 'PIPERACILLIN')#, 'PIPERACILLIN/TAZOBACTAM')
   for (abx in ampS_expS) {
      w <- which(df$BUG == 'Enterococcus faecalis' & is.na(df[[abx]]) & df$AMPICILLIN == 0L)
      df[[abx]][w] <- 0L
   }
   rm(ampS_expS, abx, w)
   
   # Enterococcus spp. (including E. faecium), susceptibility to these agents (above) is uncommon and 
   # isolates resistant to ampicillin should not be reported susceptible to either amoxicillin or piperacillin (with and without inhibitor)
   # E faecium ONLY!!! - if AMP-R, report TZP-R !!
   ampR_expR <- c('AMOXICILLIN', 'PIPERACILLIN')#, 'PIPERACILLIN/TAZOBACTAM')
   for (abx in ampR_expR) {
      w <- which(df$BUG == 'Enterococcus faecium' & is.na(df[[abx]]) & df$AMPICILLIN == 1L)
      df[[abx]][w] <- 1L
   }
   rm(ampR_expR, abx, w)
   
   # df[unlist(wmiss_list),] %>%
   #    filter(BUG %in% enterococci) %>%
   #    count(`PIPERACILLIN/TAZOBACTAM`, BUG)
   
   # doxycycline appears to be assumed resistant if not measured
   w <- which(df$BUG %in% enterococci & is.na(df$DOXYCYCLINE))
   df$DOXYCYCLINE[w] <- 1L
   
   # still need cipro and levo
   cat('\tEnterococci\n')
   
   return(distinct(df))
}

























