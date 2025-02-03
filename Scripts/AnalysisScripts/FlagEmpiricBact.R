# How do we decide what to predict? Whatever was prescribed empirically!
# What is the baseline performance? Whatever was prescribed empirically!
# Heuristic: if a physician prescribed it empirically, they predict Susceptible.
# Can we get the physician score?
library(dplyr)
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_clean_ASTs.Rdata')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_CLEANED_2015_2024_AbxAdmin.Rdata')
source(file = '~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/CreateNamedAbxAbbreviations.R')

# handle ABX administration data
# remove antifungals
abxDF <- abxDF %>%
   select(-START_DAY, -END_DATE) %>%
   filter(!grepl(',', ABX)) %>%
   filter(!ABX %in% c('CASPOFUNGIN', 'FLUCONAZOLE', 'METRONIDAZOLE', 'VORICONAZOLE', 'POSACONAZOLE', 'CLOTRIMAZOLE', 
                      'MICAFUNGIN', 'KETOCONAZOLE', 'AMPHOTERICIN', 'TERBINAFINE', 'ITRACONAZOLE', 'FLUCYTOSINE')) %>%
   distinct()
w <- which(abxDF$ABX == 'PENICILLIN G')
abxDF$ABX[w] <- 'BENZYLPENICILLIN'
rm(w)
gc()


# handle ASTs
# take only blood with certain bugs and time ranges
# make sure that each index is actaully an index with respect to:
#     other blood cultures AND non-blood cultures by 30 days
astDF_non_blood <- astDF %>% filter(!BLOOD) %>% select(PERSON_ID, RESULT_DAY) %>% distinct()
astDF <- astDF %>%
   filter(BLOOD) %>% 
   filter(lubridate::year(ORDER_DAY) %in% 2015:2024) %>%
   filter(!(grepl('Staph', BUG) & !grepl('Staphylococcus (aureus|lugdunensis)', BUG))) %>%
   filter(!(grepl('Streptococc', BUG) & !grepl('Streptococcus (agalactiae|pyogenes|pneumoniae)', BUG))) %>%
   filter(!grepl('Cryptococcus|Aspergillus|Candida', BUG)) %>%
   filter(BUG != 'Did not match')
astDF_non_blood <- astDF_non_blood %>% filter(PERSON_ID %in% unique(astDF$PERSON_ID))
common_bugs <- sort(table(astDF$BUG), decreasing = TRUE)
common_bugs <- names(common_bugs)[common_bugs >= 200L]
common_bugs <- c(
   "Staphylococcus aureus", "Staphylococcus lugdunensis", 
   
   "Enterococcus faecalis", "Enterococcus faecium",
   
   "Streptococcus agalactiae", "Streptococcus pneumoniae", "Streptococcus pyogenes",
   
   "Escherichia coli", 
   "Enterobacter aerogenes", "Enterobacter cloacae",
   "Citrobacter freundii", 
   "Klebsiella aerogenes", "Klebsiella oxytoca", "Klebsiella pneumoniae", "Klebsiella variicola", 
   "Morganella morganii",
   "Proteus mirabilis",
   "Serratia marcescens", 
   
   "Pseudomonas aeruginosa", "Stenotrophomonas maltophilia", "Acinetobacter baumannii"
)
astDF <- astDF %>% filter(BUG %in% common_bugs)
ast <- astDF %>%
   select(PERSON_ID, ORDER_DATE, ORDER_DAY) %>%
   distinct() %>%
   group_by(PERSON_ID, ORDER_DAY) %>%
   slice_min(ORDER_DATE) %>%
   ungroup()
# index blood culture by 30 days
ast <- ast %>%
   group_by(PERSON_ID) %>%
   filter(is.na(lag(ORDER_DAY)) | as.integer(ORDER_DAY - lag(ORDER_DAY)) >= 14L) %>%
   ungroup()
# index any culture by 30 days
# ast <- ast %>%
#    mutate(D14_BF_ORDER_DAY = ORDER_DAY - 1L) %>%
#    left_join(
#       x = .,
#       y = astDF_non_blood,
#       by = join_by(
#          PERSON_ID,
#          ORDER_DAY >= RESULT_DAY,
#          D14_BF_ORDER_DAY <= RESULT_DAY
#       )
#    ) %>%
#    filter(is.na(RESULT_DAY)) %>%
#    select(-D14_BF_ORDER_DAY, -RESULT_DAY)





# join AST + ABX to get just empiric abx
df <- ast %>%
   mutate( # empiric prescription time frame -48 hours to +12 hours
      JOIN_START = ORDER_DATE - 86400 * 2,
      JOIN_END = ORDER_DATE + 86400 * 0.5
   ) %>%
   left_join(
      x = .,
      y = abxDF,
      by = join_by(
         PERSON_ID,
         between(y$START_DATE, x$JOIN_START, x$JOIN_END)
      )
   ) %>%
   select(-JOIN_START, -JOIN_END, -START_DATE, -ORDER_DATE) %>% 
   distinct() %>%
   filter(!is.na(ABX))
df <- df %>%
   left_join(
      x = .,
      y = astDF %>%
         select(PERSON_ID, ORDER_DAY, BUG, CEFEPIME:ESBL) %>%
         distinct(),
      relationship = 'many-to-many',
      by = join_by(
         PERSON_ID,
         ORDER_DAY
      )
   )
# df_og <- df


# combine antibiotics into classes:
abx_class <- list(
   penicillin = c('AMPICILLIN', 'OXACILLIN', 'NAFCILLIN', 'AMOXICILLIN', 'PENICILLIN V', 'BENZYLPENICILLIN', 'DICLOXACILLIN'),
   penicillin_blinhibitor = c('AMOXICILLIN/CLAVULANATE', 'AMPICILLIN/SULBACTAM', 'PIPERACILLIN/TAZOBACTAM'),
   monobactam = c('AZTREONAM'),
   cephalosporin_1g = c('CEFAZOLIN', 'CEPHALEXIN', 'CEFADROXIL'),
   cephalosporin_2g = c('CEFUROXIME', 'CEFOXITIN', 'CEFOTETAN'),
   cephalosporin_3g = c('CEFTRIAXONE', 'CEFTAZIDIME', 'CEFOTAXIME', 'CEFPODOXIME', 'CEFDINIR', 'CEFIDEROCOL'),
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
   tetracycline = c('DOXYCYCLINE', 'TIGECYCLINE', 'TETRACYCLINE', 'MINOCYCLINE'),
   fluoroquinolone = c('CIPROFLOXACIN', 'LEVOFLOXACIN', 'MOXIFLOXACIN', 'GATIFLOXACIN', 'OFLOXACIN'),
   glycopeptide = c('VANCOMYCIN', 'TELAVANCIN', 'DALBAVANCIN'), 
   oxazolidinone = c('LINEZOLID'), 
   lipopeptide = c('DAPTOMYCIN'), 
   polypeptide = c('BACITRACIN'), 
   nitrofuran = c('NITROFURANTOIN'),
   polymyxin = c('POLYMYXIN B', 'COLISTIN'),
   phosphonic = c('FOSFOMYCIN'),
   ansamycin = c('RIFAMPIN'),
   tiacumicin = c('FIDAXOMICIN'),
   antiTB = c('ETHAMBUTOL', 'ISONIAZID', 'RIFAMPIN')
)
abx_class_named <- setNames(gsub(pattern = '[0-9]+$', 
                                 replacement = '', 
                                 x = names(unlist(abx_class))),
                            unname(unlist(abx_class)))
bug_class <- list(
   Enterobacterales = enterobacterales,
   Staphylococci = grep('Staph', unique(df$BUG), value=TRUE),
   Streptococci = grep('Strep', unique(df$BUG), value=TRUE),
   Enterococci = grep('Enterococcus', unique(df$BUG), value=TRUE),
   NonFermGN = c('Acinetobacter baumannii', 'Pseudomonas aeruginosa', 'Stenotrophomonas maltophilia')
)
bug_class_named <- setNames(
   gsub('[0-9]+$', '', names(unlist(bug_class))),
   unname(unlist(bug_class))
)

df <- df_og

df <- df %>% 
   mutate(
      ABXC = unname(abx_class_named[ABX]),
      BUGC = unname(bug_class_named[BUG])
   ) %>%
   relocate(ABXC, .after=ABX) %>%
   relocate(BUGC, .after=BUG)





# narrow down to cultures that were treated empirically with relatively common abx
# - at least 100 infections were treated empirically with this
common_abx <- df %>% select(PERSON_ID, ORDER_DAY, ABX) %>% distinct() %>% count(ABX, sort=TRUE) %>% filter(n >= 100L) %>% pull(ABX)
common_abx <- common_abx[!common_abx %in% c('BACITRACIN', 'POLYMYXIN B')]
df <- df %>%
   #group_by(PERSON_ID, ORDER_DAY) %>% # 37,743 events
   filter(ABX %in% common_abx)# %>% # 36,954 events
#ungroup()

# join with ASTs again to get BUGs

df[common_abx[!common_abx %in% names(astDF)]] <- NA_real_

# first, for staph aureus missing oxacillin AST, remove
w <- which(df$BUG == 'Staphylococcus aureus' & is.na(df$OXACILLIN)) # 20
df <- df[-w,]
rm(w)

# antibiotics for which there are no AST results for any isolate
# if they don't appear in the common_abx set, remove them as columns
# no_ast_abx <- names(which(sapply(df %>% select(!c(PERSON_ID, ORDER_DAY, BUG)), function(x) sum(!is.na(x)) == 0L)))
# no_ast_abx <- no_ast_abx[!no_ast_abx %in% common_abx]
# df <- df %>% select(-!!no_ast_abx)
# rm(no_ast_abx)



##### BEGIN JOINT DISCORDANCE + IMPUTATION PROCESS #####

# # we're going to compress to unique ABX - BUG combos for imputation
# # but don't want to lose track of the original data, including ID and DATE
# group_rows_list <- df %>% group_by_all() %>% ungroup(PERSON_ID, ORDER_DAY) %>% group_rows()
# group_rows_idx <- rep(NA, nrow(df))
# for (i in seq_along(group_rows_list)) group_rows_idx[group_rows_list[[i]]] <- i
# rm(i)
# df$GROUP <- group_rows_idx
# df <- df %>% relocate(GROUP, .after=ORDER_DAY)




# before we bring in AST results to flag as concordant/discordant
# for which we must impute,
# how many unique ABX - BUG combos are there to worry about?
# what about unique ABX - BUG - AST results combos?
df %>% count(BUG, ABX, sort=TRUE) # 518

df$FLAG <- NA


# how many that need to be tested are missing?
# each row has an ABX prescribed, wcols contains the corresponding ABX AST result column
wcols <- match(df$ABX, names(df))
# how to access this smattering of cells in the matrix quickly?
# go column by column, using tapply?
wrows_list <- tapply(X = seq_along(wcols),
                     INDEX = wcols,
                     FUN = list)
rm(wcols)
names(wrows_list) <- unname(sapply(wrows_list, function(x) df$ABX[x[1]]))
wrows_list <- wrows_list[c('AMPICILLIN', # 'OXACILLIN', # penicillin
                           'AMOXICILLIN/CLAVULANATE', 'AMPICILLIN/SULBACTAM', 'PIPERACILLIN/TAZOBACTAM', # penicillin + beta-lactamase inhibitor
                           'AZTREONAM', # monobactams
                           'CEFAZOLIN', 'CEPHALEXIN', # 1st gen
                           'CEFUROXIME', 'CEFOXITIN', # 2nd gen 
                           'CEFTRIAXONE',  # 3rd gen
                           'CEFEPIME', # 4th gen
                           'CEFTAROLINE', # 5th gen
                           'ERTAPENEM',  'MEROPENEM', # carbapenem
                           'TRIMETHOPRIM/SULFAMETHOXAZOLE', 
                           'AZITHROMYCIN', 'ERYTHROMYCIN', # macrolide
                           'CLINDAMYCIN', # lincosamide
                           'GENTAMICIN', 'TOBRAMYCIN', # aminoglycoside
                           'DOXYCYCLINE', # tetraycline
                           'CIPROFLOXACIN', 'LEVOFLOXACIN', # fluoroquinolone
                           'VANCOMYCIN', 'LINEZOLID', 'DAPTOMYCIN')]
#'BACITRACIN', # 'NITROFURANTOIN', 
# 'POLYMYXIN B')]
# lapply(wrows_list, function(x) df$ABX[x[1]]) # just to confirm

makeFlag <- function(vec) {
   stopifnot(length(vec) > 0L)
   vec_copy <- as.character(vec)
   if (any(is.na(vec)))
      vec_copy[is.na(vec)] <- 'Unknown'
   if (any(vec_copy == '0'))
      vec_copy[vec_copy == '0'] <- 'Concordant'
   if (any(vec_copy == '1'))
      vec_copy[vec_copy == '1'] <- 'Discordant'
   return(vec_copy)
}



# calculate how many for each drug are missing
wmiss_list <- setNames(vector('list', length(wrows_list)),
                       names(wrows_list))
for (i in seq_along(wrows_list)) {
   abx <- names(wrows_list)[i]
   w <- wrows_list[[i]]
   df$FLAG[w] <- makeFlag(df[[abx]][w])
   wmiss_list[[i]] <- w[which(df$FLAG[w] == 'Unknown')]
}
rm(i, abx)
length(unlist(wmiss_list)) # 39,158

par(mar=c(5,4,1,1))
x <- sort(lengths(wmiss_list) / lengths(wrows_list))
b <- barplot(x, horiz=TRUE, names.arg = NA, xlim=c(0,1))
text(x=-0.01, y=b, adj=c(1, 0.5), labels=abbr[names(x)], xpd=NA)
rm(b, x)



# IMPUTE using EUCAST expected phenotypes 
# to get the easy ones out of the way
source('~/Desktop/EHR-mining/UsefulDataForCleaning/ASTimputation/EUCAST_expected_phenotypes.R')
rs <- rs %>% 
   filter(BUG %in% unique(df$BUG),
          ABX %in% unique(df$ABX))
# loop over each abx
# find which rows in df are missing for that abx
# get the imputation rules, apply them
for (i in seq_along(wmiss_list)) {
   rows <- wmiss_list[[i]]
   abx <- names(wmiss_list)[i]
   cat(i, abx, length(rows), '\n')
   
   rules <- rs %>% filter(ABX == abx) %>% select(-ABX)
   if (nrow(rules) == 0L) {
      next
   }
   rules <- setNames(rules$value, rules$BUG)
   df[[abx]][rows] <- unname(rules[df$BUG[rows]])
}
rm(i, rows, abx, rules, rs)


## UPDATE MISSING ##
wmiss_list <- setNames(vector('list', length(wrows_list)),
                       names(wrows_list))
for (i in seq_along(wrows_list)) {
   abx <- names(wrows_list)[i]
   w <- wrows_list[[i]]
   df$FLAG[w] <- makeFlag(df[[abx]][w])
   wmiss_list[[i]] <- w[which(df$FLAG[w] == 'Unknown')]
}
rm(i, abx)
length(unlist(wmiss_list)) # 21,579

par(mar=c(5,4,1,1))
x <- sort(lengths(wmiss_list) / lengths(wrows_list))
b <- barplot(x, horiz=TRUE, names.arg = NA, xlim=c(0,1))
text(x=-0.01, y=b, adj=c(1, 0.5), labels=abbr[names(x)], xpd=NA)
rm(b, x)




### START IMPUTING USING EXPERT RULES ###

### Staphylococcus aureus ###

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




## UPDATE missing ##
wmiss_list <- setNames(vector('list', length(wrows_list)),
                       names(wrows_list))
for (i in seq_along(wrows_list)) {
   abx <- names(wrows_list)[i]
   w <- wrows_list[[i]]
   df$FLAG[w] <- makeFlag(df[[abx]][w])
   wmiss_list[[i]] <- w[which(df$FLAG[w] == 'Unknown')]
}
rm(i, abx)
length(unlist(wmiss_list)) # 11,231
par(mar=c(5,4,1,1))
x <- sort(lengths(wmiss_list) / lengths(wrows_list))
b <- barplot(x, horiz=TRUE, names.arg = NA, xlim=c(0,1))
text(x=-0.01, y=b, adj=c(1, 0.5), labels=abbr[names(x)], xpd=NA)
rm(b, x)




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
   df %>% filter(grepl('Staphylococcus', BUG), ABX %in% c('LEVOFLOXACIN', 'CIPROFLOXACIN')) %>% count(NORFLOXACIN) # all missing
   df %>% filter(grepl('Staphylococcus', BUG), ABX == 'LEVOFLOXACIN') %>% count(LEVOFLOXACIN)   # 144 missing susceptibility status
   df %>% filter(grepl('Staphylococcus', BUG), ABX == 'CIPROFLOXACIN') %>% count(CIPROFLOXACIN) # 81 missing
   
   # IF susceptible to erythromycin and clindamycin,
   # THEN report as susceptible to all macrolides and lincosamides
   # affects (AZITH)
   w <- which(grepl('Staphylococcus', df$BUG) & df$CLINDAMYCIN == 0L & df$ERYTHROMYCIN == 0L) # 4,152
   df$AZITHROMYCIN[w] <- 0L # all were missing (except 4 were already S)
   df %>% filter(grepl('Staphylococcus', BUG), ABX == 'AZITHROMYCIN') %>% count(AZITHROMYCIN) # 512 missing
   df %>% filter(grepl('Staphylococcus', BUG), ABX == 'AZITHROMYCIN') %>% filter(is.na(AZITHROMYCIN)) %>% count(ERYTHROMYCIN, CLINDAMYCIN, sort=TRUE) # the rest that are missing are mostly R to either ERY or CLI
   
   # what about empiric clindamycin or erythromycin and missing AST?
   df %>% filter(grepl('Staphylococcus', BUG), ABX == 'CLINDAMYCIN') %>% count(CLINDAMYCIN, ERYTHROMYCIN)  # 7 missing (4 are ERY-R, 1 is ERY-S)
   df %>% filter(grepl('Staphylococcus', BUG), ABX == 'ERYTHROMYCIN') %>% count(ERYTHROMYCIN, CLINDAMYCIN) # 9 missing (3 are CLI-R, 5 are CLI-S)
   
   
   # DOXY 
   # if S-TET, then S-DOX, S-MINO, S-TIGE
   # if R-TET, then either R-DOX, R-MINO    OR    report DOX and MINO individually
   w <- which(grepl('Staphylococcus', df$BUG) & is.na(df$DOXYCYCLINE) & df$TETRACYCLINE == 0L) # 10,813 missing (53 need it (received DOXY empirically))
   df$DOXYCYCLINE[w] <- 0L
   df %>% filter(grepl('Staphylococcus', BUG)) %>% count(TETRACYCLINE, DOXYCYCLINE)
   w <- which(grepl('Staphylococcus', df$BUG) & df$TETRACYCLINE == 1L) # 873 (table(df$DOXYCYCLINE[w])) - some DOX-S
   df$DOXYCYCLINE[w] <- 1L
   
   # GENTAMICIN, TOBRAMYCIN
   # don't have any other aminoglycoside susceptibility data...
   df %>% filter(grepl('Staphylococcus', BUG), ABX == 'GENTAMICIN') %>% count(GENTAMICIN) # 9 missing
   df %>% filter(grepl('Staphylococcus', BUG), ABX == 'GENTAMICIN', is.na(GENTAMICIN)) %>% count(STREPTOMYCIN_SYNERGY, STREPTOMYCIN_HIGH_LEVEL, STREPTOMYCIN, NEOMYCIN, KANAMYCIN, AMIKACIN, TOBRAMYCIN, GENTAMICIN, GENTAMICIN_HIGH_LEVEL, GENTAMICIN_SYNERGY)
   df %>% filter(grepl('Staphylococcus', BUG), ABX == 'TOBRAMYCIN') %>% count(TOBRAMYCIN, GENTAMICIN, STREPTOMYCIN_SYNERGY, STREPTOMYCIN_HIGH_LEVEL, STREPTOMYCIN, NEOMYCIN, KANAMYCIN, AMIKACIN, TOBRAMYCIN, GENTAMICIN, GENTAMICIN_HIGH_LEVEL, GENTAMICIN_SYNERGY) # all 14 missing
   
   # NITROFURANTION
   df %>% filter(grepl('Staphylococcus', BUG), ABX == 'NITROFURANTOIN') %>% count(NITROFURANTOIN) # all 8 missing
   
   # BACITRACIN
   df %>% filter(grepl('Staphylococcus', BUG), ABX == 'BACITRACIN') %>% count(BACITRACIN) # all 47 missing
   
   
   # STAPHYLOCOCCUS LUGDUNENSIS
   # Beta-lactams + AZITHRO + CIPRO + DOXY + ERTAPENEM + ERYTHRO
   df[unlist(wmiss_list),] %>% filter(BUG == 'Staphylococcus lugdunensis') %>% count(ABX, sort=TRUE)
   df %>% filter(BUG == 'Staphylococcus lugdunensis') %>% count(OXACILLIN) # 35 (out of 250 have mecA?)
   
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
      w <- which(is.na(df[[abx]]) & df$BUG == 'Staphylococcus lugdunensis' & df$OXACILLIN == 0L & df$ABX == abx)
      df[[abx]][w] <- 0L
   }
   rm(abx, w, mrsl_R_abx)
}



## UPDATE missing ##
wmiss_list <- setNames(vector('list', length(wrows_list)),
                       names(wrows_list))
for (i in seq_along(wrows_list)) {
   abx <- names(wrows_list)[i]
   w <- wrows_list[[i]]
   df$FLAG[w] <- makeFlag(df[[abx]][w])
   wmiss_list[[i]] <- w[which(df$FLAG[w] == 'Unknown')]
}
rm(i, abx)
length(unlist(wmiss_list)) # 10,658
par(mar=c(5,4,1,1))
x <- sort(lengths(wmiss_list) / lengths(wrows_list))
b <- barplot(x, horiz=TRUE, names.arg = NA, xlim=c(0,1))
text(x=-0.01, y=b, adj=c(1, 0.5), labels=abbr[names(x)], xpd=NA)
rm(b, x)
# big problems seen below

miss_bugs <- matrix(0,
                    nrow=length(unique(df$BUG)), 
                    ncol=length(unique(df$ABX)),
                    dimnames = list(common_bugs, 
                                    names(wmiss_list)))
for (i in seq_along(wmiss_list)) {
   rows <- wmiss_list[[i]]
   abx <- names(wmiss_list)[i]
   bugs_table <- table(df$BUG[rows])
   miss_bugs[names(bugs_table), abx] <- bugs_table
}
rm(i, abx, rows, bugs_table)
miss_bugs[grep('Staphylococcus', rownames(miss_bugs)), ]
# still missing for SA: AZITHRO, CLINDA, ERYTHRO, GENTA, TOBRA, DOX, CIPRO, LEVO, BACITRACIN
# still missing for SL: AZITHRO, ERYTHRY, CIPRO, BACITRACIN





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
df[unlist(wmiss_list),] %>% filter(BUG == 'Streptococcus pneumoniae', ABX == 'CIPROFLOXACIN') %>% count(CIPROFLOXACIN, LEVOFLOXACIN, MOXIFLOXACIN) # 15 rows

# DOX-S if TET-S
w <- which(df$BUG == 'Streptococcus pneumoniae' & df$TETRACYCLINE == 0L) # 1093
df$DOXYCYCLINE[w] <- 0L
w <- which(df$BUG == 'Streptococcus pneumoniae' & df$TETRACYCLINE == 1L) # 217
df$DOXYCYCLINE[w] <- 1L

## UPDATE missing ##
wmiss_list <- setNames(vector('list', length(wrows_list)),
                       names(wrows_list))
for (i in seq_along(wrows_list)) {
   abx <- names(wrows_list)[i]
   w <- wrows_list[[i]]
   df$FLAG[w] <- makeFlag(df[[abx]][w])
   wmiss_list[[i]] <- w[which(df$FLAG[w] == 'Unknown')]
}
rm(i, abx)
length(unlist(wmiss_list)) # 9,304
df %>% summarise(COVERED = any(FLAG == 'Concordant'), .by=c(PERSON_ID, ORDER_DAY)) %>% count(COVERED)


# STILL MISSING: AZTREONAM, AMP, AMP/SUL, AMOX/CLAV, PIP/TAZ, CEFEPIME, CEFTRIAXONE, CEFAZOLIN, CEPHALEXIN, AZITH, CIPRO, DOXY, GENT
df[unlist(wmiss_list),] %>% 
   filter(BUG %in% c('Streptococcus pneumoniae', 'Streptococcus agalactiae', 'Streptococcus pyogenes'), 
          PENICILLIN == 1L) %>% # 192
   count(ABX, sort=TRUE)


par(mar=c(5,4,1,1))
x <- sort(lengths(wmiss_list) / lengths(wrows_list))
b <- barplot(x, horiz=TRUE, names.arg = NA, xlim=c(0,1))
text(x=-0.01, y=b, adj=c(1, 0.5), labels=abbr[names(x)], xpd=NA)
rm(b, x)
miss_bugs <- matrix(0,
                    nrow=length(unique(df$BUG)), 
                    ncol=length(unique(df$ABX)),
                    dimnames = list(common_bugs, 
                                    names(wmiss_list)))
for (i in seq_along(wmiss_list)) {
   rows <- wmiss_list[[i]]
   abx <- names(wmiss_list)[i]
   bugs_table <- table(df$BUG[rows])
   miss_bugs[names(bugs_table), abx] <- bugs_table
}
rm(i, abx, rows)
miss_bugs[grep('Strep', rownames(miss_bugs)), ]








##### ENTEROBACTERIACEAE #####
enterobacterales <- c('Escherichia coli', 'Klebsiella pneumoniae', 'Proteus mirabilis', 
                      'Enterobacter cloacae', 'Klebsiella oxytoca', 'Serratia marcescens', 
                      'Enterobacter aerogenes',
                      'Klebsiella aerogenes', 'Klebsiella variicola', 'Citrobacter freundii', 'Morganella morganii')
miss_bugs[enterobacterales, colSums(miss_bugs[enterobacterales,]) != 0]

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





## UPDATE missing ##
wmiss_list <- setNames(vector('list', length(wrows_list)),
                       names(wrows_list))
for (i in seq_along(wrows_list)) {
   abx <- names(wrows_list)[i]
   w <- wrows_list[[i]]
   df$FLAG[w] <- makeFlag(df[[abx]][w])
   wmiss_list[[i]] <- w[which(df$FLAG[w] == 'Unknown')]
}
rm(i, abx)
length(unlist(wmiss_list)) # 4,439
df %>% summarise(COVERED = any(FLAG == 'Concordant'), .by=c(PERSON_ID, ORDER_DAY)) %>% count(COVERED)

par(mar=c(5,4,1,1))
x <- sort(lengths(wmiss_list) / lengths(wrows_list))
b <- barplot(x, horiz=TRUE, names.arg = NA, xlim=c(0,1))
text(x=-0.01, y=b, adj=c(1, 0.5), labels=abbr[names(x)], xpd=NA)
rm(b, x)
miss_bugs <- matrix(0,
                    nrow=length(unique(df$BUG)), 
                    ncol=length(unique(df$ABX)),
                    dimnames = list(common_bugs,
                                    names(wmiss_list)))
for (i in seq_along(wmiss_list)) {
   rows <- wmiss_list[[i]]
   abx <- names(wmiss_list)[i]
   bugs_table <- table(df$BUG[rows])
   miss_bugs[names(bugs_table), abx] <- bugs_table
}
rm(i, abx, rows)
miss_bugs[enterobacterales,]
rm(esbl_R_set)
# missing: tobra, doxy, cipro, levo, dapto, a few polymixin B, bacitracin




##### ENTEROCOCCUS #####
enterococci <- c('Enterococcus faecium', 'Enterococcus faecalis')
miss_bugs[enterococci,] # missing amp/sul, pip/tazo, carbapenems, doxy, fluoros, bacitracin

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

df[unlist(wmiss_list),] %>%
   filter(BUG %in% enterococci) %>%
   count(`PIPERACILLIN/TAZOBACTAM`, BUG)

# doxycycline appears to be assumed resistant if not measured
w <- which(df$BUG %in% enterococci & is.na(df$DOXYCYCLINE))
df$DOXYCYCLINE[w] <- 1L

# still need cipro and levo



## UPDATE missing ##
wmiss_list <- setNames(vector('list', length(wrows_list)),
                       names(wrows_list))
for (i in seq_along(wrows_list)) {
   abx <- names(wrows_list)[i]
   w <- wrows_list[[i]]
   df$FLAG[w] <- makeFlag(df[[abx]][w])
   wmiss_list[[i]] <- w[which(df$FLAG[w] == 'Unknown')]
}
rm(i, abx)
length(unlist(wmiss_list)) # 2,936
df %>% summarise(COVERED = any(FLAG == 'Concordant'), .by=c(PERSON_ID, ORDER_DAY)) %>% count(COVERED)

par(mar=c(5,4,1,1))
x <- sort(lengths(wmiss_list) / lengths(wrows_list))
b <- barplot(x, horiz=TRUE, names.arg = NA, xlim=c(0,1))
text(x=-0.01, y=b, adj=c(1, 0.5), labels=abbr[names(x)], xpd=NA)
rm(b, x)
miss_bugs <- matrix(0,
                    nrow=length(unique(df$BUG)), 
                    ncol=length(unique(df$ABX)),
                    dimnames = list(common_bugs,
                                    names(wmiss_list)))
for (i in seq_along(wmiss_list)) {
   rows <- wmiss_list[[i]]
   abx <- names(wmiss_list)[i]
   bugs_table <- table(df$BUG[rows])
   miss_bugs[names(bugs_table), abx] <- bugs_table
}
rm(i, abx, rows)




# WHAT IS STILL MISSING ???
sort(rowSums(miss_bugs)) # S aures, S pneum, E coli, S agal, E faec, K pneum, S pyogen
full_join(
   df %>% filter(FLAG == 'Unknown') %>% count(BUG, sort=TRUE) %>% rename(missing = n),
   df %>% count(BUG, sort=TRUE),
   by = join_by(BUG)
) %>%
   mutate(fm = round(missing / n * 100, 2)) %>%
   arrange(desc(fm))
# S pneum, S. malto, A baumannii, S agalactiae,


full_join(
   df %>% filter(FLAG == 'Unknown') %>% count(ABX, sort=TRUE) %>% rename(missing = n),
   df %>% count(ABX, sort=TRUE),
   by = join_by(ABX)
) %>%
   mutate(fm = round(missing / n * 100, 2)) %>%
   arrange(desc(missing))
# doxy, levo, tobra, azithro, cipro, genta

df %>%
   filter(FLAG == 'Unknown') %>%
   count(BUG, ABX, sort=TRUE)


full_join(
   df %>% filter(FLAG == 'Discordant') %>% count(ABX, sort=TRUE) %>% rename(missing = n),
   df %>% count(ABX, sort=TRUE),
   by = join_by(ABX)
) %>%
   mutate(fm = round(missing / n * 100, 2)) %>%
   arrange(desc(fm))
full_join(
   df %>% filter(FLAG == 'Discordant') %>% count(BUG, sort=TRUE) %>% rename(discordant = n),
   df %>% count(BUG, sort=TRUE),
   by = join_by(BUG)
) %>%
   mutate(fd = round(discordant / n * 100, 2)) %>%
   arrange(desc(fd))
# doxy, levo, tobra, azithro, cipro, genta









# missing by bug and class:
m <- df %>% 
   count(BUGC, ABXC) %>% 
   tidyr::pivot_wider(
      names_from = BUGC, 
      values_from = n, 
      values_fill = 0
   )
n <- m$ABXC
m <- as.matrix(m %>% select(-ABXC))
rownames(m) <- n
m <- m[c('glycopeptide', 'lipopeptide', 'oxazolidinone',
         'penicillin', 'penicillin_blinhibitor', 'monobactam',
         'cephalosporin_1g', 'cephalosporin_2g', 'cephalosporin_3g', 'cephalosporin_4g', 'cephalosporin_5g',
         'carbapenem', 
         'sulfonamide_dihydrofolate', 'aminoglycoside', 'fluoroquinolone', 'macrolide', 'lincosamide', 'tetracycline'),
       c('NonFermGN', 'Enterobacterales', 'Enterococci', 'Staphylococci', 'Streptococci')]
rm(n)


# draw heatmap
drawHeatmap <- function(m) {
   col_vec <- colorRampPalette(c('white', 'black'))(max(m))
   xpos <- seq_len(ncol(m))
   ypos <- seq_len(nrow(m))
   width <- 0.5
   plot(NA, xaxs='i', yaxs='i',
        ylim=c(max(ypos)+width, min(ypos)-width),
        xlim=c(min(xpos)-width, max(xpos)+width))
   text(x = rep(xpos, each=length(ypos)),
        y = rep(ypos, times=length(xpos)),
        labels = m)
   rect(xleft = rep(xpos, each=length(ypos)) - width,
        xright = rep(xpos, each=length(ypos)) + width,
        ybottom = rep(ypos, times=length(xpos)) - width,
        ytop = rep(ypos, times=length(xpos)) + width,
        border = NA,
        col = col_vec[m])
}

par(mfrow=c(1,2))
drawHeatmap(round(apply(m, 2, function(x) x / sum(x)) * 100))
drawHeatmap(t(round(apply(m, 1, function(x) x / sum(x)) * 100)))
rm(m)








#### WHAT ARE THE TREATMENTS WE CARE ABOUT PREDICTING??? ####
# 1) whether or not a carbapenem will be required
#        if they prescribed it, they predicted it was necessary, can we do better?
#        ground truth: bug requires carbapenem (ESBL enterbacterales)
t <- df %>% 
   filter(BUG %in% enterobacterales) %>% 
   select(PERSON_ID, ORDER_DAY, BUG, ESBL) %>%
   distinct() %>% 
   mutate(year = lubridate::year(ORDER_DAY)) %>%
   summarise(
      ESBL = as.integer(any(ESBL == 1L)),
      .by=c(PERSON_ID, ORDER_DAY, year)
   ) %>% # count(ESBL)                                # 4,187 positive examples
   select(ESBL, year) %>% table()  
t['1',] / colSums(t) # not really increasing in prevalence over time
rm(t)
df %>% select(PERSON_ID, ORDER_DAY) %>% distinct()    # 37,705 negative examples


calcMetrics <- function(t) {
   sens <- t[1,1] / sum(t[1,])
   spec <- t[2,2] / sum(t[2,])
   prec <- t[1,1] / sum(t[,1])
   ngpv <- t[2,2] / sum(t[,2])
   r <- function(x) round(x * 100, 1)
   cat('Sens:', r(sens), '\n')
   cat('Prec:', r(prec), '\n')
   cat('Spec:', r(spec), '\n')
   cat('Ngpv:', r(ngpv), '\n')
}

flags <- c('ESBL', 'VANCOMYCIN')
t <- df %>%
   select(PERSON_ID, ORDER_DAY, !!flag, ABX) %>%
   rename(X = !!flag) %>%
   summarise(
      X = as.integer(any(X == 1L)),
      PEN = as.integer(any(grepl('PENEM$', ABX))),
      .by = c(PERSON_ID, ORDER_DAY)
   ) %>%
   mutate(X = ifelse(is.na(X), 0, X)) %>% # count(ESBL)
   select(X, PEN) %>%
   table()
names(dimnames(t))[1] <- flag
t <- t[2:1, 2:1]
calcMetrics(t)


save(df, file = '~/Desktop/EHR/EHR work/RdataFiles/R01/bact_empiric_flag_ast_abx.Rdata')



















