# How do we decide what to predict? Whatever was prescribed empirically!
# What is the baseline performance? Whatever was prescribed empirically!
# Heuristic: if a physician prescribed it empirically, they predict Susceptible.
# Can we get the physician score?
library(dplyr)
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_clean_ASTs.Rdata')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_CLEANED_2015_2024_AbxAdmin.Rdata') # 12,373,076 rows
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
ast <- ast %>%
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
# first, for staph aureus missing oxacillin AST, remove
w <- which(ast$BUG == 'Staphylococcus aureus' & is.na(ast$OXACILLIN)) # 20
ast <- ast[-w,]
rm(w)

ast$CEPHALEXIN <- NA
sum(sapply(ast %>% select(CEFEPIME:CEPHALEXIN), function(x) sum(!is.na(x)))) # 646,078
# ast_og <- ast

# IMPUTE MISSING ASTs
source('~/Desktop/EHR-mining/UsefulDataForCleaning/ASTimputation/ASTimputation.R')
ast <- imputeASTs(ast)
sum(sapply(ast %>% select(CEFEPIME:CEPHALEXIN), function(x) sum(!is.na(x)))) # 646,105 --> 1,982,857



abxDF <- abxDF %>% filter(PERSON_ID %in% unique(ast$PERSON_ID))

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
   filter(!is.na(ABX)) %>% 
   relocate(ABX, .after=BUG)

source('~/Desktop/EHR-mining/UsefulDataForCleaning/getAbxBugClassFxn.R')
df <- getBugClass(df)
df <- getAbxClass(df)

# df_og <- df
df <- df_og




# narrow down to cultures that were treated empirically with relatively common abx
# - at least 100 infections were treated empirically with this
common_abx <- df %>% select(PERSON_ID, ORDER_DAY, ABX) %>% distinct() %>% count(ABX, sort=TRUE) %>% filter(n >= 100L) %>% pull(ABX)
common_abx <- common_abx[!common_abx %in% c('BACITRACIN', 'POLYMYXIN B')]
df <- df %>% filter(ABX %in% common_abx)# %>% # 39,979 --> 39,630 events


# get discordance / concordance flag:
df$FLAG <- NA
df <- df %>% relocate(FLAG, .after=ABXC)


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


## Assess missingness and discordance flags ##
wmiss_list <- setNames(vector('list', length(wrows_list)),
                       names(wrows_list))
for (i in seq_along(wrows_list)) {
   abx <- names(wrows_list)[i]
   w <- wrows_list[[i]]
   df$FLAG[w] <- makeFlag(df[[abx]][w])
   wmiss_list[[i]] <- w[which(df$FLAG[w] == 'Unknown')]
}
rm(i, abx)
length(unlist(wmiss_list)) # 3,055 (out of 76K)
par(mar=c(5,4,1,1))
x <- sort(lengths(wmiss_list) / lengths(wrows_list))
b <- barplot(x, horiz=TRUE, names.arg = NA, xlim=c(0,1))
text(x=-0.01, y=b, adj=c(1, 0.5), labels=abbr[names(x)], xpd=NA)
rm(b, x)

# which bugs/drugs are missing?
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


df %>% count(FLAG)
df %>%
   summarise(COVERED = any(FLAG == 'Concordant'), 
             .by=c(PERSON_ID, ORDER_DAY)) %>% 
   count(COVERED)



save(df, file = '~/Desktop/EHR/EHR work/RdataFiles/R01/bact_empiric_flag_ast_abx.Rdata')







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




















