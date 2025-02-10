library(dplyr)
load(file = '~/Desktop/EHR/EHR work/RdataFiles/R01/bact_empiric_flag_ast_abx.Rdata')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_clean_ASTs.Rdata')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_CLEANED_2015_2024_AbxAdmin.Rdata')
source('~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/CreateNamedAbxAbbreviations.R')
source('~/Desktop/EHR-mining/UsefulDataForCleaning/getAbxBugClassFxn.R')
source('~/Desktop/EHR-mining/UsefulDataForCleaning/ASTimputation/ASTimputation.R')


# We want to find signals of association between prior patient/facility/time-specific
# antibiotic use + resistance + pathogens --> modern empiric Rx + resistance + pathogens

# To do so, we need: 
#     historical data on prescriptions + some way to encode
#     historical data on resistance + some way to encode


##### START data preprocessing #####
# empiric flag by visit
# add visit flags for MRSA, VRE, ESBL, MDR-Pseudomonas
df <- df %>%
   mutate(APP_R = case_when(`PIPERACILLIN/TAZOBACTAM` == 1L ~ 1L, .default = 0L),
          CEF_R = case_when(CEFTAZIDIME == 1L | CEFEPIME == 1L ~ 1L, .default = 0L),
          FLQ_R = case_when(CIPROFLOXACIN == 1L | LEVOFLOXACIN == 1L ~ 1L, .default = 0L),
          AMI_R = case_when(TOBRAMYCIN == 1L | GENTAMICIN == 1L | AMIKACIN == 1L ~ 1L, .default = 0L),
          CAR_R = case_when(MEROPENEM == 1L | IMIPENEM == 1L ~ 1L, .default = 0L)) %>% 
   mutate(MDR_SCORE = APP_R + CEF_R + FLQ_R + AMI_R + CAR_R) %>%
   group_by(PERSON_ID, ORDER_DAY) %>%
   mutate(
      DISCORDANT = as.integer(!any(FLAG == 'Concordant')),
      MRSA = as.integer(any(BUG == 'Staphylococcus aureus' & OXACILLIN == 1L)),
      VRE = as.integer(any(grepl('Enterococcus', BUG) & VANCOMYCIN == 1L)),
      ESBL = as.integer(any(!is.na(ESBL) & ESBL == 1L)),
      MDR_Pseud = as.integer(any(BUG == 'Pseudomonas aeruginosa' & MDR_SCORE >= 3L)),
   ) %>%
   select(!c(APP_R, CEF_R, FLQ_R, AMI_R, CAR_R, MDR_SCORE, FLAG)) %>%
   ungroup()

# get resistance flags if R to any antibiotic in each class
df <- getAbxClassResistanceFlags(df)

# check missingness of resistance phenotypes
x <- df %>%
   select(!c(ABX, ABXC, BUG, BUGC)) %>%
   distinct() %>%
   select(!!unique(df$ABX))
x <- sapply(x, FUN = function(x) sum(!is.na(x))) / nrow(x)
data.frame(completeness = round(sort(x), 3))
keep_phenotypes_abx <- names(x[x > 0.9])
rm(x)
df <- df %>% 
   select(PERSON_ID, ORDER_DAY, ABX, ABXC, BUG, BUGC, DISCORDANT, MRSA, VRE, ESBL, MDR_Pseud, !!unique(df$ABX), contains('ResC_')) %>% 
   distinct()
df$ABX <- unname(abbr[df$ABX])
all_bugs <- unique(df$BUG)
df$BUG <- gsub('^([A-Z])[a-z]+ ([a-z]+)$', '\\1_\\2', df$BUG)

# resistance phenotypes present in current infection
res1DF <- df %>% 
   select(PERSON_ID, ORDER_DAY, !!keep_phenotypes_abx) %>% # by abx
   distinct() %>%
   mutate(
      across(.cols = !c(PERSON_ID, ORDER_DAY),
             .fns = ~ as.integer(ifelse(is.na(.), 0, .)))
   ) %>%
   summarise(
      across(.cols = everything(), 
             .fns = max),
      .by = c(PERSON_ID, ORDER_DAY)
   )
names(res1DF)[!names(res1DF) %in% c('PERSON_ID', 'ORDER_DAY')] <- paste0('Res1_', abbr[names(res1DF)[!names(res1DF) %in% c('PERSON_ID', 'ORDER_DAY')]])

resCDF <- df %>% 
   select(PERSON_ID, ORDER_DAY, contains('ResC_')) %>%
   summarise(
      across(.cols = contains('ResC_'),
             .fns = max),
      .by = c(PERSON_ID, ORDER_DAY)
   )
##### END data preprocessing #####


getDF <- function(df, var, prefix='') {
   df <- df %>%
      select(PERSON_ID, ORDER_DAY, !!var) %>%
      distinct() %>%
      mutate(X = 1L) %>%
      tidyr::pivot_wider(
         id_cols = c(PERSON_ID, ORDER_DAY),
         values_from = X,
         names_from = !!var,
         values_fill = 0,
         names_prefix = prefix
      )
   return(df)
}


##### START outcomes #####
outcomesDF <- inner_join(
   x = inner_join(
      x = inner_join(
         # prescribed empiric abx for current infection
         x = getDF(df, 'ABX', 'Rx1_'), # by abx
         y = getDF(df, 'ABXC', 'RxC_'), # by abx class
         by = join_by(PERSON_ID, ORDER_DAY)
      ),
      y = inner_join(
         # pathogens isolated in current infection
         x = getDF(df, 'BUG', 'Bug1_'), # by bug
         y = getDF(df, 'BUGC', 'BugC_'), # by bug group
         by = join_by(PERSON_ID, ORDER_DAY)
      ),
      by = join_by(PERSON_ID, ORDER_DAY)
   ),
   y = inner_join(
      # resistance phenotypes in current infection
      x = res1DF, # by abx_class
      y = resCDF,
      by = join_by(PERSON_ID, ORDER_DAY)
   ),
   by = join_by(PERSON_ID, ORDER_DAY)
)
outcomesDF <- full_join(
   x = df %>% 
      select(PERSON_ID, ORDER_DAY, MRSA, VRE, ESBL, MDR_Pseud, DISCORDANT) %>% 
      distinct(),
   y = outcomesDF,
   by = join_by(PERSON_ID, ORDER_DAY)
)
outcomesDF %>% names
rm(res1DF, resCDF, keep_phenotypes_abx)

save(outcomesDF, file = '~/Desktop/EHR/EHR work/RdataFiles/R01/bact_outcomes_variables.Rdata')
##### END outcomes #####


##### START historical data featurization #####
## ASTs ##
# take all astDF: 
#     -- keep only these bugs
#     -- add bug class
#     -- keep only relevant abx
#     -- impute all missing
#     -- abbreviate abx names
#     -- add abx class - check if any in class are R
#     -- only some patients
#     -- remove outcomesDF ASTs
#     abbreviate bug names
# astDF_og <- astDF
astDF <- astDF %>%
   filter(
      PERSON_ID %in% unique(outcomesDF$PERSON_ID),
      BUG %in% all_bugs
   )
w <- which(astDF$BUG == 'Staphylococcus aureus' & is.na(astDF$OXACILLIN)) # 688
astDF <- astDF[-w,]
rm(w)
# remove the ASTs we are looking at??
astDF <- astDF %>%
   left_join(
      x = .,
      y = outcomesDF %>% select(PERSON_ID, ORDER_DAY) %>% mutate(X = 1L),
      by = join_by(PERSON_ID, ORDER_DAY)
   ) %>%
   filter(is.na(X)) %>% # this step removes any AST data that is in our outcomes dataframe
   select(-X) %>%
   distinct()

# impute missing ASTs
sum(sapply(astDF %>% select(CEFEPIME:CEPHALEXIN), function(x) sum(!is.na(x))))
astDF <- imputeASTs(astDF)

# get visit flags for MRSA, VRE, ESBL, MDR-Pseudomonas
astDF <- astDF %>%
   mutate(APP_R = case_when(`PIPERACILLIN/TAZOBACTAM` == 1L ~ 1L, .default = 0L),
          CEF_R = case_when(CEFTAZIDIME == 1L | CEFEPIME == 1L ~ 1L, .default = 0L),
          FLQ_R = case_when(CIPROFLOXACIN == 1L | LEVOFLOXACIN == 1L ~ 1L, .default = 0L),
          AMI_R = case_when(TOBRAMYCIN == 1L | GENTAMICIN == 1L | AMIKACIN == 1L ~ 1L, .default = 0L),
          CAR_R = case_when(MEROPENEM == 1L | IMIPENEM == 1L ~ 1L, .default = 0L)) %>% 
   mutate(MDR_SCORE = APP_R + CEF_R + FLQ_R + AMI_R + CAR_R) %>%
   group_by(PERSON_ID, ORDER_DAY) %>%
   mutate(
      MRSA = as.integer(any(BUG == 'Staphylococcus aureus' & OXACILLIN == 1L)),
      VRE = as.integer(any(grepl('Enterococcus', BUG) & VANCOMYCIN == 1L)),
      ESBL = as.integer(any(!is.na(ESBL) & ESBL == 1L)),
      MDR_Pseud = as.integer(any(BUG == 'Pseudomonas aeruginosa' & MDR_SCORE >= 3L)),
   ) %>%
   select(!c(APP_R, CEF_R, FLQ_R, AMI_R, CAR_R, MDR_SCORE)) %>%
   ungroup()

astDF <- getBugClass(astDF)
astDF$BUG <- gsub('^([A-Z])[a-z]+ ([a-z]+)$', '\\1_\\2', astDF$BUG)
astDF <- getAbxClassResistanceFlags(astDF)

# which antibiotics to keep individual and class resistance flags for?
#     -- all classes
#     -- individuals
keep_abx <- gsub('Res1_', '', grep('Res1', names(outcomesDF), value=TRUE))
keep_abx <- unname(setNames(names(abbr), abbr)[keep_abx])
astDF <- astDF %>%
   select(PERSON_ID, ORDER_DAY, BUG, BUGC, MRSA, VRE, ESBL, MDR_Pseud, !!keep_abx, contains('ResC')) %>%
   distinct()
names(astDF)[names(astDF) %in% keep_abx] <- paste0('Res1_', 
                                                   abbr[names(astDF)[names(astDF) %in% keep_abx]])
# missing ASTs will be imputed as 0 for now
astDF <- astDF %>%
   mutate(across(.cols = contains('Res1_'),
                 .fns = ~ as.integer(ifelse(is.na(.), 0L, .))))


astDF <- astDF %>%
   rename(PREV_ORDER_DAY = ORDER_DAY) %>%
   inner_join(
      x = .,
      y = outcomesDF %>%
         select(PERSON_ID, ORDER_DAY) %>%
         mutate(
            X = 1L,
            JOIN_END = ORDER_DAY - 7L
         ),
      by = join_by(
         PERSON_ID,
         PREV_ORDER_DAY <= JOIN_END
      )
   ) %>%
   distinct() %>%
   mutate(DAYS_BEFORE_ORDER = as.integer(ORDER_DAY - PREV_ORDER_DAY)) %>% 
   select(!c(PREV_ORDER_DAY, X, JOIN_END)) %>%
   relocate(ORDER_DAY, DAYS_BEFORE_ORDER, .after=PERSON_ID)

names(astDF)[names(astDF) %in% c('MRSA', 'VRE', 'ESBL', 'MDR_Pseud')] <- paste0('ResX_',
                                                                                names(astDF)[names(astDF) %in% c('MRSA', 'VRE', 'ESBL', 'MDR_Pseud')])


 
### ABXs ###
# tasks:
#     -- only patients we care about
#     -- remove admins from infection
#     -- class-wise antibiotic Rx
#     only specific antibiotics

abxDF <- abxDF %>% filter(PERSON_ID %in% unique(outcomesDF$PERSON_ID))
w <- which(abxDF$ABX == 'PENICILLIN G')
abxDF$ABX[w] <- 'BENZYLPENICILLIN'
rm(w)

abxDF <- abxDF %>%
   inner_join(
      x = .,
      y = outcomesDF %>% 
         select(PERSON_ID, ORDER_DAY) %>%
         mutate(
            X = 1L,
            JOIN_END = ORDER_DAY - 7L
         ),
      by = join_by(
         PERSON_ID,
         START_DAY <= JOIN_END
      )
   ) %>%
   select(PERSON_ID, ORDER_DAY, START_DAY, ABX) %>%
   distinct() %>%
   mutate(DAYS_BEFORE_ORDER = as.integer(ORDER_DAY - START_DAY)) %>% 
   select(-START_DAY)

abx_class_named <- setNames(gsub(pattern = '[0-9]+$', 
                                 replacement = '', 
                                 x = names(unlist(abx_class))),
                            unname(unlist(abx_class)))
abxDF <- abxDF %>% 
   mutate(
      ABXC = unname(abx_class_named[ABX]),
      ABX = unname(abbr[ABX])
   ) %>%
   relocate(ABXC, .after=ABX) %>%
   filter(!is.na(ABXC))




# FEATURIZE ABXs AND ASTs by time since index culture
# ever, past 2 years, past 180 days, past 60 days
getTimeFrameFeatures <- function(df, var, days) {
   flag2 <- ifelse(grepl('(BUG|ABX)C$', var), 'C', '1')
   flag1 <- ifelse(grepl('BUGC?', var), '_Bug', '_Rx')
   prefix <- paste0('D_', days, flag1, flag2, '_')
   if (days == 'ever') days <- 50000
   if (var == 'ABX') df <- df %>% filter(ABX %in% unname(abbr[keep_abx]))
   df %>%
      filter(DAYS_BEFORE_ORDER <= days) %>%
      select(PERSON_ID, ORDER_DAY, !!var) %>%
      distinct() %>%
      mutate(X = 1L) %>%
      tidyr::pivot_wider(
         id_cols = c(PERSON_ID, ORDER_DAY),
         names_from = !!var,
         values_from = X,
         values_fill = 0L,
         names_prefix = prefix
      )
}

abxFeaturesDF <- full_join(
   x = full_join(
      x = full_join(
         x = getTimeFrameFeatures(abxDF, 'ABXC', 'ever'),
         y = getTimeFrameFeatures(abxDF, 'ABX', 'ever'),
         by = join_by(PERSON_ID, ORDER_DAY)
      ),
      y = full_join(
         x = getTimeFrameFeatures(abxDF, 'ABXC', 730),
         y = getTimeFrameFeatures(abxDF, 'ABX', 730),
         by = join_by(PERSON_ID, ORDER_DAY)
      ),
      by = join_by(PERSON_ID, ORDER_DAY)
   ),
   y = full_join(
      x = full_join(
         x = getTimeFrameFeatures(abxDF, 'ABXC', 365),
         y = getTimeFrameFeatures(abxDF, 'ABX', 365),
         by = join_by(PERSON_ID, ORDER_DAY)
      ),
      y = full_join(
         x = getTimeFrameFeatures(abxDF, 'ABXC', 90),
         y = getTimeFrameFeatures(abxDF, 'ABX', 90),
         by = join_by(PERSON_ID, ORDER_DAY)
      ),
      by = join_by(PERSON_ID, ORDER_DAY)
   ),
   by = join_by(PERSON_ID, ORDER_DAY)
)

astFeaturesDF <- full_join(
   x = full_join(
      x = full_join(
         x = getTimeFrameFeatures(astDF, 'BUGC', 'ever'),
         y = getTimeFrameFeatures(astDF, 'BUG', 'ever'),
         by = join_by(PERSON_ID, ORDER_DAY)
      ),
      y = full_join(
         x = getTimeFrameFeatures(astDF, 'BUGC', 730),
         y = getTimeFrameFeatures(astDF, 'BUG', 730),
         by = join_by(PERSON_ID, ORDER_DAY)
      ),
      by = join_by(PERSON_ID, ORDER_DAY)
   ),
   y = full_join(
      x = full_join(
         x = getTimeFrameFeatures(astDF, 'BUGC', 365),
         y = getTimeFrameFeatures(astDF, 'BUG', 365),
         by = join_by(PERSON_ID, ORDER_DAY)
      ),
      y = full_join(
         x = getTimeFrameFeatures(astDF, 'BUGC', 90),
         y = getTimeFrameFeatures(astDF, 'BUG', 90),
         by = join_by(PERSON_ID, ORDER_DAY)
      ),
      by = join_by(PERSON_ID, ORDER_DAY)
   ),
   by = join_by(PERSON_ID, ORDER_DAY)
)

# add resistance flags to astFeaturesDF
# ever, past 2 years, past 180 days, past 60 days
getResTimeFrameFeatures <- function(df, days) {
   prefix <- paste0('D_', days, '_')
   if (days == 'ever') days <- 50000
   df <- df %>%
      filter(DAYS_BEFORE_ORDER <= days) %>%
      select(-BUG, -BUGC, -DAYS_BEFORE_ORDER) %>%
      distinct() %>%
      summarise(
         across(.cols = matches('^Res(1|C|X)_'),
                .fns = max),
         .by = c(PERSON_ID, ORDER_DAY)
      )
   names(df)[!names(df) %in% c('PERSON_ID', 'ORDER_DAY')] <- paste0(prefix,
                                                                    names(df)[!names(df) %in% c('PERSON_ID', 'ORDER_DAY')])
   return(df)
}

resFeaturesDF <- full_join(
   x = full_join(
      x = getResTimeFrameFeatures(astDF, 'ever'),
      y = getResTimeFrameFeatures(astDF, 730),
      by = join_by(PERSON_ID, ORDER_DAY)
   ),
   y = full_join(
      x = getResTimeFrameFeatures(astDF, 365),
      y = getResTimeFrameFeatures(astDF, 90),
      by = join_by(PERSON_ID, ORDER_DAY)
   ),
   by = join_by(PERSON_ID, ORDER_DAY)
)

astFeaturesDF <- full_join(
   x = astFeaturesDF,
   y = resFeaturesDF,
   by = join_by(PERSON_ID, ORDER_DAY)
)

featuresDF <- full_join(
   x = abxFeaturesDF,
   y = astFeaturesDF,
   by = join_by(PERSON_ID, ORDER_DAY)
)

# %>%
#    mutate(
#       across(.cols = !c(PERSON_ID, ORDER_DAY),
#              .fns = ~ ifelse(is.na(.), 0, .))
#    )


save(featuresDF, file = '~/Desktop/EHR/EHR work/RdataFiles/R01/bact_features_variables.Rdata')













