library(dplyr)
load(file = '~/Desktop/EHR-mining/UsefulDataForCleaning/bug_groups.Rdata')
source('~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/CreateNamedAbxAbbreviations.R')
abbr <- c(abbr, 'PENICILLIN G' = 'PEN')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_clean_ASTs.Rdata')
years <- 2016:2024

Rdata_file_path <- '~/Desktop/EHR/EHR work/RdataFiles/BuildCohortDatasets/RawData/'
load(file = paste0(Rdata_file_path, 'abxDF_clean_2015_2024_bacteremia.Rdata'))
load(file = paste0(Rdata_file_path, 'encsDF_raw_2015_2024_bacteremia.Rdata'))

fxn_file_path <- '~/Desktop/EHR-mining/BuildCohortDatasets/FeaturizationScripts/'
source(paste0(fxn_file_path, 'featurizeEncounters.R'))
source(paste0(fxn_file_path, 'featurizeDemographics.R'))


astDF_og <- astDF
# Keep only 'index' blood cultures - first blood culture in record or first in 60 days
astDF <- astDF_og %>%
   filter(
      SITE == 'blood',
      lubridate::year(ORDER_DAY) %in% years,
      BUG %in% unname(unlist(bug_groups))
   )
astDF$DAYS_SINCE_PRV <- astDF$ORDER_DAY - lag(astDF$ORDER_DAY)
astDF$DAYS_SINCE_PRV[c(1, which(astDF$PERSON_ID != lag(astDF$PERSON_ID)))] <- NA_integer_
astDF <- astDF %>%
   filter(is.na(DAYS_SINCE_PRV) | DAYS_SINCE_PRV > 60L) %>%
   select(-DAYS_SINCE_PRV)


# abbreviate bug names
# astDF <- astDF %>%
#    mutate(BUG = gsub('^([A-Z])[a-z]+ ([a-z]+)$', '\\1_\\2', BUG))
# 
# entero_abbr <- gsub('^([A-Z])[a-z]+ ([a-z]+)$', '\\1_\\2', bug_groups$Enterobacterales)
astDF <- astDF %>%
   mutate(BUG = case_when(
      ESBL == 1L & BUG %in% unlist(bug_groups$Enterobacterales) ~ paste0(BUG, '_ESBL'),
      OXACILLIN == 1L & BUG == 'Staphylococcus aureus' ~ paste0(BUG, '_MRSA'),
      VANCOMYCIN == 1L & BUG %in% c('Enterococcus faecalis', 'Enterococcus faecium') ~ paste0(BUG, '_VRE'),
      .default = BUG
   )) %>%
   select(PERSON_ID, ORDER_DAY, RESULT_DAY, ORDER_DATE, RESULT_DATE, BUG)
df <- astDF


# Join Demographics.
df <- featurizeDemographics(df=df, var_location='BUG', index_var='ORDER_DAY')

# Join Encounters.
df <- featurizeEncounters(df=df, encsDF=encsDF %>% filter(PERSON_ID %in% unique(df$PERSON_ID)))

# Join Antibiotics Administration.
# Goal is to assess antibiotic treatment heterogeneity on a bug-resistance basis
# So let's only look at AbxAdmin after AST results are returned
df_og <- df
df <- left_join(
   x = df_og %>% 
      mutate(
         JOIN_START = ORDER_DATE,
         JOIN_END = ORDER_DATE + 4 * 86400
      ),
   y = abxDF,
   by = join_by(
      PERSON_ID,
      between(y$START_DATE, x$JOIN_START, x$JOIN_END)
   )
) %>%
   mutate(
      ABX_DAY = as.integer(START_DAY - ORDER_DAY),
      ABX_TIME = as.numeric(lubridate::as.duration(START_DATE - ORDER_DATE)) / 86400
   )

df <- df %>%
   group_by(PERSON_ID, ORDER_DAY) %>%
   filter(
      any(ABX_TIME < 2)
   ) %>%
   ungroup()
df_og2 <- df




# Which time frame to keep?
# How sensitive are the cohort/treatment group sample sizes to this choice?
getDisJointness <- function(bug, abx1, abx2, index_var='ORDER', window_start=1, window_end=4) {
   # stopifnot('Need >1 abx' = length(abx) > 1L)
   
   cat('trt1:', abx1, '\n')
   cat('trt2:', abx2, '\n')
   
   x <- df_og2 %>%
      filter(
         BUG %in% bug, 
         ABX %in% c(abx1, abx2),
         ABX_TIME > window_start & ABX_TIME < window_end
      ) %>%
      mutate(ABX_DAY = floor(ABX_TIME)) %>%
      select(PERSON_ID, ORDER_DAY, ABX, ABX_DAY) %>%
      distinct() %>%
      # group_by(PERSON_ID, ORDER_DAY, ABX_DAY) %>%
      mutate(
         ABX = case_when(
            ABX %in% abx1 ~ paste(unname(abbr[abx1]), collapse='/'),
            ABX %in% abx2 ~ paste(unname(abbr[abx2]), collapse='/'),
            .default = ABX
         )
         #ABX = unname(abbr[ABX]),
      ) %>%
      distinct()
   abx <- unique(x$ABX)
   x <- x %>%
      mutate(ABX_DAY = 1L) %>%
      tidyr::pivot_wider(
         values_from = ABX_DAY,
         names_from = ABX,
         values_fn = sum,
         values_fill = 0L
      ) %>%
      select(-PERSON_ID, -ORDER_DAY)
   t <- sort(table(apply(x, 1, function(r) paste(abx[r > 0L], collapse='+'))))
   t <- setNames(data.frame(t), c('ABX', 'Count'))
   return(t)
   
   # t <- t[!grepl('+', t$ABX, fixed=TRUE), ]
   # t <- setNames(t$Count, t$ABX)
   # barplot(t)
}




top_bugs <- df %>% select(PERSON_ID, ORDER_DAY, BUG) %>% distinct() %>% count(BUG, sort=TRUE)
window_start <- 1
window_end <- 4
index_var <- 'ORDER'
bug_list <- list()
idx <- 1L

while (TRUE) {
   print(top_bugs, n=25L)
   b <- as.integer(strsplit(readline('Choose bug row: '), split=',')[[1]])
   
   bug <- top_bugs$BUG[b]
   print(bug)
   
   top_abxs <- df_og2 %>% 
      filter(
         BUG %in% bug,
         ABX_TIME > window_start & ABX_TIME < window_end
      ) %>%
      mutate(ABX_DAY = floor(ABX_TIME)) %>%
      select(PERSON_ID, ORDER_DAY, ABX, ABX_DAY) %>%
      distinct() %>%
      count(ABX, sort=TRUE)
   print(top_abxs, n=20)
   trt1 <- as.integer(strsplit(readline(prompt='Choose group 1: '), split=',')[[1]])
   trt2 <- as.integer(strsplit(readline(prompt='Choose group 2: '), split=',')[[1]])
   
   dj <- getDisJointness(bug=bug, abx1=top_abxs$ABX[trt1], abx2=top_abxs$ABX[trt2], index_var=index_var, window_start=window_start, window_end=window_end)
   print(dj)
   x <- readline('Good? [y, n]: ')
   if (x == 'y') {
      bug_list[[idx]] <- dj
      names(bug_list)[idx] <- bug
      idx <- idx + 1L
   }
}



# The 11 bacteremia pathogen-abx cohorts that the AMP team has roughly agreed are worth investigating.
cohort_info <- list(
   list('bug' = 'E_coli',        'abx' = c('CEFTRIAXONE', 'PIPERACILLIN/TAZOBACTAM')),
   list('bug' = 'K_pneumoniae',  'abx' = c('CEFTRIAXONE', 'PIPERACILLIN/TAZOBACTAM')),
   list('bug' = 'P_mirabilis',   'abx' = c('CEFTRIAXONE', 'PIPERACILLIN/TAZOBACTAM')),
   list('bug' = 'P_aeruginosa',  'abx' = c('CEFEPIME',    'PIPERACILLIN/TAZOBACTAM')),
   list('bug' = c('E_cloacae', 'C_freundii', 'K_aerogenes', 'S_marcescens', 'P_rettgeri',
                  'P_stuartii', 'P_alcalifaciens', 'M_morganii'), 
        'abx' = c('CEFEPIME', 'PIPERACILLIN/TAZOBACTAM')),
   list('bug' = 'S_aureus',      'abx' = c('OXACILLIN+NAFCILLIN',   'CEFAZOLIN')),
   list('bug' = 'S_aureus_MRSA', 'abx' = c('VANCOMYCIN',  'DAPTOMYCIN'), switch=TRUE),
   list('bug' = 'S_aureus_MRSA', 'abx' = c('VANCOMYCIN', 'DAPTOMYCIN')),
   list('bug' = 'E_faecalis',    'abx' = c('VANCOMYCIN',  'AMPICILLIN'), switch=TRUE),
   list('bug' = 'E_faecalis',    'abx' = c('VANCOMYCIN',  'AMPICILLIN')),
   list('bug' = 'E_faecium_VRE', 'abx' = c('DAPTOMYCIN',  'LINEZOLID'))
)

# Add some extra info about each.
for (cohort in seq_along(cohort_info)) {
   cohort_bug <- cohort_info[[cohort]]$bug
   cohort_abx <- cohort_info[[cohort]]$abx
   cohort_name <- paste0(ifelse(length(cohort_bug) > 1L, 'ampCprod', cohort_bug),
                         ifelse('switch' %in% names(cohort_info[[cohort]]), '_switch', ''), '_',
                         paste0(sapply(strsplit(cohort_abx, split='+', fixed=TRUE), function(abx) paste(abbr[abx], collapse='or')), collapse='_'), 
                         '_treatment')
   cohort_info[[cohort]]$cohort_name <- cohort_name
   
   if (length(cohort_bug) == 1L) {
      clean_bug_name <- gsub('([A-Z])_([a-z]+)(_[A-Z]+)?', '\\1. \\2', cohort_bug)
      if (any(grepl('([A-Z])_([a-z]+)_([A-Z]+)$', cohort_bug))) 
         clean_bug_name <- paste0(clean_bug_name, gsub('([A-Z])_([a-z]+)_?([A-Z]+)?', ' (\\3)', cohort_bug))
      clean_bug_name <- paste(clean_bug_name, collapse=', ')
      if ('switch' %in% names(cohort_info[[cohort]])) clean_bug_name <- paste0(clean_bug_name, ' switch')   
   } else {
      clean_bug_name <- 'ampC Producers'
   }
   cohort_info[[cohort]]$clean_bug_name <- clean_bug_name
}

# Save it.
save(cohort_info, file = '~/Desktop/EHR-mining/BuildCohortDatasets/CohortInfo/treatment/cohort_info.Rdata')









# This is old...but coul
#### RESISTANCE COHORT SAMPLES SIZES ####
blood_astDF <- astDF_og %>%
   filter(
      SITE == 'blood',
      lubridate::year(ORDER_DAY) %in% years,
      BUG %in% unname(unlist(bug_groups))
   )
blood_astDF$DAYS_SINCE_PRV <- blood_astDF$ORDER_DAY - lag(blood_astDF$ORDER_DAY)
blood_astDF$DAYS_SINCE_PRV[c(1, which(blood_astDF$PERSON_ID != lag(blood_astDF$PERSON_ID)))] <- NA_integer_
blood_astDF <- blood_astDF %>%
   filter(is.na(DAYS_SINCE_PRV) | DAYS_SINCE_PRV > 60L) %>%
   select(-DAYS_SINCE_PRV, -ORDER_DATE, -RESULT_DATE, -RESULT_DAY, -SITE) %>%
   mutate(across(CEFEPIME:ESBL, ~ ifelse(is.na(.), -1L, .)))

top_bugs <- blood_astDF %>% count(BUG, sort=TRUE)

while (TRUE) {
   print(top_bugs, n=15)
   b <- as.integer(strsplit(readline('Choose bug row: '), split=',')[[1]])
   print(top_bugs[b,])
   bug <- top_bugs$BUG[b]
   
   t <- sapply(X = blood_astDF %>% 
             filter(BUG %in% bug) %>% 
             select(CEFEPIME:ESBL), 
          FUN = function(x) {
             return(c(
                'U' = sum(x == -1L),
                'R' = sum(x == 1L),
                'S' = sum(x == 0L)
             ))
          })
   t <- data.frame(t(t))
   t <- t[order(t$U), ]
   t$ABX <- rownames(t)
   t <- tibble(t)
   t <- t %>%
      relocate(ABX, .before=U) %>%
      mutate(
         tested = R + S,
         pct_tested = tested / (U + R + S),
         pct_susceptible = ifelse(tested > 0, S / tested, NA),
         balance = ifelse(tested > 0, 1 - abs((S - R) / tested), NA),
         # Weighted score: adjust weights as needed
         rank_score = pct_tested * 0.4 + pct_susceptible * 0.3 + balance * 0.5
      ) %>%
      arrange(desc(rank_score))
   print(t)
   x <- readline('Next?')
}




cohort_info <- list(
   list(bug = 'E_coli',       abx = 'TRIMETHOPRIM/SULFAMETHOXAZOLE', estimand = 'ATT'),
   list(bug = 'E_coli',       abx = 'AMPICILLIN',                    estimand = 'ATE'), 
   list(bug = 'E_coli',       abx = 'AMPICILLIN/SULBACTAM',          estimand = 'ATE'), 
   list(bug = 'E_coli',       abx = 'CIPROFLOXACIN',                 estimand = 'ATT'), 
   list(bug = 'E_coli',       abx = 'CEFAZOLIN',                     estimand = 'ATT'), 
   list(bug = 'K_pneumoniae', abx = 'TRIMETHOPRIM/SULFAMETHOXAZOLE', estimand = 'ATT'),
   list(bug = 'P_aeruginosa', abx = 'CIPROFLOXACIN',                 estimand = 'ATT'), 
   list(bug = 'P_aeruginosa', abx = 'AZTREONAM',                     estimand = 'ATT'), 
   list(bug = 'P_mirabilis',  abx = 'CIPROFLOXACIN',                 estimand = 'ATE'),
   list(bug = 'P_mirabilis',  abx = 'TRIMETHOPRIM/SULFAMETHOXAZOLE', estimand = 'ATE'),
   list(bug = 'S_aureus',     abx = 'OXACILLIN',                     estimand = 'ATE'),
   list(bug = 'S_aureus',     abx = 'CLINDAMYCIN',                   estimand = 'ATE'), 
   list(bug = 'E_faecalis',   abx = 'RIFAMPIN',                      estimand = 'ATE'), 
   list(bug = 'E_faecium',    abx = 'VANCOMYCIN',                    estimand = 'ATE')
   #list(bug = 'E_faecium',    abx = 'AMPICILLIN',                    estimand = 'ATE')
)

save(cohort_info, file='~/Desktop/EHR-mining/Scripts/AnalysisScripts/BuildCohortDatasets/info/resistance/cohort_info.Rdata')








x <- df_og2 %>%
   filter(
      BUG == 'S_agalactiae',
      ABX_TIME > 0 & ABX_TIME <= 5
   ) %>%
   mutate(ABX_DAY = floor(ABX_TIME)) %>%
   select(PERSON_ID, ORDER_DAY, ABX, ABX_DAY) %>%
   distinct()

x %>% count(ABX, sort=TRUE)

x %>%
   filter(ABX %in% c('VANCOMYCIN', 'CEFTRIAXONE', 'PIPERACILLIN/TAZOBACTAM', 'CEFAZOLIN', 'CEFEPIME', 'PENICILLIN G')) %>%
   mutate(X = 1L) %>%
   summarise(
      n = n(),
      .by = c(ABX, ABX_DAY)
   ) %>%
   tidyr::pivot_wider(
      names_from = ABX_DAY,
      values_from = n
   )









