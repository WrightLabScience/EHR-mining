library(dplyr)
load(file = '~/Desktop/EHR/EHR work/RdataFiles/R01/bact_empiric_flag_ast_abx.Rdata')
source('~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/CreateNamedAbxAbbreviations.R')

# We want to find signals of association between prior patient/facility/time-specific
# antibiotic use + resistance + pathogens --> modern empiric Rx + resistance + pathogens

# To do so, we need: 
#     historical data on prescriptions + some way to encode
#     historical data on resistance + some way to encode


##### data preprocessing #####
df <- df %>% select(PERSON_ID, ORDER_DAY, ABX, ABXC, BUG, BUGC, !!unique(df$ABX)) %>% distinct()
df$ABX <- unname(abbr[df$ABX])
df$BUG <- gsub('^([A-Z])[a-z]+ ([a-z]+)$', '\\1_\\2', df$BUG)
# flag by visit?
##### END #####


df %>%
   select(PERSON_ID, ORDER_DAY, ABX) %>%
   distinct() %>%
   summarise(
      PENEM = as.integer(any(grepl('PENEM$', ABX))),
      .by = c(PERSON_ID, ORDER_DAY)
   )

getDF <- function(var, prefix='') {
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

##### OUTCOMES #####
# outDF
getDF('ABX', 'Rx_')
getDF('ABXC', 'Rx_')
getDF('BUG')
getDF('BUGC')
sapply(df %>% select(!c(ABX, ABXC, BUG, BUGC, PERSON_ID, ORDER_DAY)),
       FUN = function(x) )





tail(sort(table(df$ABX)), n=12)
tail(sort(table(df$ABXC)), n=12)


























