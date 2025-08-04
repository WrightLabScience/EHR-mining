library(dplyr)
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_CLEANED_2015_2024_AbxAdmin.Rdata')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_clean_ASTs.Rdata')
load(file = '~/Desktop/EHR-mining/UsefulDataForCleaning/bug_groups.Rdata')

df <- astDF %>% filter(BUG %in% bug_groups$Enterobacterales, BLOOD)
df <- df %>% 
   group_by(PERSON_ID, ORDER_DAY, ORDER_DATE) %>%
   mutate(DAYS_SINCE_PRV = as.integer(ORDER_DAY - lag(ORDER_DAY))) %>%
   filter(is.na(DAYS_SINCE_PRV) | DAYS_SINCE_PRV > 42L) %>%
   ungroup() %>%
   select(-DAYS_SINCE_PRV)
# df_og <- df
abxDF <- abxDF %>%
   filter(PERSON_ID %in% unique(df$PERSON_ID))


# Join ABXs + ASTs
df <- df %>%
   mutate(
      JOIN_START = ORDER_DATE - 3600 * 48,
      JOIN_END = ORDER_DATE + 3600 * 12
   ) %>%
   left_join(
      x = .,
      y = abxDF,
      by = join_by(
         PERSON_ID,
         JOIN_START <= START_DATE,
         JOIN_END >= START_DATE
      )
   ) %>% 
   filter(!is.na(ABX)) %>%
   select(-JOIN_START, -JOIN_END) %>%
   mutate(
      ABX_DAY = as.integer(START_DAY - ORDER_DAY),
      ABX_TIME = as.numeric(lubridate::as.duration(START_DATE - ORDER_DATE))
   ) %>%
   relocate(START_DATE:START_DAY, .after=RESULT_DAY)

df %>%
   select(PERSON_ID, ORDER_DAY, ABX) %>%
   count(ABX, sort=TRUE)






















