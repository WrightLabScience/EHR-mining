#### JOIN DEATH DATA WITH CONCORDANCE ANALYSIS DATA ####
library(dplyr)
load(file = '~/Desktop/EHR/EHR work/ALL/Rdata/ASTs_AbxAdmin_blood_2017_imputed_flagged.Rdata')
load(file = '~/Desktop/EHR/EHR work/Bacteremia2017/data/CleanedDemographics.Rdata')

length(unique(empDF$PERSON_ID)) # 4,007
length(unique(dths$PERSON_ID)) # 4,546
length(intersect(empDF$PERSON_ID, dths$PERSON_ID)) # 2,660

empDF %>%
   #mutate(CoNS = grepl('Staphylococcus|Coag', BUG) & !grepl('aureus', BUG)) %>%
   mutate(NoEmpGiven = lengths(EMPIRIC0) == 0L) %>%
   count(BUG, NoEmpGiven) %>%
   filter(n > 10L) %>%
   tidyr::pivot_wider(values_from = n, names_from = NoEmpGiven) %>%
   na.omit() %>%
   mutate(prop = `TRUE` / (`TRUE` + `FALSE`)) %>%
   arrange(desc(prop))

empDF <- empDF %>%
   left_join(x = .,
             y = dths,
             by = join_by(PERSON_ID),
             relationship = 'many-to-one') %>%
   mutate(AGE = as.integer(ORDER_DATE - DOB) / 365) %>%
   mutate(SURV_TIME = as.integer(DEATH_DATE - ORDER_DATE)) %>%
   filter(!is.na(GENDER))

x <- empDF$AGE
hist(x, breaks=diff(range(x)), xlab='Years', main = 'Age')

x <- empDF$SURV_TIME
table(is.na(x)) # about half
table(x <= 30) #  ~9.6%
table(x <= 60) # ~12.7%
table(x <= 365)
hist(x, breaks=diff(range(x, na.rm=T)), xlim=c(0, 365), xlab='Days', main='Survival time after order')


empDF %>%
   filter(AGE < 65) %>%
   filter(SURV_TIME <= 730) %>%
   filter(FLAG01 %in% c('CONCORDANT', 'DISCORDANT', 'No empiric therapy given')) %>%
   summarise(n=n(),
             median = median(SURV_TIME, na.rm=T),
             mean = mean(SURV_TIME, na.rm=T),
             sd = sd(SURV_TIME, na.rm=T),
             .by = FLAG01)
empDF %>%
   filter(AGE < 65) %>%
   filter(SURV_TIME <= 730) %>%
   filter(FLAG0 %in% c('CONCORDANT', 'DISCORDANT', 'No empiric therapy given')) %>%
   summarise(n = n(),
             median = median(SURV_TIME, na.rm=T),
             mean = mean(SURV_TIME, na.rm=T),
             sd = sd(SURV_TIME, na.rm=T),
             .by = FLAG0)
empDF %>%
   filter(AGE < 65) %>%
   filter(SURV_TIME <= 730) %>%
   filter(FLAG1 %in% c('CONCORDANT', 'DISCORDANT', 'No empiric therapy given')) %>%
   summarise(n = n(),
             mean = mean(SURV_TIME, na.rm=T),
             median = median(SURV_TIME, na.rm=T),
             sd = sd(SURV_TIME, na.rm=T),
             .by = FLAG1)
empDF %>%
   filter(AGE < 65) %>%
   filter(SURV_TIME <= 365) %>%
   filter(TFLAG %in% c('CONCORDANT', 'DISCORDANT', 'No targeted therapy given')) %>%
   summarise(n = n(),
             median = median(SURV_TIME, na.rm=T),
             mean = mean(SURV_TIME, na.rm=T),
             sd = sd(SURV_TIME, na.rm=T),
             .by = TFLAG)


empDF %>%
   filter(TFLAG %in% c('CONCORDANT', 'DISCORDANT', 'No targeted therapy given')) %>%
   summarise(n = n(),
             median = median(AGE, na.rm=T),
             mean = mean(AGE, na.rm=T),
             sd = sd(AGE, na.rm=T),
             .by = TFLAG)
empDF %>%
   filter(FLAG1 %in% c('CONCORDANT', 'DISCORDANT', 'No targeted therapy given')) %>%
   summarise(n = n(),
             median = median(AGE, na.rm=T),
             mean = mean(AGE, na.rm=T),
             sd = sd(AGE, na.rm=T),
             .by = FLAG1)
empDF %>%
   filter(FLAG0 %in% c('CONCORDANT', 'DISCORDANT', 'No targeted therapy given')) %>%
   summarise(n = n(),
             median = median(AGE, na.rm=T),
             mean = mean(AGE, na.rm=T),
             sd = sd(AGE, na.rm=T),
             .by = FLAG0)
empDF %>%
   filter(FLAG01 %in% c('CONCORDANT', 'DISCORDANT', 'No targeted therapy given')) %>%
   summarise(n = n(),
             median = median(AGE, na.rm=T),
             mean = mean(AGE, na.rm=T),
             sd = sd(AGE, na.rm=T),
             .by = FLAG01)




# patients not given empiric treatments are generally slightly younger
empDF %>%
   mutate(not_given = grepl('given', FLAG01)) %>%
   summarise(median = median(AGE, na.rm=T),
             sd = sd(AGE, na.rm=T),
             .by = not_given)
# patients not given targeted therapies are basically same age
empDF %>%
   mutate(not_given = grepl('given', TFLAG)) %>%
   summarise(median = median(AGE, na.rm=T),
             sd = sd(AGE, na.rm=T),
             .by = not_given)




empDF %>%
   filter(SURV_TIME <= 365) %>%
   mutate(CoNS = grepl('Staphylococcus|Coag', BUG) & !grepl('aureus', BUG)) %>%
   summarise(median = median(SURV_TIME, na.rm=T),
             mean = mean(SURV_TIME, na.rm=T),
             .by = CoNS)

empDF %>%
   filter(SURV_TIME <= 365) %>%
   filter(grepl('Staphylococcus|Coag', BUG) & !grepl('aureus', BUG)) %>%
   mutate(LONG_DELAY = DELAY > 5) %>%
   summarise(n = n(),
             median = median(SURV_TIME, na.rm=T),
             mean = mean(SURV_TIME, na.rm=T),
             .by = LONG_DELAY)


# DID more resistance = higher rates of discordance?
   

empDF %>%
   mutate(CoNS = grepl('Staphylococcus|Coag', BUG) & !grepl('aureus', BUG),
          SA = BUG == 'Staphylococcus aureus') %>%
   summarise(median = median(DELAY),
             .by = c(CoNS, SA))




