library(dplyr)
library(survival)
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_clean_ASTs.Rdata')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_DEMO.Rdata')

data.frame(sort(sapply(astDF %>%
   filter(BUG == 'Klebsiella pneumoniae') %>%
   select(CEFEPIME:ESBL), function(x) {
      n <- sum(!is.na(x))
      n
      #c(sum(x, na.rm=T) / n, n)
   })))

bug <- 'Escherichia coli'
abx <- c('GENTAMICIN', 'PIPERACILLIN/TAZOBACTAM', 'TRIMETHOPRIM/SULFAMETHOXAZOLE', 'CEFTRIAXONE', 'ESBL')[5]

astDF %>%
   filter(BUG == bug) %>%
   select(PERSON_ID, ORDER_DAY, SITE, !!abx) %>%
   rename(abx = !!abx) %>%
   inner_join(
      x = .,
      y = dth,
      by = join_by(PERSON_ID)
   ) %>%
   mutate(
      female = GENDER == 'FEMALE',
      age = as.integer(ORDER_DAY - DOB) / 365,
      time = as.integer(DEATH_DATE - ORDER_DAY),
      status = as.integer(!is.na(time) & time < 30L)
   ) %>% # summarise(n=n(), sum(time < 30L, na.rm=T) / n, .by=abx)
   coxph(Surv(time, status) ~ abx + SITE + age + female,
         data = .) %>%
   summary()












