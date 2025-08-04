library(dplyr)
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/data_path_name.Rdata')
load(file = paste0(data_path_name, 'ASTs_AbxAdmin_blood_2017_2023_imputed_flagged_WIDE_survival.Rdata'))
load(file = paste0(data_path_name, 'LABS_cleaned_blood_2017subset.Rdata'))
load(file = paste0(data_path_name, 'ENCOUNTERS_cleaned_blood_2017subset.Rdata'))
infs <- empDF %>% 
   filter(substr(ORDER_DAY, 1, 4) == '2017') %>% # 5,755 rows
   filter(AGE >= 18) %>%
   # select(PERSON_ID, AGE, DEATH_DATE, GENDER, ORDER_DATE, ORDER_DAY, FLAGe, ABX_TAR, BUG, RESISTANT, SUSCEPTIBLE) %>%
   select(!c(CEFEPIME:DELAFLOXACIN)) %>%
   filter(n() == 1L, .by = c(PERSON_ID, ORDER_DAY)) %>%
   distinct() # 4,562

# take patient intersect
patient_id_intersection <- intersect(labs$PERSON_ID, intersect(infs$PERSON_ID, encs$PERSON_ID)) # 2,861
infs <- infs %>% filter(PERSON_ID %in% patient_id_intersection) # 3,131 rows
labs <- labs %>% filter(PERSON_ID %in% patient_id_intersection) # 2.87 million
encs <- encs %>% filter(PERSON_ID %in% patient_id_intersection) # 6,799
rm(patient_id_intersection)


# join encounters and ASTs
df <- infs %>%
   left_join(y = encs %>%
                   mutate(ADMIT_DAY = as.Date(substr(ADMIT_DATE,1,10)),
                          DISCHARGE_DAY = as.Date(substr(DISCHARGE_DATE,1,10)),
                          JOIN_START = ADMIT_DAY - 5,
                          JOIN_END = DISCHARGE_DAY + 5),
             by = join_by(PERSON_ID,
                          ORDER_DAY >= JOIN_START,
                          ORDER_DAY <= JOIN_END)) %>%
   select(-JOIN_START, -JOIN_END)

# there are some ASTs being joined to multiple close encounters
# df %>%
#    filter(n() > 1L, .by=c(PERSON_ID, ORDER_DAY)) %>% 
#    mutate(SINCE = ifelse(test = PERSON_ID == lag(PERSON_ID), 
#                          yes = as.integer(ADMIT_DAY - lag(DISCHARGE_DAY)), 
#                          no = NA)) %>%
#    select(PERSON_ID, ORDER_DAY, ADMIT_DAY, DISCHARGE_DAY)

# combine these encounters
df <- rbind(df %>% 
               filter(n() > 1L, .by = c(PERSON_ID, ORDER_DAY)) %>% 
               filter(ORDER_DAY >= ADMIT_DAY & ORDER_DAY <= DISCHARGE_DAY),
            df %>%
               filter(n() == 1L, .by = c(PERSON_ID, ORDER_DAY))) %>%
   arrange(PERSON_ID, ORDER_DAY) %>%
   filter(!is.na(ADMIT_DATE)) %>% # lost ~100
   filter(DEATH_DATE >= DISCHARGE_DAY | is.na(DEATH_DATE)) # lose 1 mistake

df$LENGTH_OF_STAY <- as.numeric(lubridate::as.duration(df$DISCHARGE_DATE - df$ADMIT_DATE)) / 86400

hist(df$LENGTH_OF_STAY, breaks=diff(range(df$LENGTH_OF_STAY, na.rm=T)), xlim=c(0, 30))
median(df$LENGTH_OF_STAY, na.rm=T) # 8 days
mean(df$LENGTH_OF_STAY, na.rm=T) # 13 days

sum(df$DEATH_DATE == df$DISCHARGE_DAY, na.rm=T) / nrow(df)      # 8%
sum(df$DEATH_DATE <= df$ORDER_DAY + 30, na.rm=T) / nrow(df)     # 11.6%
sum(df$DEATH_DATE <= df$DISCHARGE_DAY + 30, na.rm=T) / nrow(df) # 13.8%



########## JOIN WITH LAB RESULTS?? ############
labs <- labs %>% 
   mutate(ORDER_DAY = as.Date(substr(ORDER_DATE,1,10))) %>%
   rename(LAB_ORDER_DAY = ORDER_DAY) %>%
   select(-ORDER_DATE, -UNITS) %>%
   distinct()


# # visualize GLUCOSE measurements in one person
# lab <- labs %>% filter(LAB == 'GLUCOSE')
# pid <- '1000000541'
# glu <- labs %>% rename(LAB_ORDER_DAY = ORDER_DAY) %>% filter(LAB == 'GLUCOSE', PERSON_ID == pid) %>% arrange(ORDER_DATE)
# window_size <- 6
# x <- lubridate::decimal_date(glu$ORDER_DATE[w])
# y <- glu$RESULT_VALUE[w]
# plot(x=x, y=y, pch=16, cex=0.5, xaxt='n', ylab='Glucose (mg/dL)')
# lines(x=x, y=stats::filter(y, rep(1/window_size, window_size)))
# axis(side=1, at=lubridate::decimal_date(unique(glu$ORDER_DAY)), labels=unique(glu$ORDER_DAY))
# points(x=lubridate::decimal_date(c(df$ORDER_DAY[df$PERSON_ID == pid], df$RESULT_DAY[df$PERSON_ID == pid])), y=c(80, 80), pch=16, col='red')
# abline(v = lubridate::decimal_date(c(df$ADMIT_DATE[df$PERSON_ID == pid], df$DISCHARGE_DATE[df$PERSON_ID == pid])))
# rm(x, y, w, window_size)



# Join LABs
df <- df %>%
   mutate(JOIN_START = ADMIT_DAY - 4,
          JOIN_END = DISCHARGE_DAY + 4) %>%
   left_join(y = labs,
             by = join_by(PERSON_ID,
                          JOIN_START <= LAB_ORDER_DAY,
                          JOIN_END >= LAB_ORDER_DAY)) %>%
   select(-ADMIT_DAY, -DISCHARGE_DAY, -JOIN_START, -JOIN_END, -LAB_ORDER_DAY)

start <- Sys.time()
df <- df %>%
   tidyr::pivot_wider(
      values_from = RESULT_VALUE,
      names_from = LAB,
      values_fn = median
   )
df <- df %>% select(-`NA`)
print(Sys.time() - start) # 1.36 million rows - 12 seconds
rm(start)



# join with vitals
source(file = '~/Desktop/EHR/EHR work/config_file.R')
vits <- tbl(conn, in_schema('AMB_ETL', 'BACT_VITALS_VW')) %>%
   collect() %>%
   mutate(CONTACT_DATE = strptime(CONTACT_DATE, format='%m/%d/%Y %T')) %>%
   arrange(PERSON_ID, CONTACT_DATE) %>%
   distinct() %>%
   mutate(NUM_HEIGHT = as.numeric(HEIGHT))

vits <- vits_og
# vits <- vits %>%
#    left_join(y = df %>% select(PERSON_ID, DOB) %>% distinct(),
#              relationship = 'many-to-one',
#              by = join_by(PERSON_ID)) %>%
#    filter(!is.na(DOB)) %>%
#    mutate(AGE_AT_MEASURE = as.numeric(lubridate::as.duration(as.Date(substr(CONTACT_DATE,1,10)) - DOB)) / 86400 / 365)

# fix height
vits <- vits %>%
   mutate(CALC_HEIGHT = case_when(
      NUM_HEIGHT > 999 ~ NA,
      NUM_HEIGHT < 6 & NUM_HEIGHT > 3 ~ NUM_HEIGHT * 12,
      NUM_HEIGHT < 3 ~ NA,
      grepl('^[0-9]{2,3}(\\.)?', HEIGHT) & NUM_HEIGHT < 90 ~ NUM_HEIGHT,
      grepl('^[0-9]{2,3}(\\.)?', HEIGHT) & NUM_HEIGHT >= 90 ~ NUM_HEIGHT / 2.54
   ))
w <- which(is.na(vits$CALC_HEIGHT) & grepl("^[0-9]' [0-9]*(\\.[0-9]*)?\"", vits$HEIGHT))
vits$CALC_HEIGHT[w] <- as.integer(gsub("^([0-9])' ([0-9]*(\\.[0-9]*)?)\"", '\\1', vits$HEIGHT[w])) * 12 + as.numeric(gsub("^([0-9])' ([0-9]*(\\.[0-9]*)?)\"", '\\2', vits$HEIGHT[w]))

# fix weight
vits$WEIGHT <- as.numeric(vits$WEIGHT)
fact1 <- 2.7
fact2 <- 35.4
vits <- vits %>%
   mutate(CALC_WEIGHT = case_when(
      WEIGHT >= 280 & WEIGHT < 1050 ~ WEIGHT / fact1,
      WEIGHT >= 1050 ~ WEIGHT / fact2,
      .default = WEIGHT
   ))

# pivot longer
vits <- vits %>% select(-WEIGHT, -HEIGHT, -NUM_HEIGHT)
vits <- vits %>%
   #group_by(PERSON_ID) %>%
   #fill(CALC_HEIGHT) %>%
   #fill(CALC_WEIGHT) %>%
   #ungroup() %>%
   tidyr::pivot_longer(
      cols = c(CALC_HEIGHT, CALC_WEIGHT),
      names_to = 'MEASURE',
      values_to = 'VALUE'
   ) %>%
   filter(!is.na(VALUE)) %>%
   mutate(CONTACT_DAY = as.Date(substr(CONTACT_DATE,1,10)))


# JOIN HEIGHT
dfH <- df %>%
   mutate(JOIN_START = as.Date(substr(ADMIT_DATE,1,10)) - 730,
          JOIN_END = as.Date(substr(ORDER_DAY,1,10)) + 3) %>%
   left_join(vits %>% filter(MEASURE == 'CALC_HEIGHT'),
             by=join_by(PERSON_ID,
                        JOIN_START <= CONTACT_DAY,
                        JOIN_END >= CONTACT_DAY), 
             relationship='many-to-many') %>% 
   mutate(TIME_DIFF = as.numeric(lubridate::as.duration(ADMIT_DATE - CONTACT_DATE)) / 86400) %>% # filter(all(is.na(VALUE)), .by=c(PERSON_ID, ORDER_DAY)) # 436 missing all (2 years prior, 3 days after)
   select(-CONTACT_DATE, -CONTACT_DAY) %>%
   rename(HEIGHT = VALUE) %>% select(-MEASURE) %>%
   group_by(PERSON_ID, ORDER_DAY) %>%
   slice_min(TIME_DIFF, with_ties = FALSE) %>%
   ungroup() %>%
   select(-TIME_DIFF)

# JOIN WEIGHT
dfHW <- dfH %>%
   left_join(vits %>% filter(MEASURE == 'CALC_WEIGHT'),
             by = join_by(PERSON_ID,
                          JOIN_START <= CONTACT_DAY,
                          JOIN_END >= CONTACT_DAY),
             relationship='many-to-many') %>%
   select(-JOIN_START, -JOIN_END) %>%
   mutate(TIME_DIFF = as.numeric(lubridate::as.duration(ADMIT_DATE - CONTACT_DATE)) / 86400) %>%
   select(-CONTACT_DATE, -CONTACT_DAY) %>%
   rename(WEIGHT = VALUE) %>% select(-MEASURE) %>%
   group_by(PERSON_ID, ORDER_DAY) %>%
   slice_min(TIME_DIFF, with_ties = FALSE) %>%
   ungroup() %>%
   select(-TIME_DIFF)

dfHW %>% count(is.na(HEIGHT), is.na(WEIGHT)) # 436 missing both height and weight
df <- dfHW


save(df, file = paste0(data_path_name, 'Blood_2017subset_ASTsAbxEncsLabsVits.Rdata'))


# how many visits are missing labs?
lab_miss_person <- apply(df %>% select(ABS_BASOPHILS:PCT_MONOCYTES), 1, function(x) sum(is.na(x)))
barplot(table(lab_miss_person), xlab='Number of labs missing', ylab='Number of patients')


pcnt_each_lab_miss <- sort(sapply(df %>% select(ABS_BASOPHILS:PCT_MONOCYTES), 
                                  function(x) sum(is.na(x))) / nrow(df)) * 100
lab_names <- names(pcnt_each_lab_miss)
par(mar=c(4, 10, 2, 2))
b <- barplot(pcnt_each_lab_miss, horiz=TRUE, names.arg=rep('', length(pcnt_each_lab_miss)), 
             xlim=c(0, 100), xlab='% missing')
text(x=-1, y=b, adj=1, cex=0.7, xpd=NA,
     labels=gsub('_', ' ', stringr::str_to_sentence(lab_names)))


rm(b, pcnt_each_lab_miss, lab_names, lab_miss_person)


# association between age and lab values?
age <- df$AGE
for (i in 12:length(df1)) {
   vals <- df1[[i]]
   l <- lm(vals ~ age)
   plot(age, vals, cex=0.5, pch=16)
   abline(l)
   r2 <- round(summary(l)$adj.r.squared, 4)
   cat(names(df1)[i], 'R^2', r2, '\n')
}



















