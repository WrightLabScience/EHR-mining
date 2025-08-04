source(file = '~/Desktop/EHR/EHR work/config_file.R')
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/data_path_name.Rdata')

vits <- tbl(conn, in_schema('AMB_ETL', 'BACT_VITALS_VW')) %>%
   collect()

vits <- vits %>%
   mutate(CONTACT_DATE = strptime(CONTACT_DATE, format='%m/%d/%Y %T')) %>%
   arrange(PERSON_ID, CONTACT_DATE) %>%
   distinct()


# HEIGHT
x <- unique(vits$HEIGHT[!is.na(vits$HEIGHT)])
w1 <- grep('^[0-9]{2,3}$', x)
w2 <- grep('^[0-9]+\\.', x)
w3 <- grep("^[0-9]' [0-9]*(\\.[0-9]*)?", x)
length(intersect(w1, w2))
length(intersect(w3, w2))
length(intersect(w1, w3))
w <- union(w1, union(w3, w2))
x[-w]

vits <- vits %>%
   mutate(NUM_HEIGHT = as.numeric(HEIGHT))
# vits_og <- vits

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

x <- vits$CALC_HEIGHT
range(x, na.rm=T)
hist(x, breaks=diff(range(x, na.rm=T)))



# WEIGHT
vits$WEIGHT <- as.numeric(vits$WEIGHT)

table(vits$WEIGHT > 300) # > half?
table(vits$WEIGHT > 1000) # > half?

hist(vits$WEIGHT, breaks=diff(range(vits$WEIGHT, na.rm=T)), xlim=c(0, 280))
hist(vits$WEIGHT, breaks=diff(range(vits$WEIGHT, na.rm=T)), xlim=c(280, 1050), ylim=c(0,40))
hist(vits$WEIGHT, breaks=diff(range(vits$WEIGHT, na.rm=T)), xlim=c(1050, 5000))
hist(vits$WEIGHT, breaks=diff(range(vits$WEIGHT, na.rm=T)), xlim=c(1900, 4000))
# abline(v=c(1999, 2079, 2159, 2239)) # WTF!!!

t <- vits %>%
   reframe(sm = median(WEIGHT[WEIGHT < 280], na.rm=T),
           me = median(WEIGHT[WEIGHT >= 280 & WEIGHT < 1050], na.rm=T),
           la = median(WEIGHT[WEIGHT >= 1050], na.rm=T),
           .by=PERSON_ID) %>%
   mutate(r1 = me / sm,
          r2 = la / me,
          r3 = la / sm)

hist(t$r1, breaks=diff(range(t$r1, na.rm=T)))
hist(t$r2, breaks=diff(range(t$r2, na.rm=T)))
hist(c(t$r1, t$r2), breaks=diff(range(c(t$r1, t$r2), na.rm=T)))
hist(t$r3, breaks=diff(range(t$r3, na.rm=T)))

sum(!is.na(t$r1)) # 124
sum(!is.na(t$r2)) # 55
sum(!is.na(t$r3)) # 2972

median(t$r1, na.rm=T) # 2.8
median(t$r2, na.rm=T) # 2.6
median(c(t$r1, t$r2), na.rm=T) # 2.8
median(t$r3, na.rm=T) # 35.4

# determine each person's bin-wise median
vits <- vits %>%
   group_by(PERSON_ID) %>%
   mutate(sm = median(WEIGHT[WEIGHT < 280], na.rm=T),
          me = median(WEIGHT[WEIGHT >= 280 & WEIGHT < 1050], na.rm=T),
          la = median(WEIGHT[WEIGHT >= 1050], na.rm=T)) %>%
   ungroup()

fact1 <- median(c(t$r1, t$r2), na.rm=T)
fact2 <- median(t$r3, na.rm=T)

vits <- vits %>%
   mutate(CALC_WEIGHT = case_when(
      WEIGHT >= 280 & WEIGHT < 1050 ~ WEIGHT / fact1,
      WEIGHT >= 1050 ~ WEIGHT / fact2,
      .default = WEIGHT
   ))

hist(vits$CALC_WEIGHT, breaks=diff(range(vits$CALC_WEIGHT, na.rm=T)))

x <- vits %>%
   reframe(WEIGHT = median(WEIGHT[WEIGHT < 280]),
           CWEIGHT = median(CALC_WEIGHT),
           .by = PERSON_ID) %>%
   mutate(D = CWEIGHT - WEIGHT)

hist(x$D, breaks=diff(range(x$D, na.rm=T)))




vits <- vits %>%
   filter(CALC_WEIGHT < 500) %>%
   filter(CALC_HEIGHT > )
select(!c(sm, me, la, WEIGHT, HEIGHT))


vits %>% count(is.na(CALC_WEIGHT))
vits %>% count(is.na(CALC_HEIGHT)) # 80,435 missing




















