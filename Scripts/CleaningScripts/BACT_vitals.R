source(file = '~/Desktop/EHR/EHR work/config_file.R')
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/data_path_name.Rdata')
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/plots_path_name.Rdata')

vits <- tbl(conn, in_schema('AMB_ETL', 'BACT_VITALS_VW')) %>%
   collect()

n <- vits %>% count(PERSON_ID) %>% select(n) %>% unlist()
length(n[n == 1L]) / length(n) # 3%
length(n[n > 1L & n < 30L]) / length(n) # 17%
hist(n, breaks=diff(range(n)), xlim=c(0,100))
plot(cumsum(n) / sum(n))


# HEIGHT
vits %>%
   group_by(PERSON_ID) %>%
   filter(any(!is.na(HEIGHT))) # 4,050 / 4,300 have a height

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
   mutate(CALC_HEIGHT = case_when(
      grepl('^[0-9]{4,}$', HEIGHT) ~ NA,
      HEIGHT %in% as.character(3:6) ~ as.integer(HEIGHT) * 12,
      grepl('^[0-9]\\.[0-9]$', HEIGHT) ~ NA,
      grepl('^[0-9]{2,3}(\\.)?', HEIGHT) & as.numeric(HEIGHT) < 90 ~ as.numeric(HEIGHT),
      grepl('^[0-9]{2,3}(\\.)?', HEIGHT) & as.numeric(HEIGHT) >= 90 ~ as.numeric(HEIGHT) / 2.54,
      grepl("^[0-9]' [0-9]*(\\.[0-9]*)?\"", HEIGHT) ~ as.integer(gsub("^([0-9])' ([0-9]*(\\.[0-9]*)?)\"", '\\1', HEIGHT)) * 12 + as.numeric(gsub("^([0-9])' ([0-9]*(\\.[0-9]*)?)\"", '\\2', HEIGHT))
   ))

x <- vits$CALC_HEIGHT
range(x, na.rm=T)
hist(x, breaks=diff(range(x, na.rm=T)))

grep('^[0-9]{2,3}(\\.)?', x, value=T)

vits %>%
   summarise(n = n(),
             min = min(HEIGHT),
             median = median(HEIGHT),
             max = max(HEIGHT),
             .by = PERSON_ID)
