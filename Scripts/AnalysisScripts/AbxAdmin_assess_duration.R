### What is the typical duration of antibiotic therapies??

library(dplyr)
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/data_path_name.Rdata')
load(file = paste0(data_path_name, 'ASTs_blood_2015_2023.Rdata'))
load(file = paste0(data_path_name, 'ALL_2017_AbxAdmin.Rdata'))

# FOR NOW ONLY 2017!!!
empDF <- empDF %>% filter(substr(ORDER_DATE, 1, 4) == '2017')     # 57,824 --> 5,908
abxDF <- abxDF %>% filter(PERSON_ID %in% unique(empDF$PERSON_ID)) # 1,288,305 --> 290,855

# abxDF <- abxDF %>%
#    group_by(PERSON_ID, ABX) %>%
#    mutate(SINCE = as.numeric(ADMIN_START_DATE - lag(ADMIN_END_DATE)) / 86400) %>%
#    ungroup()
# x <- abxDF$SINCE
# x <- x[!is.na(x)]
# hist(x, breaks=diff(range(x)), xlim=c(0, 365))

# Figure out which rows are to be combined with adjacent rows
abxDF <- abxDF %>%
   mutate(across(contains('DATE'), ~ as.Date(substr(., 1, 10)))) %>%
   distinct() %>%  # 290,855 --> 154,875
   mutate(PERSON_ID = as.character(PERSON_ID)) %>%
   rename(START_DATE = ADMIN_START_DATE,
          END_DATE = ADMIN_END_DATE) %>% 
   arrange(PERSON_ID, ABX, START_DATE) %>%
   mutate(END_DATE = case_when(
      is.na(END_DATE) ~ START_DATE,
      .default = END_DATE
   ))

GAP_SIZE_DAYS <- 2
abxDF <- abxDF %>%
   mutate(SINCE = as.integer(START_DATE - lag(END_DATE)),
          UNTIL = as.integer(lead(START_DATE) - END_DATE)) %>%
   mutate(COMB_NEXT = UNTIL <= GAP_SIZE_DAYS & PERSON_ID == lead(PERSON_ID) & ABX == lead(ABX),
          COMB_PREV = SINCE <= GAP_SIZE_DAYS & PERSON_ID == lag(PERSON_ID) & ABX == lag(ABX))

# use a run length encoding procedure to determine where the runs of TRUE are
# TRUE meaning the current row ought to be combined with the next row
r <- rle(abxDF$COMB_NEXT)
lengthsOfRuns <- r$lengths[which(r$values)]
startsOfRuns <- cumsum(r$lengths)[which(r$values)] - lengthsOfRuns + 1
endsOfRuns <- startsOfRuns + lengthsOfRuns - 1
endsOfRuns <- endsOfRuns + abxDF$COMB_PREV[endsOfRuns + 1]
runs <- mapply(FUN = ':', startsOfRuns, endsOfRuns, SIMPLIFY = FALSE)

start <- Sys.time()
for (i in seq_along(runs)) {
   block <- runs[[i]]
   abxDF$START_DATE[block] <- min(abxDF$START_DATE[block])
   abxDF$END_DATE[block] <- max(abxDF$END_DATE[block])
}
print(Sys.time() - start) 
# 1.0 minutes to combine 26,198 (gap=2) runs
# 1.9 minutes for 32,185 runs

# Step 2: check the dates after runs to see if they fall inside the previous slot
# rows at the ends of runs where the next row needs to be combined with that run
w <- endsOfRuns[which(as.integer(abxDF$START_DATE[endsOfRuns + 1] - abxDF$END_DATE[endsOfRuns]) <= 2 
                      & abxDF$PERSON_ID[endsOfRuns + 1] == abxDF$PERSON_ID[endsOfRuns]
                      & abxDF$ABX[endsOfRuns + 1] == abxDF$ABX[endsOfRuns])] + 1
if (length(w) > 0) {
   count <- 1
   while(length(w) > 0) {
      print(count)
      abxDF$START_DATE[w] <- abxDF$START_DATE[w - 1]
      abxDF$END_DATE[w] <- abxDF$END_DATE[w - 1]
      w <- w[which(as.integer(abxDF$START_DATE[w+1] - abxDF$END_DATE[w]) <= 1
                   & abxDF$PERSON_ID[w+1] == abxDF$PERSON_ID[w] 
                   & abxDF$ABX[endsOfRuns + 1] == abxDF$ABX[endsOfRuns])] + 1
      count <- count + 1
   } # 2 loops (gap=1), 0 if gap=2
}

# combine overlapping (now identical) rows
abxDF <- abxDF %>%
   select(-SINCE, -UNTIL, -COMB_NEXT, -COMB_PREV) %>%
   distinct() %>% # 194,034 --> 45,446
   arrange(PERSON_ID, ABX, START_DATE)

# 0 rows - GOOD!
abxDF %>%
   mutate(SINCE = as.integer(START_DATE - lag(END_DATE)), UNTIL = as.integer(lead(START_DATE) - END_DATE)) %>%
   mutate(COMB_NEXT = UNTIL <= GAP_SIZE_DAYS & PERSON_ID == lead(PERSON_ID) & ABX == lead(ABX), COMB_PREV = SINCE <= GAP_SIZE_DAYS & PERSON_ID == lag(PERSON_ID) & ABX == lag(ABX)) %>%
   filter(COMB_NEXT | COMB_PREV)

rm(start, i, block, count, w, runs, endsOfRuns, startsOfRuns, lengthsOfRuns, r, GAP_SIZE_DAYS)



# How long are antibiotic therapy courses?
abxDF <- abxDF %>%mutate(DURATION = as.integer(END_DATE - START_DATE))
x <- abxDF$DURATION
median(x) # 2
mean(x) # 3.3
quantile(x, c(0.1, 0.5, 0.75, 0.9, 0.99))
# 10% 50% 75% 90% 99% 
#   0   2   4   8  25
#x[x >= 42] <- 42
t <- table(x)
barplot(t[1:30], xlab='Days', main='Length of consecutive\nadministration days (by antibiotic)')

barplot(table(abxDF$DURATION[abxDF$ABX == 'VANCOMYCIN'])[1:30], xlab='Days', main='VAN')
barplot(table(abxDF$DURATION[abxDF$ABX == 'PIPERACILLIN/TAZOBACTAM'])[1:30], xlab='Days', main='TZP')
barplot(table(abxDF$DURATION[abxDF$ABX == 'CEFTRIAXONE'])[1:30], xlab='Days', main='CRO')
barplot(table(abxDF$DURATION[abxDF$ABX == 'CEFAZOLIN'])[1:30], xlab='Days', main='CFZ')

abxDF %>%
   filter(n() > 1000L, .by=ABX) %>%
   summarise(n = n(),
             median = median(DURATION),
             mean = mean(DURATION),
             #q99 = quantile(DURATION, 0.99),
             .by = ABX) %>%
   arrange(desc(n))

abxDF <- abxDF %>% select(-DURATION)

abxDF <- abxDF %>%
   group_by(PERSON_ID, ABX) %>%
   mutate(SINCE = as.integer(START_DATE - lag(END_DATE))) %>%
   ungroup()

x <- abxDF$SINCE
table(is.na(x)) # more than half are NA
x <- x[!is.na(x)]
range(x) # 3, 385
median(x) # 25
mean(x) # 48.8
hist(x, breaks=diff(range(x)), xlim=c(0, 50))






