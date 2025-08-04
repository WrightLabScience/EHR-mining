library(dplyr)
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/data_path_name.Rdata')
load(file = paste0(data_path_name, 'ASTs_AbxAdmin_blood_2017_2023_imputed_flagged_WIDE_survival.Rdata'))
load(file = paste0(data_path_name, 'LABS_cleaned_blood_2017subset.Rdata'))
inf <- empDF %>% 
   filter(substr(ORDER_DAY, 1, 4) == '2017') %>% # 5,755 rows
   filter(AGE >= 18) %>%
   summarise(RESULT_DAY = min(RESULT_DAY), ###################################### CHANGE THE REST FROM HERE
             .by = c(PERSON_ID, ORDER_DAY, DEATH_DATE, AGE)) %>%
   relocate(DEATH_DATE, AGE, .after=PERSON_ID) # 4,930 rows

patient_id_intersection <- intersect(labs$PERSON_ID, empDF$PERSON_ID) # 3,324
inf <- inf %>% filter(PERSON_ID %in% patient_id_intersection) # 3,518 rows
patient_id_intersection <- intersect(labs$PERSON_ID, inf$PERSON_ID) # 3,118
labs <- labs %>% filter(PERSON_ID %in% patient_id_intersection) # 3.2 million


source('~/Desktop/EHR/EHR work/config_file.R')

chunks <- list(1:1000, 1001:2000, 2001:3000, 3001:length(patient_id_intersection))
encs <- tibble()
for (i in seq_along(chunks)) {
   print(i)
   encs <- rbind(encs,
                 tbl(conn, in_schema('AMB_ETL', 'SENS_ENCOUNTER_VW')) %>%
                    filter(PERSON_ID %in% local(patient_id_intersection[chunks[[i]]]),
                           substr(ADMIT_DATE, 7, 10) == '2017') %>%
                    collect())
}
rm(patient_id_intersection, i, chunks)

nrow(encs) # 7,859 rows
length(unique(encs$PERSON_ID)) # 3,080 people

# basic cleaning
encs <- encs %>% 
   mutate(across(contains('DATE'), ~ strptime(., format='%m/%d/%Y %T'))) %>%
   arrange(PERSON_ID, ADMIT_DATE) %>%
   mutate(ADMIT_DAY = as.Date(substr(ADMIT_DATE,1,10)),
          DISCHARGE_DAY = as.Date(substr(DISCHARGE_DATE,1,10))) %>%
   distinct() # 7,851 rows
encs_og <- encs

## COMBINE OVERLAPPING / ABUTTING ENCOUNTERS
encs$SINCE <- as.integer(encs$ADMIT_DAY - lag(encs$DISCHARGE_DAY))
encs$SINCE[encs$PERSON_ID != lag(encs$PERSON_ID)] <- NA
encs$UNTIL <- lead(encs$SINCE)
encs$UNTIL[encs$PERSON_ID != lead(encs$PERSON_ID)] <- NA
encs$COMB_NEXT <- encs$UNTIL <= 2
encs$COMB_PREV <- encs$SINCE <= 2

# use a run length encoding procedure to determine where the runs of TRUE are
# TRUE meaning the current row ought to be combined with the next row
r <- rle(encs$COMB_NEXT)
lengthsOfRuns <- r$lengths[which(r$values)]
startsOfRuns <- cumsum(r$lengths)[which(r$values)] - lengthsOfRuns + 1
endsOfRuns <- startsOfRuns + lengthsOfRuns - 1
endsOfRuns <- endsOfRuns + encs$COMB_PREV[endsOfRuns + 1]
runs <- mapply(FUN = ':', startsOfRuns, endsOfRuns, SIMPLIFY = FALSE) # 505

start <- Sys.time()
for (i in seq_along(runs)) {
   block <- runs[[i]]
   encs$ADMIT_DATE[block] <- min(encs$ADMIT_DATE[block])
   encs$DISCHARGE_DATE[block] <- max(encs$DISCHARGE_DATE[block])
}
print(Sys.time() - start) # < 1 second

# Step 2: check the dates after runs to see if they fall inside the previous slot
encs$ADMIT_DAY <- as.Date(substr(encs$ADMIT_DATE, 1, 10))
encs$DISCHARGE_DAY <- as.Date(substr(encs$DISCHARGE_DATE, 1, 10))

# rows at the ends of runs where the next row needs to be combined with that run
w <- endsOfRuns[which(as.integer(encs$ADMIT_DAY[endsOfRuns + 1] - encs$DISCHARGE_DAY[endsOfRuns]) <= 2 & encs$PERSON_ID[endsOfRuns + 1] == encs$PERSON_ID[endsOfRuns])] + 1
if (length(w) > 0) {
   count <- 1
   while(length(w) > 0) {
      print(count)
      encs$ADMIT_DATE[w] <- encs$ADMIT_DATE[w - 1]
      encs$DISCH_DATE[w] <- encs$DISCH_DATE[w - 1]
      encs$ADMIT_DAY <- as.Date(substr(encs$ADMIT_DATE, 1, 10))
      encs$DISCH_DAY <- as.Date(substr(encs$DISCH_DATE, 1, 10))
      w <- w[which(as.integer(encs$ADMIT_DAY[w+1] - encs$DISCH_DAY[w]) <= 1 & encs$PERSON_ID[w+1] == encs$PERSON_ID[w])] + 1
      count <- count + 1
   }
} # w is length 0!


# combine rows with same ADMISSION and DISCHARGE DATES
encs <- encs %>%
   select(PERSON_ID, ADMIT_DATE, DISCHARGE_DATE, ENCOUNTER_TYPE, FACILITY, ADMIT_SOURCE) %>%
   distinct() %>%
   summarise(ENCOUNTER_TYPE = paste(sort(unique(ENCOUNTER_TYPE)), collapse=' + '),
             FACILITY = paste(sort(unique(FACILITY)), collapse=' + '),
             ADMIT_SOURCE = paste(sort(unique(ADMIT_SOURCE)), collapse=' + '),
             .by = c(PERSON_ID, ADMIT_DATE, DISCHARGE_DATE)) # 7,295 rows

save(encs, file = paste0('ENCOUNTERS_cleaned_blood_2017subset.Rdata'))












