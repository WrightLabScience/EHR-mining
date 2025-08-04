combineOverlappingEncounters <- function(encs, please_arrange=FALSE) {
   if (!any(names(encs) == 'ADMIT_DAY'))
      encs$ADMIT_DAY <- lubridate::as_date(encs$ADMIT_DATE)
   if (!any(names(encs) == 'DISCHARGE_DAY'))
      encs$DISCHARGE_DAY <- lubridate::as_date(encs$DISCHARGE_DATE)
   
   if (please_arrange)
      encs <- encs %>% arrange(PERSON_ID, ADMIT_DATE)
   
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
   runs <- mapply(FUN = ':', startsOfRuns, endsOfRuns, SIMPLIFY = FALSE) # ~7K
   cat('Runs: ', length(runs), '\n')
   
   encs$COMBINED_ADJ <- FALSE
   start <- Sys.time()
   for (i in seq_along(runs)) {
      if (i %% 100 == 0) print(i)
      block <- runs[[i]]
      encs$COMBINED_ADJ[block] <- TRUE
      encs$ADMIT_DATE[block] <- min(encs$ADMIT_DATE[block])
      encs$DISCHARGE_DATE[block] <- max(encs$DISCHARGE_DATE[block])
   }
   print(Sys.time() - start) # ~2 minutes for ~7K runs
   rm(i, block)
   
   # Step 2: check the dates after runs to see if they fall inside the previous slot
   encs$ADMIT_DAY <- as.Date(substr(encs$ADMIT_DATE, 1, 10))
   encs$DISCHARGE_DAY <- as.Date(substr(encs$DISCHARGE_DATE, 1, 10))
   
   # rows at the ends of runs where the next row needs to be combined with that run
   w <- endsOfRuns[which(as.integer(encs$ADMIT_DAY[endsOfRuns + 1] - encs$DISCHARGE_DAY[endsOfRuns]) <= 2 & encs$PERSON_ID[endsOfRuns + 1] == encs$PERSON_ID[endsOfRuns])] + 1
   print(length(w)) # 5
   count <- 1
   while(length(w) > 0) {
      print(count)
      encs$ADMIT_DATE[w] <- encs$ADMIT_DATE[w - 1]
      encs$DISCHARGE_DATE[w] <- encs$DISCHARGE_DATE[w - 1]
      encs$ADMIT_DAY[w] <- as.Date(substr(encs$ADMIT_DATE[w], 1, 10))
      encs$DISCHARGE_DAY[w] <- as.Date(substr(encs$DISCHARGE_DATE[w], 1, 10))
      w <- w[which(as.integer(encs$ADMIT_DAY[w+1] - encs$DISCHARGE_DAY[w]) <= 1 & encs$PERSON_ID[w+1] == encs$PERSON_ID[w])] + 1
      count <- count + 1
   } # 4 rounds
   rm(w, count, r, lengthsOfRuns, startsOfRuns, endsOfRuns, runs, start)
   
   
   # combine rows with same ADMISSION and DISCHARGE DATES
   encs <- encs %>%
      select(PERSON_ID, ADMIT_DATE, DISCHARGE_DATE, ADMIT_DAY, DISCHARGE_DAY, ENCOUNTER_TYPE, FACILITY, ADMIT_SOURCE, ADMIT_TYPE, COMBINED_ADJ) %>%
      distinct() %>%
      summarise(
         ENCOUNTER_TYPE = paste(ENCOUNTER_TYPE, collapse=', '),
         FACILITY = paste(FACILITY, collapse=', '),
         ADMIT_SOURCE = paste(ADMIT_SOURCE, collapse=', '),
         ADMIT_TYPE = paste(ADMIT_TYPE, collapse=', '),
         .by = c(PERSON_ID, ADMIT_DATE, DISCHARGE_DATE, ADMIT_DAY, DISCHARGE_DAY, COMBINED_ADJ)) %>% # ~134K rows
      mutate(LENGTH_OF_STAY = as.numeric(lubridate::as.duration(DISCHARGE_DATE - ADMIT_DATE)) / 86400)
   
   return(encs)
}
