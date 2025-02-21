featurizeEncounters <- function(df, encsDF, time_frame=c(-4L, 4L)) {
   cat('Preprocessing data.\n')
   df <- df %>%
      select(PERSON_ID, ORDER_DAY, ORDER_DATE) %>%
      distinct() %>%
      slice_min(ORDER_DATE, by=c(PERSON_ID, ORDER_DAY))
   
   cat('Joining encounters.\n')
   df <- df %>%
      mutate(
         BEFORE_ORDER_DAY = ORDER_DAY + time_frame[1],
         AFTER_ORDER_DAY = ORDER_DAY + time_frame[2]
      ) %>%
      left_join(
         x = .,
         y = encsDF,
         by = join_by(
            PERSON_ID,
            overlaps(x_lower = ORDER_DAY,
                   x_upper = AFTER_ORDER_DAY,
                   y_lower = ADMIT_DAY,
                   y_upper = DISCHARGE_DAY)
         )
      )
   
   cat('Computing features.\n')
   df <- df %>%
      mutate(
         TIME_BW_ADMIT_ORDER = as.numeric(lubridate::as.duration(ORDER_DATE - ADMIT_DATE)) / 86400,
         DAYS_BW_ADMIT_ORDER = as.integer(ORDER_DAY - ADMIT_DAY)
      ) %>%
      group_by(PERSON_ID, ORDER_DAY) %>%
      mutate(
         HOSP_ACQ = as.integer(any(TIME_BW_ADMIT_ORDER >= 2)),
         NURSING_HOME_HOSPICE = as.integer(any(grepl('Nursing|hospice', ADMIT_SOURCE))),
         EMERGENCY_DEPT = as.integer(any(grepl('ED|Emergency|Trauma|Urgent', ADMIT_TYPE))),
         MATERNAL = as.integer(any(grepl('Baby|Maternal', ADMIT_SOURCE) | grepl('Newborn', ADMIT_TYPE))),
         TRANSFER = as.integer(any(grepl('Trans', ADMIT_SOURCE))),
         FACILITIES = list(unique(FACILITY))
      ) %>%
      ungroup() %>%
      select(PERSON_ID, ORDER_DAY, HOSP_ACQ, NURSING_HOME_HOSPICE, EMERGENCY_DEPT, MATERNAL, TRANSFER, FACILITIES) %>%
      distinct()
   
   cat('Final processing.\n')
   # handle missing
   w <- which(is.na(df$HOSP_ACQ))
   df$MISSING_ENCOUNTER_DATA <- FALSE
   df$MISSING_ENCOUNTER_DATA[w] <- TRUE
   # df$HOSP_ACQ[w] <- 0L
   df$NURSING_HOME_HOSPICE[w] <- NA_integer_
   df$EMERGENCY_DEPT[w] <- NA_integer_
   df$MATERNAL[w] <- NA_integer_
   df$TRANSFER[w] <- NA_integer_
   
   # get the admit and transfer facilities
   df$FACILITY1 <- sapply(df$FACILITIES, '[', 1L)
   df$FACILITY2 <- sapply(df$FACILITIES, '[', 2L)
   w <- which(is.na(df$FACILITY2) & !is.na(df$FACILITY1))
   df$FACILITY2[w] <- df$FACILITY1[w]
   # translate to readable facility names
   source('~/Desktop/EHR-mining/Scripts/CleaningScripts/CleanFacilityNames.R')
   df <- df %>%
      select(-FACILITIES) %>%
      mutate(
         across(.cols = c(FACILITY1, FACILITY2),
                .fns = ~ getFacilityNames(.))
      )
   
   names(df)[!names(df) %in% c('PERSON_ID', 'ORDER_DAY')] <- paste0('ENC_', names(df)[!names(df) %in% c('PERSON_ID', 'ORDER_DAY')])
   
   return(df)
}
