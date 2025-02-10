##### PREP RAW DATA #####
library(dplyr)
source(file = '~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/CreateNamedAbxAbbreviations.R')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_clean_ASTs.Rdata')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_CLEANED_2015_2024_AbxAdmin.Rdata')

# get blood cultures
astDF_non_blood <- astDF %>% filter(!BLOOD) %>% select(PERSON_ID, RESULT_DAY) %>% distinct()
astDF <- astDF %>%
   filter(BLOOD) %>% 
   filter(lubridate::year(ORDER_DAY) %in% 2015:2024) %>%
   filter(!(grepl('Staph', BUG) & !grepl('Staphylococcus (aureus|lugdunensis)', BUG))) %>%
   filter(!(grepl('Streptococc', BUG) & !grepl('Streptococcus (agalactiae|pyogenes|pneumoniae)', BUG))) %>%
   filter(!grepl('Cryptococcus|Aspergillus|Candida', BUG)) %>%
   filter(BUG != 'Did not match')
astDF_non_blood <- astDF_non_blood %>% filter(PERSON_ID %in% unique(astDF$PERSON_ID))
common_bugs <- sort(table(astDF$BUG), decreasing = TRUE)
common_bugs <- names(common_bugs)[common_bugs >= 200L]
common_bugs <- c(
   "Staphylococcus aureus", "Staphylococcus lugdunensis", 
   
   "Enterococcus faecalis", "Enterococcus faecium",
   
   "Streptococcus agalactiae", "Streptococcus pneumoniae", "Streptococcus pyogenes",
   
   "Escherichia coli", 
   "Enterobacter aerogenes", "Enterobacter cloacae",
   "Citrobacter freundii", 
   "Klebsiella aerogenes", "Klebsiella oxytoca", "Klebsiella pneumoniae", "Klebsiella variicola", 
   "Morganella morganii",
   "Proteus mirabilis",
   "Serratia marcescens", 
   
   "Pseudomonas aeruginosa", "Stenotrophomonas maltophilia", "Acinetobacter baumannii"
)
astDF <- astDF %>% 
   filter(BUG %in% common_bugs)

astDF <- astDF %>%
   mutate(RESULT_DELAY = as.numeric(lubridate::as.duration(RESULT_DATE - ORDER_DATE)) / 86400)
astDF <- astDF %>% filter(RESULT_DELAY <= 8)
# first, for staph aureus missing oxacillin AST, remove
w <- which(astDF$BUG == 'Staphylococcus aureus' & is.na(astDF$OXACILLIN)) # 20
astDF <- astDF[-w,]
rm(w)

# IMPUTE MISSING ASTs
source('~/Desktop/EHR-mining/UsefulDataForCleaning/ASTimputation/ASTimputation.R')
astDF <- imputeASTs(astDF)
rm(imputeASTs, rs)

# FLAG as concordant or discordants
# TODO:

astDF$BUG <- gsub('^([A-Z])[a-z]+ ([a-z]+)$', '\\1_\\2', astDF$BUG)

# take only index blood culture by 42 days
ast <- astDF %>%
   select(-BLOOD, -RESPIRATORY) %>%
   group_by(PERSON_ID) %>%
   filter(is.na(lag(ORDER_DAY)) | as.integer(ORDER_DAY - lag(ORDER_DAY)) >= 42L | as.integer(ORDER_DAY - lag(ORDER_DAY)) == 0L) %>%
   ungroup()
names(ast)[grep('^(RESULT|ORDER)_DA(TE|Y)', names(ast))] <- paste0('INDEX_', names(ast)[grep('^(RESULT|ORDER)_DA(TE|Y)', names(ast))])
# now, each visit day is defined by a person_id and an INDEX_ORDER_DATE
# there are gonna be other positive blood cultures during those visits,
# let's add them and indicate the ORDER_DATE
secondary_asts <- ast %>%
   select(PERSON_ID, INDEX_ORDER_DATE, INDEX_RESULT_DATE, INDEX_ORDER_DAY) %>%
   distinct() %>%
   mutate(JOIN_END = INDEX_ORDER_DAY + 21) %>%
   left_join(
      x = .,
      y = astDF %>% 
         select(PERSON_ID, ORDER_DATE, RESULT_DATE, ORDER_DAY, BUG) %>% 
         distinct() %>%
         rename(OTHER_BUG = BUG),
      by = join_by(
         PERSON_ID,
         INDEX_ORDER_DAY < ORDER_DAY,
         JOIN_END >= ORDER_DAY
      )
   ) %>%
   mutate(
      ORDER_DATE = as.numeric(lubridate::as.duration(ORDER_DATE - INDEX_ORDER_DATE)) / 86400,
      RESULT_DATE = as.numeric(lubridate::as.duration(RESULT_DATE - INDEX_ORDER_DATE)) / 86400,
   ) %>%
   filter(!is.na(ORDER_DATE)) %>%
   summarise(
      SECONDARY_BUGS = list(OTHER_BUG),
      SECONDARY_ORDER_TIME = list(ORDER_DATE),
      SECONDARY_RESULT_TIME = list(RESULT_DATE),
      .by = c(PERSON_ID, INDEX_ORDER_DAY)
   )

ast <- ast %>%
   left_join(
      x = .,
      y = secondary_asts,
      by = join_by(PERSON_ID, INDEX_ORDER_DAY)
   ) %>%
   relocate(
      SECONDARY_BUGS, SECONDARY_ORDER_TIME, SECONDARY_RESULT_TIME,
      .after = BUG
   )

names(ast)[grep('^INDEX_', names(ast))] <- gsub('INDEX_', '', grep('^INDEX_', names(ast), value=TRUE))



# handle ABX administration data
# remove antifungals
abxDF <- abxDF %>%
   filter(PERSON_ID %in% unique(ast$PERSON_ID)) %>%
   filter(!ABX %in% c('CASPOFUNGIN', 'FLUCONAZOLE', 'METRONIDAZOLE', 'VORICONAZOLE', 'POSACONAZOLE', 'CLOTRIMAZOLE', 
                      'MICAFUNGIN', 'KETOCONAZOLE', 'AMPHOTERICIN', 'TERBINAFINE', 'ITRACONAZOLE', 'FLUCYTOSINE')) %>%
   select(-START_DAY, -END_DATE) %>%
   distinct()
w <- which(abxDF$ABX == 'PENICILLIN G')
abxDF$ABX[w] <- 'BENZYLPENICILLIN'
rm(w)
abx <- abxDF

# save so don't have to recompute
save(ast, file = '~/Desktop/EHR/EHR work/RdataFiles/R01/blood_asts.Rdata')
save(abx, file = '~/Desktop/EHR/EHR work/RdataFiles/R01/blood_abxs.Rdata')

# load previously computed data
library(dplyr)
load(file = '~/Desktop/EHR/EHR work/RdataFiles/R01/blood_asts.Rdata')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/R01/blood_abxs.Rdata')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_DEMO.Rdata')

# join AST + ABX to get all proximal administrations
df <- ast %>%
   select(PERSON_ID, ORDER_DAY, ORDER_DATE, RESULT_DATE, BUG, contains('SECONDARY')) %>%
   distinct() %>%
   left_join(
      x = .,
      y = dth,
      by = 'PERSON_ID'
   ) %>%
   mutate(
      AGE = as.integer(ORDER_DAY - DOB) / 365,
      time = as.integer(DEATH_DATE - ORDER_DAY)
   ) %>%
   mutate( # empiric prescription time frame -48 hours to +12 hours
      JOIN_START = ORDER_DATE - 86400 * 2,
      JOIN_END = ORDER_DATE + 86400 * 21
   ) %>%
   left_join(
      x = .,
      y = abx,
      by = join_by(
         PERSON_ID,
         between(y$START_DATE, x$JOIN_START, x$JOIN_END)
      )
   ) %>%
   select(!c(DOB, PATIENT_STATUS, DEATH_DATE, JOIN_START, JOIN_END)) %>%
   distinct() %>% 
   relocate(ABX, .after=ORDER_DATE) %>%
   filter(!is.na(ABX)) %>%
   mutate(
      TIME_DIFF_HOURS = as.numeric(lubridate::as.duration(START_DATE - ORDER_DATE)) / 86400,
      RESULT_TIME = as.numeric(lubridate::as.duration(RESULT_DATE - ORDER_DATE)) / 86400
   ) %>%
   mutate(
      TIME_DIFF_12hours = round(TIME_DIFF_HOURS * 2) / 2,
      TIME_DIFF_DAYS = floor(TIME_DIFF_HOURS),
      RESULT_TIME_12hours = round(RESULT_TIME * 2) / 2
   ) %>%
   select(-RESULT_DATE, -START_DATE)


# save for future uses
save(df, file = '~/Desktop/EHR/EHR work/RdataFiles/R01/blood_asts_abxs.Rdata')
##### END #####


library(dplyr)
source(file = '~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/CreateNamedAbxAbbreviations.R')
source(file = '~/Desktop/WrightLab/UsefulRfuncs/newBarplotFxn.R')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/R01/blood_asts_abxs.Rdata')

bug_groups <- list(
   Staphylococci = c("Staphylococcus aureus", "Staphylococcus lugdunensis"),
   
   Enterococci = c("Enterococcus faecalis", "Enterococcus faecium"),
   
   Streptococci = c("Streptococcus agalactiae", "Streptococcus pneumoniae", "Streptococcus pyogenes"),
   
   Enterobacterales = c("Escherichia coli", 
                        "Enterobacter aerogenes", "Enterobacter cloacae",
                        "Citrobacter freundii", 
                        "Klebsiella aerogenes", "Klebsiella oxytoca", "Klebsiella pneumoniae", "Klebsiella variicola", 
                        "Morganella morganii",
                        "Proteus mirabilis",
                        "Serratia marcescens"),
   
   NonFermGN = c("Pseudomonas aeruginosa", "Stenotrophomonas maltophilia", "Acinetobacter baumannii")
)


# when do switches happen?
# plot the total number of unique ABXs that a patient has received up to that point
df_og <- df

df <- df_og %>%
   select(PERSON_ID, ORDER_DAY, RESULT_TIME_12hours, BUG, ABX, TIME_DIFF_12hours) %>%
   distinct() %>% 
   summarise(
      x = min(TIME_DIFF_12hours), 
      .by=c(PERSON_ID, ORDER_DAY, RESULT_TIME_12hours, BUG, ABX)
   ) %>% 
   group_by(PERSON_ID, ORDER_DAY, RESULT_TIME_12hours, BUG) %>%
   count(x) %>%
   mutate(n = cumsum(n)) %>%
   ungroup() %>%
   tidyr::pivot_wider(
      names_from = x,
      values_from = n,
      names_prefix = 'd'
   )

names_order <- c('PERSON_ID', 'ORDER_DAY', 'RESULT_TIME_12hours', 'BUG',
                 gsub('\\.', '_', gsub('-', 'm', paste0('d', as.character(seq(-2, 21, by=0.5))))))
names(df) <- gsub('\\.', '_', gsub('-', 'm', names(df)))
df <- df[names_order]
rm(names_order)



# do every row in parallel, column by column
w <- which(is.na(df$dm2))
df$dm2[w] <- 0L
for (i in (which(names(df) == 'dm2') + 1):length(df)) {
   w <- which(is.na(df[[i]]))
   if (length(w) == 0L)
      next
   # for this column, set the NAs equal to the previous column
   df[[i]][w] <- df[[i-1]][w]
}
any(sapply(df, function(x) any(is.na(x)))) # it worked.
rm(w, i)




# PLOT SWTICHES
plotAbxSwitchPlot <- function(df, bug_group='all', pdf_flag=FALSE) {
   if (bug_group != 'all') {
      df <- df %>% filter(BUG %in% bug_groups[[bug_group]])
   }
   cat(bug_group, nrow(df), '\n')
   
   # get frequency of each transition at each day
   tranDF <- tibble()
   for (i in which(names(df) == 'd0'):(length(df)-1)) {
      # cat(names(df)[i], '')
      x <- df %>% select(!!names(df)[i], !!names(df)[i+1])
      names(x) <- c('first', 'second')
      x <- x %>% count(first, second)
      x$from_day <- names(df)[i]
      x$to_day <- names(df)[i+1]
      tranDF <- rbind(tranDF, x)
   }
   day_names <- setNames(seq(0, 21, 0.5),
                         names(df)[which(names(df) == 'd0'):length(df)])
   
   tranDF <- tranDF %>% 
      mutate(
         from_day = unname(day_names[from_day]),
         to_day = unname(day_names[to_day])
      )
   
   # first and second tell the y values
   # for from_day and to_day (x values), respectively
   # n is the frequency of that change
   
   # this isn't super helpful...other than to say plenty of patient continue to switch
   # at each day, what is the mean/median number of antibiotics given up to that point?
   means <- tranDF %>%
      summarise(
         mean = weighted.mean(x=first, w=n), 
         median = median(rep(first, times=n)),
         q90 = quantile(rep(first, times=n), 0.9),
         .by=from_day
      )
   
   
   # I want to plot a point at each RESULT_TIME (x-value)
   # at the height of the patients' corresponding # unique (y-value)
   # to get that y-value, I need find which column it is contained in
   # use match() to get columns
   res_times_list <- lapply(unique(df$RESULT_TIME_12hours),
                            FUN = function(res_time_i) {
                               w_col <- match(paste0('d', gsub('\\.', '_', res_time_i)), names(df))
                               w_rows <- which(df$RESULT_TIME_12hours == res_time_i)
                               df[[w_col]][w_rows]
                            })
   names(res_times_list) <- unique(df$RESULT_TIME_12hours)
   
   res_times_list <- res_times_list[order(as.numeric(names(res_times_list)))]
   
   tables <- sapply(res_times_list, table)
   
   if (pdf_flag)
      pdf(file = paste0('~/Desktop/EHR-mining/Scripts/AnalysisScripts/R01/switch_plots/totalAbxByDay_', bug_group, '.pdf'),
          width = 8)
   plot(NA, yaxt='n',
        xlim=c(0, 21), ylim=c(0, 12), 
        xlab='Days', ylab='# unique abx given', 
        main=paste0(bug_group, ' (n = ', prettyNum(nrow(df), big.mark=','), ')'))
   axis(side=2, las=1)
   segments(x0 = tranDF$from_day, 
            x1 = tranDF$to_day,
            y0 = tranDF$first,
            y1 = tranDF$second,
            col = '#00000055',
            lwd = (log(tranDF$n)/2)+0.1)
   
   for (i in seq_along(tables)) {
      points(x = rep(as.numeric(names(tables[i])), length(tables[[i]])),
             y = as.integer(names(tables[[i]])),
             cex = (log(tables[[i]])/4)+0.1,
             pch = 16)
   }
   
   lines(x=means$from_day, y=means$mean, col='red', lwd=6)
   lines(x=means$from_day, y=means$median, col='blue', lwd=6)
   lines(x=means$from_day, y=means$q90, col='purple', lwd=6)
   legend('topleft', inset=c(0.025, 0), bty='n',
          legend = c('AST results available', 'real patients', 'mean', 'median', '90th percentile'), 
          col = c('black', 'black', 'red', 'blue', 'purple'), 
          lwd = rep(c(NA, 2, 6), times=c(1, 1, 3)), 
          pch = rep(c(16, NA), times=c(1, 4)),
          pt.cex = 2)
   if (pdf_flag)
      dev.off()
}



# plot all bugs
{
   plotAbxSwitchPlot(df, pdf_flag=TRUE)
   plotAbxSwitchPlot(df, bug_group='Staphylococci', pdf_flag=TRUE) # slightly earlier switches, higher mean number of abx by end
   plotAbxSwitchPlot(df, bug_group='Streptococci', pdf_flag=TRUE) # almost no early treatment, much earlier switches than average
   plotAbxSwitchPlot(df, bug_group='Enterococci', pdf_flag=TRUE) # earlier switches, higher starting
   plotAbxSwitchPlot(df, bug_group='Enterobacterales', pdf_flag=TRUE) # later switches
   plotAbxSwitchPlot(df, bug_group='NonFermGN', pdf_flag=TRUE) # earlier switches
   
   
   pdf(file = '~/Desktop/EHR-mining/Scripts/AnalysisScripts/R01/switch_plots/totalAbxByDay_all_separate.pdf', width=12, height=8)
   par(mfrow=c(2,3), mar=c(4,4,2.5,1.5), tck=-0.015, mgp=c(2, 0.5, 0))
   for (bg in names(bug_groups)) {
      plotAbxSwitchPlot(df, bug_group=bg)
   }
   dev.off()
   rm(bg)
}





pdf(file = '~/Desktop/EHR-mining/Scripts/AnalysisScripts/R01/switch_plots/num_unique_abx_per_visit.pdf')
df_og %>%
   select(PERSON_ID, ORDER_DAY, ABX) %>%
   distinct() %>%
   count(PERSON_ID, ORDER_DAY) %>%
   pull(n) %>%
   plotBarplot(main = 'Number of unique antibiotics administered during stay', 
               xlab = '# unique abx')
dev.off()






# above captures when a new antibiotic is added
# but not when an antibiotic that hasn't been seen in a few days gets re-introduced


uids <- df_og %>% 
   filter(lubridate::year(ORDER_DAY) >= 2023L) %>%
   select(PERSON_ID, ORDER_DAY) %>% 
   distinct() %>%
   slice_sample(prop = 1L, replace=FALSE)

b <- 30L
N <- 20L * b
chunks <- mapply(seq(1, N, b), seq(b, N, b), FUN = ':')

for (n in seq_len(ncol(chunks))) {
   cat(n, '')
   # setup pdf and plotting area for this chunk of 20 patient timelines
   pdf(file = paste0('~/Desktop/EHR-mining/Scripts/AnalysisScripts/R01/switch_plots/patient_abx_sequence_plots/chunk_', n, '.pdf'), 
       width=12)
   par(mfrow=c(5,6), mar=c(0.5,2,1,0.5), mgp=c(1,0.2,0), tck=-0.015, cex.main=1, cex.lab=1, cex.axis=1)
   
   for (i in chunks[,n]) {
      pid <- uids$PERSON_ID[i]
      ord <- uids$ORDER_DAY[i]
      x <- df_og %>%
         filter(PERSON_ID == pid, ORDER_DAY == ord) %>%
         select(PERSON_ID, ORDER_DAY, time, RESULT_TIME_12hours, BUG, ABX, contains('SECONDARY'), TIME_DIFF_12hours) %>%
         distinct() %>%
         mutate(X = 1L) %>%
         tidyr::pivot_wider(
            values_from = X,
            names_from = ABX,
            values_fill = 0
         )
      yr <- lubridate::year(x$ORDER_DAY[1])
      other_bugs <- x$SECONDARY_BUGS[[1]]
      other_order_times <- x$SECONDARY_ORDER_TIME[[1]]
      other_result_times <- x$SECONDARY_RESULT_TIME[[1]]
      
      w <- which(names(x) == 'TIME_DIFF_12hours')
      plot(NA, xlim=c(0, 14), ylim=c(0, length(x)-w+1), main=paste(c(unique(x$BUG), other_bugs), collapse=','), xlab='', yaxt='n', ylab='', xaxt='n')
      legend('topright', title=yr, legend=NA, bty='n', title.font=2)
      axis(side=1, at=seq(0, 21, 3), tck=0.015, labels=rep('', length(seq(0, 21, 3))))
      axis(side=2, at=1:(length(x)-w), las=1, labels=unname(abbr[names(x)[(w+1):length(x)]]))
      abline(v = c(0, other_order_times), lty=2, col='blue')
      abline(v = c(unique(x$RESULT_TIME_12hours), other_result_times), lty=3)
      abline(v = unique(x$time), lty=2, lwd=2, col='red')
      for (i in 1:(length(x)-w)) {
         col <- i + w
         xvals <- x$TIME_DIFF_12hours[x[[col]] == 1L]
         yvals <- rep(i, sum(x[[col]] == 1L))
         segments(x0=xvals-0.25, x1=xvals+0.25, y0=yvals, lwd=5)
      }
   }
   dev.off()
}









