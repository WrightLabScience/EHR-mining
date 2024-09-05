library(dplyr)

load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_clean_ASTs.Rdata')
print(object.size(astDF), units='Mb') # 1,975.3 Mb

# consider adding some missing AST antibiotics as columns to impute later!
# but first get only the drugs that are actually prescribed in the context of blood cultures
# will be like 20 drugs! (not 120+)
# 
# ast_abx <- names(astDF)[!names(astDF) %in% c('PERSON_ID', 'ORDER_PROC_ID', 'ORDER_DATE', 'RESULT_DATE', 'BLOOD', 'BUG')]
# 
# load(file = '~/Desktop/EHR/EHR work/ALL/Rdata/abxDF_all.Rdata')
# med_abx <- unique(abxDF$ABX)
# 
# length(ast_abx) # 96
# length(med_abx) # 105
# length(intersect(ast_abx, med_abx)) # 73
# length(union(ast_abx, med_abx)) # 128
# missing_ast_abx <- setdiff(med_abx, ast_abx)
# length(missing_ast_abx) # 32


# sometimes the same exact same bug is repeated but with slightly different times (same date) and therefore different order_proc_ids
empDF <- astDF %>%
   # filter(BLOOD) %>% # 151,258
   # filter(!grepl('Aspergillus|Candida|Cryptococcus', BUG)) %>%
   # select(-BLOOD) %>%
   select(-ORDER_PROC_ID) %>%
   mutate(ORDER_DAY =  as.Date(substr(ORDER_DATE,1,10)),
          RESULT_DAY =  as.Date(substr(RESULT_DATE,1,10))) %>%
   relocate(ORDER_DAY, RESULT_DAY, .after=RESULT_DATE) %>%
   distinct() %>% # 150,716
   arrange(PERSON_ID, ORDER_DATE, RESULT_DATE, BUG)
rm(astDF)
gc()


# how often is the same bug isolated on same day with same eventual AST profile, but with different result_date?
# fairly often, combine them and take the minimum result_date
# WHEN ALL IS EQUAL, TAKE MIN ORDER_DATE
# start <- Sys.time()
# empDF <- empDF %>%
#    group_by_all() %>% ungroup(ORDER_DATE) %>%
#    slice_min(ORDER_DATE) %>% # 150,716 --> 148,846
#    ungroup()
# print(Sys.time() - start) # ~30 seconds

#### groups have same order DAY and result DAY, but not necessarily time
# take the minimum of both (often its the same row, occasionally )
empDF <- empDF %>%
   group_by_all() %>%
   ungroup(RESULT_DATE, ORDER_DATE) %>%
   summarise(ORDER_DATE = min(ORDER_DATE),
             RESULT_DATE = min(RESULT_DATE)) %>% # 150,716 --> 144,824
   relocate(ORDER_DATE, RESULT_DATE, .before=ORDER_DAY) %>%
   ungroup()
print(nrow(empDF)) # 1,837,358 --> 1,808,910

empDF <- empDF %>%
   select(-ORDER_DAY, -RESULT_DAY) %>%
   group_by_all() %>% ungroup(RESULT_DATE) %>% # 144,824 --> 143,972
   slice_min(RESULT_DATE) %>%
   ungroup()
print(nrow(empDF)) # 1,808,910 --> 1,804,376




# PLOT - number of blood culture per month
if (FALSE) {
   num_blood_cults_month <- table(substr(empDF$ORDER_DATE, 1, 7))
   getYear <- function(x) as.integer(substr(x,1,4)) + (as.integer(substr(x,6,7)) - 1) / 12
   names(num_blood_cults_month) <- getYear(names(num_blood_cults_month))
   plot(num_blood_cults_month, type='l', xaxt='n')
   axis(side=1, at=seq(2004, 2024, 1))
   abline(v = c(2015, 2023))
   abline(h = c(700, 1000), lty=2)
   rm(num_blood_cults_month, getYear)
}


# only take 2015 - 2023 - 143,972 --> 98,717
# (including non-blood) 1,804,376 --> 1,240,090
empDF <- empDF %>% filter(substr(ORDER_DATE, 1, 4) %in% as.character(2015:2023))
print(nrow(empDF))
length(unique(empDF$PERSON_ID)) # 55,726 (446,636)


# plotting bugs vs. year
if (FALSE) {
   # PLOT PREVALENCE OF BUGS/GENUS PER YEAR
   # BY GENUS
   bugs_year <- empDF %>%
      mutate(Genus = case_when(
         BUG == 'Coagulase Negative Staph' ~ 'Staphylococcus',
         grepl('Group [A-Z] Streptococci', BUG) ~ 'Streptococcus',
         grepl('Streptococc', BUG, ignore.case=TRUE) ~ 'Streptococcus',
         lengths(strsplit(BUG, ' ')) == 2L ~ gsub('^([A-Z][a-z]+) [a-z]+$', '\\1', BUG),
         grepl('.+ complex$', BUG) ~ gsub('^([A-Z][a-z]+) .+', '\\1', BUG)
      )) %>%
      filter(n() > 1000L, .by=Genus) %>%
      mutate(year = substr(ORDER_DATE,1,4)) %>%
      select(Genus, year) %>%
      table()
   colors <- c('#ffe119', '#4363d8', '#f58231', '#dcbeff', '#800000', '#000075', '#333333', '#888888', '#000000', '#911eb4')
   col_vec <- setNames(colors, rownames(bugs_year)); rm(colors)
   {
      pdf(file = '~/Desktop/EHR/EHR work/ALL/plots/GenusByYear.pdf', width = 10, height = 7.5)
      par(mar = c(2.5, 4, 1.5, 12), mgp = c(3, 0.4, 0), tck=-0.01)
      plot(NA, xlim=c(2014.8, 2023.2), ylim=c(100, 5000), log='y', yaxt='n', xaxt='n', xaxs='i',
           ylab = 'Number of isolates', xlab='',
           main = paste0('Prevalence of genus in blood cultures'))
      axis(side = 2, at = c(100, 200, 300, 500, 1000, 2000, 5000), las=1)
      axis(side = 2, at = seq(100, 4700, 200), labels=rep('', length(seq(100, 4700, 200))), tck=-0.005)
      axis(side = 1, at = 2015:2023)
      title(xlab = 'Year', line=1.5)
      for (b in seq_len(nrow(bugs_year))) {
         bug <- rownames(bugs_year)[b]
         bug_counts <- bugs_year[b,]
         lines(x=2015:2023, bug_counts, lwd=1.6, col=col_vec[bug])
      }
      text(x = 2023.3,
           y = bugs_year[,ncol(bugs_year)] + c('Candida' = 0, 'Enterobacter' = 0, 'Enterococcus' = 0, 'Escherichia' = 0, 
                                               'Klebsiella' = 50, 'Proteus' = 0, 'Pseudomonas' = 0, 'Serratia' = 0,
                                               'Staphylococcus' = 0, 'Streptococcus' = -30),
           labels = paste0(rownames(bugs_year), ' (n = ', format(rowSums(bugs_year), big.mark=','), ')'),
           xpd = NA,
           adj = 0,
           col = col_vec[rownames(bugs_year)])
      dev.off()
   }
   rm(b, bug, bug_counts, bugs_year, col_vec)
   
   # BY BUG
   topN <- 10
   common_bugs <- names(head(sort(table(empDF$BUG), decreasing = TRUE), n = topN))
   bug_by_year <- empDF %>%
      filter(BUG %in% common_bugs) %>%
      mutate(year = as.integer(substr(ORDER_DATE, 1, 4))) %>%
      filter(year %in% 2015:2023) %>%
      select(year, BUG) %>%
      table()
   years <- as.integer(rownames(bug_by_year))
   colnames(bug_by_year)[1] <- 'Coag. Neg. Staph.'
   colnames(bug_by_year)[2:ncol(bug_by_year)] <- gsub('([A-Z])[a-z]+ ([a-z]+)', '\\1. \\2', colnames(bug_by_year)[2:ncol(bug_by_year)])
   colors <- c('#ffe119', '#4363d8', '#f58231', '#dcbeff', '#800000', '#000075', '#a9a9a9', '#000000', '#fabed4', '#911eb4')
   col_vec <- setNames(colors, colnames(bug_by_year))
   {
      pdf(file = '~/Desktop/EHR/EHR work/ALL/plots/BugByYear.pdf', width = 10, height = 7.5)
      par(mar = c(2.5, 4, 1.5, 12), mgp = c(3, 0.4, 0), tck=-0.01)
      plot(NA, xlim=c(2014.8, 2023.2), ylim=c(140, 2500), log='y', yaxt='n', xaxt='n', xaxs='i',
           ylab = 'Number of isolates', xlab='',
           main = paste0('Prevalence of top ', topN, ' most common bugs in blood cultures'))
      axis(side = 2, at = c(200, 300, 400, 500, 1000, 2000), las = 1)
      axis(side = 2, at = seq(200, 2500, 100), labels=rep('', length(seq(200, 2500, 100))), tck=-0.005)
      axis(side = 1, at = seq(years[1],years[length(years)],1))
      title(xlab = 'Year', line=1.5)
      for (b in seq_len(ncol(bug_by_year))) {
         bug <- colnames(bug_by_year)[b]
         bug_counts <- bug_by_year[,b]
         lines(x = years, y = bug_counts, col=col_vec[bug], lwd=1.6)
      }
      text(x = 2023.3, y = bug_by_year[nrow(bug_by_year),]
           + c('CNS' = 22, 'faecalis' = 0, 'faecium' = -10, 'EC' = 0, 'KP' = 0, 'Pmir' = -28, 'PA' = 10, 'SA' = 0, 'Sepi' = 0, 'Shom' = -20),
           labels = paste0(names(col_vec), ' (n = ', format(colSums(bug_by_year), big.mark=','), ')'), 
           cex=1, xpd=NA, col=col_vec, adj = 0)
      dev.off()
   }
   rm(topN, common_bugs, years, b, bug, bug_counts, col_vec, colors, bug_by_year)
}


# determine proximity of subsequent cultures - do they cluster into encounters?
# since we do not have encounters data, this will have to do
if (FALSE) {
   empDF %>%
      group_by(PERSON_ID) %>% 
      tally() %>% 
      count(n > 1L)
   
   x <- empDF %>%
      group_by(PERSON_ID) %>%
      mutate(across(contains('DATE'), ~ as.Date(substr(., 1, 10)))) %>%
      mutate(DAYS_TO_NEXT_CULTURE = as.integer(lead(ORDER_DATE) - ORDER_DATE)) %>%
      ungroup() %>%
      select(DAYS_TO_NEXT_CULTURE) %>%
      unlist()
   x <- x[!is.na(x)]
   mean(x) # 126 days LOL
   median(x) # 2 days LOL
   range(x) # 0 - 3,255 (~9 years)
   table(x == 0)[2] / length(x) # 35% is same day as first
   pdf(file = '~/Desktop/EHR/EHR work/ALL/plots/TimeBWSubsequentCultures.pdf', height=9)
   par(mfrow = c(2, 1), mgp=c(2, 0.5, 0), tck=-0.015)
   hist(x, breaks=diff(range(x)), ylim=c(0, 500), xlim=c(0, 365))
   abline(v = 7)
   barplot(table(x)[1:30])
   dev.off()
}




# REMOVE NON-INDEX CULTURES - FLAG OTHER CULTURES
empDF <- empDF %>%
   mutate(ORDER_DAY =  as.Date(substr(ORDER_DATE,1,10)),
          RESULT_DAY =  as.Date(substr(RESULT_DATE,1,10))) %>%
   relocate(ORDER_DAY, RESULT_DAY, .after=RESULT_DATE) %>%
   mutate(FIRST_RECORD = PERSON_ID != lag(PERSON_ID)) %>%
   mutate(MULT_ISO = (ORDER_DAY == lag(ORDER_DAY) & PERSON_ID == lag(PERSON_ID))
          | (ORDER_DAY == lead(ORDER_DAY) & PERSON_ID == lead(PERSON_ID))) %>%
   relocate(FIRST_RECORD, MULT_ISO, .after=RESULT_DAY)
empDF$FIRST_RECORD[1] <- TRUE
empDF$MULT_ISO[1] <- (empDF$ORDER_DAY[1] == lead(empDF$ORDER_DAY)[1] & empDF$PERSON_ID[1] == lead(empDF$PERSON_ID)[1])
empDF$MULT_ISO[nrow(empDF)] <- (empDF$ORDER_DAY[nrow(empDF)] == lag(empDF$ORDER_DAY)[nrow(empDF)]  & empDF$PERSON_ID[nrow(empDF)] == lag(empDF$PERSON_ID)[nrow(empDF)])
empDF %>% count(MULT_ISO)
empDF %>% count(FIRST_RECORD)
empDF %>% select(MULT_ISO, FIRST_RECORD) %>% table()

empDF$INDEX_RECORD <- empDF$FIRST_RECORD | (empDF$PERSON_ID == lag(empDF$PERSON_ID) & empDF$ORDER_DAY - lag(empDF$ORDER_DAY) >= 180)
empDF %>% select(INDEX_RECORD, FIRST_RECORD) %>% table()
empDF <- empDF %>% 
   relocate(INDEX_RECORD, .after=MULT_ISO) %>% 
   select(-FIRST_RECORD)


# WHEN INDEX RECORD HAD MULTIPLE ISOLATES ON THE SAME ORDER DAY
# MARK THOSE OTHER ISOLATES AS INDEX AS WELL
w <- which(empDF$INDEX_RECORD & empDF$MULT_ISO)
length(w) # 89,742
shift <- 1
while(length(w) > 0) {
   cat(shift, length(w), '\n')
   empDF$INDEX_RECORD[w + shift] <- TRUE
   shift <- shift + 1
   w <- w[which(empDF$ORDER_DAY[w + shift] == empDF$ORDER_DAY[w] & empDF$MULT_ISO[w + shift])]
}
rm(w, shift)

## DETERMINE IF INDEX CULTURES ARE "CONTAMINATED" BY SUBSEQUENT CULTURES
# THAT IS, WERE SUBSEQUENT CULTURES ORDERED OR HAD RESULTS RETURNED IN THE EMPIRIC WINDOW?
# INDEX CULTURES
empDF_idx <- empDF %>%
   filter(INDEX_RECORD) %>%
   group_by(PERSON_ID, ORDER_DAY) %>%
   slice_min(RESULT_DATE, with_ties = FALSE) %>%
   ungroup() %>%
   select(PERSON_ID, ORDER_DAY, RESULT_DAY) %>% 
   distinct()

# NON-INDEX CULTURES - POTENTIAL "CONTAMINATORS"
empDF_all <- empDF %>%
   # filter(!INDEX_RECORD) %>%
   select(PERSON_ID, ORDER_DAY, RESULT_DAY, RESULT_DATE) %>%
   group_by(PERSON_ID, ORDER_DAY) %>%
   slice_min(RESULT_DATE, with_ties = FALSE) %>%
   ungroup() %>%
   select(-RESULT_DATE) %>%
   distinct()



empDF_idx_joined <- left_join(x = empDF_idx,
          y = empDF_all,
          by = join_by(PERSON_ID,
                       x$ORDER_DAY < y$ORDER_DAY)) %>% # 1,319,430
   # filter(RESULT_DAY.x > RESULT_DAY.y) # 3,090 - cool!
   filter(!is.na(ORDER_DAY.y)) %>% # 970,822
   arrange(PERSON_ID, ORDER_DAY.x) %>%
   group_by(PERSON_ID, ORDER_DAY.x) %>%
   slice_min(ORDER_DAY.y) %>%
   ungroup() %>% 
   rename(ORDER_DAY = ORDER_DAY.x, RESULT_DAY = RESULT_DAY.x, 
          NEXT_ORDER_DAY = ORDER_DAY.y, NEXT_RESULT_DAY = RESULT_DAY.y)

length(unique(empDF_idx_joined$PERSON_ID)) # 180,148  (292,853 index infections)
length(unique(empDF_idx$PERSON_ID))        # 446,636  (631,461 index infections)

empDF     # all - with all features
empDF_all # all - just person and dates
empDF_idx # only index cultures
empDF_idx_joined # only index cultures in which THERE IS a next culture, along with the date of that next culture

empDF_index <- empDF %>%
   filter(INDEX_RECORD) %>% # 762,226
   left_join(x = .,
             y = empDF_idx_joined,
             by = join_by(PERSON_ID, ORDER_DAY, RESULT_DAY)) %>% # 762,226
   relocate(NEXT_ORDER_DAY, NEXT_RESULT_DAY, .before=BUG)
   

### COULD SAVE ENTIRE EMPIRIC DATASET RIGHT HERE
# save(empDF, file = '')

## ONLY KEEP BLOODSTREAM INFECTIONS
empDF <- empDF_index %>% filter(BLOOD)

empDF %>%
   filter(!is.na(NEXT_ORDER_DAY)) %>%      # 23,674 (out of 55,989) have a next order
   filter(NEXT_ORDER_DAY < RESULT_DAY) %>% # 8,501 - next order day occurs before index result day
   filter(NEXT_RESULT_DAY < RESULT_DAY)    # 1,946 get results for the second culture before the index

x <- empDF$NEXT_ORDER_DAY - empDF$ORDER_DAY
x <- x[!is.na(x)]
median(x) # 12 days
table(x <= 7) # 47%
table(x == 1) # 23.8%
barplot(table(x), xlim=c(0,30))


x <- empDF %>%
   filter(!is.na(NEXT_ORDER_DAY)) %>%      # 23,674 (out of 55,989) have a next order
   filter(NEXT_ORDER_DAY < RESULT_DAY) %>%
   mutate(x = NEXT_ORDER_DAY - ORDER_DAY) %>%
   select(x) %>%
   unlist(x)
table(x == 1) # 66.3% occur on day after
table(x == 2) # 22.3% occur on day after
barplot(table(x), xlim=c(0,30))

save(empDF, file = '~/Desktop/EHR/EHR work/RdataFiles/ASTs_blood_2015_2023.Rdata')









