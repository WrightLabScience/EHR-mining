library(dplyr)

load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_clean_ASTs.Rdata')
print(object.size(astDF), units='Mb') # 2,028.7 Mb

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
astDF_og <- astDF
astDF <- astDF %>%
   select(-ORDER_PROC_ID, -RESULT_DAY, -PATH_NAME) %>%
   distinct() %>%
   mutate(ORDER_DAY =  as.Date(substr(ORDER_DATE,1,10))) %>%
   relocate(ORDER_DAY, .after=RESULT_DATE)
gc()




# how often is the same bug isolated on same day with same eventual AST profile, but with different result_date?
# fairly often, combine them and take the minimum result_date
#### groups have same order DAY and result DAY, but not necessarily time
# take the minimum of both (often its the same row, occasionally )
astDFm <- astDF %>% group_by_all() %>% ungroup(RESULT_DATE, ORDER_DATE) %>% filter(n() > 1L)
astDFm <- astDFm %>%
   summarise(ORDER_DATE = min(ORDER_DATE),
             RESULT_DATE = min(RESULT_DATE)) %>%
   ungroup()

astDF <- rbind(astDFm, 
               astDF %>% 
                  group_by_all() %>% 
                  ungroup(RESULT_DATE, ORDER_DATE) %>% 
                  filter(n() == 1L))
rm(astDFm)

astDF <- astDF %>%
   mutate(RESULT_DAY = as.Date(substr(RESULT_DATE,1,10))) %>%
   relocate(ORDER_DATE, RESULT_DATE, RESULT_DAY, .after=ORDER_DAY) %>%
   arrange(PERSON_ID, ORDER_DAY, RESULT_DAY)





# PLOT - number of blood culture per month
if (FALSE) {
   num_blood_cults_month <- table(substr(astDF$ORDER_DATE, 1, 7))
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
astDF <- astDF %>% filter(substr(ORDER_DATE, 1, 4) %in% as.character(2015:2023))
length(unique(astDF$PERSON_ID)) # 55,726 (446,601)


# plotting bugs vs. year
if (FALSE) {
   # PLOT PREVALENCE OF BUGS/GENUS PER YEAR
   # BY GENUS
   bugs_year <- astDF %>%
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
   common_bugs <- names(head(sort(table(astDF$BUG), decreasing = TRUE), n = topN))
   bug_by_year <- astDF %>%
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
   astDF %>%
      group_by(PERSON_ID) %>% 
      tally() %>% 
      count(n > 1L)
   
   x <- astDF %>%
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
astDF <- astDF %>%
   mutate(FIRST_RECORD = PERSON_ID != lag(PERSON_ID)) %>%
   mutate(MULT_ISO = (ORDER_DAY == lag(ORDER_DAY) & PERSON_ID == lag(PERSON_ID)) | (ORDER_DAY == lead(ORDER_DAY) & PERSON_ID == lead(PERSON_ID))) %>%
   relocate(FIRST_RECORD, MULT_ISO, .after=RESULT_DAY)
astDF$FIRST_RECORD[1] <- TRUE
astDF$MULT_ISO[1] <- (astDF$ORDER_DAY[1] == lead(astDF$ORDER_DAY)[1] & astDF$PERSON_ID[1] == lead(astDF$PERSON_ID)[1])
astDF$MULT_ISO[nrow(astDF)] <- (astDF$ORDER_DAY[nrow(astDF)] == lag(astDF$ORDER_DAY)[nrow(astDF)]  & astDF$PERSON_ID[nrow(astDF)] == lag(astDF$PERSON_ID)[nrow(astDF)])
astDF %>% count(MULT_ISO)       # nearly 1/3
astDF %>% count(FIRST_RECORD)   # > 1/3
astDF %>% select(MULT_ISO, FIRST_RECORD) %>% table()

astDF$INDEX_RECORD <- astDF$FIRST_RECORD | (astDF$PERSON_ID == lag(astDF$PERSON_ID) & astDF$ORDER_DAY - lag(astDF$ORDER_DAY) >= 180)
astDF %>% select(INDEX_RECORD, FIRST_RECORD) %>% table()
astDF <- astDF %>% 
   relocate(INDEX_RECORD, .after=MULT_ISO) %>% 
   select(-FIRST_RECORD)


# WHEN INDEX RECORD HAD MULTIPLE ISOLATES ON THE SAME ORDER DAY
# MARK THOSE OTHER ISOLATES AS INDEX AS WELL
w <- which(astDF$INDEX_RECORD & astDF$MULT_ISO)
length(w) # 89,742
shift <- 1
while(length(w) > 0) {
   cat(shift, length(w), '\n')
   astDF$INDEX_RECORD[w + shift] <- TRUE
   shift <- shift + 1
   w <- w[which(astDF$ORDER_DAY[w + shift] == astDF$ORDER_DAY[w] & astDF$MULT_ISO[w + shift])]
}
rm(w, shift)

## DETERMINE IF INDEX CULTURES ARE "CONTAMINATED" BY SUBSEQUENT CULTURES
# THAT IS, WERE SUBSEQUENT CULTURES ORDERED OR HAD RESULTS RETURNED IN THE EMPIRIC WINDOW?
# INDEX CULTURES
astDF_idx <- astDF %>%
   filter(INDEX_RECORD) %>%
   group_by(PERSON_ID, ORDER_DAY) %>%
   slice_min(RESULT_DATE, with_ties = FALSE) %>%
   ungroup() %>%
   select(PERSON_ID, ORDER_DAY, RESULT_DAY) %>% 
   distinct()

# NON-INDEX CULTURES - POTENTIAL "CONTAMINATORS"
astDF_all <- astDF %>%
   # filter(!INDEX_RECORD) %>%
   select(PERSON_ID, ORDER_DAY, RESULT_DAY, RESULT_DATE) %>%
   group_by(PERSON_ID, ORDER_DAY) %>%
   slice_min(RESULT_DATE, with_ties = FALSE) %>%
   ungroup() %>%
   select(-RESULT_DATE) %>%
   distinct()



astDF_idx_joined <- left_join(x = astDF_idx,
                              y = astDF_all,
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

length(unique(astDF_idx_joined$PERSON_ID)) # 180,148  (292,853 index infections)
length(unique(astDF_idx$PERSON_ID))        # 446,636  (631,461 index infections)

astDF     # all - with all variables
astDF_all # all - just person and dates
astDF_idx # only index cultures
astDF_idx_joined # only index cultures in which THERE IS a next culture, along with the date of that next culture

astDF_index <- astDF %>%
   filter(INDEX_RECORD) %>% # 762,226
   left_join(x = .,
             y = astDF_idx_joined,
             by = join_by(PERSON_ID, ORDER_DAY, RESULT_DAY)) %>% # 762,226
   relocate(NEXT_ORDER_DAY, NEXT_RESULT_DAY, .before=BUG)


### COULD SAVE ENTIRE INDEX DATASET RIGHT HERE
# save(astDF, file = '')

## ONLY KEEP BLOODSTREAM INFECTIONS
astDF %>% count(BLOOD)
astDF_index %>% count(BLOOD)
astDF <- astDF_index %>% filter(BLOOD)
rm(astDF_all, astDF_idx, astDF_idx_joined, astDF_index)

astDF %>%
   #filter(!is.na(NEXT_ORDER_DAY)) %>%      # 23,671 have a next order
   #filter(NEXT_ORDER_DAY < RESULT_DAY) %>% # 8,500 - next order day occurs before index result day
   filter(NEXT_RESULT_DAY < RESULT_DAY)    # 1,948 get results for the second culture before the index

x <- astDF$NEXT_ORDER_DAY - astDF$ORDER_DAY
x <- x[!is.na(x)]
median(x) # 12 days
length(x[x <= 7]) / length(x) # 47%
length(x[x == 1]) / length(x) # 24%
barplot(table(x), xlim=c(0,500), ylim=c(0,100))


x <- astDF %>%
   filter(!is.na(NEXT_ORDER_DAY)) %>%      # 23,674 (out of 55,989) have a next order
   filter(NEXT_ORDER_DAY < RESULT_DAY) %>%
   mutate(x = NEXT_ORDER_DAY - ORDER_DAY) %>%
   select(x) %>%
   unlist(x)
table(x == 1) # 66.3% occur on day after
table(x == 2) # 22.3% occur on day after
barplot(table(x), xlim=c(0,30))

save(astDF, file = '~/Desktop/EHR/EHR work/RdataFiles/ASTs_blood_2015_2023.Rdata')









