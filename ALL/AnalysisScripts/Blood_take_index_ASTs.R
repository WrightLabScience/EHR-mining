library(dplyr)
library(tidyr)

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
   filter(BLOOD) %>% # 151,258
   #filter(!grepl('Aspergillus|Candida|Cryptococcus', BUG)) %>%
   select(-BLOOD, -ORDER_PROC_ID) %>%
   mutate(across(c(ORDER_DATE, RESULT_DATE), ~ as.Date(substr(.,1,10)))) %>%
   distinct() %>% # 144,824
   arrange(PERSON_ID, ORDER_DATE, RESULT_DATE, BUG)
rm(astDF)
gc()


# how often is the same bug isolated on same day with same eventual AST profile, but with different result_date?
# fairly often, combine them and take the minimum result_date
empDF <- empDF %>%
   group_by_all() %>%
   ungroup(RESULT_DATE) %>%
   slice_min(RESULT_DATE) %>% # 144,824 --> 141,726
   ungroup()


num_blood_cults_month <- table(substr(empDF$ORDER_DATE, 1, 7))
getYear <- function(x) as.integer(substr(x,1,4)) + (as.integer(substr(x,6,7)) - 1) / 12
names(num_blood_cults_month) <- getYear(names(num_blood_cults_month))
plot(num_blood_cults_month, type='l', xaxt='n')
axis(side=1, at=seq(2004, 2024, 1))
abline(v = c(2015, 2023))
abline(h = c(700, 1000), lty=2)
rm(num_blood_cults_month, getYear)


# only take 2015 - 2023 - 97,202
empDF <- empDF %>% filter(substr(ORDER_DATE, 1, 4) %in% as.character(2015:2023))
length(unique(empDF$PERSON_ID)) # 55,726



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



# I want to take only the first culture in >= 30 days
#  how many of these have >1 isolate? (exclude them for now)
#  how to do this?

empDF$WITHIN_PRV <- FALSE
empDF$WITHIN_PRV[empDF$ORDER_DATE <= lag(empDF$RESULT_DATE) & empDF$ORDER_DATE > lag(empDF$ORDER_DATE)] <- TRUE
empDF$WITHIN_PRV[empDF$PERSON_ID != lag(empDF$PERSON_ID)] <- FALSE

empDF$DAYS_SINCE_PRV <- as.integer(empDF$ORDER_DATE - lag(empDF$RESULT_DATE))
empDF$DAYS_SINCE_PRV[empDF$PERSON_ID != lag(empDF$PERSON_ID)] <- NA
#empDF$DAYS_SINCE_PRV[empDF$WITHIN_PRV] <- NA

empDF %>% count(WITHIN_PRV, DAYS_SINCE_PRV >= 60)
# 60
#   WITHIN_PRV `DAYS_SINCE_PRV >= 30`      n
# 1 FALSE      FALSE                  22,097  -  subsequent cultures that occurred during the same "visit"
# 2 FALSE      TRUE                   10,190  -  subsequent cultures that represent the next "visit"
# 3 FALSE      NA                     55,726  -  each person's first culture record
# 4 TRUE       NA                      9,189  -  cultures that occurred in the window of the previous culture

# 30
#   WITHIN_PRV `DAYS_SINCE_PRV >= 30`      n
# 1 FALSE      FALSE                   4,740  -  subsequent cultures that occurred during the same "visit"
# 2 FALSE      TRUE                   12,211  -  subsequent cultures that represent the next "visit"
# 3 FALSE      NA                     55,213  -  each person's first culture record
# 4 TRUE       NA                     23,208  -  cultures that occurred in the window of the previous culture

empDF$MULT_ISOLATES <- FALSE
empDF$MULT_ISOLATES[empDF$ORDER_DATE == lead(empDF$ORDER_DATE)] <- TRUE
empDF$MULT_ISOLATES[which(empDF$MULT_ISOLATES) + 1] <- TRUE
empDF$DAYS_SINCE_PRV[empDF$MULT_ISOLATES] <- NA

empDF <- empDF %>%
   relocate(WITHIN_PRV, DAYS_SINCE_PRV, MULT_ISOLATES, .before=BUG)

empDF %>% count(WITHIN_PRV, MULT_ISOLATES, DAYS_SINCE_PRV >= 60)


# 97,202 --> 82,042
empDF <- empDF %>% filter(!WITHIN_PRV, is.na(DAYS_SINCE_PRV) | DAYS_SINCE_PRV > 60)



# plot num visits / time bw visits
if (FALSE) {
   x <- empDF %>% group_by(PERSON_ID) %>% tally() %>% count(n)
   x <- x$nn
   barplot(x, main='Number of "index" blood cultures per patient', names.arg = 1:length(x), xlab='"Index" cultures')
   
   x <- empDF$DAYS_SINCE_PRV
   x <- x[!is.na(x)]
   hist(x, breaks=diff(range(x)), xlab='Days', main='Days between "visits"')
   
   x <- as.integer(empDF$RESULT_DATE - empDF$ORDER_DATE)
   hist(x, breaks=diff(range(x)), xlim=c(0,12))
   
   rm(x)
}


empDF <- empDF %>% select(-WITHIN_PRV, -DAYS_SINCE_PRV, -MULT_ISOLATES)

save(empDF, file = '~/Desktop/EHR/EHR work/RdataFiles/ASTs_blood_2015_2023.Rdata')





