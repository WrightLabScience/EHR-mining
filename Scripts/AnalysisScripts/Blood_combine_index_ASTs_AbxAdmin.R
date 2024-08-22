library(dplyr)

load(file = '~/Desktop/EHR/EHR work/RdataFiles/ASTs_blood_2015_2023.Rdata')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_2017_AbxAdmin.Rdata')

# FOR NOW ONLY 2017!!!
empDF <- empDF %>% filter(substr(ORDER_DATE, 1, 4) == '2017')     # 82,042 --> 9,452
abxDF <- abxDF %>% filter(PERSON_ID %in% unique(empDF$PERSON_ID)) # 1,288,305 --> 369,272

length(unique(empDF$PERSON_ID)) # 6,719
length(unique(abxDF$PERSON_ID)) # 6,304
length(intersect(empDF$PERSON_ID, abxDF$PERSON_ID)) # 6,304 - all! (better than medication dispensing data!!)



# add in a column describing the date of the most recent antibiotic prescription
abxDF <- abxDF %>%
   mutate(across(contains('DATE'), ~ as.Date(substr(., 1, 10)))) %>%
   distinct() %>%  # 290,855 --> 154,875
   mutate(PERSON_ID = as.character(PERSON_ID)) %>%
   rename(START_DATE = ADMIN_START_DATE,
          END_DATE = ADMIN_END_DATE) %>%
   arrange(PERSON_ID, START_DATE) %>%
   group_by(PERSON_ID) %>%
   mutate(DATE_OF_MOST_RECENT_ABX = case_when(
      lag(END_DATE) < START_DATE ~ lag(END_DATE)
   )) %>%
   ungroup()

# mark but remove "visits" with >1 index isolate
empDF <- empDF %>% 
   group_by(PERSON_ID, ORDER_DATE) %>% 
   mutate(MULT_ISOLATES = ifelse(n() > 1L, TRUE, FALSE)) %>% 
   slice_min(RESULT_DATE, with_ties = FALSE) %>% 
   ungroup()


# JOIN
empDF <- empDF %>%
   mutate(JOIN_START = ORDER_DATE - 30,
          JOIN_END = RESULT_DATE + 30) %>%
   left_join(x = .,
             y = abxDF,
             by = join_by(PERSON_ID,
                          JOIN_START <= START_DATE,
                          JOIN_END >= START_DATE)) %>%
   # if the most recent Rx was before empiric window, get num days
   mutate(DAYS_SINCE_PRV_ABX = ifelse(DATE_OF_MOST_RECENT_ABX < ORDER_DATE, as.integer(ORDER_DATE - DATE_OF_MOST_RECENT_ABX), NA)) %>%
   select(-JOIN_START, -JOIN_END) %>%
   group_by(PERSON_ID, ORDER_DATE) %>% # 7,521
   relocate(ABX, START_DATE, END_DATE, .before=BUG) %>%
   filter(!is.na(ABX)) %>% # 6,765 (~10%)            
   filter(!any(DAYS_SINCE_PRV_ABX <= 14, na.rm=TRUE)) %>% # 6,295 (~%) keep for now
   select(-DATE_OF_MOST_RECENT_ABX, -DAYS_SINCE_PRV_ABX) %>%
   ungroup()

# 79,009 (5,192 "visits")
empDF %>% group_by(PERSON_ID, ORDER_DATE)

# which are antibiotics are given
t <- empDF %>% 
   select(PERSON_ID, ORDER_DATE, ABX) %>%
   distinct() %>%
   select(ABX) %>% table()
t <- head(sort(t, decreasing=TRUE), n=12)
t <- t / 6765
x <- barplot(t, horiz=TRUE, plot=F)
{
   pdf(file = '~/Desktop/EHR/EHR work/ALL/plots/BloodAbxAdmin/AbxCounts.pdf', width=6, height=5)
   par(mar=c(2, 10, 1, 1))
   barplot(t, horiz=TRUE, names.arg=NA, xlim=c(0,1))
   text(x = -0.01, y = x, adj=1, labels=stringr::str_to_sentence(names(t)), xpd=NA)
   text(x = t+0.01, y=x, adj=0, labels=paste0(round(t,3)*100, '%'))
   dev.off()
}
rm(t, x)



# still contains rows where antibiotics were prescribed leading up to empiric window
empDF %>% count(RESULT_DATE - ORDER_DATE, sort=TRUE)
# 3, 6, 2, 4, 5, 7

# PLOTTING
plot_path <- '~/Desktop/EHR/EHR work/ALL/plots/BloodAbxAdmin/'


# 81,766 days of prescription
x <- empDF %>%
   select(PERSON_ID, ORDER_DATE, RESULT_DATE, START_DATE) %>%
   distinct() %>%
   mutate(X = as.integer(START_DATE - ORDER_DATE),
          D = as.integer(RESULT_DATE - ORDER_DATE))

makeBarplot <- function(main='') {
   t <- x %>%
      filter(D == as.integer(main)) %>%
      mutate(RX = 1) %>%
      tidyr::pivot_wider(id_cols = c(PERSON_ID, ORDER_DATE),
                         values_from = RX,
                         names_from = X,
                         values_fill = 0) %>%
      select(`0`, `1`, `2`, `3`, `4`, `5`, `6`, `7`, `8`, `9`, `10`)
   h <- colSums(t, na.rm=T) / nrow(t)
   b <- barplot(h, plot=F)
   barplot(h, xlab='Days since blood culture order', ylim=c(0, 1), main=paste0(main, ' day result delay\n(n = ', prettyNum(nrow(t),big.mark=','), ')'), ylab='% patients receiving abx')
   abline(v = b[names(h) == main], lwd=1.1, lty=2)
   text(x = b[names(h) == main]+0.1, y=0.95, adj=0, labels='result day', font=3)
}
{
   pdf(file = '~/Desktop/EHR/EHR work/ALL/plots/BloodAbxAdmin/TimingOfAbxAdmin.pdf', width = 10, height = 5.5)
   par(mfrow=c(2,3), mgp=c(1.6, 0.4, 0), mar=c(3,3,3.5,1), tck=-0.015)
   makeBarplot(main='2')
   makeBarplot(main='3')
   makeBarplot(main='4')
   makeBarplot(main='5')
   makeBarplot(main='6')
   dev.off()
}
rm(makeBarplot)


# how many people have been "covered" by each day?
h <- vector('list', length(2:7))
for (d in seq_along(h)) {
   t <- x %>%
      filter(D == as.integer(d+1)) %>%
      mutate(RX = 1) %>%
      tidyr::pivot_wider(id_cols = c(PERSON_ID, ORDER_DATE),
                         values_from = RX,
                         names_from = X,
                         values_fill = 0) %>%
      select(`0`, `1`, `2`, `3`, `4`, `5`, `6`, `7`, `8`, `9`, `10`)
   for (i in seq_len(nrow(t))) {
      print(i)
      w <- which(t[i,] == 1L)[1]
      if (is.na(w)) next
      t[i, w:ncol(t)] <- 1L
   }
   h[[d]] <- sapply(t, sum) / nrow(t)
}
names(h) <- 2:7
plot(NA, xlim=c(0, 10), ylim=c(0, 1))
for (d in 1:6) {
   lines(x = 0:10, y = h[[d]])
}
x %>%
   mutate(RX = 1) %>%
   tidyr::pivot_wider(id_cols = c(PERSON_ID, ORDER_DATE),
                      values_from = RX,
                      names_from = X,
                      values_fill = 0) %>%
   select(`0`, `1`, `2`, `3`, `4`, `5`, `6`, `7`, `8`, `9`, `10`) %>%
   count(`0`, `1`, sort=TRUE)
rm(h, t, x, d, i, w)


# how often are antibiotics administered before result day? - ~always
x <- empDF %>%
   slice_min(START_DATE, by = c(PERSON_ID, ORDER_DATE), with_ties = FALSE) %>%
   mutate(BEFORE_RESULT = START_DATE < RESULT_DATE)
x %>% count(BEFORE_RESULT)





empDF <- empDF %>%
   reframe(EMPIRIC0 = list(sort(unique(ABX[START_DATE == ORDER_DATE]))),
           EMPIRIC1 = list(sort(unique(ABX[START_DATE == (ORDER_DATE+1)]))),
           EMPIRIC2 = list(sort(unique(ABX[START_DATE == (ORDER_DATE+2) & START_DATE < RESULT_DATE]))),
           EMPIRIC3 = list(sort(unique(ABX[START_DATE == (ORDER_DATE+3) & START_DATE < RESULT_DATE]))),
           EMPIRIC4 = list(sort(unique(ABX[START_DATE == (ORDER_DATE+4) & START_DATE < RESULT_DATE]))),
           EMPIRIC5 = list(sort(unique(ABX[START_DATE == (ORDER_DATE+5) & START_DATE < RESULT_DATE]))),
           TARGETED = list(sort(unique(ABX[START_DATE == RESULT_DATE]))),
           .by = c(PERSON_ID, ORDER_DATE, RESULT_DATE, BUG, MULT_ISOLATES, CEFEPIME:DELAFLOXACIN)) %>%
   mutate(DELAY = as.integer(RESULT_DATE - ORDER_DATE)) %>%
   relocate(DELAY, EMPIRIC0, EMPIRIC1, EMPIRIC2, EMPIRIC3, EMPIRIC4, EMPIRIC5, TARGETED, .before=BUG)

#empDF <- empDF %>% filter(!MULT_ISOLATES)


save(empDF, file = '~/Desktop/EHR/EHR work/RdataFiles/ASTs_AbxAdmin_blood_2015_2023.Rdata')


# empiric and targeted course sizes
makePlot <- function(d, main='', targ=FALSE) {
   x <- empDF[[paste0('EMPIRIC', d)]][empDF$DELAY > d]
   if (targ) x <- empDF$TARGETED
   n <- sum(empDF$DELAY > d)
   x <- table(lengths(x)) / n
   x <- c(x[names(x) %in% as.character(0:3)], sum(x[!names(x) %in% as.character(0:3)]))
   names(x)[names(x) == ''] <- '>=3'
   b <- barplot(x, plot=F)
   barplot(x, ylim=c(0, 1), xlab='Number of antibiotics', main=paste0(main, ' (n = ', n, ')'))
   text(x = b, y = x+0.05, labels=paste0(round(x, 3) * 100, '%'))
}
{
   pdf(file = paste0(plot_path, 'CourseSize.pdf'), height=6.5, width=9)
   par(mfrow=c(3,3), mar=c(4, 3, 2, 1))
   makePlot(0, 'day of order')
   makePlot(1, '1 day after order')
   makePlot(2, '2 days after order')
   makePlot(3, '3 days after order')
   makePlot(4, '4 days after order')
   makePlot(5, '5 days after order')
   makePlot(0, 'Day of result', targ=TRUE)
   dev.off()
}
rm(makePlot)




# DID THERAPIES NARROW IN THE TARGETED WINDOW?
# 41% narrowed
empDF %>% filter(lengths(TARGETED) > 0) # 3,503
empDF %>% filter(lengths(TARGETED) > 0) %>% count(lengths(TARGETED) < lengths(EMPIRIC0))
empDF %>% filter(lengths(TARGETED) > 0) %>% count(lengths(TARGETED) == lengths(EMPIRIC1))


t <- empDF %>%
   select(EMPIRIC0, EMPIRIC1, EMPIRIC2, EMPIRIC3, TARGETED) %>%
   mutate(across(everything(), ~ lengths(.)))
m <- sapply(t, mean)

plot(x = 0:4,
     y = m,
     ylim=c(0,10),
     type = 'l',
     lwd = 2)
apply(t, 1, function(y) lines(x=0:6, y=y, lwd=0.1, col='#44444411'))
rm(t, m)



#########################################################
#### WHAT IS PRESCRIBED EMPIRICALLY vs. TARGETED??? #####
#########################################################
{
   abbr <- tibble(read.table(file = '~/Desktop/EHR/EHR work/data/ABX_ABBR.txt', header = TRUE)) %>% 
      select(-Class) %>%
      mutate(Antibiotic_Name = gsub('-', '/', Antibiotic_Name),
             Antibiotic_Name = gsub('_', ' ', Antibiotic_Name),
             Antibiotic_Name = toupper(Antibiotic_Name))
   abbr <- setNames(abbr$Abbreviation, abbr$Antibiotic_Name)
   
   abbr2 <- tibble(read.table(file = '~/Desktop/EHR/EHR work/data/ABX_ABBR.txt', header = TRUE)) %>% 
      select(-Class) %>%
      mutate(Antibiotic_Name = gsub('-', '/', Antibiotic_Name),
             Antibiotic_Name = gsub('_', ' ', Antibiotic_Name),
             Antibiotic_Name = toupper(Antibiotic_Name)) %>%
      filter(Abbreviation %in% c('CIP', 'CFZ', 'SAM', 'MET', 'FEP', 'AZM', 'CRO', 'TZP', 'VAN', 'LVX')) %>% arrange(Abbreviation)
   
   # EMPIRIC
   emp0 <- sapply(empDF$EMPIRIC0[empDF$DELAY > 0L], function(x) paste(abbr[x], collapse='+'))
   emp1 <- sapply(empDF$EMPIRIC1[empDF$DELAY > 1L], function(x) paste(abbr[x], collapse='+'))
   emp2 <- sapply(empDF$EMPIRIC2[empDF$DELAY > 2L], function(x) paste(abbr[x], collapse='+'))
   emp3 <- sapply(empDF$EMPIRIC3[empDF$DELAY > 3L], function(x) paste(abbr[x], collapse='+'))
   emp4 <- sapply(empDF$EMPIRIC4[empDF$DELAY > 4L], function(x) paste(abbr[x], collapse='+'))
   emp5 <- sapply(empDF$EMPIRIC5[empDF$DELAY > 5L], function(x) paste(abbr[x], collapse='+'))
   empT <- sapply(empDF$TARGETED, function(x) paste(abbr[x], collapse='+'))
   
   n <- 8
   b <- barplot(rep(1, n), horize=T, plot=F)
   makeBarplot <- function(d, main='', targ=FALSE) {
      x <- empDF[[paste0('EMPIRIC', d)]][empDF$DELAY > d]
      if (targ) x <- empDF$TARGETED
      #x <- x[lengths(x) > 0]
      x <- sapply(x, function(x) paste(abbr[x], collapse='+'))
      h <- head(sort(table(unlist(x)), decreasing = TRUE) / length(x), n)
      names(h)[names(h) == ''] <- 'No abx'
      barplot(h, horiz=TRUE, names.arg=NA, xlim=c(0,1), xpd=NA)
      axis(side=2, at=b, labels=names(h), las=1, tick=F)
      title(main=paste0(main, '\n(n = ', length(x), ')'), line=0.25)
      text(x = h + 0.01, y = b, adj = 0, xpd=NA,
           labels = paste0(round(h,3)*100, '%'))
   }
   makeLegend <- function() {
      par(mar=c(0,0,0,0))
      plot(NA, xlim=c(1, 2), ylim=c(0,nrow(abbr)+1.5), axes=F, ann=F)
      text(x = 1.25, y = 1:nrow(abbr), adj=0, xpd=NA,
           labels = paste0(abbr$Abbreviation, ' - ', abbr$Antibiotic_Name))
   }
   {
      pdf(file = paste0(plot_path, '/WhichAntibioticsEmpiric.pdf'), height=4.5, width=8)
      par(mfrow=c(2,3), mar=c(2.5, 5.5, 2.5, 2), mgp=c(1.25,0.4,0), tck=-0.015)
      makeBarplot(0, main='Empiric = day of order')
      makeBarplot(1, main='Empiric = 1 day after order')
      makeBarplot(2, main='Empiric = 2 days after order')
      makeBarplot(3, main='Empiric = 3 days after order')
      makeBarplot(4, main='Empiric = 4 days after order')
      makeBarplot(0, main='Targeted = day of result', targ=TRUE)
      dev.off()
   }
   rm(makeBarplot, makeLegend, n, b, abbr, abbr2)
}



# DOES THE EMPIRIC THERAPY DEPEND UPON BUG?
head(empDF$EMPIRIC_ABBR)
empDF$EMPIRIC_ABBR <- sapply(empDF$EMPIRIC_ABBR, paste, collapse='+')

n <- 8
b <- barplot(rep(1, n), horize=T, plot=F)
makeBarplot <- function(bug = '') {
   if (bug == 'Gram negatives') {
      x <- empDF$EMPIRIC_ABBR[grep('Klebsiella|Escherichia|Acinetobacter|Citrobacter|Enterobacter|Morganella|Serratia|Proteus|Pseudomonas', empDF$BUG)]
   } else if (bug == 'Gram positives') {
      x <- empDF$EMPIRIC_ABBR[grep('Enterococ|Streptoc|Staph', empDF$BUG)]
   } else {
      x <- empDF$EMPIRIC_ABBR[which(empDF$BUG == bug)]
   }
   h <- sort(table(unlist(x)), decreasing = TRUE) / length(x)
   h <- head(h, n)
   if (any(names(h) == '')) names(h)[names(h) == ''] <- 'No abx'
   barplot(h, horiz=TRUE, names.arg=NA, xlim=c(0,0.5))
   axis(side=2, at=b, labels=names(h), las=1, tick=F)
   title(main=paste0(bug, ' (n = ', length(x), ')'), line=0.25, xpd=NA)
   text(x = h + 0.01, y = b, adj = 0, xpd=NA,
        labels = paste0(round(h,3)*100, '%'))
}

{
   pdf(file = '~/Desktop/EHR/EHR work/ALL/plots/BloodAbxAdmin/EmpiricByBug.pdf', height=9, width=7)
   par(mfcol=c(3, 2), mar=c(2.5, 7, 1.5, 2), mgp=c(1.25,0.4,0), tck=-0.015)
   makeBarplot('Staphylococcus epidermidis')
   makeBarplot('Staphylococcus aureus')
   makeBarplot('Gram positives')
   makeBarplot('Escherichia coli')
   makeBarplot('Klebsiella pneumoniae')
   makeBarplot('Gram negatives')
   dev.off()
}
bugs <- empDF %>% count(BUG, sort=TRUE)
rm(makeBarplot, bugs, b, n)



head(empDF$TARGETED_ABBR)
empDF$TARGETED_ABBR <- sapply(empDF$TARGETED_ABBR, paste, collapse='+')

n <- 8
b <- barplot(rep(1, n), horize=T, plot=F)
makeBarplot <- function(bug = '') {
   if (bug == 'Gram negatives') {
      x <- empDF$TARGETED_ABBR[grep('Klebsiella|Escherichia|Acinetobacter|Citrobacter|Enterobacter|Morganella|Serratia|Proteus|Pseudomonas', empDF$BUG)]
   } else if (bug == 'Gram positives') {
      x <- empDF$TARGETED_ABBR[grep('Enterococ|Streptoc|Staph', empDF$BUG)]
   } else {
      x <- empDF$TARGETED_ABBR[which(empDF$BUG == bug)]
   }
   h <- sort(table(unlist(x)), decreasing = TRUE) / length(x)
   h <- head(h, n)
   if (any(names(h) == '')) names(h)[names(h) == ''] <- 'No abx'
   barplot(h, horiz=TRUE, names.arg=NA, xlim=c(0,0.5))
   axis(side=2, at=b, labels=names(h), las=1, tick=F)
   title(main=paste0(bug, ' (n = ', length(x), ')'), line=0.25, xpd=NA)
   text(x = h + 0.01, y = b, adj = 0, xpd=NA,
        labels = paste0(round(h,3)*100, '%'))
}

{
   pdf(file = '~/Desktop/EHR/EHR work/ALL/plots/BloodAbxAdmin/TargetedByBug.pdf', height=9, width=7)
   par(mfcol=c(3, 2), mar=c(2.5, 7, 1.5, 2), mgp=c(1.25,0.4,0), tck=-0.015)
   makeBarplot('Staphylococcus epidermidis')
   makeBarplot('Staphylococcus aureus')
   makeBarplot('Gram positives')
   makeBarplot('Escherichia coli')
   makeBarplot('Klebsiella pneumoniae')
   makeBarplot('Gram negatives')
   dev.off()
}
bugs <- empDF %>% count(BUG, sort=TRUE)
rm(makeBarplot, bugs, b, n)





rm(plot_path)

empDF <- empDF %>%
   select(-EMPIRIC_ABBR, -TARGETED_ABBR)















