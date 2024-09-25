library(dplyr)
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/data_path_name.Rdata')
# load(file = '~/Desktop/EHR/EHR work/RdataFiles/ASTs_AbxAdmin_blood_2017_2023_flagged.Rdata')
load(file = paste0(data_path_name, 'ASTs_AbxAdmin_blood_2017_2023_imputed_flagged.Rdata'))
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/antibiotic_names/abbr_named.Rdata')
abbr <- c(abbr, setNames('No abx', 'No abx'))








# Merge flag0 for multiple blood isolate cases
{
   # x <- empDF %>%
   #    # filter(MULT_BLOOD_ISO, FLAG0 != 'No empiric therapy given') %>%
   #    summarise(FLAG0 = paste(sort(unique(FLAG0)), collapse='+'),
   #              .by = c(PERSON_ID, ORDER_DAY)) %>%
   #    mutate(NEW_FLAG0 = case_when(
   #       FLAG0 == 'CONCORDANT+NOT TESTED' ~ 'CONCORDANT',
   #       FLAG0 == 'CONCORDANT+DISCORDANT' ~ 'DISCORDANT',
   #       FLAG0 == 'DISCORDANT+NOT TESTED' ~ 'DISCORDANT',
   #       FLAG0 == 'CONCORDANT+DISCORDANT+NOT TESTED' ~ 'DISCORDANT',
   #       .default = FLAG0
   #    ))
   # x %>% count(FLAG0, sort=TRUE)
   # x %>% count(NEW_FLAG0, sort=TRUE) # 10.1% of given antibiotics (7% of infections)
   # rm(x)
   # empDF <- rbind(
   #    empDF %>%
   #       group_by(PERSON_ID, ORDER_DAY) %>%
   #       filter(n() > 1L) %>% 
   #       mutate(mFLAG0 = case_when(
   #          all(FLAG0 == 'No empiric therapy given') ~ 'No empiric therapy given',
   #          any(FLAG0 == 'DISCORDANT') ~ 'DISCORDANT',
   #          any(FLAG0 == 'CONCORDANT') ~ 'CONCORDANT',
   #          all(FLAG0 == 'NOT TESTED') ~ 'NOT TESTED'
   #       )) %>%
   #       ungroup(),
   #    empDF %>% 
   #       filter(n() == 1L, .by = c(PERSON_ID, ORDER_DAY)) %>% 
   #       mutate(mFLAG0 = NA)
   # ) %>% relocate(mFLAG0, .before=FLAG0)
   # empDF %>% select(mFLAG0, FLAG0) %>% table()
}






empDF %>% filter(FLAG0 != 'No empiric therapy given') %>% count(FLAG0, sort=TRUE) %>% mutate(r = n / sum(n) * 100)
empDF %>%
   filter(grepl('aureus', BUG) | !grepl('Staph|Candida', BUG), FLAG0 != 'No empiric therapy given') %>%
   count(FLAG0, sort=TRUE) %>% 
   mutate(r = n / sum(n) * 100)

empDF %>%
   filter(FLAG0 != 'No empiric therapy given') %>%
   count(FLAG0, sort=TRUE) %>% 
   mutate(r = n / sum(n) * 100)

empDF %>% count(FLAGT, sort=TRUE) %>% mutate(r = n / nrow(empDF) * 100)



# for which bugs is not tested overrepresented?
t <- empDF %>%
   filter(grepl('aureus', BUG) | !grepl('Staph|Candida', BUG)) %>%
   mutate(NT = FLAG0 == 'NOT TESTED') %>%
   select(BUG, NT) %>%
   table()
t <- data.frame(NOT_TESTED = t[,'TRUE'],
                TESTED = t[,'FALSE'],
                row.names = rownames(t))
t <- t[rowSums(t) > 50L, ]
t$rate <- round(t$NOT_TESTED / rowSums(t), 2)
t <- t[order(t$rate), ]
t
rm(t)




empDF %>%
   summarise(num_disc = sum(FLAG0 == 'DISCORDANT'),
             total = sum(FLAG0 != 'No empiric therapy given'),
             disc_rate = round(num_disc / total * 100, 1),
             .by = BUG) %>%
   filter(total > 100L) %>%
   arrange(desc(disc_rate)) %>%
   View()
# E. faecium ≈ 69%
# E. faecalis ≈ 33%
# P. aeruginosa ≈ 33%
# E. coli ≈ 8.1%
# S. aureus ≈ 8.8%

empDF %>%
   summarise(num_disc = sum(FLAGT == 'DISCORDANT'),
             total = sum(FLAGT %in% c('DISCORDANT', 'CONCORDANT')),
             disc_rate = round(num_disc / total * 100, 1),
             .by = BUG) %>%
   filter(total > 100L) %>%
   arrange(desc(disc_rate)) %>%
   View()
# E. faecium ~ 10.6%




makeBarplot <- function(x, main) {
   n <- 8
   b <- barplot(rep(1,n), plot=F)
   h <- head(sort(table(unlist(x)), decreasing = TRUE), n)
   names(h)[names(h) == ''] <- 'No abx'
   barplot(h, horiz=TRUE, names.arg=NA, xpd=NA, xaxt='n', main=paste0('n = ', length(x)))
   axis(side=2, at=b, labels=names(h), las=1, tick=F)
   text(x = h + 1, y = b, adj = 0, xpd=NA, labels = h)
}




bug <- 'Enterobacter cloacae'
t <- sapply(empDFo %>% filter(BUG == bug) %>% select(CEFEPIME:DELAFLOXACIN), function(x) c(TOTAL = sum(!is.na(x)), R = sum(x,na.rm=T)))
t[, order(colSums(t))]

empDF %>% filter(BUG == bug) %>% select(CEFTRIAXONE, `PIPERACILLIN/TAZOBACTAM`) %>% table()
empDF %>% filter(BUG == bug, FLAG0 != 'No empiric therapy given') %>% summarise(n(), .by=FLAG0)
empDF %>% filter(BUG == bug, FLAG0 != 'No empiric therapy given') %>% 
   summarise(n(), disc = sum(FLAG0 == 'DISCORDANT') / n(), .by = `PIPERACILLIN/TAZOBACTAM`)






x %>% count(EMP0_ABBR, sort=TRUE)
data.frame(head(sort(table(unlist(x$EMPIRIC0)), decreasing = TRUE), n=8))

x <- empDF %>% filter(BUG == 'Enterococcus faecium', FLAG0 == 'DISCORDANT')
data.frame(head(sort(table(sapply(x$EMPIRIC0, function(x) paste(abbr[x], collapse='+'))), decreasing = TRUE), n=10))
data.frame(head(sort(table(unlist(x$EMPIRIC0)), decreasing = TRUE), n=8))




# VRE has discordant empiric rate of 71% vs. 23% for vancomycin-susceptible
e <- empDF %>% filter(BUG == 'Enterococcus faecium') %>% mutate(across(CEFEPIME:DELAFLOXACIN, ~ ifelse(is.na(.), -1, .)))
t <- e %>% 
   mutate(DISC = FLAG0 == c('DISCORDANT')) %>%
   select(VANCOMYCIN, DISC) %>%
   filter(VANCOMYCIN >= 0) %>%
   table()
t <- round(t(apply(t, 1, function(x) x / sum(x))) * 100, 1)
t

e <- e %>% filter(FLAG0 == c('DISCORDANT'))
e %>% count(EMPIRIC0, sort=TRUE) %>% View()
e %>% filter(grepl('PIPERACILLIN/TAZOBACTAM', EMPIRIC0)) %>% select(`PIPERACILLIN/TAZOBACTAM`, AMPICILLIN) %>% table()
e %>% filter(grepl('VANCOMYCIN', EMPIRIC0)) %>% select(VANCOMYCIN, `PIPERACILLIN/TAZOBACTAM`) %>% table()



# MRSA leads to discordant empiric therapy 18% of time vs. <1% for MSSA
e <- empDF %>% filter(grepl('Staphylococcus aureus', BUG)) %>% mutate(across(CEFEPIME:DELAFLOXACIN, ~ ifelse(is.na(.), -1, .)))
t <- e %>%
   mutate(DISC = FLAG0 == c('DISCORDANT')) %>%
   select(OXACILLIN, DISC) %>%
   filter(OXACILLIN != -1) %>%
   table()
t <- round(t(apply(t, 1, function(x) x / sum(x))) * 100, 1)


e <- empDF %>% filter(grepl('Staphylococcus aureus', BUG)) %>% mutate(across(CEFEPIME:DELAFLOXACIN, ~ ifelse(is.na(.), -1, .)))
t <- e %>%
   mutate(DISC = FLAG0 == c('DISCORDANT')) %>%
   select(OXACILLIN, DISC) %>%
   filter(OXACILLIN != -1) %>%
   table()
t <- round(t(apply(t, 1, function(x) x / sum(x))) * 100, 1)


# Extended-spectrum cephalosporin=resistant E. coli 37% vs. 3%!
e <- empDF %>% filter(BUG == 'Escherichia coli')# %>% mutate(across(CEFEPIME:DELAFLOXACIN, ~ ifelse(is.na(.), -1, .)))
a <- sapply(e[24:length(e)], function(x) sum(x >= 0, na.rm=T))
a <- a[a > 0]

e <- e %>% 
   mutate(ESCR = as.integer(CEFTAZIDIME | CEFTRIAXONE | CEFOTAXIME | CEFEPIME)) %>%
   mutate(DISC = FLAG0 == c('DISCORDANT'))
t <- e %>% 
   select(ESCR, DISC) %>% 
   table()
t <- round(t(apply(t, 1, function(x) x / sum(x))) * 100, 1)
t







empDF %>%
   filter(FLAG0 == 'NOT TESTED', BUG == 'Staphylococcus aureus') %>% #select(OXACILLIN, AMPICILLIN, PENICILLIN) %>% table()
   count(EMPIRIC0, sort=TRUE) %>%
   View()

empDF %>%
   filter(FLAG0 == 'NOT TESTED', BUG == 'Staphylococcus aureus') %>%
   filter(grepl('PIPERACILLIN/TAZOBACTAM', EMPIRIC0)) %>%
   mutate(across(CEFEPIME:DELAFLOXACIN, ~ ifelse(is.na(.), -1, .))) %>%
   select(OXACILLIN, `PIPERACILLIN/TAZOBACTAM`) %>%
   table()



e <- empDF %>% filter(FLAG0 == 'NOT TESTED', BUG == 'Escherichia coli') %>% mutate(across(CEFEPIME:DELAFLOXACIN, ~ ifelse(is.na(.), -1, .)))
e %>%
   count(EMPIRIC0, sort=TRUE) %>%
   View()
e <- e %>% mutate(CEFAZ_CEFOT = as.integer(CEFAZOLIN | CEFOTAXIME))
e %>% select(CEFTRIAXONE, CEFAZ_CEFOT) %>% table()
rm(e)




# # HOW?? s
# tdisc <- empDF %>% filter(FLAGT == 'DISCORDANT') # only 40 (out of >4,000)
# tdisc$TARGETED <- sapply(tdisc$TARGETED, paste, collapse='+')
# tdisc %>% count(BUG) %>% left_join(bugs, by=join_by(BUG)) %>% mutate(f = n.x / n.y) %>% arrange(desc(f))
# tdisc %>% count(TARGETED, sort=TRUE)
# tdisc %>% count(BUG, TARGETED, sort=TRUE)
# rm(tdisc)
# 
# 
# edisc <- empDF %>% filter(FLAG0 == 'DISCORDANT') # 250 (out of >4,000)
# edisc$EMPIRIC0 <- sapply(edisc$EMPIRIC0, paste, collapse='+')
# edisc %>% count(BUG, sort=TRUE)
# edisc %>% count(TARGETED, sort=TRUE)
# edisc %>% count(BUG, TARGETED, sort=TRUE)
# rm(tdisc)


# previous stuff from dispense data - update this soon
{
   t <- empDF %>%
      mutate(FLAGT = case_when(
         grepl('^R', FLAGT) ~ 'maybeDISCORDANT',
         FLAGT %in% c('NOT TESTED', 'NT + nev', 'NEVER TESTED') ~ 'NOT/NEVER TESTED',
         .default = FLAGT
      ),
      EFLAG = case_when(
         grepl('^R', EFLAG) ~ 'maybeDISCORDANT',
         EFLAG %in% c('NOT TESTED', 'NT + nev', 'NEVER TESTED') ~ 'NOT/NEVER TESTED',
         .default = EFLAG
      )) %>%
      select(EFLAG, FLAGT) %>%
      table()
   t <- t / rowSums(t)
   t
   # mostly therapies staying concordant or going from discordant --> concordant
   # lots of concordant empiric therapies were entirely discontinuted (no targeted therapy)
   # surprising amounts of concordant empiric prescriptions changing to discordant targeted therapies (3.4% of concordant empiric therapies)
   
   
   
   bugs <- head(sort(table(empDF$BUG[!empDF$EFLAG %in% c('CONCORDANT', 'DISCORDANT', 'No empiric therapy given')]), decreasing=TRUE), n=12)
   
   for (i in seq_along(bugs)) {
      w <- which(!empDF$EFLAG %in% c('CONCORDANT', 'DISCORDANT', 'No empiric therapy given') & empDF$BUG == names(bugs)[i])
      e <- empDF$EMPIRIC[w]
      s <- data.frame(round(sort(table(unlist(e)), decreasing=TRUE) / length(w) * 100))
      s$pcnt_NA <- round(sapply(s$Var1, function(x) sum(is.na(empDF[[x]][w]))) / length(w) * 100)
      s <- s[s$pcnt_NA > 0, ]
      s <- s[s$Freq > 5, ]
      print(bugs[i])
      print(s)
      print('________________________________________')
   }
   
   # Staphylococcus epidermidis: n = 1174
   #                      Var1 Freq pcnt_NA
   # 1             CEFTRIAXONE   38      68
   # 3 PIPERACILLIN/TAZOBACTAM   22      67
   # 4                CEFEPIME   13      69
   # 5            LEVOFLOXACIN   13      53
   # 7           CIPROFLOXACIN    7       3
   
   # Coagulase Negative Staph: n = 745 
   #                      Var1 Freq pcnt_NA
   # 1             CEFTRIAXONE   54      98
   # 2            AZITHROMYCIN   29      99
   # 3 PIPERACILLIN/TAZOBACTAM   14      99
   # 4    AMPICILLIN/SULBACTAM    6      90
   # 5           CIPROFLOXACIN    6      69
   
   # Escherichia coli: n = 587 
   #                      Var1 Freq pcnt_NA
   # 1             CEFTRIAXONE   38      68
   # 3 PIPERACILLIN/TAZOBACTAM   22      67
   # 4                CEFEPIME   13      69
   # 5            LEVOFLOXACIN   13      53
   # 7           CIPROFLOXACIN    7       3
   
   # Staphylococcus hominis: 376 
   #                      Var1 Freq pcnt_NA
   # 1             CEFTRIAXONE   52      96
   # 2            AZITHROMYCIN   25     100
   # 3 PIPERACILLIN/TAZOBACTAM   16     100
}









t <- empDF %>% filter(BUG == 'Staphylococcus aureus') %>% 
   group_by(PERSON_ID, ORDER_DAY, BUG) %>% 
   select(PERSON_ID, ORDER_DAY, START_DATE) %>% 
   mutate(START_DATE = as.Date(substr(START_DATE,1,10))) %>% 
   distinct() %>% 
   summarise(min = min(START_DATE)) %>% 
   mutate(first = as.integer(min - ORDER_DAY))
t <- table(t$first)

barplot(t)
sum(t[names(t) %in% as.character(-30:-1)])
t[names(t) == '0']
t[names(t) == '1']
sum(t[as.integer(names(t)) > 0])





load(file = paste0(data_path_name, 'ASTs_AbxAdmin_blood_2017_2023.Rdata'))





# Empiric course size and contents - depend up on gram stain?
{
   n <- 12
   b <- barplot(rep(1, n), horize=T, plot=F)
   makeBarplot <- function(bug = '') {
      if (bug == 'Gram negatives') x <- empDF$EMPIRIC0[grep('Klebsiella|Escherichia|Acinetobacter|Citrobacter|Enterobacter|Morganella|Serratia|Proteus|Pseudomonas', empDF$BUG)]
      if (bug == 'Gram positives') x <- empDF$EMPIRIC0[grep('Enterococ|Streptoc|Staph', empDF$BUG)]
      x <- sapply(x, function(a) paste(abbr[a], collapse='+'))
      h <- sort(table(unlist(x)), decreasing = TRUE) / length(x)
      h <- head(h, n)
      if (any(names(h) == '')) names(h)[names(h) == ''] <- 'No abx'
      barplot(h, horiz=TRUE, names.arg=NA, xlim=c(0,0.5))
      axis(side=2, at=b, labels=names(h), las=1, tick=F)
      title(main=paste0(bug, ' (n = ', prettyNum(length(x), big.mark=','), ')'), line=0.25, xpd=NA)
      text(x = h + 0.01, y = b, adj = 0, xpd=NA,
           labels = paste0(round(h,3)*100, '%'))
   }
   # depend up on gram stain?
   {
      pdf(file = paste0(plots_path_name, 'EmpiricByGram.pdf'), height=3.8, width=8)
      par(mfrow=c(1, 2), mar=c(2.5, 6.5, 1.5, 1), mgp=c(1.25,0.4,0), tck=-0.015)
      # makeBarplot('Staphylococcus epidermidis')
      # makeBarplot('Staphylococcus aureus')
      # makeBarplot('Enterococcus faecalis')
      makeBarplot('Gram positives')
      # makeBarplot('Escherichia coli')
      # makeBarplot('Klebsiella pneumoniae')
      # makeBarplot('Pseudomonas aeruginosa')
      makeBarplot('Gram negatives')
      dev.off()
   }
   
   rm(makeBarplot, makePlot, n, b, abbr)
}

# DID THERAPIES NARROW IN THE TARGETED WINDOW?
sum(lengths(empDF$TARGETED) > 0) / nrow(empDF) # 66%
sum(lengths(empDF$TARGETED) < lengths(empDF$EMPIRIC0)) / nrow(empDF) # 33% narrowed
sum(lengths(empDF$TARGETED) > lengths(empDF$EMPIRIC0)) / nrow(empDF) # 29.9% broadened
sum(lengths(empDF$TARGETED) == lengths(empDF$EMPIRIC0)) / nrow(empDF) # 36.9% stayed same



# WHAT IS PRESCRIBED EMPIRICALLY vs. TARGETED?
{
   abbr <- tibble(read.table(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/antibiotic_names/ABX_ABBR.txt', header = TRUE)) %>% 
      select(-Class) %>%
      mutate(Antibiotic_Name = gsub('-', '/', Antibiotic_Name),
             Antibiotic_Name = gsub('_', ' ', Antibiotic_Name),
             Antibiotic_Name = toupper(Antibiotic_Name))
   abbr <- setNames(abbr$Abbreviation, abbr$Antibiotic_Name)
   
   abbr2 <- tibble(read.table(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/antibiotic_names/ABX_ABBR.txt', header = TRUE)) %>% 
      select(-Class) %>%
      mutate(Antibiotic_Name = gsub('-', '/', Antibiotic_Name),
             Antibiotic_Name = gsub('_', ' ', Antibiotic_Name),
             Antibiotic_Name = toupper(Antibiotic_Name)) %>%
      filter(Abbreviation %in% c('CIP', 'CFZ', 'SAM', 'MET', 'FEP', 'AZM', 'CRO', 'TZP', 'VAN', 'LVX')) %>% arrange(Abbreviation)
   
   
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
      pdf(file = paste0(plots_path_name, 'WhichAntibioticsEmpiric.pdf'), height=4.5, width=8)
      par(mfrow=c(2,3), mar=c(2.5, 5.5, 2.5, 2), mgp=c(1.25,0.4,0), tck=-0.015)
      makeBarplot(0, main='Empiric = day of order')
      makeBarplot(1, main='Empiric = 1 day after order')
      makeBarplot(2, main='Empiric = 2 days after order')
      makeBarplot(3, main='Empiric = 3 days after order')
      makeBarplot(4, main='Empiric = 4 days after order')
      makeBarplot(0, main='Targeted = day of result', targ=TRUE)
      dev.off()
   }
   
   makeBarplot <- function(year='') {
      x <- empDF$EMPIRIC0[substr(empDF$ORDER_DAY,1,4) == year]
      x <- sapply(x, function(x) paste(abbr[x], collapse='+'))
      h <- head(sort(table(unlist(x)), decreasing = TRUE) / length(x), n)
      names(h)[names(h) == ''] <- 'No abx'
      barplot(h, horiz=TRUE, names.arg=NA, xlim=c(0,1), xpd=NA)
      axis(side=2, at=b, labels=names(h), las=1, tick=F)
      title(main=paste0('Empiric - ', year, '\n(n = ', length(x), ')'), line=0.25)
      text(x = h + 0.01, y = b, adj = 0, xpd=NA,
           labels = paste0(round(h,3)*100, '%'))
   }
   {
      pdf(file = paste0(plots_path_name, 'WhichAntibioticsEmpiricByYear.pdf'), height=6, width=8)
      par(mfrow=c(3,3), mar=c(2.5, 5.5, 2.5, 2), mgp=c(1.25,0.4,0), tck=-0.015)
      makeBarplot(year='2017')
      makeBarplot(year='2018')
      makeBarplot(year='2019')
      makeBarplot(year='2020')
      makeBarplot(year='2021')
      makeBarplot(year='2022')
      makeBarplot(year='2022')
      makeBarplot(year='2023')
      dev.off()
   }
   
   makeBarplot <- function(year='') {
      x <- empDF$TARGETED[substr(empDF$ORDER_DAY,1,4) == year]
      x <- sapply(x, function(x) paste(abbr[x], collapse='+'))
      h <- head(sort(table(unlist(x)), decreasing = TRUE) / length(x), n)
      names(h)[names(h) == ''] <- 'No abx'
      barplot(h, horiz=TRUE, names.arg=NA, xlim=c(0,1), xpd=NA)
      axis(side=2, at=b, labels=names(h), las=1, tick=F)
      title(main=paste0('Targeted - ', year, '\n(n = ', length(x), ')'), line=0.25)
      text(x = h + 0.01, y = b, adj = 0, xpd=NA,
           labels = paste0(round(h,3)*100, '%'))
   }
   {
      pdf(file = paste0(plots_path_name, 'WhichAntibioticsTargetedByYear.pdf'), height=6, width=8)
      par(mfrow=c(3,3), mar=c(2.5, 5.5, 2.5, 2), mgp=c(1.25,0.4,0), tck=-0.015)
      makeBarplot(year='2017')
      makeBarplot(year='2018')
      makeBarplot(year='2019')
      makeBarplot(year='2020')
      makeBarplot(year='2021')
      makeBarplot(year='2022')
      makeBarplot(year='2022')
      makeBarplot(year='2023')
      dev.off()
   }
   rm(makeBarplot, makeLegend, n, b, abbr, abbr2)
}

# DOES THE EMPIRIC THERAPY DEPEND UPON BUG?
{
   data.frame(head(sort(table(unlist(empDF$BUG)), decreasing=TRUE), n=14))
   abbr <- tibble(read.table(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/antibiotic_names/ABX_ABBR.txt', header = TRUE)) %>% 
      select(-Class) %>%
      mutate(Antibiotic_Name = gsub('-', '/', Antibiotic_Name),
             Antibiotic_Name = gsub('_', ' ', Antibiotic_Name),
             Antibiotic_Name = toupper(Antibiotic_Name))
   abbr <- setNames(abbr$Abbreviation, abbr$Antibiotic_Name)
   n <- 8
   b <- barplot(rep(1, n), horize=T, plot=F)
   
   makeBarplot <- function(bug = '', targ=FALSE) {
      if (targ) e <- empDF$TARGETED
      if (!targ) e <- empDF$EMPIRIC0
      if (bug == 'Gram negatives') {
         x <- e[grep('Klebsiella|Escherichia|Acinetobacter|Citrobacter|Enterobacter|Morganella|Serratia|Proteus|Pseudomonas', empDF$BUG)]
      } else if (bug == 'Gram positives') {
         x <- e[grep('Enterococ|Streptoc|Staph', empDF$BUG)]
      } else {
         x <- e[grep(bug, empDF$BUG)]
      }
      x <- sapply(x, function(a) paste(abbr[a], collapse='+'))
      h <- sort(table(unlist(x)), decreasing = TRUE) / length(x)
      h <- head(h, n)
      if (any(names(h) == '')) names(h)[names(h) == ''] <- 'No abx'
      barplot(h, horiz=TRUE, names.arg=NA, xlim=c(0,0.5), col=ifelse(names(h) == 'No abx', 'gray', '#00000066'))
      axis(side=2, at=b, labels=names(h), las=1, tick=F)
      title(main=paste0(bug, ' (n = ', prettyNum(length(x), big.mark=','), ')'), line=0.25, xpd=NA)
      text(x = h + 0.01, y = b, adj = 0, xpd=NA,
           labels = paste0(round(h,3)*100, '%'))
      text(x=0.5, adj=1, y=b[length(b)]+0.5, font=4, xpd=NA,
           labels=ifelse(targ, 'Targeted', 'Empiric'))
   }
   {
      pdf(file = paste0(plots_path_name, 'EmpiricByBug.pdf'), height=6.5, width=12)
      par(mfrow=c(2, 4), mar=c(2.5, 5, 1.5, 1), mgp=c(1.25,0.4,0), tck=-0.015)
      makeBarplot('Staphylococcus epidermidis')
      makeBarplot('Staphylococcus aureus')
      makeBarplot('Enterococcus faecalis')
      makeBarplot('Gram positives')
      makeBarplot('Escherichia coli')
      makeBarplot('Klebsiella pneumoniae')
      makeBarplot('Pseudomonas aeruginosa')
      makeBarplot('Gram negatives')
      dev.off()
   }
   {
      pdf(file = paste0(plots_path_name, 'TargetedByBug.pdf'), height=6.5, width=12)
      par(mfrow=c(2, 4), mar=c(2.5, 5, 1.5, 1), mgp=c(1.25,0.4,0), tck=-0.015)
      makeBarplot('Staphylococcus epidermidis', targ=T)
      makeBarplot('Staphylococcus aureus', targ=T)
      makeBarplot('Enterococcus faecalis', targ=T)
      makeBarplot('Gram positives', targ=T)
      makeBarplot('Escherichia coli', targ=T)
      makeBarplot('Klebsiella pneumoniae', targ=T)
      makeBarplot('Pseudomonas aeruginosa', targ=T)
      makeBarplot('Gram negatives', targ=T)
      dev.off()
   }
   rm(makeBarplot, b, n, abbr)
}





empDF %>% select(MULT_BLOOD_ISO, MULT_ISO) %>% table() # 61.3%
empDF %>% count(substr(ORDER_DAY,1,4)) # ~6,000 per year
empDF %>% filter(!MULT_BLOOD_ISO, !MULT_ISO) %>% count(substr(ORDER_DAY,1,4)) # ~3,300 per year

# let's look at S aureus...
# what is treated with?
x <- empDF %>%
   #filter(!MULT_BLOOD_ISO) %>%
   #filter(!MULT_ISO) %>% 
   filter(grepl('Staphylococcus aureus', BUG))
x %>% count(MULT_BLOOD_ISO)
x %>% count(MULT_ISO)

# number of single isolate S aureus cases ~doubled in 2020 (COVID?)
# number of cases in which multiple isolates were detected has declined 
empDF %>% filter(grepl('Staphylococcus aureus', BUG)) %>% mutate(year = substr(ORDER_DAY,1,4)) %>% select(year, MULT_BLOOD_ISO) %>% table()
empDF %>% filter(grepl('Staphylococcus aureus', BUG)) %>% mutate(year = substr(ORDER_DAY,1,4)) %>% select(year, MULT_ISO) %>% table()
t <- x %>% count(substr(ORDER_DAY,1,7), MULT_ISO) %>% tidyr::pivot_wider(names_from=MULT_ISO, values_from=n)
s <- seq(1, nrow(t), 3)
plot(NA, xlim=c(1, nrow(t)), ylim=c(0, 110), type='l', xaxt='n', xlab='', ylab='Cases', main=expression('Number of'~italic('S. aureus')~'infections over time'))
points(t$`FALSE`, type='l', lwd=2)
points(t$`TRUE`, type='l', lwd=2, col='red')
axis(side=1, at=seq_len(nrow(t))[s], labels=t$`substr(ORDER_DAY, 1, 7)`[s], las=3)
abline(v = seq_len(nrow(t))[s], lty=2, lwd=0.25)
legend('topright',
       legend=c('1', '>1'), title='No. isolates', col=c('black', 'red'), lty=1, lwd=2)
rm(t, s)


# not generally the case for all bugs
t <- empDF %>% mutate(year = substr(ORDER_DAY,1,4)) %>% select(year, MULT_ISO) %>% table()
t <- apply(t, 1, function(x) x['TRUE'] / sum(x))
t
t <- empDF %>% mutate(year = substr(ORDER_DAY,1,4)) %>% select(year, MULT_BLOOD_ISO) %>% table()
t <- apply(t, 1, function(x) x['TRUE'] / sum(x))
t
rm(t)




# BRING BACK AST RESULTS TO ASSESS STRAIN-SPECIFIC TARGETED THERAPIES
empDF_ <- empDF
load(file = paste0(data_path_name, 'ASTs_blood_2015_2023.Rdata'))
astDF <- empDF
empDF <- empDF_
rm(empDF_)
gc()
astDF <- astDF %>%
   select(PERSON_ID, ORDER_DAY, BUG, CEFEPIME:DELAFLOXACIN) %>%
   tidyr::pivot_longer(cols = CEFEPIME:DELAFLOXACIN,
                       values_to = 'STATUS',
                       names_to = 'ANTIBIOTIC',
                       values_drop_na = TRUE) %>%
   tidyr::pivot_wider(id_cols = c(PERSON_ID, ORDER_DAY, BUG),
                      values_from = STATUS,
                      names_from = ANTIBIOTIC,
                      values_fn = max)

# GET ANTIBIOTIC ABBREVIATIONS
abbr <- tibble(read.table(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/antibiotic_names/ABX_ABBR.txt', header = TRUE)) %>% 
   select(-Class) %>%
   mutate(Antibiotic_Name = gsub('-', '/', Antibiotic_Name),
          Antibiotic_Name = gsub('_', ' ', Antibiotic_Name),
          Antibiotic_Name = toupper(Antibiotic_Name))
abbr <- setNames(abbr$Abbreviation, abbr$Antibiotic_Name)
empDF$TARG_ABBR <- sapply(empDF$TARGETED, function(x) {
   a <- abbr[x]
   w <- which(is.na(a))
   if (length(w) > 0L)
      a[w] <- x[w]
   return(a)
})
names(empDF$TARG_ABBR) <- NULL
empDF$TARG_ABBR_COMB <- sapply(empDF$TARG_ABBR, paste, collapse='+')
rm(abbr)

# NARROW DOWN TO SINGLE BUGS, JOIN WITH ASTs
empDF_ <- empDF %>% filter(!MULT_ISO)
empDF_$BUG <- unlist(empDF_$BUG)
empDF_ <- left_join(x=empDF_, y=astDF, by=join_by(PERSON_ID, BUG, ORDER_DAY))
empDF_ <- empDF_ %>% filter(!CONTAM %in% c(1, 2))



makeBarplot <- function(bug, phen, main, it=F, separate=F, leg=F) {
   par(mgp=c(1.5, 0.5, 0), mar=c(2,5,2,1), cex.main=1)
   w <- grep(bug, empDF_$BUG)
   if (length(phen) == 1L) w <- w[which(empDF_[[names(phen)]][w] == phen)]
   if (length(phen) > 1L)  w <- Reduce(intersect, sapply(seq_along(phen), function(i) w[which(empDF_[[names(phen)[i]]][w] == phen[i])]))
   x <- empDF_$TARG_ABBR_COMB[w]
   n <- 8
   b <- barplot(rep(1, n), horize=T, plot=F)
   if (separate) {
      t <- strsplit(x, '+', fixed=TRUE)
      t[lengths(t) == 0L] <- ''
      t <- sort(table(unlist(t)), decreasing = TRUE)
      t <- head(t, n)
      abx_names <- names(t)
      abx_names[abx_names == ''] <- 'No abx'
      t <- sapply(seq_along(t), function(i){
         if (names(t)[i] == '') return(c('FALSE'=t[i], 'TRUE'=0L))
         tab <- table(grepl('+', x[grep(names(t)[i], x)], fixed=TRUE))
         if (length(tab) == 2L) return(tab)
         if (length(tab) == 1L & names(tab) == 'TRUE') return(c('FALSE'=0L, tab))
         if (length(tab) == 1L & names(tab) == 'FALSE') return(c(tab, 'TRUE'=0L))
      })
      t <- t[2:1,]
      colnames(t) <- abx_names
      rownames(t)[rownames(t) == 'FALSE'] <- 'Monotherapy'
      rownames(t)[rownames(t) == 'TRUE'] <- 'Combination'
      h <- t / length(x)
      names_h <- colnames(h)
      quan <- colSums(t)
      h_vals <- colSums(h)
   } else {
      t <- sort(table(x), decreasing = TRUE)
      t <- head(t, n)
      names(t)[names(t) == ''] <- 'No abx'
      h <- t / length(x)
      names_h <- names(h)
      quan <- t
      h_vals <- h
   }
   barplot(h, horiz=TRUE, names.arg=rep('',n), xlim=c(0,1), xpd=NA, legend=leg)#, args.legend=list(bty='n'))#, col=ifelse(names_h == 'No abx', 'white', '#666666'))
   text(x = h_vals + 0.01, y = b, adj = 0, xpd=NA, labels = quan)
   axis(side=2, at=b, labels=names_h, las=1, tick=F)
   if (it) title(main = substitute(bolditalic(main)), line=1)
   if (!it) title(main = main, line=1)
   title(main=paste0('\n(n = ', prettyNum(length(x), big.mark=','), ') - ', ifelse(separate, 'single', 'whole')), line=0)
}

# Staphylococcus aureus - MSSA = oxacillin or cefazolin, MRSA = vancomycin
{
   bug <- 'Staphylococcus aureus'
   pdf(file = paste0(plots_path_name, 'Targeted/', gsub(' ', '_', bug), '.pdf'), height=6, width=8)
   par(mfrow=c(2,2))
   makeBarplot(bug, c('OXACILLIN'=0), 'MSSA')
   makeBarplot(bug, c('OXACILLIN'=1), 'MRSA')
   makeBarplot(bug, c('OXACILLIN'=0), 'MSSA', separate=T)
   makeBarplot(bug, c('OXACILLIN'=1), 'MRSA', separate=T, leg=T)
   dev.off()
   rm(bug)  
}

# CoNS - TODO!

# Enterococcus - amp-S = ampicillin, van-R/amp-R - linezolid, daptomycin, tigecycline
{
   bug <- 'Enterococcus'
   pdf(file = paste0(plots_path_name, 'Targeted/', gsub(' ', '_', bug), '.pdf'), height=5, width=11)
   par(mfrow=c(2,4))
   makeBarplot(bug, c('VANCOMYCIN'=1, 'AMPICILLIN'=1), 'VAN-R, AMP-R')
   makeBarplot(bug, c('VANCOMYCIN'=1, 'AMPICILLIN'=0), 'VAN-R, AMP-S')
   makeBarplot(bug, c('VANCOMYCIN'=0, 'AMPICILLIN'=1), 'VAN-S, AMP-R')
   makeBarplot(bug, c('VANCOMYCIN'=0, 'AMPICILLIN'=0), 'VAN-S, AMP-S')
   makeBarplot(bug, c('VANCOMYCIN'=1, 'AMPICILLIN'=1), 'VAN-R, AMP-R', separate=T)
   makeBarplot(bug, c('VANCOMYCIN'=1, 'AMPICILLIN'=0), 'VAN-R, AMP-S', separate=T)
   makeBarplot(bug, c('VANCOMYCIN'=0, 'AMPICILLIN'=1), 'VAN-S, AMP-R', separate=T)
   makeBarplot(bug, c('VANCOMYCIN'=0, 'AMPICILLIN'=0), 'VAN-S, AMP-S', separate=T, leg=T)
   dev.off()
   rm(bug)
}

# Pseudomonas aeruginosa - cefepime or pip/tazo
{
   bug <- 'Pseudomonas aeruginosa'
   pdf(file = paste0(plots_path_name, 'Targeted/', gsub(' ', '_', bug), '.pdf'), height=6, width=4)
   par(mfrow = c(2,1))
   makeBarplot(bug, character(0L), 'P. aeruginosa', it=T)
   makeBarplot(bug, character(0L), 'P. aeruginosa', it=T, separate=T, leg=T)
   dev.off()
   rm(bug)
}

# Group A Strep - penicillin
{
   bug <- 'Streptococcus pyogenes'
   pdf(file = paste0(plots_path_name, 'Targeted/', gsub(' ', '_', bug), '.pdf'), height=6, width=4)
   par(mfrow = c(2,1))
   makeBarplot(bug, character(0L), 'Group A Strep')
   makeBarplot(bug, character(0L), 'Group A Strep', separate=T, leg=T)
   dev.off()
   rm(bug)
}

# Stenotrophomonas maltophilia - trim/sulfa
{
   bug <- 'Stenotrophomonas maltophilia'
   pdf(file = paste0(plots_path_name, 'Targeted/', gsub(' ', '_', bug), '.pdf'), height=6, width=4)
   par(mfrow = c(2,1))
   makeBarplot(bug, character(0L), 'S. maltophilia', it=T)
   makeBarplot(bug, character(0L), 'S. maltophilia', it=T, separate=T, leg=T)
   dev.off()
   rm(bug)
}

rm(makeBarplot)



#### STAPHYLOCOCCUS AUREUS treatment changes over time  + number of treatments are equal ####
{
   x <- empDF_ %>% filter(BUG == 'Staphylococcus aureus') # 2,493
   
   # Plot treatment %s over time for MRSA vs. MSSA
   {
      tr <- x %>%
         filter(OXACILLIN == 1L) %>%
         mutate(year = substr(ORDER_DAY,1,4)) %>%
         summarise(NO_TRT = sum(lengths(TARGETED) == 0L) / n(),
                   VAN = sum(grepl('VANCOMYCIN', TARGETED)) / n(),
                   DAP = sum(grepl('DAPTOMYCIN', TARGETED)) / n(),
                   TZP = sum(grepl('PIPERACILLIN/TAZOBACTAM', TARGETED)) / n(),    
                   FEP = sum(grepl('CEFEPIME', TARGETED)) / n(),  
                   CPT = sum(grepl('CEFTAROLINE', TARGETED)) / n(),  
                   .by = year) %>%
         arrange(year)
      ts <- x %>%
         filter(OXACILLIN == 0L) %>%
         mutate(year = substr(ORDER_DAY,1,4)) %>%
         summarise(NO_TRT = sum(lengths(TARGETED) == 0L) / n(),
                   VAN = sum(grepl('VANCOMYCIN', TARGETED)) / n(),
                   CFZ = sum(grepl('CEFAZOLIN', TARGETED)) / n(),
                   TZP = sum(grepl('PIPERACILLIN/TAZOBACTAM', TARGETED)) / n(),    
                   OXA = sum(grepl('OXACILLIN', TARGETED)) / n(),  
                   FEP = sum(grepl('CEFEPIME', TARGETED)) / n(),  
                   .by = year) %>%
         arrange(year)
      col_vec <- c('NO_TRT' = 'black', 'VAN' = 'blue', 'OXA' = 'red', 'CFZ' = 'darkorange', 'DAP' = 'darkgreen', 'TZP' = 'lightblue3', 'FEP' = 'tan1', 'CPT' = 'purple3')
      makePlot <- function(t, main, l) {
         plot(NA, yaxt='n', xaxs='i',
              xlim=c(2016.9, 2023.1), ylim=c(0, 0.62), 
              ylab='', xlab='Year', 
              main=paste0(main, ' (n = ', prettyNum(sum(l, na.rm=T), big.mark=','), ')'))
         abline(h = seq(0, 0.7, 0.1), lty=2, lwd=0.25)
         axis(side=2, at=seq(0, 0.7, 0.1), labels=paste0(seq(0, 0.7, 0.1) * 100, '%'), las=1)
         for (i in 2:length(t)) lines(x = t$year, y = t[[i]], col=col_vec[names(t)[i]], lwd=1.5)
         text(x = 2023.15, y = t[nrow(t),-1], adj=0, labels=names(t)[-1], xpd=NA, col=col_vec[names(t)[-1]])
      }
      {
         pdf(file = paste0(plots_path_name, 'Targeted/Staphylococcus_aureus_by_year.pdf'), width = 9, height=4)
         par(mfrow=c(1,2), mar=c(2.5,2.5,1.2,4), mgp=c(1.5, 0.3, 0), tck=-0.015, cex.main=1)
         makePlot(tr, 'MRSA', x$OXACILLIN == 1L)
         makePlot(ts, 'MSSA', x$OXACILLIN == 0L)
         dev.off()
      }
      rm(col_vec, tr, ts, makePlot)
   }
   
   # how many treatments?
   mrsa <- table(lengths(x$TARGETED[x$OXACILLIN == 1L]))
   mrsa <- mrsa / sum(mrsa)
   mssa <- table(lengths(x$TARGETED[x$OXACILLIN == 0L]))
   mssa <- mssa / sum(mssa)
   lens <- rbind(MRSA = mrsa[1:5], 
                 MSSA = mssa[1:5])
   par(mfrow=c(1,1))
   barplot(lens, beside=TRUE, legend=TRUE, xlab='No. of antibiotics administered')
   rm(x, lens, mrsa, mssa)
}







x <- empDF_ %>% filter(BUG == 'Escherichia coli') # 2,837
x <- x[sapply(x, function(x) sum(!is.na(x)) > 0L)]

barplot(table(lengths(x$TARGETED)))
x$TARG_ABBR[lengths(x$TARG_ABBR) == 0L] <- 'No abx'
x$TARG_ABBR_COMB <- sapply(x$TARG_ABBR, paste, collapse='+')
data.frame(head(sort(round(table(x$TARG_ABBR_COMB) / nrow(x) * 100, 1), decreasing = TRUE), n=12))
data.frame(head(sort(round(table(unlist(x$TARG_ABBR)) / nrow(x) * 100, 1), decreasing = TRUE), n=12))

sapply(x[18:length(x)], function(x) round(sum(x, na.rm=T) / length(x) * 100, 1))

# cross-resistances
x %>% select(CIPROFLOXACIN, LEVOFLOXACIN) %>% table()
x %>% select(AMPICILLIN, TETRACYCLINE) %>% table()
x %>% select(`TRIMETHOPRIM/SULFAMETHOXAZOLE`, TETRACYCLINE) %>% table()

t <- apply(x[18:length(x)], 1, sum, na.rm=TRUE)
table(t == 0) # half
table(t > 6)
barplot(table(t))


bugs <- c('Staphylococcus aureus', 'Enterococcus', 'Pseudomonas aeruginosa', 'Klebsiella pneumoniae', 'Escherichia coli', 'Streptococcus pyogenes')
par(mfrow=c(2,3))
for (b in seq_along(bugs)) {
   x <- empDF_ %>% filter(grepl(bugs[b], BUG))
   sapply(x[18:length(x)], function(x) sum(!is.na(x)))
   t <- apply(x[18:length(x)], 1, sum, na.rm=TRUE)
   barplot(table(t), main=paste0(bugs[b], ' (n = ', prettyNum(nrow(x), big.mark=','), ')'))
}
rm(b, x, t, bugs)









# CAN we find specific patterns of empiric --> targeted treatment?
# e.g., Staph aureus (VAN) --> (OXA) after learning its MSSA


# COOL!
empDF_ %>% filter(BUG == 'Staphylococcus aureus') %>%
   mutate(OXACILLIN = ifelse(is.na(OXACILLIN), 0, OXACILLIN)) %>%
   filter(grepl('VANCOMYCIN', EMPIRIC0) | grepl('VANCOMYCIN', EMPIRIC1)) %>% # most S aureus were treated "empirically" with VAN
   summarise(VAN_targeted = sum(grepl('VANCOMYCIN', TARGETED)),
             OXA_targeted = sum(grepl('OXACILLIN', TARGETED)),
             .by = OXACILLIN)



empDF_ %>% filter(grepl('Enterococcus', BUG)) %>%
   mutate(across(c(VANCOMYCIN, AMPICILLIN, PENICILLIN), ~ ifelse(is.na(.), 0, .))) %>%
   mutate(BETA_LACTAM = as.integer(AMPICILLIN == 1L | PENICILLIN == 1L)) %>%
   summarise(VAN_targ = sum(grepl('VANCOMYCIN', TARGETED)),
             DAP_targ = sum(grepl('DAPTOMYCIN', TARGETED)),
             LZD_targ = sum(grepl('LINEZOLID', TARGETED)),
             FEP_targ = sum(grepl('CEFEPIME', TARGETED)),
             CRO_targ = sum(grepl('CEFTRIAXONE', TARGETED)),
             TZP_targ = sum(grepl('PIPERACILLIN/TAZOBACTAM', TARGETED)),
             BL_targ = sum(grepl('PIPERACILLIN/TAZOBACTAM|CEFEPIME|CEFTRIAXONE|AZTREONAM|AMPICILLIN', TARGETED)),
             .by = c(BETA_LACTAM, VANCOMYCIN)) %>%
   arrange(VANCOMYCIN, BETA_LACTAM)


















