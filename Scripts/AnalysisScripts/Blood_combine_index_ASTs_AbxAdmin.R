library(dplyr)
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/data_path_name.Rdata')
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/plots_path_name.Rdata')
plots_path_name <- paste0(plots_path_name, 'BloodAbxAdmin/')
load(file = paste0(data_path_name, 'ASTs_blood_2015_2023.Rdata'))
load(file = paste0(data_path_name, 'ALL_CLEANED_2017_2023_AbxAdmin.Rdata'))

# only have medication administration data from 2017 - 2023
astDF <- empDF %>% filter(substr(ORDER_DATE, 1, 4) %in% as.character(2017:2023)) #    55,959 -->    44,376
abxDF <- abxDF %>% filter(PERSON_ID %in% unique(astDF$PERSON_ID))                # 9,681,748 --> 2,588,417
rm(empDF)

length(unique(astDF$PERSON_ID)) # 35,672
length(unique(abxDF$PERSON_ID)) # 33,642
length(intersect(astDF$PERSON_ID, abxDF$PERSON_ID)) # 33,642


# add in a column describing the date of the most recent antibiotic prescription
abxDF <- abxDF %>%
   mutate(PERSON_ID = as.character(PERSON_ID)) %>%
   rename(START_DATE = ADMIN_START_DATE,
          END_DATE = ADMIN_END_DATE) %>%
   arrange(PERSON_ID, START_DATE)# %>%
# mutate(DATE_OF_MOST_RECENT_ABX = case_when(
#    START_DATE > lag(END_DATE) & PERSON_ID == lag(PERSON_ID) ~ lag(END_DATE)
# ))

# mark but remove "visits" with >1 index isolate
# astDF <- astDF %>%
#    mutate(CONTAM = case_when(
#       NEXT_RESULT_DAY <= RESULT_DAY ~ 2,
#       NEXT_ORDER_DAY < RESULT_DAY ~ 1
#    )) %>%
#    select(-NEXT_ORDER_DAY, -NEXT_RESULT_DAY) %>%
#    relocate(CONTAM, .before=BUG)

# astDF %>%
#    filter(grepl('Enterococcus', BUG)) %>%
#    mutate(BUG2 = case_when(
#       BUG == 'Staphylococcus aureus' & OXACILLIN == 1 ~ 'MRSA',
#       grepl('Enterococcus', BUG) & VANCOMYCIN == 1 ~ 'VRE',
#       grepl('Enterococcus', BUG) & VANCOMYCIN == 0 ~ 'Ent',
#    )) %>%
#    select(BUG2, AMPICILLIN) %>% table()

astDF <- astDF %>%
   select(!CEFEPIME:DELAFLOXACIN) %>%
   summarise(ORDER_DATE = min(ORDER_DATE),
             RESULT_DATE = min(RESULT_DATE),
             BUG = list(sort(unique(BUG))),
             NEXT_ORDER_DAY = min(NEXT_ORDER_DAY),
             NEXT_RESULT_DAY = min(NEXT_RESULT_DAY),
             .by = c(PERSON_ID, ORDER_DAY, MULT_ISO)) %>%
   mutate(ORDER_DAY = as.Date(substr(ORDER_DATE,1,10)),
          RESULT_DAY = as.Date(substr(RESULT_DATE,1,10))) %>%
   mutate(CONTAM = case_when(
      NEXT_RESULT_DAY <= RESULT_DAY ~ 2,
      NEXT_ORDER_DAY < RESULT_DAY ~ 1
   )) %>%
   mutate(MULT_BLOOD_ISO = lengths(BUG) > 1L) %>%
   select(-NEXT_ORDER_DAY, -NEXT_RESULT_DAY) %>%
   relocate(RESULT_DAY, .after=ORDER_DAY) %>%
   relocate(ORDER_DATE, RESULT_DATE, .before=ORDER_DAY) %>%
   relocate(BUG, .after=CONTAM)
# group_by(PERSON_ID, ORDER_DAY) %>%
# mutate(MULT_BLOOD_ISO = n() > 1L) %>%
# slice_min(RESULT_DATE, with_ties = FALSE) %>%
# ungroup() %>%
# relocate(MULT_BLOOD_ISO, .before=BUG) # ~1/2 of the multiple isolate infections had multiple blood isolates


## DIFFERENCES BETWEEN ORDER_TIMES for multiple isolates
# d <- astDF %>%
#    group_by(PERSON_ID, ORDER_DAY) %>%
#    filter(n() > 1L) %>%
#    filter(length(unique(ORDER_DATE)) > 1) %>%
#    mutate(D = as.numeric(diff(range(ORDER_DATE))) / 60) %>%
#    ungroup() %>%
#    select(D) %>%
#    unlist()
# d <- unname(d)
# hist(d)
# rm(d)
# 
# astDF1 <- astDF %>%
#    summarise(ORDER_DATE = min(ORDER_DATE),
#              RESULT_DATE = min(RESULT_DATE),
#              RESULT_DAY = min(RESULT_DAY),
#              BUG = list(sort(unique(BUG))),
#              .by = c(PERSON_ID, ORDER_DAY))
# barplot(table(lengths(astDF1$BUG)), main='Number of pathogens isolated')
# 
# astDF2 <- astDF %>%
#    group_by(PERSON_ID, ORDER_DAY) %>%
#    slice_min(RESULT_DATE, with_ties = FALSE) %>%
#    ungroup()
# 
# sapply(names(astDF1), function(x) all(astDF1[[x]] == astDF2[[x]]))
# w <- which(astDF1$ORDER_DATE != astDF2$ORDER_DATE)
# # table(astDF1$NEXT_ORDER_DAY != astDF2$NEXT_ORDER_DAY)
# # table(astDF1$NEXT_RESULT_DAY != astDF2$NEXT_RESULT_DAY)
# length(w) # 468 (out of 37,983)
# # how often is the min(ORDER_DATE) earlier than the one that we get when choosing the min(RESULT_DATE)
# # answer: EVERY SINGLE TIME
# table(astDF1$ORDER_DATE[w] < astDF2$ORDER_DATE[w])
# rm(astDF2, w)



# JOIN
empDF <- astDF %>% # 37,983 "infections"
   mutate(JOIN_START = ORDER_DAY - 30,
          JOIN_END = RESULT_DAY + 30) %>%
   left_join(x = .,
             y = abxDF,
             by = join_by(PERSON_ID,
                          JOIN_START <= START_DATE,
                          JOIN_END >= START_DATE)) %>% # 1,187,737
   # if the most recent Rx was before empiric window, get num days
   # mutate(DAYS_SINCE_PRV_ABX = ifelse(test = DATE_OF_MOST_RECENT_ABX < ORDER_DATE, 
   #                                    yes = as.numeric(lubridate::as.duration(ORDER_DATE - DATE_OF_MOST_RECENT_ABX))/86400, 
   #                                    no = NA)) %>%
   select(-JOIN_START, -JOIN_END) %>%
   group_by(PERSON_ID, ORDER_DAY) %>% # 37,983
   relocate(ABX, START_DATE, END_DATE, .before=BUG) %>%
   ungroup()


# determine the first 
df <- empDF %>%
   select(PERSON_ID, ORDER_DAY, RESULT_DAY, START_DATE, BUG) %>%
   distinct() %>%
   group_by(PERSON_ID, ORDER_DAY) %>%
   slice_min(START_DATE, with_ties=FALSE) %>%
   ungroup() %>%
   mutate(START_DAY = as.Date(substr(START_DATE,1,10))) %>%
   mutate(ABX_PROX_ORDER = as.integer(START_DAY - ORDER_DAY))

sum(is.na(df$START_DATE)) / nrow(df) # 10.9% had no abx in the 30 days beofre or after ORDER-->RESULT
sum(df$ABX_PROX_ORDER < 0, na.rm=T) / nrow(df)  # 18% had a prescription in the 30 days leading up to the ORDER_DAY


### PLOT PROPORTION OF TOP BUGS THAT WERE UNTREATED
{
   # bugs <- list(Treated = unlist(df$BUG[!is.na(df$START_DATE)]),
   #              Untreated = unlist(df$BUG[is.na(df$START_DATE)]))
   bugs <- tapply(df$BUG, is.na(df$START_DATE), unlist)
   names(bugs) <- c('Treated', 'Untreated')
   lengths(bugs) # 36,564 treated, 4,458 untreated
   x <- sort(table(bugs$Treated), decreasing=TRUE)
   t <- data.frame(unname(x), row.names=names(x))
   t <- t[names(t) == 'Freq']
   names(t) <- 'Treated'
   t$Untreated <- sort(table(bugs$Untreated), decreasing=TRUE)[rownames(t)]
   t$Untreated[is.na(t$Untreated)] <- 0
   t$Total <- t$Treated + t$Untreated
   t$FracUntr <- t$Untreated / t$Total
   t <- t[t$Total >= 100, ]
   t <- t[order(t$FracUntr, decreasing=TRUE), ]
   t <- t[t$FracUntr > 0.07, ]
   t <- t[which(t[, 'Untreated'] > 0.07), ]
   w <- which(grepl('^([A-Z])[a-z]+ ([a-z]+)$', rownames(t)) & !grepl('.+ species$', rownames(t)))
   rownames(t)[w] <- gsub('^([A-Z])[a-z]+ ([a-z]+)$', '\\1. \\2', rownames(t)[w])
   rownames(t)[rownames(t) == 'Coagulase Negative Staph'] <- 'Coag. Neg. Staph'
   rownames(t)[rownames(t) == 'Viridans Streptococci'] <- 'Viridans Strep'
   rownames(t)[rownames(t) == 'Streptococcus mitis/oralis group'] <- 'S. mitis/oralis group'
   rownames(t)[rownames(t) == 'Did not match'] <- 'unknown/rare pathogen'
   t <- t(as.matrix(t))
   b <- barplot(t[1:2,], horiz=TRUE, plot=FALSE)
   {
      pdf(file = paste0(plots_path_name, 'BugsTreatedVsUntreated.pdf'), height=7, width=8.5)
      par(mar=c(3, 12.1, 3.25, 4))
      barplot(t['FracUntr',], xlim=c(0,1), col='#555555', horiz=TRUE, names.arg=rep('', ncol(t)), main='Proportion of bloodstream infections where no antibiotic given')
      text(x=-0.025, y=b, adj=1, labels=colnames(t), xpd=NA)
      text(x=1, y=b, adj=1, labels=prettyNum(t['Untreated',], big.mark=','), xpd=NA)
      text(x=1, y=b[length(b)] + unique(round(diff(b),1)), adj=1, labels='No. infections', font=2, xpd=NA)
      text(x=t['FracUntr',]+0.005, y=b, labels=round(t['FracUntr',], 2), adj=0, xpd=NA)
      dev.off()
   }
   rm(b, t, x, bugs, w)
}

# WHAT WERE THE BUGS FOR THOSE THAT WERE NOT GIVEN ABX UNTIL RESULTS WERE KNOWN?
{
   bugs <- tapply(df$BUG, df$START_DAY >= df$RESULT_DAY, unlist)
   lengths(bugs) # 420 true
   names(bugs) <- c('Emp', 'Targ')
   x <- sort(table(bugs$Emp), decreasing=TRUE)
   t <- data.frame(unname(x), row.names=names(x))
   t <- t[names(t) == 'Freq']
   names(t) <- 'Emp'
   t$Targ <- sort(table(bugs$Targ), decreasing=TRUE)[rownames(t)]
   t$Targ[is.na(t$Targ)] <- 0
   t$Total <- t$Emp + t$Targ
   t$FracUntr <- t$Targ / t$Total
   t <- t[t$Targ > 4, ]
   t <- t[order(t$FracUntr, decreasing=TRUE), ]
   w <- which(grepl('^([A-Z])[a-z]+ ([a-z]+)$', rownames(t)) & !grepl('.+ species$', rownames(t)))
   rownames(t)[w] <- gsub('^([A-Z])[a-z]+ ([a-z]+)$', '\\1. \\2', rownames(t)[w])
   rownames(t)[rownames(t) == 'Coagulase Negative Staph'] <- 'Coag. Neg. Staph'
   #rownames(t)[rownames(t) == 'Viridans Streptococci'] <- 'Viridans Strep'
   #rownames(t)[rownames(t) == 'Streptococcus mitis/oralis group'] <- 'S. mitis/oralis group'
   #rownames(t)[rownames(t) == 'Did not match'] <- 'unknown/rare pathogen'
   t <- t(as.matrix(t))
   b <- barplot(t[1:2,], horiz=TRUE, plot=FALSE)
   {
      pdf(file = paste0(plots_path_name, 'BugsNotTreatedUntilTargeted.pdf'), height=5, width=7)
      par(mar=c(3, 12.1, 3.25, 4))
      barplot(t['FracUntr',], xlim=c(0,1), col='#555555', horiz=TRUE, names.arg=rep('', ncol(t)), main='Proportion of bloodstream infections with no empiric therapy')
      text(x=-0.025, y=b, adj=1, labels=colnames(t), xpd=NA)
      text(x=1, y=b, adj=1, labels=prettyNum(t['Targ',], big.mark=','), xpd=NA)
      text(x=1, y=b[length(b)] + unique(round(diff(b),1)), adj=1, labels='No. infections', font=2, xpd=NA)
      text(x=t['FracUntr',]+0.005, y=b, labels=round(t['FracUntr',], 2), adj=0, xpd=NA)
      dev.off()
   }
   rm(b, t, x, bugs, w)
}

# WHAT WERE THE BUGS FOR THOSE THAT WERE GIVEN ABX BEFORE ORDER DATE?
{
   bugs <- tapply(df$BUG, df$START_DAY < df$ORDER_DAY, unlist)
   lengths(bugs) # 7,424 true
   names(bugs) <- c('Emp', 'Prior')
   x <- sort(table(bugs$Emp), decreasing=TRUE)
   t <- data.frame(unname(x), row.names=names(x))
   t <- t[names(t) == 'Freq']
   names(t) <- 'Emp'
   t$Prior <- sort(table(bugs$Prior), decreasing=TRUE)[rownames(t)]
   t$Prior[is.na(t$Prior)] <- 0
   t$Total <- t$Emp + t$Prior
   t$FracPrior <- t$Prior / t$Total
   t <- t[order(t$FracPrior, decreasing=TRUE), ]
   t <- t[t$Prior > 4, ]
   t <- t[t$Total >= 100, ]
   t <- t[t$FracPrior >= 0.2, ]
   w <- which(grepl('^([A-Z])[a-z]+ ([a-z]+)$', rownames(t)) & !grepl('.+ species$', rownames(t)))
   rownames(t)[w] <- gsub('^([A-Z])[a-z]+ ([a-z]+)$', '\\1. \\2', rownames(t)[w])
   #rownames(t)[rownames(t) == 'Coagulase Negative Staph'] <- 'Coag. Neg. Staph'
   #rownames(t)[rownames(t) == 'Viridans Streptococci'] <- 'Viridans Strep'
   #rownames(t)[rownames(t) == 'Streptococcus mitis/oralis group'] <- 'S. mitis/oralis group'
   rownames(t)[rownames(t) == 'Did not match'] <- 'unknown/rare pathogen'
   t <- t(as.matrix(t))
   b <- barplot(t[1:2,], horiz=TRUE, plot=FALSE)
   {
      pdf(file = paste0(plots_path_name, 'BugsTreatedPriorToEmpiric.pdf'), height=5.5, width=8.5)
      par(mar=c(3, 12.1, 3.25, 4))
      barplot(t['FracPrior',], xlim=c(0,1), col='#555555', horiz=TRUE, names.arg=rep('', ncol(t)), main='Proportion of bloodstream infections with pre-empiric treatment')
      text(x=-0.025, y=b, adj=1, labels=colnames(t), xpd=NA)
      text(x=1, y=b, adj=1, labels=prettyNum(t['Prior',], big.mark=','), xpd=NA)
      text(x=1, y=b[length(b)] + unique(round(diff(b),1)), adj=1, labels='No. infections', font=2, xpd=NA)
      text(x=t['FracPrior',]+0.005, y=b, labels=round(t['FracPrior',], 2), adj=0, xpd=NA)
      dev.off()
   }
   rm(b, t, x, bugs, w)
}

# WHAT WERE THE AST ORDER NAMES FOR THOSE NOT TREATED?
{
   # load(file = '~/Desktop/EHR/EHR work/RdataFiles/lab_micro_result_all_vw.Rdata')
   # astoDF <- astoDF %>%
   #    filter(PERSON_ID %in% unique(df$PERSON_ID[is.na(df$START_DATE)])) %>%
   #    select(PERSON_ID, ORDER_NAME, ORDER_DATE, RESULT_DATE, RESULT_VALUE, SPECIMEN_TYPE) %>%
   #    mutate(across(contains('DATE'), ~ strptime(., format='%m/%d/%Y %T')))
   # x <- df %>% 
   #    filter(is.na(START_DATE)) %>% # 4,131
   #    left_join(y = astoDF,
   #              by = join_by(PERSON_ID, ORDER_DATE, RESULT_DATE)) # 35,997
   # x %>% filter(RESULT_VALUE == 'NOT DETECTED')
}


df <- df %>% filter(!is.na(START_DATE), ABX_PROX_ORDER >= 0) # 27,005
sum(df$START_DAY >= df$RESULT_DAY, na.rm=T) / nrow(df) # 1.6%
sum(df$ABX_PROX_ORDER == 0, na.rm=T) / nrow(df) # 78.2% had first abx on order day
sum(df$ABX_PROX_ORDER == 1, na.rm=T) / nrow(df) # 17.4% had first abx on next day
sum(df$ABX_PROX_ORDER == 2, na.rm=T) / nrow(df) #  2.3% had first abx on 3rd day


# PLOT WHEN ANTIBIOTICS ARE ADMINISTERED RELATIVE TO AST ORDER DAY
{
   t <- table(df$ABX_PROX_ORDER)
   any(diff(as.integer(names(t))) > 1)
   t <- t / sum(t)
   t <- c(t[1:4], sum(t[5:length(t)]))
   names(t)[length(t)] <- '>3'
   t <- round(t*100, 1)
   b <- barplot(t, plot=FALSE)
   {
      pdf(file = paste0(plots_path_name, 'DateOfFirstAbx.pdf'), width=6, height=5)
      par(tck=-0.015, mgp=c(1.75, 0.5, 0), mar=c(3,1,2,1))
      barplot(t, ylim=c(0,85), yaxt='n', names.arg=NA, main='When was the first antibiotic administered?', xlab='Days since AST order')
      text(x = b, y = -3, labels=names(t), xpd=NA)
      text(x = b, y = t + 2.5, labels=paste0(t, '%'))
      dev.off()
   }
   rm(t, b)
}

# PLOT WHICH ANTIBIOTICS ARE GIVEN ON A TREATED INFECTION BASIS
{
   num_infections <- nrow(empDF %>% filter(!is.na(START_DATE)) %>% group_by(PERSON_ID, ORDER_DATE) %>% group_keys())
   t <- empDF %>% 
      select(PERSON_ID, ORDER_DAY, ABX) %>%
      distinct() %>%
      select(ABX) %>% table()
   t <- head(sort(t, decreasing=TRUE), n=12)
   t <- t / num_infections
   x <- barplot(t, horiz=TRUE, plot=F)
   {
      pdf(file = paste0(plots_path_name, 'AbxCounts.pdf'), width=6, height=5)
      par(mar=c(2, 10, 1, 1))
      barplot(t, horiz=TRUE, names.arg=NA, xlim=c(0,1), xaxt='n', main='Percentage of infections treated')
      text(x = -0.01, y = x, adj=1, labels=stringr::str_to_sentence(names(t)), xpd=NA)
      text(x = t+0.01, y=x, adj=0, labels=paste0(round(t,3)*100, '%'))
      dev.off()
   }
   rm(t, x, num_infections)
}


# still contains rows where antibiotics were prescribed leading up to empiric window
df <- empDF %>%
   filter(!is.na(START_DATE)) %>%
   filter(START_DATE >= ORDER_DATE) %>%
   mutate(across(contains('DATE'), ~ as.Date(substr(.,1,10)))) %>%
   select(PERSON_ID, ORDER_DATE, RESULT_DATE, START_DATE) %>%
   distinct() %>%
   mutate(X = as.integer(START_DATE - ORDER_DATE),
          D = as.integer(RESULT_DATE - ORDER_DATE))

# PLOT NUMBER OF PATIENTS RECEIVING ANTIBIOTICS ON EACH DAY - BROKEN DOWN BY RESULT DELAY
{
   df %>% count(D, sort=TRUE)
   # 3, 6, 2, 4, 5, 7
   makeBarplot <- function(main='') {
      x <- df %>% filter(D == as.integer(main))
      x <- x %>% mutate(RX=1) %>% 
         tidyr::pivot_wider(id_cols = c(PERSON_ID, ORDER_DATE), names_from=X, values_from=RX, values_fill=0) %>% 
         select(`0`, `1`, `2`, `3`, `4`, `5`, `6`, `7`, `8`, `9`, `10`)
      num_infections <- nrow(x)
      t <- colSums(x) / num_infections
      x <- x %>%
         mutate(`1` = ifelse(`0` == 1, 1, `1`)) %>%
         mutate(`2` = ifelse(`1` == 1, 1, `2`)) %>%
         mutate(`3` = ifelse(`2` == 1, 1, `3`)) %>%
         mutate(`4` = ifelse(`3` == 1, 1, `4`)) %>%
         mutate(`5` = ifelse(`4` == 1, 1, `5`)) %>%
         mutate(`6` = ifelse(`5` == 1, 1, `6`)) %>%
         mutate(`7` = ifelse(`6` == 1, 1, `7`)) %>%
         mutate(`8` = ifelse(`7` == 1, 1, `8`)) %>%
         mutate(`9` = ifelse(`8` == 1, 1, `9`)) %>%
         mutate(`10` = ifelse(`9` == 1, 1, `10`))
      tt <- colSums(x) / num_infections
      b <- barplot(t, plot=F)
      barplot(t, xlab='Days since blood culture order', ylim=c(0, 1), ylab='% patients receiving abx',
              main=paste0(main, ' day result delay\n(n = ', prettyNum(num_infections,big.mark=','), ')'))
      points(x = b, y = tt, pch=16, xpd=NA)
      lines(x = b, y = tt, xpd=NA)
      abline(v = b[names(t) == main], lwd=1.1, lty=2)
      text(x = b[names(t) == main]+0.1, y=0.95, adj=0, labels='result day', font=3)
   }
   {
      pdf(file = paste0(plots_path_name, 'TimingOfAbxAdmin.pdf'), width=9, height=6)
      par(mfrow=c(2,3), mgp=c(1.6, 0.4, 0), mar=c(3,3,3.5,1), tck=-0.015)
      makeBarplot(main='2')
      makeBarplot(main='3')
      makeBarplot(main='4')
      makeBarplot(main='5')
      makeBarplot(main='6')
      makeBarplot(main='7')
      dev.off()
   }
   rm(makeBarplot)
}

rm(df, abxDF)




empDF <- empDF %>%
   mutate(START_DAY = as.Date(substr(START_DATE,1,10))) %>%
   #select(-ORDER_DATE, -RESULT_DATE, -BLOOD) %>%
   reframe(EMPIRIC0 = list(sort(unique(ABX[START_DAY == ORDER_DAY]))),
           EMPIRIC1 = list(sort(unique(ABX[START_DAY == (ORDER_DAY+1)]))),
           EMPIRIC2 = list(sort(unique(ABX[START_DAY == (ORDER_DAY+2) & START_DAY < RESULT_DAY]))),
           EMPIRIC3 = list(sort(unique(ABX[START_DAY == (ORDER_DAY+3) & START_DAY < RESULT_DAY]))),
           EMPIRIC4 = list(sort(unique(ABX[START_DAY == (ORDER_DAY+4) & START_DAY < RESULT_DAY]))),
           EMPIRIC5 = list(sort(unique(ABX[START_DAY == (ORDER_DAY+5) & START_DAY < RESULT_DAY]))),
           TARGETED = list(sort(unique(ABX[START_DAY == RESULT_DAY]))),
           .by = c(PERSON_ID, ORDER_DAY, RESULT_DAY, BUG, MULT_BLOOD_ISO, CONTAM, MULT_ISO)) %>% #, INDEX_RECORD, CEFEPIME:DELAFLOXACIN)) %>%
   mutate(DELAY = as.integer(RESULT_DAY - ORDER_DAY)) %>%
   relocate(DELAY, EMPIRIC0, EMPIRIC1, EMPIRIC2, EMPIRIC3, EMPIRIC4, EMPIRIC5, TARGETED, .before=BUG)


save(empDF, file = paste0(data_path_name, 'ASTs_AbxAdmin_blood_2017_2023.Rdata'))



# do multi-bug infections get treated with more stuff? empirically? targeted?
{
   t <- empDF %>%
      mutate(num_bugs = lengths(BUG)) %>%
      filter(num_bugs < 6L) %>%
      summarise(n(),
                meanT = mean(lengths(TARGETED)),
                meanE = mean(lengths(EMPIRIC0)),
                .by = num_bugs)
   
   {
      pdf(file = paste0(plots_path_name, 'NumAbxVsNumBugs.pdf'))
      plot(t$meanT, pch=16, type='b', ylim=c(1, 2), xlab='No. of infecting bugs', ylab='Mean no. of antibiotic administered')
      points(t$meanE, pch=16, type='b', col='gray')
      legend('topright', legend=c('targeted', 'empiric'), col=c('black', 'gray'), pch=16, lty=1)
      dev.off()
   }
   rm(t)
}

# empiric and targeted course sizes
{
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
      pdf(file = paste0(plots_path_name, 'CourseSize.pdf'), height=6.5, width=9)
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
}



# DID THERAPIES NARROW IN THE TARGETED WINDOW?
sum(lengths(empDF$TARGETED) > 0) / nrow(empDF) # 64.6%
sum(lengths(empDF$TARGETED) < lengths(empDF$EMPIRIC0)) / nrow(empDF) # 35.7% narrowed
sum(lengths(empDF$TARGETED) > lengths(empDF$EMPIRIC0)) / nrow(empDF) # 25.5% broadened
sum(lengths(empDF$TARGETED) == lengths(empDF$EMPIRIC0)) / nrow(empDF) # 38.8% stayed same



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
      pdf(file = paste0(plots_path_name, 'TargetedByBug_.pdf'), height=6.5, width=12)
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
empDF %>% count(substr(ORDER_DAY,1,4)) # ~5,400 per year
empDF %>% filter(!MULT_BLOOD_ISO, !MULT_ISO) %>% count(substr(ORDER_DAY,1,4)) %>% summarise(sum(n)) # ~3,300 per year

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








