library(dplyr)
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/data_path_name.Rdata')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ASTs_AbxAdmin_blood_2017_2023_flagged.Rdata')
load(file = paste0(data_path_name, 'ASTs_AbxAdmin_blood_2017_2023_imputed_flagged.Rdata'))
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/antibiotic_names/abbr_named.Rdata')
abbr <- c(abbr, setNames('No abx', 'No abx'))






############################
library(dplyr)
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/data_path_name.Rdata')
load(file = paste0(data_path_name, 'ASTs_AbxAdmin_blood_2017_2023_imputed_flagged_LONG.Rdata'))
chours <- 3600

empDF_ <- empDF
# remove abx admin times that are missing time-stamp
# any missing time-stamp is still present, but the START_DATE is now NA
empDF_$START_DATE[substr(empDF_$START_DATE, 12, 12) == ''] <- NA

empDF <- empDF %>%
   mutate(ABX_PROX_ORDER = as.numeric(lubridate::as.duration(START_DATE - ORDER_DATE)) / 3600,
          pEMPIRIC = START_DATE < (ORDER_DATE - 48*chours),
          EMPIRIC = START_DATE >= (ORDER_DATE - 48*chours) & START_DATE < (ORDER_DATE + 12*chours),
          BETWEEN = START_DATE >= (ORDER_DATE + 12*chours) & START_DATE < RESULT_DATE,
          TARGETED = START_DATE >= RESULT_DATE & START_DATE <= (RESULT_DATE + 72*chours),
          lTARGETED = START_DATE > (RESULT_DATE + 72*chours)) %>%
   group_by_all() %>% 
   ungroup(ABX, START_DATE, END_DATE, pEMPIRIC, EMPIRIC, BETWEEN, TARGETED, lTARGETED, FLAG, ABX_PROX_ORDER) %>%
   reframe(TIME_TO_FIRST_ABX = min(ABX_PROX_ORDER),
           TIME_TO_CONC = list(ABX_PROX_ORDER[FLAG == 'CONCORDANT']),
           ABX_EMPp = list(sort(unique(ABX[pEMPIRIC]))),
           ABX_EMP = list(sort(unique(ABX[EMPIRIC]))),
           ABX_BTW = list(sort(unique(ABX[BETWEEN]))),
           ABX_TAR = list(sort(unique(ABX[TARGETED]))),
           ABX_TARl = list(sort(unique(ABX[lTARGETED]))),
           FLAGp = paste(sort(unique(FLAG[pEMPIRIC])), collapse='+'),
           FLAGe = paste(sort(unique(FLAG[EMPIRIC])), collapse='+'),
           FLAGb = paste(sort(unique(FLAG[BETWEEN])), collapse='+'),
           FLAGt = paste(sort(unique(FLAG[TARGETED])), collapse='+'),
           FLAGl = paste(sort(unique(FLAG[lTARGETED])), collapse='+')) %>%
   ungroup() %>%
   relocate(ABX_EMPp:ABX_TARl, FLAGp:FLAGl, TIME_TO_FIRST_ABX, TIME_TO_CONC, .before=BUG)
empDF_ <- empDF_ %>%
   mutate(ABX_PROX_ORDER = as.numeric(lubridate::as.duration(START_DATE - ORDER_DATE)) / 3600,
          pEMPIRIC = START_DATE < (ORDER_DATE - 48*chours),
          EMPIRIC = START_DATE >= (ORDER_DATE - 48*chours) & START_DATE < (ORDER_DATE + 12*chours),
          BETWEEN = START_DATE >= (ORDER_DATE + 12*chours) & START_DATE < RESULT_DATE,
          TARGETED = START_DATE >= RESULT_DATE & START_DATE <= (RESULT_DATE + 72*chours),
          lTARGETED = START_DATE > (RESULT_DATE + 72*chours)) %>%
   group_by_all() %>% 
   ungroup(ABX, START_DATE, END_DATE, pEMPIRIC, EMPIRIC, BETWEEN, TARGETED, lTARGETED, FLAG, ABX_PROX_ORDER) %>%
   reframe(TIME_TO_FIRST_ABX = min(ABX_PROX_ORDER),
           TIME_TO_CONC = list(ABX_PROX_ORDER[FLAG == 'CONCORDANT']),
           ABX_EMPp = list(sort(unique(ABX[pEMPIRIC]))),
           ABX_EMP = list(sort(unique(ABX[EMPIRIC]))),
           ABX_BTW = list(sort(unique(ABX[BETWEEN]))),
           ABX_TAR = list(sort(unique(ABX[TARGETED]))),
           ABX_TARl = list(sort(unique(ABX[lTARGETED]))),
           FLAGp = paste(sort(unique(FLAG[pEMPIRIC])), collapse='+'),
           FLAGe = paste(sort(unique(FLAG[EMPIRIC])), collapse='+'),
           FLAGb = paste(sort(unique(FLAG[BETWEEN])), collapse='+'),
           FLAGt = paste(sort(unique(FLAG[TARGETED])), collapse='+'),
           FLAGl = paste(sort(unique(FLAG[lTARGETED])), collapse='+')) %>%
   ungroup() %>%
   relocate(ABX_EMPp:ABX_TARl, FLAGp:FLAGl, TIME_TO_FIRST_ABX, TIME_TO_CONC, .before=BUG)

# is an antibiotic ever ONLY logged with missing time stamp?
# compare ABX_EMP before and after chopping time stamps...
df <- cbind(empDF %>% select(PERSON_ID, ORDER_DATE, BUG, ABX_EMP, ABX_BTW, ABX_TAR) %>% rename(OG_ABX_EMP = ABX_EMP, OG_ABX_BTW = ABX_BTW, OG_ABX_TAR = ABX_TAR),
            empDF_ %>% select(ABX_EMP, ABX_BTW, ABX_TAR) %>% rename(NW_ABX_EMP = ABX_EMP, NW_ABX_BTW = ABX_BTW, NW_ABX_TAR = ABX_TAR))
df <- tibble(df)

sum(lengths(df$NW_ABX_EMP) == 0) - sum(lengths(df$OG_ABX_EMP) == 0) # 1,297
sum(lengths(df$OG_ABX_EMP) > lengths(df$NW_ABX_EMP))  # 4,453

sum(lengths(df$NW_ABX_BTW) == 0) - sum(lengths(df$OG_ABX_BTW) == 0) # 66
sum(lengths(df$OG_ABX_BTW) > lengths(df$NW_ABX_BTW))  # 3,987

sum(lengths(df$NW_ABX_TAR) == 0) - sum(lengths(df$OG_ABX_TAR) == 0) # 239
sum(lengths(df$OG_ABX_TAR) > lengths(df$NW_ABX_TAR))  # 2,197
####################








empDF %>%
   summarise(num_disc = sum(FLAGe == 'DISCORDANT'),
             total = sum(FLAGe != 'No empiric therapy given'),
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





#### STAPHYLOCOCCUS AUREUS treatment changes over time  + number of treatments are equal ####
{
   x <- empDF %>% filter(BUG == 'Staphylococcus aureus') # 2,493
   
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
         par(mfrow=c(1,2), mar=c(2.5,2.5,1.2,4), mgp=c(1.5, 0.3, 0), tck=-0.015, cex.main=1)
         makePlot(tr, 'MRSA', x$OXACILLIN == 1L)
         makePlot(ts, 'MSSA', x$OXACILLIN == 0L)
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


# CAN we find specific patterns of empiric --> targeted treatment?
# e.g., Staph aureus (VAN) --> (OXA) after learning its MSSA


# COOL!
empDF %>% filter(BUG == 'Staphylococcus aureus') %>%
   filter(grepl('VANCOMYCIN', EMPIRIC) | grepl('VANCOMYCIN', TARGETED)) %>% # most S aureus were treated "empirically" with VAN
   summarise(VAN_targeted = sum(grepl('VANCOMYCIN', TARGETED)),
             OXA_targeted = sum(grepl('OXACILLIN', TARGETED)),
             .by = OXACILLIN)



















