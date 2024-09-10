library(dplyr)
abbr <- tibble(read.table(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/antibiotic_names/ABX_ABBR.txt', header = TRUE)) %>% 
   select(-Class) %>%
   mutate(Antibiotic_Name = gsub('-', '/', Antibiotic_Name),
          Antibiotic_Name = gsub('_', ' ', Antibiotic_Name),
          Antibiotic_Name = toupper(Antibiotic_Name))
abbr <- setNames(abbr$Abbreviation, abbr$Antibiotic_Name)

load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/data_path_name.Rdata')
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/plots_path_name.Rdata')
plots_path_name <- paste0(plots_path_name, 'BloodAbxAdmin/')
#load(file = paste0(data_path_name, 'ASTs_AbxAdmin_blood_2017_2023_imputed.Rdata'))
load(file = paste0(data_path_name, 'ASTs_AbxAdmin_blood_2017_2023.Rdata'))
# empDF <- empDF %>% filter(!MULT_BLOOD_ISO)

empDF$EMPIRIC01 <- mapply(empDF$EMPIRIC0, empDF$EMPIRIC1, FUN=function(x, y) list(sort(unique(c(unlist(x), unlist(y))))))
empDF <- empDF %>% relocate(EMPIRIC01, .before=EMPIRIC0)
empDF$FLAG01 <- empDF$FLAG1 <- empDF$FLAG0 <- empDF$FLAGT <- NA

getFLAG <- function(check, vals) {
   # if any AST value is 0 (susceptible), the therapy is concordant
   if (any(vals == 0, na.rm = T)) return("CONCORDANT")
   
   # all are in table
   if (all(check)) {
      if (all(is.na(vals)))                   return('NOT TESTED')
      if (all(!is.na(vals)) & all(vals == 1)) return('DISCORDANT')  # all tested and all resistant
      if (any(vals == 1) & any(is.na(vals)))  return('DISCORDANT')      # >=1 not tested, >=1 resistant
      print('if Im here, something went wrong')
   }
   
   # all are not in table
   if (all(!check)) return('NOT TESTED')
   
   # NT + never
   r <- any(vals == 1, na.rm=TRUE)
   nt <- any(is.na(vals))
   nev <- any(!check)
   
   if (r & nt & !nev) return('DISCORDANT')
   if (!r & nt & nev) return('NOT TESTED')
   if (r & !nt & nev) return('DISCORDANT')
   if (r & nt & nev)  return('DISCORDANT')
}

start <- Sys.time()
for (i in 1:nrow(empDF)) {
   print(i)
   
   if (length(empDF$EMPIRIC01[[i]]) == 0){
      empDF$FLAG01[i] <- 'No empiric therapy given'
   } else {
      check <- empDF$EMPIRIC01[[i]] %in% names(empDF)  # indicates which abx given have corresponding rows in empDF 
      vals <- empDF[i, empDF$EMPIRIC01[[i]][check]]     # values associated with abx given
      empDF$FLAG01[i] <- getFLAG(check, vals)
   }

   if (length(empDF$EMPIRIC1[[i]]) == 0){
      empDF$FLAG1[i] <- 'No empiric therapy given'
   } else {
      check <- empDF$EMPIRIC1[[i]] %in% names(empDF)  # indicates which abx given have corresponding rows in empDF 
      vals <- empDF[i, empDF$EMPIRIC1[[i]][check]]     # values associated with abx given
      empDF$FLAG1[i] <- getFLAG(check, vals)
   }
      
   if (length(empDF$EMPIRIC0[[i]]) == 0){
      empDF$FLAG0[i] <- 'No empiric therapy given'
   } else {
      check <- empDF$EMPIRIC0[[i]] %in% names(empDF)  # indicates which abx given have corresponding rows in empDF 
      vals <- empDF[i, empDF$EMPIRIC0[[i]][check]]     # values associated with abx given
      empDF$FLAG0[i] <- getFLAG(check, vals)
   }
   
   if (length(empDF$TARGETED[[i]]) == 0) {
      empDF$FLAGT[i] <- 'No targeted therapy given'
   } else {
      check <- empDF$TARGETED[[i]] %in% names(empDF)  # indicates which abx given have corresponding rows in empDF 
      vals <- empDF[i, empDF$TARGETED[[i]][check]]     # values associated with abx given
      empDF$FLAGT[i] <- getFLAG(check, vals)
   }
}
print(Sys.time() - start) # 1.18 minutes
rm(check, i, vals, getFLAG, start)


empDF <- empDF %>% relocate(FLAG0, FLAG1, FLAG01, FLAGT, .before=EMPIRIC01)
save(empDF, file = '~/Desktop/EHR/EHR work/RdataFiles/ASTs_AbxAdmin_blood_2017_2023_flagged.Rdata')




###########################################################################################################
library(dplyr)
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ASTs_AbxAdmin_blood_2017_2023_flagged.Rdata')
empDFo <- empDF

load(file = '~/Desktop/EHR/EHR work/RdataFiles/ASTs_AbxAdmin_blood_2017_2023_imputed_flagged.Rdata')


# Merge flag0 for multiple blood isolate cases
empDF %>% count(FLAG0, sort=TRUE)
x <- empDF %>%
   # filter(MULT_BLOOD_ISO, FLAG0 != 'No empiric therapy given') %>%
   summarise(FLAG0 = paste(sort(unique(FLAG0)), collapse='+'),
             .by = c(PERSON_ID, ORDER_DAY)) %>%
   mutate(NEW_FLAG0 = case_when(
      FLAG0 == 'CONCORDANT+NOT TESTED' ~ 'CONCORDANT',
      FLAG0 == 'CONCORDANT+DISCORDANT' ~ 'DISCORDANT',
      FLAG0 == 'DISCORDANT+NOT TESTED' ~ 'DISCORDANT',
      FLAG0 == 'CONCORDANT+DISCORDANT+NOT TESTED' ~ 'DISCORDANT',
      .default = FLAG0
   ))
x %>% count(FLAG0, sort=TRUE)
x %>% count(NEW_FLAG0, sort=TRUE) # 10.1% of given antibiotics (7% of infections)
rm(x)
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


empDF %>% count(FLAG0, sort=TRUE) %>% mutate(r = n / nrow(empDF) * 100)
empDF %>% count(FLAG01, sort=TRUE) %>% mutate(r = n / nrow(empDF) * 100)
empDF %>% count(FLAG1, sort=TRUE) %>% mutate(r = n / nrow(empDF) * 100)
empDF %>% count(FLAGT, sort=TRUE) %>% mutate(r = n / nrow(empDF) * 100)


# Do mult-bug infections have higher rates of discordance?
# They have slightly higher rates of being treated either concordantly or discordantly
t <- empDF %>% select(MULT_BLOOD_ISO, FLAG0) %>% table()
t <- apply(t, 1, function(x) x / sum(x))
t
rm(t)

t <- table(empDF$FLAG0)
t <- t / sum(t) # 11.5% not tested (down from 25%)
t



### ONLY IF USING NON-IMPUTED DATA ###
# does the number of resistances have higher rates of discordance?
{
   # getDF2 <- function(df) {
   #    # t <- apply(df[24:length(df)], 1, function(x) paste(sum(x, na.rm=T), sum(!is.na(x), na.rm=T), sep='/'))
   # }
   getDF <- function(df, prop=FALSE) {
      if (prop) {
         t <- apply(df[24:length(df)], 1, function(x) sum(x, na.rm=T) / sum(!is.na(x), na.rm=T))
         t <- round(t, 1)
         #t <- cut(t, breaks=seq(0, 1, 0.05))
      } else { 
         t <- apply(df[24:length(df)], 1, function(x) sum(x, na.rm=T))
      }
      t <- table(t, empDF$FLAG0 == 'DISCORDANT')
      t <- t(t)
      t <- apply(t, 2, function(x) {
         p <- x['TRUE'] / sum(x)
         se <- 1.96 * sqrt((p * (1 - p)) / sum(x))
         names(p) <- NULL
         names(se) <- NULL
         return(c(est = p, se = se, n = sum(x)))
      })
      t <- data.frame(t(t))
      t <- t[t$n >= 10, ]
      return(t)
   }
   ta <- getDF(empDFo)
   tp <- getDF(empDFo, prop=TRUE)
   ti <- getDF(empDF)
   
   tp <- apply(df[24:length(df)], 1, function(x) sum(x, na.rm=T) / sum(!is.na(x), na.rm=T))
   tp[is.na(tp)] <- 0
   disc <- as.integer(empDF$FLAG0 == 'DISCORDANT')
   
   m <- glm(disc ~ tp, family='binomial')
   s <- summary(m)
   
   plot(x=tp, y=, pch=16, col='#000000aa')
   
   
   
   makePlot <- function(t, xlab='', ylab=FALSE, prop=FALSE) {
      if (prop) {
         xvals <- seq(0.05, 0.95, 0.05)
      } else {
         xvals <- as.numeric(rownames(t))
      }
      plot(x = xvals,
           y = t$est,
           pch = 16, cex = 0.75,
           ylim = c(0, 1), yaxt = 'n', xaxt='n',
           xlab = paste0('Resistant phenotypes ', xlab), ylab='')
      if (ylab) title(ylab = 'Proportion discordant', line=2)
      if (prop) {
         axis(side=1, at=seq(0.05, 0.95, 0.05))
      } else {
         axis(side=1)
      }
      axis(side=2, at=seq(0,1,0.2), las=1)
      arrows(x0 = xvals, y0 = t$est - t$se, y1 = t$est + t$se,
             code=3,
             length=0.075,
             angle=90)
   }
   {
      pdf(file = paste0(plots_path_name, 'DiscordRateVsNumResistances.pdf'), width=10, height=3)
      par(mfrow=c(1,3), mar=c(3,3,1,1), tck=-0.01, mgp=c(1.7, 0.5, 0))
      makePlot(ta, 'measured (absolute number)', ylab=TRUE)
      makePlot(tp, 'measured (proportion - binned)', prop=TRUE)
      makePlot(ti, '(imputed number)')
      dev.off()
   }
}



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
   summarise(num_disc = sum(FLAG1 == 'DISCORDANT'),
             total = sum(FLAG1 %in% c('DISCORDANT', 'CONCORDANT')),
             disc_rate = round(num_disc / total * 100, 1),
             .by = BUG) %>%
   filter(total > 100L) %>%
   arrange(desc(disc_rate)) %>%
   View()
# E. faecium ~ 40%

empDF %>%
   summarise(num_disc = sum(FLAGT == 'DISCORDANT'),
             total = sum(FLAGT %in% c('DISCORDANT', 'CONCORDANT')),
             disc_rate = round(num_disc / total * 100, 1),
             .by = BUG) %>%
   filter(total > 100L) %>%
   arrange(desc(disc_rate)) %>%
   View()
# E. faecium ~ 10.6%


x <- empDF %>% filter(BUG == 'Enterococcus faecium')
#x <- x %>% filter(FLAG0 == 'DISCORDANT')
data.frame(head(sort(table(sapply(x$EMPIRIC0, function(x) paste(abbr[x], collapse='+'))), decreasing = TRUE), n=8))
data.frame(head(sort(table(unlist(x$EMPIRIC0)), decreasing = TRUE), n=8))
x %>%
   filter(FLAG0 %in% c('CONCORDANT', 'DISCORDANT', 'probDISCORDANT'))




# VRE has discordant empiric rate of 71% vs. 23% for 
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
















