library(dplyr)
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/data_path_name.Rdata')
load(file = paste0(data_path_name, 'ASTs_AbxAdmin_blood_2017_imputed.Rdata'))
# load(file = paste0(data_path_name, 'ASTs_AbxAdmin_blood_2015_2023.Rdata'))

empDF <- empDF %>% filter(!MULT_ISOLATES)

empDF$EMPIRIC01 <- mapply(empDF$EMPIRIC0, empDF$EMPIRIC1, FUN=function(x, y) list(sort(unique(c(unlist(x), unlist(y))))))
empDF <- empDF %>% relocate(EMPIRIC01, .before=EMPIRIC0)
empDF$FLAG01 <- empDF$FLAG1 <- empDF$FLAG0 <- empDF$TFLAG <- NA

getFLAG <- function(check, vals) {
   # if any AST value is 0 (susceptible), the therapy is concordant
   if (any(vals == 0, na.rm = T)) return("CONCORDANT")
   
   # all are in table
   if (all(check)) {
      if (all(is.na(vals)))                   return('NOT TESTED')
      if (all(!is.na(vals)) & all(vals == 1)) return('DISCORDANT')  # all tested and all resistant
      if (any(vals == 1) & any(is.na(vals)))  return('R + NT')      # >=1 not tested, >=1 resistant
      print('if Im here, something went wrong')
   }
   
   # all are not in table
   if (all(!check)) return('NEVER TESTED')
   
   # NT + never
   r <- any(vals == 1, na.rm=TRUE)
   nt <- any(is.na(vals))
   nev <- any(!check)
   
   if (r & nt & !nev) return('R + NT')
   if (!r & nt & nev) return('NT + nev')
   if (r & !nt & nev) return('R + nev')
   if (r & nt & nev)  return('R + NT + nev')
}

for (i in 1:nrow(empDF)){
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
      empDF$TFLAG[i] <- 'No targeted therapy given'
   } else {
      check <- empDF$TARGETED[[i]] %in% names(empDF)  # indicates which abx given have corresponding rows in empDF 
      vals <- empDF[i, empDF$TARGETED[[i]][check]]     # values associated with abx given
      empDF$TFLAG[i] <- getFLAG(check, vals)
   }
}
rm(check, i, vals, getFLAG)
empDF <- empDF %>%
   relocate(FLAG0, FLAG1, FLAG01, TFLAG, .before=EMPIRIC01)

save(empDF, file = '~/Desktop/EHR/EHR work/RdataFiles/ASTs_AbxAdmin_blood_2017_imputed_flagged.Rdata')

# send to Megan
# empDF <- empDF %>% select(!c(PERSON_ID, ORDER_DATE, RESULT_DATE, MULT_ISOLATES))
# save(empDF, file = '~/Desktop/EHR/EHR work/MegansProject/deID_2017_ASTsAbxAdmin_imputed_flagged.Rdata')

empDF %>% 
   mutate(FLAG0 = case_when(
      grepl('^R', FLAG0) ~ 'maybeDISCORDANT',
      FLAG0 %in% c('NOT TESTED', 'NT + nev', 'NEVER TESTED') ~ 'NOT/NEVER TESTED',
      .default = FLAG0
   )) %>% count(FLAG0, sort=TRUE)
empDF %>% 
   mutate(FLAG01 = case_when(
      grepl('^R', FLAG01) ~ 'maybeDISCORDANT',
      FLAG01 %in% c('NOT TESTED', 'NT + nev', 'NEVER TESTED') ~ 'NOT/NEVER TESTED',
      .default = FLAG01
   )) %>% count(FLAG01, sort=TRUE)
empDF %>% 
   mutate(TFLAG = case_when(
      grepl('^R', TFLAG) ~ 'maybeDISCORDANT',
      TFLAG %in% c('NOT TESTED', 'NT + nev', 'NEVER TESTED') ~ 'NOT/NEVER TESTED',
      .default = TFLAG
   )) %>% count(TFLAG, sort=TRUE)

# MEGAN IS WORKING ON THIS
# bugs <- empDF %>% count(BUG, sort=TRUE)
# 
# # HOW?? s
# tdisc <- empDF %>% filter(TFLAG == 'DISCORDANT') # only 40 (out of >4,000)
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
      mutate(TFLAG = case_when(
         grepl('^R', TFLAG) ~ 'maybeDISCORDANT',
         TFLAG %in% c('NOT TESTED', 'NT + nev', 'NEVER TESTED') ~ 'NOT/NEVER TESTED',
         .default = TFLAG
      ),
      EFLAG = case_when(
         grepl('^R', EFLAG) ~ 'maybeDISCORDANT',
         EFLAG %in% c('NOT TESTED', 'NT + nev', 'NEVER TESTED') ~ 'NOT/NEVER TESTED',
         .default = EFLAG
      )) %>%
      select(EFLAG, TFLAG) %>%
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
















