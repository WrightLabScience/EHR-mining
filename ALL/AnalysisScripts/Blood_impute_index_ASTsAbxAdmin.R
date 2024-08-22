### AST IMPUTATION
# 2017 - blood cultures - with med_admin

library(dplyr)
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ASTs_AbxAdmin_blood_2017.Rdata')

# From Kadri et al, 2021
source(file = '~/Desktop/EHR/EHR work/ASTimputation/ImputationRulesClean.R')

abx_names <- names(empDF)[14:length(empDF)]
imp_rules$Antibiotic[!unique(imp_rules$Antibiotic) %in% abx_names] # CEFADROXIL, CEPHALEXIN, DICLOXACILLIN
abx_names[!abx_names %in% imp_rules$Antibiotic] # lots
rm(abx_names)

Viridans <- c("Viridans Streptococci",
              "Alpha Hemolytic Streptococci",
              #"Streptococcus pneumoniae",
              "Streptococcus mitis",
              "Streptococcus mutans",
              "Streptococcus sanguinis",
              "Streptococcus mitis/oralis group",
              "Streptococcus salivarius")
BetaHemolytic <- c("Beta Hemolytic Streptococci",
                   "Group A Streptococci",
                   "Group B Streptococci",
                   "Group C Streptococci",
                   "Group G Streptococci",
                   "Streptococcus dysgalactiae",
                   "Streptococcus pyogenes",
                   "Streptococcus agalactiae",
                   "Streptococcus equi",
                   "Streptococcus gallolyticus",
                   "Streptococcus bovis",
                   "Streptococcus equinis",
                   "Streptococcus suis")

# find which rows of empDF correspond to which bugs
bug_rows <- setNames(vector('list', length(imp_rules) - 1), names(imp_rules)[-1])
for (i in seq_along(bug_rows)) {
   print(i)
   bug_name <- names(bug_rows)[i]
   if (grepl('Beta hemolytic', bug_name)) {
      bug_rows[[i]] <- which(empDF$BUG %in% BetaHemolytic)
      next
   }
   
   if (grepl('Alpha hemolytic|viridans', bug_name)) {
      bug_rows[[i]] <- which(empDF$BUG %in% Viridans)
      next
   }
   
   if (grepl('staphylococcus aureus', bug_name)) {
      w_sa <- grep('Staphylococcus aureus', empDF$BUG)
      if (grepl('sensitive', bug_name)) bug_rows[[i]] <- w_sa[which(empDF$OXACILLIN[w_sa] == 0)]
      if (grepl('resistant', bug_name)) bug_rows[[i]] <- w_sa[which(empDF$OXACILLIN[w_sa] == 1)]
      next
   }
   
   bug_rows[[i]] <- grep(bug_name, empDF$BUG)
}
lengths(bug_rows)
sum(lengths(bug_rows))
u_bug_rows <- sort(unlist(bug_rows))
names(u_bug_rows) <- NULL
any(diff(u_bug_rows) == 0) # no overlap!
head(sort(table(empDF$BUG[-u_bug_rows]), decreasing = TRUE), n=12) # frequent bugs for which we don't have imputation rules
rm(u_bug_rows, i, bug_name, w_sa)



empCOPY <- empDF

#input: a value from imp_rules and an ast row number
#output: 0, 1, or NA
op <- function(val, r) {
   if (val == 'S' || val == 'R') return(ifelse(val == 'S', 0L, 1L))  #if srW val is already "0" or "1", output integer version
   if (val == 'NA') return(NA)                            #if srW val is already "NA", output NA
   
   # E indicates that the ast score of the named antibiotic should be the output, unless that score is also missing (NA)
   if (grepl("^E", val)) {
      rest <- gsub("^(E )([A-Z/]+)(; )(.+)", "\\4", val)                  #stores everything after semicolon
      abx <- toupper(gsub("^(E )([A-Z/]+)(;.+)", "\\2", val))             #extracts abx name
      col <- which(colnames(empDF) == abx)                                          #empDF BC empCOPY VALS COULD ALREADY BE CHANGED?
      if (length(col) == 0) return(op(rest, r))                     #if abx name doesn't match, run op on rest of imp_rules val
      newval <- empDF[[col]][r]                                                     #empDF BC empCOPY VALS COULD ALREADY BE CHANGED?
      if (is.na(newval)) return(op(rest, r))                                        #if abx not tested, run op on rest of imp_rules val
      return(newval)                                             #output: abx val if 1 or 0
   }
   
   # S indicates that the ast score of the named antibiotic should be the output only if the score is 0
   if (grepl("^S IF", val)) {
      rest <- gsub("^(S IF )([A-Z/]+)(; )(.+)", "\\4", val)
      abx <- toupper(gsub("^(S IF )([A-Z/]+)(;.+)", "\\2", val))
      col <- which(colnames(empDF) == abx)
      if (length(col) == 0) return(op(rest, r))                     #if abx name doesn't match, run op on rest of imp_rules val
      newval <- empDF[[col]][r]
      if (is.na(newval) | newval == 1) { 
         return(op(rest, r))                                        #if abx test score is missing or 1, run op on rest of imp_rules val
      } else {
         return(0)                                                  #output: abx val if 0
      }
   }
   
   # output: 4 if no case matches - shouldn't happen(?)
   print('ERROR')
   return(4)
}


for (b in 2:ncol(imp_rules)) { #row by row in imp_rules (bug/species by bug/species)
   print('_________________________')
   print(names(imp_rules[b]))
   rws <- bug_rows[[names(imp_rules[b])]]               #ast rows related to bug/species
   if (length(rws) == 0) next
   
   #for (a in 1:25) {
   for (a in seq_along(imp_rules$Antibiotic)) { # drug by drug
      #print(a)
      abx <- imp_rules$Antibiotic[a]
      val <- imp_rules[[b]][a]
      if (val == "NA") next                                # checking if val is "straightforward" (NA)
      
      cl <- which(colnames(empCOPY) == abx)   # find drug col in ast 
      if (length(cl) == 0) next                            # if drug not found, move on
      
      w <- rws[which(is.na(empCOPY[[cl]][rws]))]  # find which rows of rws (this bug) are missing
      if (length(w) == 0) next                    # if no AST results are missing, move on
      
      if (val == 'S' | val == 'R') {                       # checking if val is "straightforward" (0 or 1)
         empCOPY[[cl]][w] <- ifelse(val == 'S', 0L, 1L)
         next
      }
      
      for (r in w) { #rows that pertain to bug in ast    # if val is not "straightforward", need to go row by row for relevant bugs in ast
         empCOPY[r, cl] <- op(val, r)                      # run value and ast row nums through op function
      }
   }
}
rm(a, abx, b, cl, r, rws, val, w, op)

which(sapply(empCOPY[14:length(empCOPY)], function(x) 4 %in% x))






# EUCAST EXPECTED PHENOTYPES:
library(readxl)
path <- '~/Desktop/EHR/EHR work/ASTimputation/EUCAST_expected_ASTs/ImputationRuleTables/'
R1 <- read_xlsx(paste0(path, 'Table1_exp_res.xlsx'))
names(R1)[1] <- "ORGANISMS"
R2 <- read_xlsx(paste0(path, 'Table2_exp_res.xlsx'))
names(R2)[1] <- "ORGANISMS"
R3 <- read_xlsx(paste0(path, 'Table3_exp_res.xlsx'))
names(R3)[1] <- "ORGANISMS"
R4 <- read_xlsx(paste0(path, 'Table4_exp_res.xlsx'))
names(R4)[1] <- "ORGANISMS"
R5 <- read_xlsx(paste0(path, 'Table5_exp_res.xlsx'))
names(R5)[1] <- "ORGANISMS"

# combined resistance DF
r <- full_join(full_join(full_join(full_join(R1, R2), R3), R4), R5)

S1 <- read_xlsx(paste0(path, 'Table1_exp_susc.xlsx'))
names(S1)[1] <- "ORGANISMS"
S2 <- read_xlsx(paste0(path, 'Table2_exp_susc.xlsx'))
names(S2)[1] <- "ORGANISMS"
S3 <- read_xlsx(paste0(path, 'Table3_exp_susc.xlsx'))
names(S3)[1] <- "ORGANISMS"

# combined susceptibility DF
s <- full_join(full_join(S1, S2), S3)

# combined resistance/susceptibility DF
rs <- full_join(r, s)
rs <- rs %>% 
   tidyr::pivot_longer(cols = -ORGANISMS,
                names_to = 'ABX',
                values_to = 'RULE') %>%
   arrange(ORGANISMS, ABX)
rs <- rbind(rs %>% 
               filter(all(is.na(RULE)), 
                      .by=c(ORGANISMS, ABX)) %>% 
               distinct(), # both are NA
            rs %>% 
               filter(any(!is.na(RULE)), 
                      .by=c(ORGANISMS, ABX)) %>% 
               filter(!is.na(RULE))) # one is not NA -> keep only not NA rules
rs <- rs %>%
   tidyr::pivot_wider(names_from = ABX,
               values_from = RULE)
rm(r, s, R1, R2, R3, R4, R5, S1, S2, S3, path)


names(rs) <- toupper(names(rs))

rs$ORGANISMS <- sub("(.+)( species)", "\\1", rs$ORGANISMS)
rs$ORGANISMS <- sub("Coagulase Negative Staphylococcus", "Coagulase Negative Staph", rs$ORGANISMS)
rs$ORGANISMS <- sub("Group \\[ABCG] Beta-Hemolytic Streptococci", "Group \\[ABCG] Streptococci", rs$ORGANISMS)


empTEST <- empCOPY
bug_rows <- sapply(rs$ORGANISMS, function(x) grep(x, empTEST$BUG))


start <- Sys.time()
for (b in 1:nrow(rs)){ #bug by bug
   print(b)
   rws <- bug_rows[[rs$ORGANISMS[b]]]
   for (a in 2:ncol(rs)){ #abx by abx
      cl <- which(colnames(rs)[a] == colnames(empTEST)) #find drug col in empTEST 
      
      if(length(cl) == 0) next                      #if drug not found, move on
      
      val <- rs[[a]][b]                             # relevant rs value
      if(is.na(val)) next                           # if rs value is NA, next
      
      if(val == "R"){ 
         empTEST[rws, cl][is.na(empTEST[rws, cl])] <- 1  # if rs value is R, replace NAs with 1
         next
      }
      
      if(val == "S"){
         empTEST[rws, cl][is.na(empTEST[rws, cl])] <- 0  # if rs value is S, replace NAs with 0
         next
      }
   }
}
end <- Sys.time()
print(end - start) 
rm(a, b, cl, end, rws, start, val, rs, bug_rows)

df <- data.frame(old = sapply(empCOPY[14:length(empCOPY)], function(x) sum(!is.na(x))),
                 new = sapply(empTEST[14:length(empTEST)], function(x) sum(!is.na(x))))
df$diff <- df$old - df$new
View(df)
rm(df)



# EUCAST expert rules
path <- '~/Desktop/EHR/EHR work/ASTimputation/EUCAST_expected_ASTs/'
StaphRules <- read.table(file = paste0(path, 'EUCAST Expert Rules v 3.2 - Staph Rules Condensed.tsv'), sep = "\t")
names(StaphRules) <- toupper(StaphRules[1,])
StaphRules <- StaphRules[-1,]
rownames(StaphRules) <- NULL
StaphRules[1, 4] <- "R if benzylpenicillin; NA" # temporary

StrepRules <- read.table(file = paste0(path, 'EUCAST Expert Rules v 3.2 - Strep Rules Condensed.tsv'), sep = "\t")
names(StrepRules) <- toupper(StrepRules[1,])
StrepRules <- StrepRules[-1,]
rownames(StrepRules) <- NULL

rules <- full_join(StaphRules, StrepRules) # join Staph and Strep rules
rm(StaphRules, StrepRules, path)

### ABX ####

# macrolides and lincosamides
MACROLIDES <- c("AZITHROMYCIN", "CLARITHROMYCIN", "ERYTHROMYCIN") # macrolides
LINCOSAMIDES <- c("CLINDAMYCIN", "LINCOMYCIN",  "PIRLIMYCIN") # lincosamides
MACLIN <- c(MACROLIDES, LINCOSAMIDES)
# fluoroquinolones
FLUOROQUINOLONES <- c("CIPROFLOXACIN", "GATIFLOXACIN",  "GEMIFLOXACIN",  "LEVOFLOXACIN",
                      "MOXIFLOXACIN", "NORFLOXACIN", "OFLOXACIN", "TEMAFLOXACIN")
# cephalosporins (split by generations) and carbapenems
CEPHg1 <- c("CEFAZOLIN", "CEPHALEXIN", "CEFALEXIN", "CEPHALOTHIN", "CEFALOTHIN",
            "CEFALOTIN", "CEFAPIRIN", "CEFRADINE", "CEFADROXIL", "CEFATRIZINE")
CEPHg2 <- c("CEFOXITIN", "CEFUROXIME", "CEFACLOR", "CEFPROZIL", "CEFMETAZOLE",
            "CEFONICID", "CEFOTETAN", "CEFBUPERAZONE")
CEPHg3 <- c("CEFTAZIDIME", "CEFTRIAXONE", "CEFOTAXIME", "CEFDINIR", "CEFPODOXIME",
            "CEFIXIME", "CEFTIZOXIME", "CEFOPERAZONE")
CEPHg4 <- c("CEFEPIME")
CEPHg5 <- c("CEFTAROLINE", "CEFTOLOZANE", "CEFTOBIPROLE")
CEPHALOSPORINS <- c(CEPHg1, CEPHg2, CEPHg3, CEPHg4, CEPHg5)
CARBAPENEMS <- c("DORIPENEM", "ERTAPENEM", "IMIPENEM", "MEROPENEM")
CEPCAR <- c(CEPHALOSPORINS, CARBAPENEMS)
# aminopenicillins
AMINOPENICILLINS <- c("AMOXICILLIN", "AMPICILLIN")

n <- function(char){
   if(char == "ALL MACROLIDES AND LINCOSAMIDES") return(MACLIN)
   if(char == "ALL FLUOROQUINOLONES") return(FLUOROQUINOLONES)
   if(char == "CEPHALOSPORINS AND CARBAPENEMS") return(CEPCAR)
   if(char == "AMINOPENICILLINS") return(AMINOPENICILLINS)
   
   if(char == "Staphylococcus") return(grep("Staph", unique(empDF$BUG), value = T))
   if(char == "Viridans group streptococci") return(Viridans)
   if(char == "Beta-Hemolytic streptococci") return(BetaHemolytic)
   
   return(char)
}

op <- function(val, r){
   #if(val == "0" || val == "1") return(as.integer(val))          # left from old version, useful if want to combine versions later
   if(val == "NA") return(NA)
   if(grepl("^E", val)){                                          #E indicates that the ast score of the named antibiotic should be the output, unless that score is also missing (NA)
      rest <- gsub("^(E )([a-z]+-?[a-z]+)(; )(.+)", "\\4", val)                  #stores everything after semicolon
      abx <- toupper(gsub("^(E )([a-z]+-?[a-z]+)(;.+)", "\\2", val))             #extracts abx name
      col <- colnames(empTEST)[colnames(empTEST) %in% abx]                                          
      if(length(col) == 0) return(op(rest, r))                     #if abx name doesn't match, run op on rest
      newval <- empTEST[[col]][r]                                                     
      if(is.na(newval)){
         return(op(rest, r))                                        #if abx not tested, run op on rest
      }else{
         return(newval)                                             #output: abx val if 1 or 0
      }
   }
   if(grepl("^R if", val)){                                       #R if indicates that the ast score of the named antibiotic should be the output only if the score is 1
      rest <- gsub("^(R if )([a-z]+-?[a-z]+)(; )(.+)", "\\4", val)
      abx <- toupper(gsub("^(R if )([a-z]+-?[a-z]+)(;.+)", "\\2", val))
      col <- which(colnames(empTEST) %in% abx)
      if(length(col) == 0) return(op(rest, r))                     
      newval <- empTEST[[col]][r]
      if(is.na(newval) | newval == 0){ 
         return(op(rest, r))                                       
      }else{
         return(1)                                                  #output: abx val if 1
      }
   }
   if(grepl("^S if", val)){                                       #S if indicates that the ast score of the named antibiotic should be the output only if the score is 0
      rest <- gsub("^(S if )([a-z]+-?[a-z]+)(; )(.+)", "\\4", val)
      abx <- toupper(gsub("^(S if )([a-z]+-?[a-z]+)(;.+)", "\\2", val))
      col <- which(colnames(empTEST) %in% abx)
      if(length(col) == 0) return(op(rest, r))                     
      newval <- empTEST[[col]][r]
      if(is.na(newval) | newval == 1){ 
         return(op(rest, r))                                        
      }else{
         return(0)                                                  #output: abx val if 0
      }
   }
   if(grepl("^A S if", val)){                                    # A S if indicated that the ast scores of both abx must be 0 for output to be 0
      rest <- gsub("^(A S if )([a-z]+-?[a-z]+ AND [a-z]+-?[a-z]+)(; )(.+)", "\\4", val)
      abx1 <- toupper(gsub("^(A S if )([a-z]+-?[a-z]+)( AND )([a-z]+-?[a-z]+)(;.+)", "\\2", val))
      abx2 <- toupper(gsub("^(A S if )([a-z]+-?[a-z]+)( AND )([a-z]+-?[a-z]+)(;.+)", "\\4", val))
      col <- which(colnames(empTEST) %in% c(abx1, abx2))
      if(length(col) == 0) return(op(rest, r))                     
      newval1 <- empTEST[[col[1]]][r]
      newval2 <- empTEST[[col[2]]][r]
      if(any(is.na(c(newval1, newval2))) | any(c(newval1 == 1, newval2 == 1))){ 
         return(op(rest, r))                                        
      }else{
         return(0)                                                  #output: abx val if both abx scores are 0
      }
   }
   print('ERROR')
   4                                                              #output: 4 if no case matches - shouldn't happen(?)
}



empTEST2 <- empTEST
bug_rows <- sapply(rules$BUG, function(x) grep(x, empTEST2$BUG))

start <- Sys.time()
for(b in 1:nrow(rules)){ #bug/species by bug/species in rules
   print(rules$BUG[b])
   for (a in 2:ncol(rules)){ #drug by drug in rules
      print(a)
      cl <- which(colnames(empTEST2) %in% n(names(rules[a]))) #find drug col in ast ... searches col name in function n to see if larger category
      if(length(cl) == 0) next                               #if drug(s) not found, move on
      val <- rules[[a]][b]
      if(is.na(val)) next
      rws <- bug_rows[[rules$BUG[b]]]         #ast rows related to bug/species ... searches bug name in function n to see if larger category
      # if(val == "0" || val == "1"){                                     # left from old version, useful if want to combine versions later
      #    empTEST2[rws, cl][is.na(empTEST2[rws, cl])] <- as.integer(val)
      #    next
      # }
      for(r in rws){ #rows that pertain to bug in ast  # row by row for relevant bugs in ast
         empTEST2[r, cl[is.na(empTEST2[r, cl])]] <- op(val, r)     #run value and ast row nums through op function
      }
   }
}
print(Sys.time() - start)
rm(a, AMINOPENICILLINS, b, BetaHemolytic, CARBAPENEMS, CEPCAR, CEPHALOSPORINS, CEPHg1, CEPHg2, CEPHg3, CEPHg4, CEPHg5, cl, FLUOROQUINOLONES,
   LINCOSAMIDES, MACLIN, MACROLIDES, r, rws, start, val, Viridans, n, op, bug_rows, rules)


df <- data.frame(old = sapply(empTEST[14:length(empTEST)], function(x) sum(!is.na(x))),
                 new = sapply(empTEST2[14:length(empTEST2)], function(x) sum(!is.na(x))))
df$diff <- df$old - df$new
View(df)
rm(df)



df <- data.frame(old = sapply(empDF[14:length(empDF)], function(x) sum(!is.na(x))),
                 new = sapply(empTEST2[14:length(empTEST2)], function(x) sum(!is.na(x))))
df$diff <- df$old - df$new
View(df)
rm(df)


empDF <- empTEST2

save(empDF, file = '~/Desktop/EHR/EHR work/RdataFiles/ASTs_AbxAdmin_blood_2017_imputed.Rdata')






