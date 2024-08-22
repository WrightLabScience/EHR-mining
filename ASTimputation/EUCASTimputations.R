library(dplyr)
library(readxl)
library(tidyr)

#R and S tables
load("/Users/megandollard/Downloads/EmpiricDF_bloodCultures_2015_2023_imputed.Rdata") # empDF

# empTEST <- data.frame(unique(empDF$BUG))
# empTEST[,2:97] <- NA
# colnames(empTEST)[1:97] <- names(empDF[6:102])

empTEST <- empDF

R1 <- read_xlsx('~/Desktop/Antibiotic Project/ExpectedASTsEUCAST/Table1_exp_res.xlsx')
names(R1)[1] <- "ORGANISMS"
R2 <- read_xlsx('~/Desktop/Antibiotic Project/ExpectedASTsEUCAST/Table2_exp_res.xlsx')
names(R2)[1] <- "ORGANISMS"
R3 <- read_xlsx('~/Desktop/Antibiotic Project/ExpectedASTsEUCAST/Table3_exp_res.xlsx')
names(R3)[1] <- "ORGANISMS"
R4 <- read_xlsx('~/Desktop/Antibiotic Project/ExpectedASTsEUCAST/Table4_exp_res.xlsx')
names(R4)[1] <- "ORGANISMS"
R5 <- read_xlsx('~/Desktop/Antibiotic Project/ExpectedASTsEUCAST/Table5_exp_res.xlsx')
names(R5)[1] <- "ORGANISMS"

# combined resistance DF
r <- full_join(full_join(full_join(full_join(R1, R2), R3), R4), R5)

S1 <- read_xlsx('~/Desktop/Antibiotic Project/ExpectedASTsEUCAST/Table1_exp_susc.xlsx')
names(S1)[1] <- "ORGANISMS"
S2 <- read_xlsx('~/Desktop/Antibiotic Project/ExpectedASTsEUCAST/Table2_exp_susc.xlsx')
names(S2)[1] <- "ORGANISMS"
S3 <- read_xlsx('~/Desktop/Antibiotic Project/ExpectedASTsEUCAST/Table3_exp_susc.xlsx')
names(S3)[1] <- "ORGANISMS"

# combined susceptibility DF
s <- full_join(full_join(S1, S2), S3)

# combined resistance/susceptibility DF
rs <- full_join(r, s)
rs <- rs %>% 
  pivot_longer(cols = -ORGANISMS,
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
  pivot_wider(names_from = ABX,
              values_from = RULE)


names(rs) <- toupper(names(rs)) # WHAT TO DO WITH COMMA IN ABX NAMES ? "POLYMYXIN B, COLISTIN"

rs$ORGANISMS <- sub("(.+)( species)", "\\1", rs$ORGANISMS)
rs$ORGANISMS <- sub("Coagulase Negative Staphylococcus", "Coagulase Negative Staph", rs$ORGANISMS)
rs$ORGANISMS <- sub("Group \\[ABCG] Beta-Hemolytic Streptococci", "Group \\[ABCG] Streptococci", rs$ORGANISMS)
rs$ORGANISMS <- sub("Salmonella Typhi", "Salmonella typhi", rs$ORGANISMS)


start <- Sys.time()
for(b in 1:nrow(rs)){ #bug by bug
  print(b)
  for (a in 2:ncol(rs)){ #abx by abx
    cl <- which(colnames(rs)[a] == colnames(empTEST)) #find drug col in empTEST 
    
    if(length(cl) == 0) next                      #if drug not found, move on
    
    val <- rs[[a]][b]                             # relevant rs value
    if(is.na(val)) next                           # if rs value is NA, next
    
    rws <- grep(rs$ORGANISMS[b], empTEST$BUG)     # rows relevant to to bug
    
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



#Expert Rules (species by species)

# rules (Staph and Strep)
StaphRules <- read.table(file = "~/Downloads/EUCAST Expert Rules v 3.2 - Staph Rules Condensed.tsv", sep = "\t")
names(StaphRules) <- toupper(StaphRules[1,])
StaphRules <- StaphRules[-1,]
rownames(StaphRules) <- NULL
StaphRules[1, 4] <- "R if benzylpenicillin; NA" # temporary

StrepRules <- read.table(file = "~/Downloads/EUCAST Expert Rules v 3.2 - Strep Rules Condensed.tsv", sep = "\t")
names(StrepRules) <- toupper(StrepRules[1,])
StrepRules <- StrepRules[-1,]
rownames(StrepRules) <- NULL

rules <- full_join(StaphRules, StrepRules) # join Staph and Strep rules

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

### BUGS ####
Viridans <- c("Viridans Streptococci",
              "Alpha Hemolytic Streptococci",
              "Streptococcus pneumoniae",
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

# output vector of names when col/row is organized by a broader category
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
  4                                                              #output: 4 if no case matches - shouldn't happen(?)
}


start <- Sys.time()
for(b in 1:nrow(rules)){ #bug/species by bug/species in rules
  print(b)
  for (a in 2:ncol(rules)){ #drug by drug in rules
    cl <- which(colnames(empTEST) %in% n(names(rules[a]))) #find drug col in ast ... searches col name in function n to see if larger category
    if(length(cl) == 0) next                               #if drug(s) not found, move on
    val <- rules[[a]][b]
    if(is.na(val)) next
    rws <- which(empTEST$BUG %in% n(rules$BUG[b]))         #ast rows related to bug/species ... searches bug name in function n to see if larger category
    # if(val == "0" || val == "1"){                                     # left from old version, useful if want to combine versions later
    #    empTEST[rws, cl][is.na(empTEST[rws, cl])] <- as.integer(val)
    #    next
    # }
    for(r in rws){ #rows that pertain to bug in ast  # row by row for relevant bugs in ast
      empTEST[r, cl[is.na(empTEST[r, cl])]] <- op(val, r)     #run value and ast row nums through op function
    }
  }
}


end <- Sys.time()
print(end - start) 




#checks
#NAs in original data 4007838
count <- integer(ncol(empDF))
for (i in 1:ncol(empDF)){
  count[i] <- sum(is.na(empDF[[i]]))
}
sum(count) 

#NAs in new data 3788496 3731284 3727816 3714561
count2 <- integer(ncol(empTEST))
for (i in 1:ncol(empTEST)){
  count2[i] <- sum(is.na(empTEST[[i]]))
}
sum(count2) 

#number of flag values
count4 <- integer(ncol(empTEST))
for (i in 7:ncol(empTEST)){
  count4[i] <- sum(empTEST[[i]] == 4, na.rm = T)
}
sum(count4) 
