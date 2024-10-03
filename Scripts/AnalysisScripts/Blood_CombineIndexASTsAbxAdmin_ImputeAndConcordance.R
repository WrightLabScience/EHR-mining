# The purpose of this script is to combine AST and AbxAdmin data (for bloodstream infections):
#     -Comine and save the long format for understanding abx admin proximity to order
#     -Collapse to one row per isolate, save this for future analyses
#     -Impute missing AST values
#     -Compute concordance/discordance flag (both on imputed and un-imputed tables)

# load necessary data
library(dplyr)
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/data_path_name.Rdata')



################################################################################
############################### START IMPUTATION ###############################
################################################################################
start <- Sys.time()
library(dplyr)
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/data_path_name.Rdata')
#load(file = paste0(data_path_name, 'ASTs_AbxAdmin_blood_2017_2023.Rdata'))
load(file = paste0(data_path_name, 'ASTs_blood_2015_2023.Rdata'))
astDF <- astDF %>% filter(substr(ORDER_DATE, 1, 4) %in% as.character(2017:2023))
Viridans <- c("Viridans Streptococci", "Alpha Hemolytic Streptococci", "Streptococcus mitis", "Streptococcus mutans", "Streptococcus sanguinis", "Streptococcus mitis/oralis group", "Streptococcus salivarius")
BetaHemolytic <- c("Beta Hemolytic Streptococci", "Group A Streptococci", "Group B Streptococci", "Group C Streptococci", "Group G Streptococci", "Streptococcus dysgalactiae", "Streptococcus pyogenes", 
                   "Streptococcus agalactiae", "Streptococcus equi", "Streptococcus gallolyticus", "Streptococcus bovis", "Streptococcus equinis", "Streptococcus suis")

# KADRI imputation rules
astCOPY <- astDF
w_first_abx_col <- which(names(astDF) == 'CEFEPIME')
{
   source(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/ASTimputation/ImputationRulesClean.R')
   abx_names <- names(astDF)[w_first_abx_col:length(astDF)]
   imp_rules$Antibiotic[!unique(imp_rules$Antibiotic) %in% abx_names] # CEFADROXIL, CEFPODOXIME, CEPHALEXIN, DICLOXACILLIN, NAFCILLIN
   abx_names[!abx_names %in% imp_rules$Antibiotic] # lots
   rm(abx_names)
   
   # find which rows of astDF correspond to which bugs
   bug_rows <- setNames(vector('list', length(imp_rules) - 1), names(imp_rules)[-1])
   for (i in seq_along(bug_rows)) {
      print(i)
      bug_name <- names(bug_rows)[i]
      if (grepl('Beta hemolytic', bug_name)) {
         bug_rows[[i]] <- which(astDF$BUG %in% BetaHemolytic)
         next
      }
      
      if (grepl('Alpha hemolytic|viridans', bug_name)) {
         bug_rows[[i]] <- which(astDF$BUG %in% Viridans)
         next
      }
      
      if (grepl('staphylococcus aureus', bug_name)) {
         w_sa <- grep('Staphylococcus aureus', astDF$BUG)
         if (grepl('sensitive', bug_name)) bug_rows[[i]] <- w_sa[which(astDF$OXACILLIN[w_sa] == 0)]
         if (grepl('resistant', bug_name)) bug_rows[[i]] <- w_sa[which(astDF$OXACILLIN[w_sa] == 1)]
         next
      }
      
      bug_rows[[i]] <- grep(bug_name, astDF$BUG)
   }
   lengths(bug_rows)
   sum(lengths(bug_rows)) / nrow(astDF) # 66.8%
   u_bug_rows <- sort(unlist(bug_rows))
   names(u_bug_rows) <- NULL
   any(diff(u_bug_rows) == 0) # no overlap!
   head(sort(table(astDF$BUG[-u_bug_rows]), decreasing = TRUE), n=12) # frequent bugs for which we don't have imputation rules
   rm(u_bug_rows, i, bug_name, w_sa)
   
   #input: a value from imp_rules and an ast row number
   #output: 0, 1, or NA
   op <- function(val, r) {
      if (val == 'S' || val == 'R') return(ifelse(val == 'S', 0L, 1L))  #if srW val is already "0" or "1", output integer version
      if (val == 'NA') return(NA)                            #if srW val is already "NA", output NA
      
      # E indicates that the ast score of the named antibiotic should be the output, unless that score is also missing (NA)
      if (grepl("^E", val)) {
         rest <- gsub("^(E )([A-Z/]+)(; )(.+)", "\\4", val)                  #stores everything after semicolon
         abx <- toupper(gsub("^(E )([A-Z/]+)(;.+)", "\\2", val))             #extracts abx name
         col <- which(colnames(astDF) == abx)                                          #astDF BC astCOPY VALS COULD ALREADY BE CHANGED?
         if (length(col) == 0) return(op(rest, r))                     #if abx name doesn't match, run op on rest of imp_rules val
         newval <- astDF[[col]][r]                                                     #astDF BC astCOPY VALS COULD ALREADY BE CHANGED?
         if (is.na(newval)) return(op(rest, r))                                        #if abx not tested, run op on rest of imp_rules val
         return(newval)                                             #output: abx val if 1 or 0
      }
      
      # S indicates that the ast score of the named antibiotic should be the output only if the score is 0
      if (grepl("^S IF", val)) {
         rest <- gsub("^(S IF )([A-Z/]+)(; )(.+)", "\\4", val)
         abx <- toupper(gsub("^(S IF )([A-Z/]+)(;.+)", "\\2", val))
         col <- which(colnames(astDF) == abx)
         if (length(col) == 0) return(op(rest, r))                     #if abx name doesn't match, run op on rest of imp_rules val
         newval <- astDF[[col]][r]
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
   
   #start <- Sys.time()
   for (b in 2:ncol(imp_rules)) { #row by row in imp_rules (bug/species by bug/species)
      cat(b, names(imp_rules[b]), '\n')
      rws <- bug_rows[[names(imp_rules[b])]]               #ast rows related to bug/species
      if (length(rws) == 0) next
      
      for (a in seq_along(imp_rules$Antibiotic)) { # drug by drug
         #print(a)
         abx <- imp_rules$Antibiotic[a]
         val <- imp_rules[[b]][a]
         if (val == "NA") next                                # checking if val is "straightforward" (NA)
         
         cl <- which(colnames(astCOPY) == abx)   # find drug col in ast 
         if (length(cl) == 0) next                            # if drug not found, move on
         
         w <- rws[which(is.na(astCOPY[[cl]][rws]))]  # find which rows of rws (this bug) are missing
         if (length(w) == 0) next                    # if no AST results are missing, move on
         
         if (val == 'S' | val == 'R') {                       # checking if val is "straightforward" (0 or 1)
            astCOPY[[cl]][w] <- ifelse(val == 'S', 0L, 1L)
            next
         }
         
         for (r in w) { #rows that pertain to bug in ast    # if val is not "straightforward", need to go row by row for relevant bugs in ast
            astCOPY[r, cl] <- op(val, r)                      # run value and ast row nums through op function
         }
      }
   }
   #print(Sys.time() - start) # ~  minutes
   rm(a, abx, b, cl, r, rws, val, w, op, imp_rules, bug_rows)
   
   # check for errors
   which(sapply(astCOPY[w_first_abx_col:length(astCOPY)], function(x) 4 %in% x))
}

# EUCAST expected phenotypes
astTEST <- astCOPY
{
   path <- '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/ASTimputation/'
   R1 <- readxl::read_xlsx(paste0(path, 'Table1_exp_res.xlsx'))
   names(R1)[1] <- "ORGANISMS"
   R2 <- readxl::read_xlsx(paste0(path, 'Table2_exp_res.xlsx'))
   names(R2)[1] <- "ORGANISMS"
   R3 <- readxl::read_xlsx(paste0(path, 'Table3_exp_res.xlsx'))
   names(R3)[1] <- "ORGANISMS"
   R4 <- readxl::read_xlsx(paste0(path, 'Table4_exp_res.xlsx'))
   names(R4)[1] <- "ORGANISMS"
   R5 <- readxl::read_xlsx(paste0(path, 'Table5_exp_res.xlsx'))
   names(R5)[1] <- "ORGANISMS"
   
   # combined resistance DF
   r <- full_join(full_join(full_join(full_join(R1, R2), R3), R4), R5)
   
   S1 <- readxl::read_xlsx(paste0(path, 'Table1_exp_susc.xlsx'))
   names(S1)[1] <- "ORGANISMS"
   S2 <- readxl::read_xlsx(paste0(path, 'Table2_exp_susc.xlsx'))
   names(S2)[1] <- "ORGANISMS"
   S3 <- readxl::read_xlsx(paste0(path, 'Table3_exp_susc.xlsx'))
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
   rs <- rs %>% tidyr::pivot_wider(names_from = ABX, values_from = RULE)
   rm(r, s, R1, R2, R3, R4, R5, S1, S2, S3)
   
   names(rs) <- toupper(names(rs))
   
   rs$ORGANISMS <- sub("(.+)( species)", "\\1", rs$ORGANISMS)
   rs$ORGANISMS <- sub("Coagulase Negative Staphylococcus", "Coagulase Negative Staph", rs$ORGANISMS)
   rs$ORGANISMS <- sub("Group \\[ABCG] Beta-Hemolytic Streptococci", "Group \\[ABCG] Streptococci", rs$ORGANISMS)
   
   
   bug_rows <- sapply(rs$ORGANISMS, function(x) grep(x, astTEST$BUG))
   sort(lengths(bug_rows))
   sum(lengths(bug_rows)) / nrow(astDF) # 70%
   u_bug_rows <- sort(unlist(bug_rows))
   head(sort(table(astDF$BUG[-u_bug_rows]), decreasing = TRUE), n=12) # frequent bugs for which we don't have imputation rules
   rm(u_bug_rows)
   
   
   #start <- Sys.time()
   for (b in 1:nrow(rs)){ #bug by bug
      rws <- bug_rows[[rs$ORGANISMS[b]]]
      for (a in 2:ncol(rs)){ #abx by abx
         cl <- which(colnames(rs)[a] == colnames(astTEST)) #find drug col in astTEST 
         
         if(length(cl) == 0) next                      #if drug not found, move on
         
         val <- rs[[a]][b]                             # relevant rs value
         if(is.na(val)) next                           # if rs value is NA, next
         
         if(val == "R"){ 
            astTEST[rws, cl][is.na(astTEST[rws, cl])] <- 1  # if rs value is R, replace NAs with 1
            next
         }
         
         if(val == "S"){
            astTEST[rws, cl][is.na(astTEST[rws, cl])] <- 0  # if rs value is S, replace NAs with 0
            next
         }
      }
   }
   #end <- Sys.time()
   #print(end - start) # 1.3 seconds
   rm(a, b, cl, end, rws, val, rs, bug_rows)
}

# EUCAST expert rules
astTEST2 <- astTEST
{
   StaphRules <- read.table(file = paste0(path, 'EUCAST Expert Rules v 3.2 - Staph Rules Condensed.tsv'), sep = "\t")
   names(StaphRules) <- toupper(StaphRules[1,])
   StaphRules <- StaphRules[-1,]
   rownames(StaphRules) <- NULL
   StaphRules[1, 4] <- "R if benzylpenicillin; NA" # tastorary
   
   StrepRules <- read.table(file = paste0(path, 'EUCAST Expert Rules v 3.2 - Strep Rules Condensed.tsv'), sep = "\t")
   names(StrepRules) <- toupper(StrepRules[1,])
   StrepRules <- StrepRules[-1,]
   rownames(StrepRules) <- NULL
   
   rules <- full_join(StaphRules, StrepRules) # join Staph and Strep rules
   rm(StaphRules, StrepRules, path)
   rules$BUG[rules$BUG == 'Beta-Hemolytic streptococci'] <- 'Beta Hemolytic Streptococci'
   rules$BUG[rules$BUG == 'Viridans group streptococci'] <- 'Viridans Streptococci'
   
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
      
      if(char == "Staphylococcus") return(grep("Staph", unique(astDF$BUG), value = T))
      if(char == "Viridans group streptococci") return(Viridans)
      if(char == "Beta-Hemolytic streptococci") return(BetaHemolytic)
      
      return(char)
   }
   
   op <- function(val, r) {
      #if(val == "0" || val == "1") return(as.integer(val))          # left from old version, useful if want to combine versions later
      if(val == "NA") return(NA)
      if(grepl("^E", val)){                                          #E indicates that the ast score of the named antibiotic should be the output, unless that score is also missing (NA)
         rest <- gsub("^(E )([a-z]+-?[a-z]+)(; )(.+)", "\\4", val)                  #stores everything after semicolon
         abx <- toupper(gsub("^(E )([a-z]+-?[a-z]+)(;.+)", "\\2", val))             #extracts abx name
         col <- colnames(astTEST)[colnames(astTEST) %in% abx]                                          
         if(length(col) == 0) return(op(rest, r))                     #if abx name doesn't match, run op on rest
         newval <- astTEST[[col]][r]                                                     
         if(is.na(newval)){
            return(op(rest, r))                                        #if abx not tested, run op on rest
         }else{
            return(newval)                                             #output: abx val if 1 or 0
         }
      }
      if(grepl("^R if", val)){                                       #R if indicates that the ast score of the named antibiotic should be the output only if the score is 1
         rest <- gsub("^(R if )([a-z]+-?[a-z]+)(; )(.+)", "\\4", val)
         abx <- toupper(gsub("^(R if )([a-z]+-?[a-z]+)(;.+)", "\\2", val))
         col <- which(colnames(astTEST) %in% abx)
         if(length(col) == 0) return(op(rest, r))                     
         newval <- astTEST[[col]][r]
         if(is.na(newval) | newval == 0){ 
            return(op(rest, r))                                       
         }else{
            return(1)                                                  #output: abx val if 1
         }
      }
      if(grepl("^S if", val)){                                       #S if indicates that the ast score of the named antibiotic should be the output only if the score is 0
         rest <- gsub("^(S if )([a-z]+-?[a-z]+)(; )(.+)", "\\4", val)
         abx <- toupper(gsub("^(S if )([a-z]+-?[a-z]+)(;.+)", "\\2", val))
         col <- which(colnames(astTEST) %in% abx)
         if(length(col) == 0) return(op(rest, r))                     
         newval <- astTEST[[col]][r]
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
         col <- which(colnames(astTEST) %in% c(abx1, abx2))
         if(length(col) == 0) return(op(rest, r))                     
         newval1 <- astTEST[[col[1]]][r]
         newval2 <- astTEST[[col[2]]][r]
         if(any(is.na(c(newval1, newval2))) | any(c(newval1 == 1, newval2 == 1))){ 
            return(op(rest, r))                                        
         }else{
            return(0)                                                  #output: abx val if both abx scores are 0
         }
      }
      print('ERROR')
      4                                                              #output: 4 if no case matches - shouldn't happen(?)
   }
   
   bug_rows <- sapply(rules$BUG, function(x) grep(x, astTEST2$BUG))
   lengths(bug_rows)
   
   #start <- Sys.time()
   for (b in 1:nrow(rules)) { #bug/species by bug/species in rules
      cat(b, rules$BUG[b], '\n')
      for (a in 2:ncol(rules)){ #drug by drug in rules
         cl <- which(colnames(astTEST2) %in% n(names(rules[a]))) #find drug col in ast ... searches col name in function n to see if larger category
         if(length(cl) == 0) next                               #if drug(s) not found, move on
         val <- rules[[a]][b]
         if(is.na(val)) next
         rws <- bug_rows[[rules$BUG[b]]]         #ast rows related to bug/species ... searches bug name in function n to see if larger category
         # if(val == "0" || val == "1"){                                     # left from old version, useful if want to combine versions later
         #    astTEST2[rws, cl][is.na(astTEST2[rws, cl])] <- as.integer(val)
         #    next
         # }
         for(r in rws){ #rows that pertain to bug in ast  # row by row for relevant bugs in ast
            # cat(a, r, '\n')
            astTEST2[r, cl[is.na(astTEST2[r, cl])]] <- op(val, r)     #run value and ast row nums through op function
         }
      }
   }
   #print(Sys.time() - start) # ~2.5 minutes
   rm(a, AMINOPENICILLINS, b, CARBAPENEMS, CEPCAR, CEPHALOSPORINS, CEPHg1, CEPHg2, CEPHg3, CEPHg4, CEPHg5, cl, FLUOROQUINOLONES,
      LINCOSAMIDES, MACLIN, MACROLIDES, r, rws, val, n, op, bug_rows, rules)
}

# check the difference between unimputed and imputed
df <- data.frame(old = sapply(astDF,    function(x) sum(!is.na(x))), new = sapply(astTEST2, function(x) sum(!is.na(x))))
df$diff <- df$new - df$old
df <- df[order(df$diff, decreasing=TRUE), ]

astDF <- astTEST2

# prep list columns containing R, S, and Not-tested for later
# I will use these list cols to determine if each individual therapy was concordant/discordant/not tested
astDF$RESISTANT <- astDF$SUSCEPTIBLE <- vector('list', nrow(astDF))
astDF <- astDF %>% relocate(RESISTANT, SUSCEPTIBLE, .before=CEFEPIME)
abx_cols <- which(names(astDF) == 'CEFEPIME'):which(names(astDF) == 'DELAFLOXACIN')
abx_names <- colnames(astDF)[abx_cols]
start <- Sys.time()
for (i in seq_len(nrow(astDF))) {
   if (i %% 500 == 0) print(i)
   astDF$RESISTANT[i] <- list(sort(unique(abx_names[astDF[i, abx_cols] == 1])))
   astDF$SUSCEPTIBLE[i] <- list(sort(unique(abx_names[astDF[i, abx_cols] == 0])))
}
print(Sys.time() - start) # ~5.3 minutes
rm(abx_cols, abx_names, i)


save(astDF, file = '~/Desktop/EHR/EHR work/RdataFiles/ASTs_blood_2017_2023_imputed.Rdata')
#save(astDF, file = '~/Desktop/EHR/EHR work/RdataFiles/ASTs_AbxAdmin_blood_2017_2023_imputed.Rdata')
print(Sys.time() - start) # ~3 minutes
################################################################################
################################ END IMPUTATION ################################
################################################################################







################################################################################
########################### START COMBINE AST + ABX ############################
################################################################################
start <- Sys.time()
library(dplyr)
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/data_path_name.Rdata')
#load(file = paste0(data_path_name, 'ASTs_blood_2015_2023.Rdata')) # 55,556
load(file = paste0(data_path_name, 'ASTs_blood_2017_2023_imputed.Rdata'))   #    44,768
load(file = paste0(data_path_name, 'ALL_CLEANED_2017_2023_AbxAdmin.Rdata')) # 9,681,748
abxDF <- abxDF %>% filter(PERSON_ID %in% unique(astDF$PERSON_ID))           # 3,030,936
length(unique(astDF$PERSON_ID)) # 33,972
length(unique(abxDF$PERSON_ID)) # 32,429
length(intersect(astDF$PERSON_ID, abxDF$PERSON_ID)) # 32,429

# prep data
abxDF <- abxDF %>%
   mutate(PERSON_ID = as.character(PERSON_ID)) %>%
   rename(START_DATE = ADMIN_START_DATE,
          END_DATE = ADMIN_END_DATE) %>%
   arrange(PERSON_ID, START_DATE)

astDF %>% count(MULT_BLOOD_ISO) # 7,091
astDF %>% group_by(PERSON_ID, ORDER_DAY) # 40,951



# JOIN abx_admin + ASTs
empDF <- astDF %>% # 1,327,723
   mutate(JOIN_START = ORDER_DAY - 7,
          JOIN_END = RESULT_DAY + 14) %>%
   left_join(x = .,
             y = abxDF,
             by = join_by(PERSON_ID,
                          JOIN_START <= START_DATE,
                          JOIN_END >= START_DATE)) %>%
   select(-JOIN_START, -JOIN_END) %>%
   relocate(ABX, START_DATE, END_DATE, .before=BUG)

### ASSESS CONCORDANCE ###
# test
chunks <- mapply(FUN = ':',
                 seq(1, 1300001, 100000),
                 c(seq(100000, 1300000, 100000), nrow(empDF)))

START <- Sys.time()
for (c in seq_along(chunks)) {
   print(c)
   chunk <- chunks[[c]]
   df <- empDF[chunk,] %>% select(ABX, RESISTANT, SUSCEPTIBLE)
   df$FLAG <- NA
   w <- which(is.na(df$ABX))
   df$FLAG[w] <- 'No therapy given'
   
   start <- Sys.time()
   for (i in seq_len(nrow(df))[-w]) {
      abx <- df$ABX[i]
      if (any(df$SUSCEPTIBLE[[i]] == abx)) {
         df$FLAG[i] <- 'CONCORDANT'
      } else if (any(df$RESISTANT[[i]] == abx)) {
         df$FLAG[i] <- 'DISCORDANT'
      } else {
         df$FLAG[i] <- 'NOTTESTED'
      }
   }
   total <- Sys.time() - start
   print(total) # ~30 seconds
   
   empDF$FLAG[chunk] <- df$FLAG
}
print(Sys.time() - START) # ~7 minutes
rm(START, start, i, c, abx, df, chunk, chunks, total, w) 


# save long format
save(empDF, file = paste0(data_path_name, 'ASTs_AbxAdmin_blood_2017_2023_imputed_flagged_LONG.Rdata'))


# antibiotic abbreviations
abbr <- tibble(read.table(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/antibiotic_names/ABX_ABBR.txt', header = TRUE)) %>% 
   select(-Class) %>%
   mutate(Antibiotic_Name = gsub('-', '/', Antibiotic_Name),
          Antibiotic_Name = gsub('_', ' ', Antibiotic_Name),
          Antibiotic_Name = toupper(Antibiotic_Name))
abbr <- setNames(abbr$Abbreviation, abbr$Antibiotic_Name)
abad <- unique(empDF$ABX)
abad <- abad[!abad %in% names(abbr)]
abad <- abad[!is.na(abad)]
abad <- abad[-grep(',', abad)]
abbr <- c(abbr, setNames(abad, abad))
rm(abad)
save(abbr, file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/antibiotic_names/abbr_named.Rdata')
rm(abbr)


if (FALSE) {
   library(dplyr)
   load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/data_path_name.Rdata')
   load(file = paste0(data_path_name, 'ASTs_AbxAdmin_blood_2017_2023_imputed_flagged_LONG.Rdata'))
}
# collapse so one row per isolate
chours <- 3600
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


plot(NA, xlim=c(-168, 432), ylim=c(1,1000))
for (i in 1:1000) {
   conc_times <- empDF$TIME_TO_CONC[[i]]
   if (length(conc_times) == 0) next
   points(x=conc_times, y=rep(i, length(conc_times)), pch=15, cex=0.25)
}
rm(i, conc_times)

# each row should have a 'TIME_TO_FIRST_CONCORDANT_THERAPY' column
# calculate the time to first concordant therapy, 
# only if it occurred > -48 hours before AST order
time_to_conc <- numeric(nrow(empDF))
w <- which(lengths(empDF$TIME_TO_CONC) == 0L)
time_to_conc[w] <- NA
for (i in seq_len(nrow(empDF))[-w]) {
   t <- empDF$TIME_TO_CONC[[i]]
   
   # only received concordant therapy pre-empirically, don't count this!
   if (all(t < -48)) {
      time_to_conc[i] <- NA
      next
   }
   
   # take only the times >= -48 hours
   time_to_conc[i] <- min(t[t >= -48])
}

empDF$TIME_TO_CONC <- sapply(empDF$TIME_TO_CONC, function(x) {
   if (any(x < -24)) return(-24)
   return(min(x))
})
empDF$TIME_TO_CONC[is.infinite(empDF$TIME_TO_CONC)] <- NA

empDF <- empDF %>%
   mutate(across(.cols = c(FLAGp:FLAGl),
                 .fns = ~ case_when(
                    . == '' ~ 'No abx given',
                    . == 'CONCORDANT+DISCORDANT' ~ 'CONCORDANT',
                    . == 'CONCORDANT+DISCORDANT+NOTTESTED' ~ 'CONCORDANT',
                    . == 'CONCORDANT+NOTTESTED' ~ 'CONCORDANT',
                    . == 'DISCORDANT+NOTTESTED' ~ 'DISCORDANT',
                    .default = .
                 )))

# save wide format
save(empDF, file = paste0(data_path_name, 'ASTs_AbxAdmin_blood_2017_2023_imputed_flagged_WIDE.Rdata'))
print(Sys.time() - start) # < 2 minutes

rm(abxDF, astDF, start, abbr, empDF)
gc()
################################################################################
############################ END COMBINE AST + ABX #############################
################################################################################









################################################################################
########################## START CONCORDANCE ANALYSIS ##########################
################################################################################
start <- Sys.time()
library(dplyr)
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/data_path_name.Rdata')
load(file = paste0(data_path_name, 'ASTs_AbxAdmin_blood_2017_2023.Rdata')); empDFo <- empDF
load(file = paste0(data_path_name, 'ASTs_AbxAdmin_blood_2017_2023_imputed.Rdata'))

empDF$FLAGp  <- empDF$FLAGe  <- empDF$FLAGb  <- empDF$FLAGt  <- empDF$FLAGl  <- NA
empDFo$FLAGp <- empDFo$FLAGe <- empDFo$FLAGb <- empDFo$FLAGt <- empDFo$FLAGl <- NA
empDF <- empDF %>% relocate(FLAGp, FLAGe, FLAGb, FLAGt, FLAGl, .before=EMPIRIC)
empDF <- empDF %>% relocate(FLAGp, FLAGe, FLAGb, FLAGt, FLAGl, .before=EMPIRIC)
col_names <- names(empDF %>% select(CEFEPIME:DELAFLOXACIN))
ABXp <- empDF$pEMPIRIC
ABXe <- empDF$EMPIRIC
ABXb <- empDF$BETWEEN
ABXt <- empDF$TARGETED
ABXl <- empDF$lTARGETED
source('~/Desktop/EHR/EHR-mining/Scripts/AnalysisScripts/ConcordanceFlagFunction.R')

for (i in seq_len(nrow(empDF))) {
   print(i)
   empDF$FLAGp[i] <- getConcordanceFlag(empDF[i,], ABXp[[i]])
   empDF$FLAGe[i] <- getConcordanceFlag(empDF[i,], ABXe[[i]])
   empDF$FLAGb[i] <- getConcordanceFlag(empDF[i,], ABXb[[i]])
   empDF$FLAGt[i] <- getConcordanceFlag(empDF[i,], ABXt[[i]])
   empDF$FLAGl[i] <- getConcordanceFlag(empDF[i,], ABXl[[i]])
   
   empDFo$FLAGp[i] <- getConcordanceFlag(empDFo[i,], ABXp[[i]])
   empDFo$FLAGe[i] <- getConcordanceFlag(empDFo[i,], ABXe[[i]])
   empDFo$FLAGb[i] <- getConcordanceFlag(empDFo[i,], ABXb[[i]])
   empDFo$FLAGt[i] <- getConcordanceFlag(empDFo[i,], ABXt[[i]])
   empDFo$FLAGl[i] <- getConcordanceFlag(empDFo[i,], ABXl[[i]])
}
print(Sys.time() - start) # ~4 minutes

save(empDF, file = paste0(data_path_name, 'ASTs_AbxAdmin_blood_2017_2023_imputed_flagged.Rdata'))
save(empDFo, file = paste0(data_path_name, 'ASTs_AbxAdmin_blood_2017_2023_flagged.Rdata'))

rm(getConcordanceFlag, start, i, col_names, emp0, emp1, emp2, emp3, emp4, emp5, emp6, empR, targ, empDF, empDFo)
gc()
################################################################################
########################### END CONCORDANCE ANALYSIS ###########################
################################################################################


rm(data_path_name)



