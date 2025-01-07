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
load(file = paste0(data_path_name, 'ASTs_blood_2017_2023.Rdata'))
Viridans <- c("Viridans Streptococci", "Alpha Hemolytic Streptococci", "Streptococcus mitis", "Streptococcus mutans", "Streptococcus sanguinis", "Streptococcus mitis/oralis group", "Streptococcus salivarius")
BetaHemolytic <- c("Beta Hemolytic Streptococci", "Group A Streptococci", "Group B Streptococci", "Group C Streptococci", "Group G Streptococci", "Streptococcus dysgalactiae", "Streptococcus pyogenes", 
                   "Streptococcus agalactiae", "Streptococcus equi", "Streptococcus gallolyticus", "Streptococcus bovis", "Streptococcus equinis", "Streptococcus suis")

# ESBL Enterobacteraceae
w_esbl <- which(astDF$ESBL == 1L)
for (abx in c("CEFOTAXIME", "CEFTAZIDIME", "CEFTRIAXONE", "CEFEPIME")) {
   w <- w_esbl[which(is.na(astDF[[abx]][w_esbl]))]
   astDF[[abx]][w] <- 1
   # print(sum(astDF[[abx]][w_esbl] == 0L, na.rm=T) / sum(!is.na(astDF[[abx]][w_esbl])))
}

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
   
   rm(a, AMINOPENICILLINS, b, CARBAPENEMS, CEPCAR, CEPHALOSPORINS, CEPHg1, CEPHg2, CEPHg3, CEPHg4, CEPHg5, cl, FLUOROQUINOLONES,
      LINCOSAMIDES, MACLIN, MACROLIDES, r, rws, val, n, op, bug_rows, rules)
}

astDF <- astTEST2
print(Sys.time() - start)
####### finish imputation #######


############ CREATE LIST COLS FOR RESISTANT AND SUSCEPTIBLE #########
# prep list columns containing R, S, and Not-tested for later
# I will use these list cols to determine if each individual therapy was concordant/discordant/not tested
start <- Sys.time()
astDF$RESISTANT <- astDF$SUSCEPTIBLE <- vector('list', nrow(astDF))
astDF <- astDF %>% relocate(RESISTANT, SUSCEPTIBLE, .before=CEFEPIME)
abx_cols <- which(names(astDF) == 'CEFEPIME'):which(names(astDF) == 'DELAFLOXACIN')
abx_names <- colnames(astDF)[abx_cols]
for (i in seq_len(nrow(astDF))) {
   if (i %% 500 == 0) print(i)
   astDF$RESISTANT[i] <- list(sort(unique(abx_names[astDF[i, abx_cols] == 1])))
   astDF$SUSCEPTIBLE[i] <- list(sort(unique(abx_names[astDF[i, abx_cols] == 0])))
}
rm(abx_cols, abx_names, i)
print(Sys.time() - start) # ~5 minutes
#################


save(astDF, file = '~/Desktop/EHR/EHR work/RdataFiles/ASTs_blood_2017_2023_imputed.Rdata')
gc()
################################################################################
################################ END IMPUTATION ################################
################################################################################








################################################################################
############################ JOIN ASTs + ENCOUNTERS ############################
################################################################################
start <- Sys.time()
library(dplyr)
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ASTs_blood_2017_2023_imputed.Rdata')
# WAIT! MAKE SURE TO GO GET NEW ENCOUNTERS IF THE BLOOD ASTs COHORT CHANGED BY INCLUDING DIFFERENT PATIENTS!!!
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ENCOUNTERS_cleaned_blood_2017_2023.Rdata')

encs <- encs %>%
   mutate(ADMIT_DAY = as.Date(substr(ADMIT_DATE,1,10)),
          DISCHARGE_DAY = as.Date(substr(DISCHARGE_DATE,1,10)),
          B4_ADMIT = ADMIT_DAY - 7,
          A4_DISCH = DISCHARGE_DAY + 7)

df <- astDF %>%
   left_join(x = .,
             y = encs,
             relationship = 'many-to-many',
             by = join_by(PERSON_ID,
                          between(ORDER_DAY, B4_ADMIT, A4_DISCH, bounds='[]')))
df %>% count(is.na(ADMIT_DAY)) # 40,377 / 44,250 matched (some matched more than 1)

# join on result_day as well
df <- rbind(df %>% filter(!is.na(ADMIT_DAY)),
            df %>% filter(is.na(ADMIT_DAY)) %>% # 3,873
               select(!c(ADMIT_DATE, DISCHARGE_DATE, ENCOUNTER_TYPE, FACILITY, ADMIT_SOURCE, ADMIT_DAY, DISCHARGE_DAY, B4_ADMIT, A4_DISCH)) %>%
               left_join(x = .,
                         y = encs,
                         by = join_by(PERSON_ID,
                                      between(RESULT_DAY, B4_ADMIT, A4_DISCH, bounds='[]'))))
df %>% count(is.na(ADMIT_DAY)) # 40,512 / 44,251
df <- df %>% 
   select(-B4_ADMIT, -A4_DISCH) %>%
   relocate(ADMIT_DAY, DISCHARGE_DAY, .after=RESULT_DAY)


# # look at whats still missing and check for nearest encounter, just to see whats up
# x <- df %>%
#    filter(is.na(ADMIT_DAY)) %>%
#    select(PERSON_ID, ORDER_DAY) %>%
#    mutate(YEAR = substr(ORDER_DAY,1,4)) %>%
#    left_join(y = encs %>% mutate(YEAR = substr(ADMIT_DAY, 1,4)),
#              relationship = 'many-to-many',
#              by = join_by(PERSON_ID, YEAR))
# x %>% count(is.na(ADMIT_DAY)) # still 2700 missing
# x <- x %>% mutate(D = abs(as.integer(ORDER_DAY - ADMIT_DAY)))
# hist(x$D, breaks=diff(range(x$D, na.rm=T)), xlim=c(0, 50))
# rm(x)


#### GROUP BY AST CULTURE ####
df <- df %>%
   group_by_all() %>%
   ungroup(ADMIT_DATE, DISCHARGE_DATE, ENCOUNTER_TYPE, FACILITY, ADMIT_SOURCE, ADMIT_DAY, DISCHARGE_DAY) # 43,735 groups (45,236 rows) 

df1 <- df %>% filter(n() == 1L) %>% ungroup() # 42,237
df2 <- df %>% filter(n() > 1L) # 1,498 groups (2,999 rows)

# visualize proximity of order to admission
par(mfrow=c(2,1))
x <- unname(ungroup(df2) %>% mutate(X = as.integer(ORDER_DAY - ADMIT_DAY)) %>% select(X) %>% unlist())
barplot(table(x), main='ORDER - ADMIT')
x <- unname(ungroup(df2) %>% mutate(X = as.integer(RESULT_DAY - ADMIT_DAY)) %>% select(X) %>% unlist())
barplot(table(x), main='RESULT - ADMIT')

# for df2, if a culture has one match that's day before order, day of, or after,
# keep that one and throw away the other
df2 %>% filter(ADMIT_DAY >= ORDER_DAY - 1 & ADMIT_DAY <= RESULT_DAY) # 1,444 groups (1,463 rows)
df2 %>% filter(ADMIT_DAY < ORDER_DAY - 1 | ADMIT_DAY > RESULT_DAY) # 1,480 groups (1,536 rows)

df2 <- df2 %>%
   reframe(
      ADMIT_DATE = min(ADMIT_DATE),
      DISCHARGE_DATE = max(DISCHARGE_DATE),
      ENCOUNTER_TYPE = paste(sort(unique(unlist(strsplit(ENCOUNTER_TYPE, split=' + ', fixed=TRUE)))), collapse=' + '),
      FACILITY = paste(sort(unique(unlist(strsplit(FACILITY, split=' + ', fixed=TRUE)))), collapse=' + '),
      ADMIT_SOURCE = paste(sort(unique(unlist(strsplit(ADMIT_SOURCE, split=' + ', fixed=TRUE)))), collapse=' + ')
   ) %>%
   mutate(ADMIT_DAY = as.Date(substr(ADMIT_DATE, 1, 10)),
          DISCHARGE_DAY = as.Date(substr(DISCHARGE_DATE, 1, 10)))

# combine the two groups again
astDF <- rbind(df1, df2) %>% arrange(PERSON_ID, ADMIT_DATE)

# save object
save(astDF, file = '~/Desktop/EHR/EHR work/RdataFiles/ASTs_blood_2017_2023_imputed_encs.Rdata')
print(Sys.time() - start) # ~ 10 seconds
################################################################################
########################## END JOIN ASTs + ENCOUNTERS ##########################
################################################################################









################################################################################
########################### START COMBINE AST + ABX ############################
################################################################################
start <- Sys.time()
library(dplyr)
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ASTs_blood_2017_2023_imputed_encs.Rdata')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_CLEANED_2017_2023_AbxAdmin.Rdata') # 9,681,748
abxDF <- abxDF %>% filter(PERSON_ID %in% unique(astDF$PERSON_ID))           # 3,030,936
length(unique(astDF$PERSON_ID)) # 33,972
length(unique(abxDF$PERSON_ID)) # 32,429
length(intersect(astDF$PERSON_ID, abxDF$PERSON_ID)) # 32,429

# prep data
abxDF <- abxDF %>%
   mutate(PERSON_ID = as.character(PERSON_ID)) %>%
   rename(START_DATE = ADMIN_START_DATE,
          END_DATE = ADMIN_END_DATE) %>%
   mutate(START_DAY = as.Date(substr(START_DATE, 1, 10))) %>%
   arrange(PERSON_ID, START_DATE)

astDF %>% group_by(PERSON_ID, ORDER_DAY) # 40,361
astDF %>% count(FACILITY, sort=TRUE)
table(lengths(strsplit(astDF$FACILITY, ' + ', fixed=TRUE))) # ~2K had 2 facilities, 69 had 3, 1 had 4
data.frame(sort(table(unlist(strsplit(astDF$FACILITY, ' + ', fixed=TRUE))), decreasing=TRUE))


# JOIN abx_admin + ASTs
astDF$JOIN_START <- astDF$ORDER_DAY
w <- which(astDF$ADMIT_DAY < astDF$ORDER_DAY)
astDF$JOIN_START[w] <- astDF$ADMIT_DAY[w]

astDF$JOIN_END <- astDF$RESULT_DAY
w <- which(astDF$DISCHARGE_DAY > astDF$RESULT_DAY)
astDF$JOIN_END[w] <- astDF$DISCHARGE_DAY[w]

empDF <- astDF %>% # 1.7 million
   left_join(x = .,
             y = abxDF,
             by = join_by(PERSON_ID,
                          JOIN_START <= START_DAY,
                          JOIN_END >= START_DAY)) %>%
   select(-JOIN_START, -JOIN_END) %>%
   relocate(ABX, START_DATE, END_DATE, .before=BUG)

### ASSESS CONCORDANCE ###
chunks <- mapply(FUN = ':',
                   seq(1,   floor(nrow(empDF) / 1e5) * 1e5 + 1, 1e5),
                 c(seq(1e5, floor(nrow(empDF) / 1e5) * 1e5,     1e5), nrow(empDF)))
empDF$FLAG <- NA
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
   print(total) # 23-24 seconds each
   
   empDF$FLAG[chunk] <- df$FLAG
}
print(Sys.time() - START) # ~8.8 minutes
rm(START, start, i, c, abx, df, chunk, chunks, total, w, astDF, abxDF) 



# indicate which window each abx record belongs to
chours <- 3600
empDF <- empDF %>%
   mutate(pEMPIRIC = START_DATE < (ORDER_DATE - 48*chours),
          EMPIRIC = START_DATE >= (ORDER_DATE - 48*chours) & START_DATE < (ORDER_DATE + 12*chours),
          BETWEEN = START_DATE >= (ORDER_DATE + 12*chours) & START_DATE < RESULT_DATE,
          eTARGETED = START_DATE >= RESULT_DATE & START_DATE < (RESULT_DATE + 12*chours),
          TARGETED = START_DATE >= RESULT_DATE + 12*chours & START_DATE <= (RESULT_DATE + 72*chours),
          lTARGETED = START_DATE > (RESULT_DATE + 72*chours),
          MISS_TIME = substr(START_DATE, 12, 12) == '' & !substr(END_DATE, 12, 13) %in% c('00', '01', '02', '03'),
          ABX_PROX_ORDER = as.numeric(lubridate::as.duration(START_DATE - ORDER_DATE)) / 3600)

# mark which prescriptions are in empiric window that I do not want to keep
empDF <- empDF %>% #        empiric, missing time in a suspicious way, and calendar day is day after order day
   mutate(REMOVE_FROM_EMP = EMPIRIC & MISS_TIME & as.Date(substr(START_DATE,1,10)) == ORDER_DAY + 1)

# old way of dealing with missing time stamps...
{
   # # prep data
   # empDF <- empDF %>% 
   #    mutate(ABX_DAY = as.Date(substr(START_DATE,1,10)),
   #           ORDER_AFTER_NOON = substr(ORDER_DATE, 12, 13) %in% as.character(12:23),
   #           MISS_TIME = substr(START_DATE, 12, 12) == '' & !substr(END_DATE, 12, 13) %in% c('00', '01', '02', '03')) # bad if 00:00:00 and END_TIME is not really close to midnight such that the abx start time might ACTUALLY BE 00:00:00
   # empDF %>% count(MISS_TIME) # 176,681
   # # use the end time as an estimate for the start time (subtract 2 hours)
   # #          START DATE is missing time stamp         END DATE exists          END DATE is NOT missing time stamp     
   # w <- which(empDF$MISS_TIME & !is.na(empDF$END_DATE) & substr(empDF$END_DATE, 12, 12) != '') # 64
   # empDF$START_DATE[w] <- empDF$END_DATE[w] - 2 * 3600
   # 
   # empDF <- empDF %>%
   #    mutate(MISS_TIME = substr(START_DATE, 12, 12) == '' & !substr(END_DATE, 12, 13) %in% c('00', '01', '02', '03'),
   #           ABX_PROX_ORDER = as.numeric(lubridate::as.duration(START_DATE - ORDER_DATE)) / 3600)
   # empDF %>% count(MISS_TIME) # 176,617
   # 
   # # mark which ABX records we will keep
   # empDF$KEEP_ABX <- !empDF$MISS_TIME; sum(empDF$KEEP_ABX, na.rm=T) # 1,146,741
   # # AST order was after noon and the ABX was the same calendar day
   # w <- which(!empDF$KEEP_ABX) # 176,617
   # w <- w[which(empDF$ORDER_AFTER_NOON[w] & empDF$ABX_DAY[w] == empDF$ORDER_DAY[w])] # 9,903
   # empDF$KEEP_ABX[w] <- TRUE; sum(empDF$KEEP_ABX, na.rm=T) # 1,156,644
   # # Abx was administered the day before the AST order - always safe
   # w <- which(!empDF$KEEP_ABX) # 166,714
   # w <- w[which(empDF$ABX_DAY[w] == empDF$ORDER_DAY[w] - 1)] # 1,667
   # empDF$KEEP_ABX[w] <- TRUE; sum(empDF$KEEP_ABX, na.rm=T) # 1,158,311
   # rm(w)
   # empDF <- empDF %>% select(-ORDER_AFTER_NOON)
   # empDF$KEEP_ABX[is.na(empDF$KEEP_ABX)] <- TRUE
   }

# save long format
save(empDF, file = '~/Desktop/EHR/EHR work/RdataFiles/ASTs_AbxAdmin_blood_2017_2023_imputed_encs_flagged_LONG.Rdata')



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
   load(file = '~/Desktop/EHR/EHR work/RdataFiles/ASTs_AbxAdmin_blood_2017_2023_imputed_encs_flagged_LONG.Rdata')
}

##### collapse to one row per isolate ##### 
empDF <- empDF %>% # 43,735
   group_by_all() %>%
   ungroup(ABX, START_DATE, START_DAY, END_DATE, pEMPIRIC, EMPIRIC, BETWEEN, eTARGETED, TARGETED, lTARGETED, FLAG, MISS_TIME, REMOVE_FROM_EMP, ABX_PROX_ORDER)
empDF <- empDF %>%
   reframe(num_abx_admin = n(),
           num_removed_from_emp = sum(REMOVE_FROM_EMP),
           TIME_TO_FIRST_ABX = list(ABX_PROX_ORDER[!REMOVE_FROM_EMP]),
           TIME_TO_FIRST_ABX_STRICT = list(ABX_PROX_ORDER[!MISS_TIME]),
           TIME_TO_CONC = list(ABX_PROX_ORDER[FLAG == 'CONCORDANT' & !REMOVE_FROM_EMP]),
           TIME_TO_CONC_STRICT = list(ABX_PROX_ORDER[FLAG == 'CONCORDANT' & !MISS_TIME]),
           ABX_EMPp = list(sort(unique(ABX[pEMPIRIC]))),
           ABX_EMP = list(sort(unique(ABX[EMPIRIC & !REMOVE_FROM_EMP]))),
           ABX_EMPr = list(sort(unique(ABX[REMOVE_FROM_EMP]))),
           ABX_BTW = list(sort(unique(ABX[BETWEEN]))),
           ABX_TARe = list(sort(unique(ABX[eTARGETED]))),
           ABX_TAR = list(sort(unique(ABX[TARGETED]))),
           ABX_TARl = list(sort(unique(ABX[lTARGETED]))),
           FLAGp = paste(sort(unique(FLAG[pEMPIRIC])), collapse='+'),
           FLAGe = paste(sort(unique(FLAG[EMPIRIC & !REMOVE_FROM_EMP])), collapse='+'),
           FLAGer = paste(sort(unique(FLAG[REMOVE_FROM_EMP])), collapse='+'),
           FLAGb = paste(sort(unique(FLAG[BETWEEN])), collapse='+'),
           FLAGet = paste(sort(unique(FLAG[eTARGETED])), collapse='+'),
           FLAGt = paste(sort(unique(FLAG[TARGETED])), collapse='+'),
           FLAGl = paste(sort(unique(FLAG[lTARGETED])), collapse='+')) %>%
   ungroup() %>%
   relocate(ABX_EMPp:ABX_TARl, FLAGp:FLAGl, TIME_TO_FIRST_ABX, TIME_TO_CONC, TIME_TO_FIRST_ABX_STRICT, TIME_TO_CONC_STRICT, num_abx_admin, num_removed_from_emp, .before=BUG)
# simple fixes
empDF$num_abx_admin[is.na(empDF$TIME_TO_FIRST_ABX)] <- 0
sum(empDF$num_abx_admin == 1L) / nrow(empDF) # 1.4%
sum(empDF$num_removed_from_emp > 0, na.rm=T) / sum(!is.na(empDF$num_removed_from_emp)) # 27.6% of antibiotic-receiving infections had at least one abx removed from the empiric window

# remove ABX_EMPr if they appear in ABX_EMP
empr <- empDF$ABX_EMPr
emp <- empDF$ABX_EMP
w <- which(lengths(empr) > 0)
for (i in w) {
   if (i %% 100 == 0) print(i)
   empri <- empr[[i]]
   empi <- emp[[i]]
   empr[[i]] <- empri[!empri %in% empi]
}
empDF$ABX_EMPr <- empr
rm(empr, emp, empri, empi, i, w)

w <- which(lengths(empDF$ABX_EMPr) == 0)
empDF$FLAGer[w] <- ''


plot(NA, xlim=c(-168, 432), ylim=c(1,2000))
for (i in 1:2000) {
   conc_times <- empDF$TIME_TO_CONC[[i]]
   if (length(conc_times) == 0) next
   points(x=conc_times, y=rep(i, length(conc_times)), pch=15, cex=0.25)
}
rm(i, conc_times)
abline(v = c(-48, -24))

# each row should have 'TIME_TO_FIRST_CONCORDANT_THERAPY' and 'TIME_TO_FIRST_ABX_THERAPY' columns
# currently these are lists that potentially includes times before -48 hours
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
empDF$TIME_TO_CONC <- time_to_conc

time_to_conc <- numeric(nrow(empDF))
w <- which(lengths(empDF$TIME_TO_CONC_STRICT) == 0L)
time_to_conc[w] <- NA
for (i in seq_len(nrow(empDF))[-w]) {
   t <- empDF$TIME_TO_CONC_STRICT[[i]]
   # only received concordant therapy pre-empirically, don't count this!
   if (all(t < -48)) {
      time_to_conc[i] <- NA
      next
   }
   # take only the times >= -48 hours
   time_to_conc[i] <- min(t[t >= -48])
}
empDF$TIME_TO_CONC_STRICT <- time_to_conc

#
time_to_abx <- numeric(nrow(empDF))
w <- which(lengths(empDF$TIME_TO_FIRST_ABX) == 0L | is.na(empDF$TIME_TO_FIRST_ABX))
time_to_abx[w] <- NA
for (i in seq_len(nrow(empDF))[-w]) {
   t <- empDF$TIME_TO_FIRST_ABX[[i]]
   # only received antibiotics pre-empirically, don't count
   if (all(t < -48)) {
      time_to_abx[i] <- NA
      next
   }
   # take only the times >= -48 hours
   time_to_abx[i] <- min(t[t >= -48])
}
empDF$TIME_TO_FIRST_ABX <- time_to_abx

#
time_to_abx <- numeric(nrow(empDF))
w <- which(lengths(empDF$TIME_TO_FIRST_ABX_STRICT) == 0L | is.na(empDF$TIME_TO_FIRST_ABX_STRICT))
time_to_abx[w] <- NA
for (i in seq_len(nrow(empDF))[-w]) {
   t <- empDF$TIME_TO_FIRST_ABX_STRICT[[i]]
   # only received antibiotics pre-empirically, don't count
   if (all(t < -48)) {
      time_to_abx[i] <- NA
      next
   }
   # take only the times >= -48 hours
   time_to_abx[i] <- min(t[t >= -48])
}
empDF$TIME_TO_FIRST_ABX_STRICT <- time_to_abx


# collapse empiric flags
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
save(empDF, file = '~/Desktop/EHR/EHR work/RdataFiles/ASTs_AbxAdmin_blood_2017_2023_imputed_encs_flagged_WIDE.Rdata')
gc()



########### WHY ARE SOME PATIENTS MISSING ENCOUNTER AND ANTIBIOTICS DATA ??? #############
# missing encoutners data?
# most are missing antibiotics
empDF %>% count(is.na(ADMIT_DAY)) # 3,651 missing
t <- empDF %>%
   mutate(noEnc = is.na(ADMIT_DAY),
          noAbx = lengths(ABX_EMPp) == 0 & lengths(ABX_EMPr) == 0 & lengths(ABX_EMP) == 0 & lengths(ABX_BTW) == 0 & lengths(ABX_TAR) == 0 & lengths(ABX_TARl) == 0) %>%
   select(noEnc, noAbx) %>% 
   table()
print(t)
print(fisher.test(t))

bugs <- names(head(sort(table(empDF$BUG), decreasing=TRUE), n=12))
for (bug in bugs) {
   print(bug)
   t <- empDF %>%
      filter(BUG == bug) %>%
      mutate(noEnc = is.na(ADMIT_DAY),
             noAbx = lengths(ABX_EMPp) == 0 & lengths(ABX_EMPr) == 0 & lengths(ABX_EMP) == 0 & lengths(ABX_BTW) == 0 & lengths(ABX_TAR) == 0 & lengths(ABX_TARl) == 0) %>%
      select(noEnc, noAbx) %>% 
      table()
   print(t)
   print(fisher.test(t))
   print('____________________________________________________')
}

empDF %>%
   mutate(noEnc = is.na(ADMIT_DAY),
          noTar = lengths(ABX_TAR) == 0,
          noEmp = lengths(ABX_EMP) == 0 & lengths(ABX_EMPr) == 0) %>%
   summarise(n = n(),
             noEnc = round(sum(noEnc) / n * 100, 1),
             noTarAbx = round(sum(noTar) / n * 100, 1),
             noEmp = round(sum(noEmp) / n * 100, 1),
             .by = BUG) %>%
   filter(n > 300) %>%
   arrange(noEnc)
################################################################################
############################ END COMBINE AST + ABX #############################
################################################################################








################################################################################
########################## ADD SURVIVAL / DEMO INFO ############################
################################################################################
library(dplyr)
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ASTs_AbxAdmin_blood_2017_2023_imputed_encs_flagged_WIDE.Rdata')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_DEMO.Rdata')

empDF <- empDF %>%
   left_join(x = .,
             y = dth,
             by = join_by(PERSON_ID),
             relationship = 'many-to-one') %>%
   mutate(AGE = as.integer(ORDER_DAY - DOB) / 365) %>%
   mutate(SURV_TIME = as.integer(DEATH_DATE - ORDER_DAY))

empDF <- empDF %>%
   relocate(ADMIT_DATE:SURV_TIME, .before=BUG)


empDF %>% count(is.na(GENDER), is.na(DOB)) # only 43 missing
empDF %>% count(AGE >= 18)                 # ~1,800 kids
empDF %>% count(is.na(SURV_TIME))          # ~17,600 have recorded deaths
empDF %>% count(SURV_TIME <= 2)            # ~1,900 died within 48 hours of blood culture
empDF %>% count(DEATH_DATE < ORDER_DAY)    # 8 died before order_day - nonsensical
empDF %>% count(DEATH_DATE < ADMIT_DAY)    # 6 died before admit_day - nonsensical
empDF %>% count(as.integer(substr(DEATH_DATE,1,4)) < 2017) # same 6

empDF <- empDF %>%
   filter(!is.na(GENDER)) %>%
   filter(AGE > 18) %>%
   filter(is.na(SURV_TIME) | SURV_TIME > 2)

save(empDF, file = '~/Desktop/EHR/EHR work/RdataFiles/ASTs_AbxAdmin_blood_2017_2023_imputed_encs_flagged_WIDE_survival.Rdata')












