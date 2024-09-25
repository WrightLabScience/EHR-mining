# The purpose of this script is to combine AST and AbxAdmin data (for bloodstream infections):
#     -Comine and save the long format for understanding abx admin proximity to order
#     -Collapse to one row per isolate, save this for future analyses
#     -Impute missing AST values
#     -Compute concordance/discordance flag (both on imputed and un-imputed tables)

# load necessary data
library(dplyr)
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/data_path_name.Rdata')




################################################################################
########################### START COMBINE AST + ABX ############################
################################################################################
start <- Sys.time()
library(dplyr)
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/data_path_name.Rdata')
load(file = paste0(data_path_name, 'ASTs_blood_2015_2023.Rdata'))
load(file = paste0(data_path_name, 'ALL_CLEANED_2017_2023_AbxAdmin.Rdata'))
astDF <- astDF %>% filter(substr(ORDER_DATE, 1, 4) %in% as.character(2017:2023)) #    55,959 -->    43,787
abxDF <- abxDF %>% filter(PERSON_ID %in% unique(astDF$PERSON_ID))                # 9,681,748 --> 2,588,330
length(unique(astDF$PERSON_ID)) # 35,666
length(unique(abxDF$PERSON_ID)) # 33,638
length(intersect(astDF$PERSON_ID, abxDF$PERSON_ID)) # 33,638
# prep data
abxDF <- abxDF %>%
   mutate(PERSON_ID = as.character(PERSON_ID)) %>%
   rename(START_DATE = ADMIN_START_DATE,
          END_DATE = ADMIN_END_DATE) %>%
   arrange(PERSON_ID, START_DATE)
astDF <- astDF %>%
   group_by(PERSON_ID, ORDER_DAY) %>%
   mutate(MULT_BLOOD_ISO = n() > 1L) %>%
   ungroup() %>%
   relocate(MULT_BLOOD_ISO, .before=BUG) %>%
   mutate(ORDER_DAY = as.Date(substr(ORDER_DATE,1,10)),
          RESULT_DAY = as.Date(substr(RESULT_DATE,1,10))) %>%
   mutate(CONTAM = case_when(
      NEXT_RESULT_DAY <= RESULT_DAY ~ 2,
      NEXT_ORDER_DAY < RESULT_DAY ~ 1
   )) %>%
   select(-NEXT_ORDER_DAY, -NEXT_RESULT_DAY) %>%
   relocate(RESULT_DAY, .after=ORDER_DAY) %>%
   relocate(ORDER_DATE, RESULT_DATE, .before=ORDER_DAY) %>%
   relocate(BUG, .after=CONTAM)

nrow(astDF) # 43,787
astDF %>% count(MULT_BLOOD_ISO) # 10,825
astDF %>% group_by(PERSON_ID, ORDER_DAY) # 37,976

## DIFFERENCES BETWEEN ORDER_TIMES for multiple isolates
{
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
}

# JOIN
empDF <- astDF %>%
   mutate(JOIN_START = ORDER_DAY - 30,
          JOIN_END = RESULT_DAY + 90) %>%
   left_join(x = .,
             y = abxDF,
             by = join_by(PERSON_ID,
                          JOIN_START <= START_DATE,
                          JOIN_END >= START_DATE)) %>%
   select(-JOIN_START, -JOIN_END) %>%
   relocate(ABX, START_DATE, END_DATE, .before=BUG)
# save long format
save(empDF, file = paste0(data_path_name, 'ASTs_AbxAdmin_blood_2017_2023_LONG.Rdata'))

# collapse so one row per isolate
empDF <- empDF %>% 
   select(-BLOOD) %>% 
   mutate(START_DAY = as.Date(substr(START_DATE,1,10))) %>%
   distinct() %>%
   relocate(CONTAM, ABX, START_DATE, END_DATE, BUG, START_DAY, .before=CEFEPIME) %>%
   group_by_all() %>% ungroup(ABX, START_DATE, END_DATE, START_DAY) %>%
   reframe(DATE_OF_MOST_RECENT_ABX = max(START_DAY[START_DAY < ORDER_DAY]),
           EMPIRIC0 = list(sort(unique(ABX[START_DAY == ORDER_DAY]))),
           EMPIRIC1 = list(sort(unique(ABX[START_DAY == (ORDER_DAY+1)]))),
           EMPIRIC2 = list(sort(unique(ABX[START_DAY == (ORDER_DAY+2) & START_DAY < RESULT_DAY]))),
           EMPIRIC3 = list(sort(unique(ABX[START_DAY == (ORDER_DAY+3) & START_DAY < RESULT_DAY]))),
           EMPIRIC4 = list(sort(unique(ABX[START_DAY == (ORDER_DAY+4) & START_DAY < RESULT_DAY]))),
           EMPIRIC5 = list(sort(unique(ABX[START_DAY == (ORDER_DAY+5) & START_DAY < RESULT_DAY]))),
           EMPIRIC6 = list(sort(unique(ABX[START_DAY == (ORDER_DAY+6) & START_DAY < RESULT_DAY]))),
           EMPIRICR = list(sort(unique(ABX[START_DAY > (ORDER_DAY+6) & START_DAY < RESULT_DAY]))),
           TARGETED = list(sort(unique(ABX[START_DAY == RESULT_DAY | START_DAY == RESULT_DAY + 1])))) %>%
   ungroup() %>%
   mutate(DELAY = as.integer(RESULT_DAY - ORDER_DAY)) %>%
   relocate(DELAY, DATE_OF_MOST_RECENT_ABX, EMPIRIC0, EMPIRIC1, EMPIRIC2, EMPIRIC3, EMPIRIC4, EMPIRIC5, EMPIRIC6, EMPIRICR, TARGETED, .before=BUG)
empDF$DATE_OF_MOST_RECENT_ABX[is.infinite(empDF$DATE_OF_MOST_RECENT_ABX)] <- NA
empDF$DAYS_SINCE_MOST_RECENT_ABX <- as.integer(empDF$ORDER_DAY - empDF$DATE_OF_MOST_RECENT_ABX)
empDF <- empDF %>% select(-DATE_OF_MOST_RECENT_ABX)
empDF <- empDF %>% relocate(DAYS_SINCE_MOST_RECENT_ABX, .before=EMPIRIC0)

# antibiotic abbreviations
abbr <- tibble(read.table(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/antibiotic_names/ABX_ABBR.txt', header = TRUE)) %>% 
   select(-Class) %>%
   mutate(Antibiotic_Name = gsub('-', '/', Antibiotic_Name),
          Antibiotic_Name = gsub('_', ' ', Antibiotic_Name),
          Antibiotic_Name = toupper(Antibiotic_Name))
abbr <- setNames(abbr$Abbreviation, abbr$Antibiotic_Name)
abad <- unique(unlist(empDF %>% select(EMPIRIC0:TARGETED)))
abad <- abad[!abad %in% names(abbr)]
abbr <- c(abbr, setNames(abad, abad))
rm(abad)
save(abbr, file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/antibiotic_names/abbr_named.Rdata')

empDF$EMP0_ABBR <- sapply(empDF$EMPIRIC0, function(x) paste(sort(abbr[x]), collapse='+'))
empDF$TARG_ABBR <- sapply(empDF$TARGETED, function(x) paste(sort(abbr[x]), collapse='+'))
empDF <- empDF %>% relocate(EMP0_ABBR, TARG_ABBR, .before=EMPIRIC0)

# save wide format
save(empDF, file = paste0(data_path_name, 'ASTs_AbxAdmin_blood_2017_2023.Rdata'))
print(Sys.time() - start) # < 2 minutes

rm(abxDF, astDF, start, abbr, empDF)
gc()
################################################################################
############################ END COMBINE AST + ABX #############################
################################################################################






################################################################################
############################### START IMPUTATION ###############################
################################################################################
start <- Sys.time()
library(dplyr)
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/data_path_name.Rdata')
load(file = paste0(data_path_name, 'ASTs_AbxAdmin_blood_2017_2023.Rdata'))
Viridans <- c("Viridans Streptococci", "Alpha Hemolytic Streptococci", "Streptococcus mitis", "Streptococcus mutans", "Streptococcus sanguinis", "Streptococcus mitis/oralis group", "Streptococcus salivarius")
BetaHemolytic <- c("Beta Hemolytic Streptococci", "Group A Streptococci", "Group B Streptococci", "Group C Streptococci", "Group G Streptococci", "Streptococcus dysgalactiae", "Streptococcus pyogenes", 
                   "Streptococcus agalactiae", "Streptococcus equi", "Streptococcus gallolyticus", "Streptococcus bovis", "Streptococcus equinis", "Streptococcus suis")

# KADRI imputation rules
empCOPY <- empDF
{
   source(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/ASTimputation/ImputationRulesClean.R')
   abx_names <- names(empDF)[16:length(empDF)]
   imp_rules$Antibiotic[!unique(imp_rules$Antibiotic) %in% abx_names] # CEFADROXIL, CEFPODOXIME, CEPHALEXIN, DICLOXACILLIN, NAFCILLIN
   abx_names[!abx_names %in% imp_rules$Antibiotic] # lots
   rm(abx_names)
   
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
   sum(lengths(bug_rows)) / nrow(empDF) # 66.8%
   u_bug_rows <- sort(unlist(bug_rows))
   names(u_bug_rows) <- NULL
   any(diff(u_bug_rows) == 0) # no overlap!
   head(sort(table(empDF$BUG[-u_bug_rows]), decreasing = TRUE), n=12) # frequent bugs for which we don't have imputation rules
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
   #print(Sys.time() - start) # ~  minutes
   rm(a, abx, b, cl, r, rws, val, w, op, imp_rules, bug_rows)
   
   # check for errors
   which(sapply(empCOPY[16:length(empCOPY)], function(x) 4 %in% x))
}

# EUCAST expected phenotypes
empTEST <- empCOPY
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
   
   
   bug_rows <- sapply(rs$ORGANISMS, function(x) grep(x, empTEST$BUG))
   sort(lengths(bug_rows))
   sum(lengths(bug_rows)) / nrow(empDF) # 70%
   u_bug_rows <- sort(unlist(bug_rows))
   head(sort(table(empDF$BUG[-u_bug_rows]), decreasing = TRUE), n=12) # frequent bugs for which we don't have imputation rules
   rm(u_bug_rows)
   
   
   #start <- Sys.time()
   for (b in 1:nrow(rs)){ #bug by bug
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
   #end <- Sys.time()
   #print(end - start) # 1.3 seconds
   rm(a, b, cl, end, rws, val, rs, bug_rows)
}

# EUCAST expert rules
empTEST2 <- empTEST
{
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
      
      if(char == "Staphylococcus") return(grep("Staph", unique(empDF$BUG), value = T))
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
   
   bug_rows <- sapply(rules$BUG, function(x) grep(x, empTEST2$BUG))
   lengths(bug_rows)
   
   #start <- Sys.time()
   for (b in 1:nrow(rules)) { #bug/species by bug/species in rules
      cat(b, rules$BUG[b], '\n')
      for (a in 2:ncol(rules)){ #drug by drug in rules
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
            # cat(a, r, '\n')
            empTEST2[r, cl[is.na(empTEST2[r, cl])]] <- op(val, r)     #run value and ast row nums through op function
         }
      }
   }
   #print(Sys.time() - start) # ~2.5 minutes
   rm(a, AMINOPENICILLINS, b, CARBAPENEMS, CEPCAR, CEPHALOSPORINS, CEPHg1, CEPHg2, CEPHg3, CEPHg4, CEPHg5, cl, FLUOROQUINOLONES,
      LINCOSAMIDES, MACLIN, MACROLIDES, r, rws, val, n, op, bug_rows, rules)
}

# check the difference between unimputed and imputed
df <- data.frame(old = sapply(empDF,    function(x) sum(!is.na(x))), new = sapply(empCOPY, function(x) sum(!is.na(x))))
df$diff <- df$new - df$old
df <- df[order(df$diff, decreasing=TRUE), ]

empDF <- empTEST2
save(empDF, file = '~/Desktop/EHR/EHR work/RdataFiles/ASTs_AbxAdmin_blood_2017_2023_imputed.Rdata')
print(Sys.time() - start) # < 4 minutes

rm(empCOPY, empTEST, empTEST2, df, Viridans, BetaHemolytic, start, empDF)
gc()
################################################################################
################################ END IMPUTATION ################################
################################################################################






################################################################################
########################## START CONCORDANCE ANALYSIS ##########################
################################################################################
start <- Sys.time()
library(dplyr)
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/data_path_name.Rdata')
load(file = paste0(data_path_name, 'ASTs_AbxAdmin_blood_2017_2023.Rdata')); empDFo <- empDF
load(file = paste0(data_path_name, 'ASTs_AbxAdmin_blood_2017_2023_imputed.Rdata'))

empDF$FLAG0  <- empDF$FLAG1  <- empDF$FLAG2  <- empDF$FLAG3  <- empDF$FLAG4  <- empDF$FLAG5  <- empDF$FLAG6  <- empDF$FLAGR  <- empDF$FLAGT  <- NA
empDFo$FLAG0 <- empDFo$FLAG1 <- empDFo$FLAG2 <- empDFo$FLAG3 <- empDFo$FLAG4 <- empDFo$FLAG5 <- empDFo$FLAG6 <- empDFo$FLAGR <- empDFo$FLAGT <- NA
empDF <- empDF %>% relocate(FLAG0, FLAG1, FLAG2, FLAG3, FLAG4, FLAG5, FLAG6, FLAGR, FLAGT, .before=EMPIRIC0)
empDF <- empDF %>% relocate(FLAG0, FLAG1, FLAG2, FLAG3, FLAG4, FLAG5, FLAG6, FLAGR, FLAGT, .before=EMPIRIC0)
col_names <- names(empDF %>% select(CEFEPIME:DELAFLOXACIN))
emp0 <- empDF$EMPIRIC0
emp1 <- empDF$EMPIRIC1
emp2 <- empDF$EMPIRIC2
emp3 <- empDF$EMPIRIC3
emp4 <- empDF$EMPIRIC4
emp5 <- empDF$EMPIRIC5
emp6 <- empDF$EMPIRIC6
empR <- empDF$EMPIRICR
targ <- empDF$TARGETED
source('~/Desktop/EHR/EHR-mining/Scripts/AnalysisScripts/ConcordanceFlagFunction.R')

for (i in seq_len(nrow(empDF))) {
   print(i)
   empDF$FLAG0[i] <- getConcordanceFlag(empDF[i,], emp0[[i]])
   empDF$FLAG1[i] <- getConcordanceFlag(empDF[i,], emp1[[i]])
   empDF$FLAG2[i] <- getConcordanceFlag(empDF[i,], emp2[[i]])
   empDF$FLAG3[i] <- getConcordanceFlag(empDF[i,], emp3[[i]])
   empDF$FLAG4[i] <- getConcordanceFlag(empDF[i,], emp4[[i]])
   empDF$FLAG5[i] <- getConcordanceFlag(empDF[i,], emp5[[i]])
   empDF$FLAG6[i] <- getConcordanceFlag(empDF[i,], emp6[[i]])
   empDF$FLAGR[i] <- getConcordanceFlag(empDF[i,], empR[[i]])
   empDF$FLAGT[i] <- getConcordanceFlag(empDF[i,], targ[[i]])
   
   empDFo$FLAG0[i] <- getConcordanceFlag(empDFo[i,], emp0[[i]])
   empDFo$FLAG1[i] <- getConcordanceFlag(empDFo[i,], emp1[[i]])
   empDFo$FLAG2[i] <- getConcordanceFlag(empDFo[i,], emp2[[i]])
   empDFo$FLAG3[i] <- getConcordanceFlag(empDFo[i,], emp3[[i]])
   empDFo$FLAG4[i] <- getConcordanceFlag(empDFo[i,], emp4[[i]])
   empDFo$FLAG5[i] <- getConcordanceFlag(empDFo[i,], emp5[[i]])
   empDFo$FLAG6[i] <- getConcordanceFlag(empDFo[i,], emp6[[i]])
   empDFo$FLAGR[i] <- getConcordanceFlag(empDFo[i,], empR[[i]])
   empDFo$FLAGT[i] <- getConcordanceFlag(empDFo[i,], targ[[i]])
}
print(Sys.time() - start) # ~6 minutes

save(empDF, file = paste0(data_path_name, 'ASTs_AbxAdmin_blood_2017_2023_imputed_flagged.Rdata'))
save(empDFo, file = paste0(data_path_name, 'ASTs_AbxAdmin_blood_2017_2023_flagged.Rdata'))

rm(getConcordanceFlag, start, i, col_names, emp0, emp1, emp2, emp3, emp4, emp5, emp6, empR, targ, empDF, empDFo)
gc()
################################################################################
########################### END CONCORDANCE ANALYSIS ###########################
################################################################################


rm(data_path_name)



