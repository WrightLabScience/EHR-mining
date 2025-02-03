StaphRules <- read.table(file = paste0(path, 'EUCAST Expert Rules v 3.2 - Staph Rules Condensed.tsv'), sep = "\t")
names(StaphRules) <- toupper(StaphRules[1,])
StaphRules <- StaphRules[-1,]
rownames(StaphRules) <- NULL
StaphRules[1, 4] <- "R if benzylpenicillin; NA" # tastorary
StaphRules <- StaphRules %>%
   tidyr::pivot_longer(
      cols = -BUG, 
      names_to = 'ABX', 
      values_drop_na = TRUE
   )

bl_list <- list(
   'PENICILLINS' = c('AMOXICILLIN', 'AMPICILLIN', 'BENZYLPENICILLIN', 'NAFCILLIN', 'METHICILLIN', ''),
   'ISOXAZOLYL_PENICILLINS' = c('OXACILLIN', 'CLOXACILLIN', 'DICLOXACILLIN', 'FLUCLOXACILLIN'),
   'PENICILLIN_BL_INH' = c('AMPICILLIN/SULBACTAM'),
   'CEPHOLOSPORINS_1_4' = c('CEFAZOLIN', 'CEPHALEXIN', 'CEFALEXIN', 'CEPHALOTHIN', 'CEFALOTHIN', 'CEFALOTIN', 'CEFAPIRIN', 'CEFRADINE', 'CEFADROXIL', 'CEFATRIZINE',
                        'CEFOXITIN', 'CEFUROXIME', 'CEFACLOR', 'CEFPROZIL', 'CEFMETAZOLE', 'CEFONICID', 'CEFOTETAN', 'CEFBUPERAZONE',
                        'CEFTAZIDIME', 'CEFTRIAXONE', 'CEFOTAXIME', 'CEFDINIR', 'CEFPODOXIME', 'CEFIXIME', 'CEFTIZOXIME', 'CEFOPERAZONE',
                        'CEFEPIME'),
   'CEPHOLOSPORINS_5' = c('CEFTAROLINE', 'CEFTOLOZANE', 'CEFTOBIPROLE')
)

StaphRules %>% filter(ABX %in% ISOXAZOLYL_PENICILLINS)

# MRSA is R to all beta_lactams except 5th gen cepholosporins

StrepRules <- read.table(file = paste0(path, 'EUCAST Expert Rules v 3.2 - Strep Rules Condensed.tsv'), sep = "\t")
names(StrepRules) <- toupper(StrepRules[1,])
StrepRules <- StrepRules[-1,]
rownames(StrepRules) <- NULL

rules <- full_join(StaphRules, StrepRules) # join Staph and Strep rules
rm(StaphRules, StrepRules, path)
rules$BUG[rules$BUG == 'Beta-Hemolytic streptococci'] <- 'Beta Hemolytic Streptococci'
rules$BUG[rules$BUG == 'Viridans group streptococci'] <- 'Viridans Streptococci'