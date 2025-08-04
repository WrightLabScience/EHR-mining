library(dplyr)
imp_rules <- readxl::read_xlsx(path = '~/Desktop/EHR-mining/UsefulDataForCleaning/ASTimputation/ImputationRulesFromRheeEtAl2020.xlsx')
imp_rules <- tidyr::pivot_longer(data = imp_rules,
                          cols = !Antibiotic,
                          names_to = 'BUG',
                          values_to = 'RULE') %>% 
   arrange(BUG) %>%
   mutate(across(everything(), toupper)) %>%
   mutate(Antibiotic = gsub('-', '/', Antibiotic)) %>%
   mutate(RULE = gsub('-S', '', RULE)) %>%
   mutate(RULE = gsub('-', '/', RULE)) %>%
   mutate(RULE = gsub('IF MISSING THEN ', '', RULE)) %>%
   mutate(RULE = gsub('ELSE ', '', RULE)) %>%
   mutate(RULE = gsub('\\.', 'NA', RULE)) %>%
   mutate(RULE = gsub(' OR ', '; S IF ', RULE)) %>%
   mutate(BUG = ifelse(grepl(' SPECIES$', BUG),
                       yes = gsub(' SPECIES$', '', BUG),
                       no = BUG)) %>%
   mutate(BUG = case_when(
      BUG == 'ALPHA-HEMOLYTIC / VIRIDANS STREPTOCOCCUS' ~ 'ALPHA HEMOLYTIC STREPTOCOCCI|VIRIDANS STREPTOCOCCI',
      BUG == 'BETA-HEMOLYTIC STREPTOCOCCUS' ~ 'BETA HEMOLYTIC STREPTOCOCCI',
      .default = BUG
   )) %>%
   tidyr::pivot_wider(names_from = BUG, 
                      values_from = RULE) %>%
   select(-`PENICILLIN-SENSITIVE STAPHYLOCOCCUS AUREUS`)
imp_rules$Antibiotic[imp_rules$Antibiotic == 'TMP/SMX'] <- 'TRIMETHOPRIM/SULFAMETHOXAZOLE'
imp_rules$Antibiotic[imp_rules$Antibiotic == 'QUINUPRISTIN/DALFOPRISTIN'] <- 'SYNERCID'
names(imp_rules) <- stringr::str_to_sentence(names(imp_rules))
