# create expected S/R phenotype rules table from raw data

path <- '~/Desktop/EHR-mining/UsefulDataForCleaning/ASTimputation/EUCAST_expected_phenotypes.xlsx'
R1 <- readxl::read_xlsx(path, sheet='exp_res_1')
R2 <- readxl::read_xlsx(path, sheet='exp_res_2')
R3 <- readxl::read_xlsx(path, sheet='exp_res_3')
R4 <- readxl::read_xlsx(path, sheet='exp_res_4')
R5 <- readxl::read_xlsx(path, sheet='exp_res_5')

# combined resistance DF
r <- full_join(full_join(full_join(full_join(R1, R2), R3), R4), R5)

S1 <- readxl::read_xlsx(path, sheet='exp_susc_1')
S2 <- readxl::read_xlsx(path, sheet='exp_susc_2')
S3 <- readxl::read_xlsx(path, sheet='exp_susc_3')

# combined susceptibility DF
s <- full_join(full_join(S1, S2), S3)

# combined resistance/susceptibility DF
rs <- full_join(r, s)
rs <- rs %>% 
   tidyr::pivot_longer(cols = -BUG,
                       names_to = 'ABX',
                       values_to = 'RULE') %>%
   arrange(BUG, ABX)
rs <- rbind(rs %>% 
               filter(all(is.na(RULE)), 
                      .by=c(BUG, ABX)) %>% 
               distinct(), # both are NA
            rs %>% 
               filter(any(!is.na(RULE)), 
                      .by=c(BUG, ABX)) %>% 
               filter(!is.na(RULE))) # one is not NA -> keep only not NA rules
rs <- rs %>% tidyr::pivot_wider(names_from = ABX, values_from = RULE)
rm(r, s, R1, R2, R3, R4, R5, S1, S2, S3)

names(rs) <- toupper(names(rs))

rs$BUG <- sub("(.+)( species)", "\\1", rs$BUG)
rs$BUG <- sub("Coagulase Negative Staphylococcus", "Coagulase Negative Staph", rs$BUG)
rs$BUG <- sub("Group \\[ABCG] Beta-Hemolytic Streptococci", "Group \\[ABCG] Streptococci", rs$BUG)


rs <- rs %>% tidyr::pivot_longer(cols = -BUG, names_to = 'ABX', values_drop_na = TRUE) # 91 rows
rs <- rs %>% mutate(value=ifelse(value=='R', 1, 0))
enterobacterales <- c("Escherichia coli", "Klebsiella pneumoniae", "Proteus mirabilis",
                      "Enterobacter cloacae", "Klebsiella oxytoca", "Serratia marcescens",
                      'Enterobacter aerogenes',
                      "Klebsiella aerogenes", "Klebsiella variicola", 'Citrobacter freundii', 'Morganella morganii')
rs <- rbind(
   rs,
   tibble(
      BUG = enterobacterales,
      ABX = 'VANCOMYCIN',
      value = 1
   )
)
rs <- rbind(
   rs,
   tibble(
      BUG = c(enterobacterales, 'Pseudomonas aeruginosa', 'Stenotrophomonas maltophilia', 'Acinetobacter baumannii', 'Group G Streptococci',
              'Streptococcus pyogenes', 'Streptococcus agalactiae', 'Streptococcus mitis', 'Streptococcus anginosus', 'Viridans Streptococci'),
      ABX = 'LINEZOLID',
      value = 0
   )
)

rs <- rbind(
   rs,
   tibble(
      BUG = c(enterobacterales, 'Pseudomonas aeruginosa', 'Stenotrophomonas maltophilia', 'Acinetobacter baumannii', 
              'Streptococcus pneumoniae', 'Streptococcus agalactiae', 'Streptococcus pyogenes',
              'Staphylococcus aureus', 'Staphylococcus lugdunensis', 
              'Enterococcus faecium', 'Enterococcus faecalis'),
      ABX = 'DAPTOMYCIN',
      value = 0
   )
)

rm(enterobacterales)

van_S_bugs <- c('Streptococcus mitis', 'Viridans Streptococci', 'Streptococcus agalactiae', 'Streptococcus pyogenes', 'Streptococcus anginosus',
              'Streptococcus salivarius', 'Group G Streptococci', 'Streptococcus mitis/oralis group', 'Beta Hemolytic Streptococci',
              'Enterococcus faecium', 'Enterococcus faecalis')
rs <- rbind(
   rs,
   tibble(
      BUG = van_S_bugs,
      ABX = 'VANCOMYCIN',
      value = 0
   )
)

rm(van_S_bugs)

rm(path)


tmp_sxt_bugs <- c('Staphylococcus aureus', 'Staphylococcus lugdunensis', 
                  'Streptococcus mitis', 'Viridans Streptococci', 'Streptococcus agalactiae', 'Streptococcus pyogenes', 'Streptococcus anginosus',
                  'Streptococcus salivarius', 'Group G Streptococci', 'Streptococcus mitis/oralis group', 'Beta Hemolytic Streptococci', 'Streptococcus pneumoniae',
                  'Enterococcus faecium', 'Enterococcus faecalis', 'Escherichia coli')
rs <- rbind(
   rs,
   tibble(
      BUG = tmp_sxt_bugs,
      ABX = rep('TRIMETHOPRIM/SULFAMETHOXAZOLE', length(tmp_sxt_bugs)),
      value = rep(0, length(tmp_sxt_bugs))
   )
)
rm(tmp_sxt_bugs)

rs <- rbind(
   rs,
   tibble(
      BUG = 'Pseudomonas aeruginosa',
      ABX = 'TRIMETHOPRIM/SULFAMETHOXAZOLE',
      value = 1
   )
)

rs <- rbind(
   rs,
   tibble(
      BUG = rep(c('Enterococcus faecium', 'Enterococcus faecalis'), each=2),
      ABX = rep(c('MEROPENEM', 'ERTAPENEM'), times=2),
      value = 1
   )
)

rs <- distinct(rs)




