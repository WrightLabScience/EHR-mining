abbr <- read.table(file = '~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/ABX_ABBR.txt', 
                   header=TRUE, 
                   sep=' ')
abbr$Antibiotic_Name <- gsub(pattern = '-', 
                             replacement = '/', 
                             x = abbr$Antibiotic_Name)
abbr$Antibiotic_Name <- gsub(pattern = '_',
                             replacement = ' ',
                             x = abbr$Antibiotic_Name)
abbr$Antibiotic_Name <- toupper(abbr$Antibiotic_Name)
abbr <- setNames(
   abbr$Abbreviation,
   abbr$Antibiotic_Name
)
