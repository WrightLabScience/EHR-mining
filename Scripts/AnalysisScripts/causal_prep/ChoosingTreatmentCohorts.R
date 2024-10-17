library(dplyr)
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/data_path_name.Rdata')
load(file = paste0(data_path_name, 'ASTs_AbxAdmin_blood_2017_2023_imputed_flagged_WIDE_survival.Rdata'))

empDF %>% 
   filter(BUG == 'Enterococcus faecium', VANCOMYCIN == 1L) %>%
   count(grepl('DAPTOMYCIN', ABX_TAR) | grepl('DAPTOMYCIN', ABX_BTW),
         grepl('LINEZOLID', ABX_TAR) | grepl('LINEZOLID', ABX_BTW))

empDF %>% 
   filter(BUG == 'Enterococcus faecium', VANCOMYCIN == 1L) %>% 
   filter(!grepl('DAPTOMYCIN', ABX_TAR) & !grepl('DAPTOMYCIN', ABX_BTW) & !grepl('LINEZOLID', ABX_TAR) & !grepl('LINEZOLID', ABX_BTW)) %>%
   View()


empDF %>% 
   filter(BUG == 'Staphylococcus aureus', OXACILLIN == 1L) %>% 
   mutate(VAN = grepl('VANCOMYCIN', ABX_TAR) | grepl('VANCOMYCIN', ABX_BTW) | grepl('VANCOMYCIN', ABX_EMP),
          DAP = grepl('DAPTOMYCIN', ABX_TAR) | grepl('DAPTOMYCIN', ABX_BTW) | grepl('DAPTOMYCIN', ABX_EMP)) %>%
   count(VAN, DAP)


