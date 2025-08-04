library(dplyr)
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/data_path_name.Rdata')
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/plots_path_name.Rdata')
plots_path_name <- paste0(plots_path_name, 'BloodAbxAdmin/')
load(file = paste0(data_path_name, 'ASTs_blood_2017_2023_imputed.Rdata'))
astDFi <- astDF
load(file = paste0(data_path_name, 'ASTs_blood_2015_2023.Rdata'))
load(file = paste0(data_path_name, 'ALL_CLEANED_2017_2023_AbxAdmin.Rdata'))

# only have medication administration data from 2017 - 2023
# prep data for abx duration analysis
astDF <- astDF %>% filter(substr(ORDER_DATE, 1, 4) %in% as.character(2017:2023)) # 55,959 --> 43,787
astDFi <- astDFi %>% filter(substr(ORDER_DATE, 1, 4) %in% as.character(2017:2023)) # 55,959 --> 43,787
abxDF <- abxDF %>% 
   filter(PERSON_ID %in% unique(astDF$PERSON_ID)) %>%                # 9,681,748 --> 2,588,330
   mutate(PERSON_ID = as.character(PERSON_ID)) %>%
   rename(START_DATE = ADMIN_START_DATE,
          END_DATE = ADMIN_END_DATE) %>%
   arrange(PERSON_ID, START_DATE)


# LOOK AT SPECIFIC CULTURES AND THEIR SURROUNDING ABX ADMINISTRATIONS TO GET A SENSE
{
   # S. aureus bacteremia 2/5 --> 2/8 (and Proteus mirabilis elsewhere 2/5 --> 2/7)
   # received PIP/TAZO, VAN, and CEFAZOLIN empirically
   # CEFAZOLIN only once first day, then resumed on result_day, after it was revealed that it was MSSA
   astDF %>% slice(2) %>% select(ORDER_DAY, RESULT_DAY, MULT_ISO, BUG, CEFAZOLIN, `PIPERACILLIN/TAZOBACTAM`, VANCOMYCIN, OXACILLIN, CEFOXITIN, AMPICILLIN, PENICILLIN)
   x <- abxDF %>% filter(PERSON_ID == astDF$PERSON_ID[2])
   x %>%
      filter(substr(START_DATE,1,4) == '2018') %>%
      mutate(across(c(START_DATE, END_DATE), ~ as.Date(substr(.,1,10))),
             RX = 1) %>%
      slice(5:58) %>%
      distinct() %>% tidyr::pivot_wider(id_cols = START_DATE, names_from = ABX, values_from = RX)
   
   # E. coli 5/6 --> 5/10 (and something else)
   # received PIP/TAZO, VAN emperically
   # VAN stopped after 4th day (must have recevied gram-negative stain results at that time)
   # after result day, PIP/TAZO stopped and ampicillin began (amoxicillin used a few days later - might be for other infection)
   astDF %>% slice(3) %>% select(ORDER_DAY, RESULT_DAY, MULT_ISO, BUG, VANCOMYCIN, `PIPERACILLIN/TAZOBACTAM`, AMPICILLIN, AMOXICILLIN)
   x <- abxDF %>% filter(PERSON_ID == astDF$PERSON_ID[3]) %>% filter(substr(START_DATE,1,4) == substr(astDF$ORDER_DAY[3],1,4))
   x %>%
      select(-END_DATE) %>%
      mutate(START_DATE = as.Date(substr(START_DATE,1,10)),
             RX = 1) %>%
      distinct() %>% tidyr::pivot_wider(id_cols = START_DATE, names_from = ABX, values_from = RX)
   
   # P. aeruginosa 10/10 --> 10/14
   # received PIP/TAZO and VAN 1 day after empiric, VAN only second day, AMP/SUL on 3rd day
   i <- 4:5; astDF %>% slice(i) %>% select(ORDER_DAY, RESULT_DAY, MULT_ISO, BUG, VANCOMYCIN, `PIPERACILLIN/TAZOBACTAM`, `AMPICILLIN/SULBACTAM`)
   x <- abxDF %>% filter(PERSON_ID == astDF$PERSON_ID[i]) %>% filter(substr(START_DATE,1,4) == substr(astDF$ORDER_DAY[i[1]],1,4))
   x %>%
      select(-END_DATE) %>%
      mutate(START_DATE = as.Date(substr(START_DATE,1,10)),
             RX = 1) %>%
      distinct() %>% tidyr::pivot_wider(id_cols = START_DATE, names_from = ABX, values_from = RX)
   
   
   # S. aureus 1/12 --> 1/18 (and another)
   # received PIP/TAZO day before empiric, VAN empirically
   # Kept both through result day, finally learned it was MSSA, but must have already been cured?
   i <- 6; astDF %>% slice(i) %>% select(ORDER_DAY, RESULT_DAY, MULT_ISO, BUG, VANCOMYCIN, `PIPERACILLIN/TAZOBACTAM`, OXACILLIN, CEFOXITIN, AMPICILLIN, PENICILLIN)
   x <- abxDF %>% filter(PERSON_ID == astDF$PERSON_ID[i]) %>% filter(substr(START_DATE,1,4) == substr(astDF$ORDER_DAY[i],1,4))
   x %>%
      select(-END_DATE) %>%
      mutate(START_DATE = as.Date(substr(START_DATE,1,10)),
             RX = 1) %>%
      distinct() %>% tidyr::pivot_wider(id_cols = START_DATE, names_from = ABX, values_from = RX)
   
   # Coag neg staph 8/22 --> 8/27 (oxa susceptible)
   # VAN 8/24-8/26, OXA 8/27-8/31, VAN 8/31-9/1
   i <- 7; astDF %>% slice(i) %>% select(PERSON_ID, ORDER_DAY, RESULT_DAY, MULT_ISO, BUG, VANCOMYCIN, `PIPERACILLIN/TAZOBACTAM`, OXACILLIN, CEFOXITIN, AMPICILLIN, PENICILLIN)
   x <- abxDF %>% filter(PERSON_ID == astDF$PERSON_ID[i]) %>% filter(substr(START_DATE,1,4) == substr(astDF$ORDER_DAY[i],1,4))
   x %>%
      select(-END_DATE) %>%
      mutate(START_DATE = as.Date(substr(START_DATE,1,10)),
             RX = 1) %>%
      distinct() %>% tidyr::pivot_wider(id_cols = START_DATE, names_from = ABX, values_from = RX)
   
   # Group A Strep 8/24 --> 8/27 (solo)
   # empiric: FEP, VAN
   # next day: swtich VAN for DAP, continue FEP for several days
   # result_day: continute FEP, add MET and CLI
   # few days later: CRO, little VAN
   # many days later, VAN only for ~ week
   # then another culture happened 2 months later
   i <- 8; astDF %>% slice(i) %>% select(PERSON_ID, ORDER_DAY, RESULT_DAY, NEXT_ORDER_DAY, NEXT_RESULT_DAY, MULT_ISO, BUG, VANCOMYCIN, `PIPERACILLIN/TAZOBACTAM`, OXACILLIN, CEFOXITIN, AMPICILLIN, PENICILLIN)
   x <- abxDF %>% filter(PERSON_ID == astDF$PERSON_ID[i]) %>% filter(substr(START_DATE,1,4) == substr(astDF$ORDER_DAY[i],1,4))
   x %>%
      select(-END_DATE) %>%
      mutate(START_DATE = as.Date(substr(START_DATE,1,10)),
             RX = 1) %>%
      distinct() %>% tidyr::pivot_wider(id_cols = START_DATE, names_from = ABX, values_from = RX) %>% print(n=50)
   
   # E. coli 4/8 --> 4/10 (solo)
   # Just PIP/TAZO starting order_day, through 9 days (1 week after result)
   i <- 9; astDF %>% slice(i) %>% select(PERSON_ID, ORDER_DAY, RESULT_DAY, NEXT_ORDER_DAY, NEXT_RESULT_DAY, MULT_ISO, BUG, VANCOMYCIN, `PIPERACILLIN/TAZOBACTAM`, OXACILLIN, CEFOXITIN, AMPICILLIN, PENICILLIN)
   x <- abxDF %>% filter(PERSON_ID == astDF$PERSON_ID[i]) %>% filter(substr(START_DATE,1,4) == substr(astDF$ORDER_DAY[i],1,4))
   x %>%
      select(-END_DATE) %>%
      mutate(START_DATE = as.Date(substr(START_DATE,1,10)),
             RX = 1) %>%
      distinct() %>% tidyr::pivot_wider(id_cols = START_DATE, names_from = ABX, values_from = RX) %>% print(n=50)
   
   # S. aureus (MRSA) 6/25 --> 6/30
   # starting 6/27: VAN alone
   # 7/4: add PIP/TAZO
   # 7/8: add CRO and MET, stop PIP/TAZO
   # VAN continued until 8/19 (stopped for a week 7/28-8/5)
   i <- which(astDF$BUG == 'Staphylococcus aureus' & astDF$OXACILLIN == 1L)[1] # 18
   astDF %>% slice(i) %>% select(PERSON_ID, ORDER_DAY, RESULT_DAY, NEXT_ORDER_DAY, NEXT_RESULT_DAY, MULT_ISO, BUG, VANCOMYCIN, `PIPERACILLIN/TAZOBACTAM`, OXACILLIN, CEFOXITIN, AMPICILLIN, PENICILLIN)
   x <- abxDF %>% filter(PERSON_ID == astDF$PERSON_ID[i]) %>% filter(substr(START_DATE,1,4) == substr(astDF$ORDER_DAY[i],1,4))
   x %>%
      select(-END_DATE) %>%
      mutate(START_DATE = as.Date(substr(START_DATE,1,10)),
             RX = 1) %>%
      distinct() %>% tidyr::pivot_wider(id_cols = START_DATE, names_from = ABX, values_from = RX) %>% print(n=50)
   x %>% filter(ABX == 'VANCOMYCIN') %>% print(n=50)
   
   
   # E. faecium 3/10/23 --> 3/13/23 (solo)
   # on PIP/TAZO, FLUCONZAOLE, some VAN before culture sampling
   # stayed on PIP/TAZO through 1 day after results
   # given 2 days VAN between order and result
   i <- which(astDF$BUG == 'Enterococcus faecium' & astDF$VANCOMYCIN == 1L)[1] # 247
   astDF %>% slice(i) %>% select(PERSON_ID, ORDER_DAY, RESULT_DAY, NEXT_ORDER_DAY, NEXT_RESULT_DAY, MULT_ISO, BUG, VANCOMYCIN, `PIPERACILLIN/TAZOBACTAM`, OXACILLIN, CEFOXITIN, AMPICILLIN, PENICILLIN)
   astDFi%>% slice(i) %>% select(PERSON_ID, ORDER_DAY, RESULT_DAY, NEXT_ORDER_DAY, NEXT_RESULT_DAY, MULT_ISO, BUG, VANCOMYCIN, `PIPERACILLIN/TAZOBACTAM`, OXACILLIN, CEFOXITIN, AMPICILLIN, PENICILLIN)
   x <- abxDF %>% filter(PERSON_ID == astDF$PERSON_ID[i]) %>% filter(substr(START_DATE,1,4) == substr(astDF$ORDER_DAY[i],1,4))
   x %>%
      select(-END_DATE) %>%
      mutate(START_DATE = as.Date(substr(START_DATE,1,10)),
             RX = 1) %>%
      distinct() %>% tidyr::pivot_wider(id_cols = START_DATE, names_from = ABX, values_from = RX) %>% print(n=50)
   
   # E. faecium 2/24/17 --> 2/26/17 (solo)
   # on PIP/TAZO, MET, some CFZ before culture sampling
   # PIP/TAZO, AZM, VAN given empirically, stopped amin on 2/25 (before AST results!)
   i <- which(astDF$BUG == 'Enterococcus faecium' & astDF$VANCOMYCIN == 1L)[2] # 478
   astDF %>% slice(i) %>% select(PERSON_ID, ORDER_DAY, RESULT_DAY, NEXT_ORDER_DAY, NEXT_RESULT_DAY, MULT_ISO, BUG, VANCOMYCIN, `PIPERACILLIN/TAZOBACTAM`, OXACILLIN, CEFOXITIN, AMPICILLIN, PENICILLIN)
   astDFi%>% slice(i) %>% select(PERSON_ID, ORDER_DAY, RESULT_DAY, NEXT_ORDER_DAY, NEXT_RESULT_DAY, MULT_ISO, BUG, VANCOMYCIN, `PIPERACILLIN/TAZOBACTAM`, OXACILLIN, CEFOXITIN, AMPICILLIN, PENICILLIN)
   x <- abxDF %>% filter(PERSON_ID == astDF$PERSON_ID[i]) %>% filter(substr(START_DATE,1,4) == substr(astDF$ORDER_DAY[i],1,4))
   x %>%
      select(-END_DATE) %>%
      mutate(START_DATE = as.Date(substr(START_DATE,1,10)),
             RX = 1) %>%
      distinct() %>% tidyr::pivot_wider(id_cols = START_DATE, names_from = ABX, values_from = RX) %>% print(n=50)
   
   # E. faecium 12/6/22 --> 12/12/22 (solo)
   # on LVX, FLU, FEP, etc. before culture sampling
   # PIP/TAZO, FLU, DAP, CEFT, given empirically
   # MER, DAP, CEFT between result and order (PIP/TAZO stopped)
   # MER stoppe,d LVX, FLU started days after AST results
   i <- which(astDF$BUG == 'Enterococcus faecium' & astDF$VANCOMYCIN == 1L)[3] # 516
   astDF %>% slice(i) %>% select(PERSON_ID, ORDER_DAY, RESULT_DAY, NEXT_ORDER_DAY, NEXT_RESULT_DAY, MULT_ISO, BUG, VANCOMYCIN, `PIPERACILLIN/TAZOBACTAM`, OXACILLIN, CEFOXITIN, AMPICILLIN, PENICILLIN)
   astDFi%>% slice(i) %>% select(PERSON_ID, ORDER_DAY, RESULT_DAY, NEXT_ORDER_DAY, NEXT_RESULT_DAY, MULT_ISO, BUG, VANCOMYCIN, `PIPERACILLIN/TAZOBACTAM`, OXACILLIN, CEFOXITIN, AMPICILLIN, PENICILLIN)
   x <- abxDF %>% filter(PERSON_ID == astDF$PERSON_ID[i]) %>% filter(substr(START_DATE,1,4) == substr(astDF$ORDER_DAY[i],1,4))
   x %>%
      select(-END_DATE) %>%
      mutate(START_DATE = as.Date(substr(START_DATE,1,10)),
             RX = 1) %>%
      distinct() %>% tidyr::pivot_wider(id_cols = START_DATE, names_from = ABX, values_from = RX) %>% print(n=50)
   
   
   # E. faecium 11/9/23 --> 11/12/23 (solo)
   # on many things in the months/weeks leading up to culture
   # CRO, MET, VAN in 1 week leading up
   # MET, FEP, PIP/TAZO, DAP, CASPOFUNGIN given empiric/next day??
   i <- which(astDF$BUG == 'Enterococcus faecium' & astDF$VANCOMYCIN == 1L)[4] # 635
   astDF %>% slice(i) %>% select(PERSON_ID, ORDER_DAY, RESULT_DAY, NEXT_ORDER_DAY, NEXT_RESULT_DAY, MULT_ISO, BUG, VANCOMYCIN, `PIPERACILLIN/TAZOBACTAM`, OXACILLIN, CEFOXITIN, AMPICILLIN, PENICILLIN)
   astDFi%>% slice(i) %>% select(PERSON_ID, ORDER_DAY, RESULT_DAY, NEXT_ORDER_DAY, NEXT_RESULT_DAY, MULT_ISO, BUG, VANCOMYCIN, `PIPERACILLIN/TAZOBACTAM`, OXACILLIN, CEFOXITIN, AMPICILLIN, PENICILLIN)
   x <- abxDF %>% filter(PERSON_ID == astDF$PERSON_ID[i]) %>% filter(substr(START_DATE,1,4) == substr(astDF$ORDER_DAY[i],1,4))
   x %>%
      select(-END_DATE) %>%
      mutate(START_DATE = as.Date(substr(START_DATE,1,10)),
             RX = 1) %>%
      distinct() %>% tidyr::pivot_wider(id_cols = START_DATE, names_from = ABX, values_from = RX) %>% print(n=50)
   
   # lets check this out...makes sense (peripheral blood)
   # source config file
   # astr <- tbl(conn, in_schema('AMB_ETL', 'LAB_MICRO_SENS_ALL_VW')) %>%
   #    filter(PERSON_ID == '1000056972') %>%
   #    collect() %>%
   #    arrange(RESULT_DATE)
   # asto <- tbl(conn, in_schema('AMB_ETL', 'LAB_MICRO_RESULT_ALL_VW')) %>%
   #    filter(PERSON_ID == '1000056972') %>%
   #    collect() %>%
   #    arrange(RESULT_DATE)
   # rm(asto, astr)
   
   rm(i, x, astDFi)
}


# Prep for joining
astDFj <- astDF %>% 
   select(PERSON_ID, ORDER_DAY, MULT_ISO, BUG) %>%
   group_by(PERSON_ID, ORDER_DAY, MULT_ISO) %>%
   summarise(ORDER_DAY = min(ORDER_DAY),
             BUG = list(sort(BUG)),
             MULT_BLOOD_ISO = n() > 1L) %>%
   ungroup() %>%
   mutate(JOIN_START = ORDER_DAY - 1,
          JOIN_END = ORDER_DAY + 90)
abxDFj <- abxDF %>% 
   mutate(across(contains('DATE'), ~ as.Date(substr(.,1,10)))) %>% 
   select(-END_DATE) %>% 
   distinct()

# JOIN
empDF <- left_join(x = astDFj,
                   y = abxDFj,
                   by = join_by(PERSON_ID,
                                JOIN_START <= START_DATE,
                                JOIN_END >= START_DATE)) %>%
   select(-JOIN_START, -JOIN_END) %>%
   relocate(ABX, START_DATE, .before=BUG)


# remove cultures without ABX
empDF <- empDF %>% 
   # count(is.na(ABX)) # 4,914
   filter(!is.na(ABX))


# how many total and over what range of days of therapy are there?
# first ignore the antibiotc, just take day of therapy
# only take 1 of cultures if there are multiple
eDF <- empDF %>% # ~620K --> ~377K
   select(PERSON_ID, ORDER_DAY, BUG, START_DATE) %>%
   distinct()
eDF1 <- eDF %>% filter(lengths(BUG) == 1L) # ~377K --> ~318K
eDF1$BUG <- unlist(eDF1$BUG)



d <- eDF %>% #filter(lengths(BUG) == 1L) %>% 
   group_by(PERSON_ID, ORDER_DAY) %>% tally() %>% ungroup() %>% select(n) %>% unlist(); names(d) <- NULL
t <- table(d)
t <- t[as.integer(names(t)) <= 70L]
barplot(t, names.arg=names(t), xlab='Days', main='Total time spent on antibiotic therapy')
median(d) # 8 (8, 9)
mean(d)   # 11.2 (10.8, 13.2)
quantile(d, c(0.1, 0.25, 0.5, 0.75, 0.9))
# length of therapy doesn't change per year
eDF %>%
   group_by(PERSON_ID, ORDER_DAY) %>% tally() %>% ungroup() %>%
   mutate(year = substr(ORDER_DAY,1,4)) %>%
   summarise(nn = n(), mean=mean(n), median = median(n), .by = year) %>%
   arrange(year)
rm(d, t)


eDF1 %>% count(BUG, sort=TRUE)
t <- eDF1 %>%
   group_by(PERSON_ID, ORDER_DAY, BUG) %>% tally() %>% ungroup() %>%
   summarise(nn = n(), median = median(n), 
             .by = BUG) %>%
   filter(nn > 100L)

t %>% arrange(desc(nn))
t %>% arrange(desc(median))
t %>% arrange(median)

t %>% arrange(desc(mean))
t %>% arrange(mean)











