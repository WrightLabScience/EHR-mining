library(dplyr)
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/data_path_name.Rdata')
load(file = paste0(data_path_name, 'ASTs_AbxAdmin_blood_2017_2023_imputed_flagged_WIDE_survival.Rdata'))
empDF <- empDF %>% filter(substr(ORDER_DAY, 1, 4) == '2017') # 5,755 rows

# read in lab results
source('~/Desktop/EHR/EHR work/config_file.R')
labs <- tbl(conn, in_schema('AMB_ETL', 'BACT_LAB_RESULT_VW')) %>% collect() # 10.8 million rows / 4,540 patients
labs <- labs %>% 
   filter(substr(ORDER_DATE, 7, 10) == '2017', substr(RESULT_DATE, 7, 10) == '2017') %>%
   mutate(across(contains('DATE'), ~ strptime(., format='%m/%d/%Y %T')))

# intersection of patients in both tables + basic cleaning (fixing data formats, removing unnecessary, etc.)
patient_id_intersection <- intersect(labs$PERSON_ID, empDF$PERSON_ID) # 3,326
empDF <- empDF %>% filter(PERSON_ID %in% patient_id_intersection) # 4,070 rows
labs <- labs %>% filter(PERSON_ID %in% patient_id_intersection) # 8.9 million rows / 3,326 patients
rm(patient_id_intersection)


# Find and clean up most common lab tests
{
   cml <- list(GLUCOSE = c('GLUCOSE(BEDSIDE TEST)', 'GLUCOSE IN SER/PLAS (MG/DL)', 'GLUCOSE, POC', 'GLUCOSE WHOLE BLOOD', 
                           'WHOLE BLOOD GLUCOSE (POCT)', 'GLUCOSE ISTAT', 'GLUCOSE, WHOLE BLOOD ISTAT'),
               HEMOGLOBIN = c('HGB', 'HEMOGLOBIN-ARTERIAL', 'HEMOGLOBIN ISTAT', 'TOTAL HGB(VENOUS)', 'HEMOGLOBIN (POCT)'),
               HEMATOCRIT = c('HEMATOCRIT(HCT)', 'HEMATOCRIT ISTAT', 'HEMATOCRIT (POCT)'),
               POTASSIUM = c('POTASSIUM(K)', 'POTASSIUM ISTAT', 'POTASSIUM(K) WHOLE BLOOD', 'POTASSIUM(K) WHOLE BLOOD &&', 'POTASSIUM (K) (POCT)', 'POTASSIUM K WHOLE BLOOD'),
               CARBON_DIOXIDE = c('CARBON DIOXIDE(CO2)'),
               CHLORIDE = c('CHLORIDE(CL)', 'CHLORIDE, WHOLE BLOOD', 'CHLORIDE ISTAT'),
               SODIUM = c('SODIUM(NA)', 'SODIUM(NA) WHOLE BLOOD', 'SODIUM ISTAT', 'SODIUM (NA) (POCT)', 'SODIUM NA WHOLE BLOOD', 'SODIUM VENOUS POC'),
               UREA_NITROGEN = c('UREA NITROGEN', 'BLOOD UREA NITROGEN', 'UREA NITROGEN ISTAT', 'UREA NITROGEN, WHOLE BLOOD IST'),
               CREATININE = c('CREATININE', 'CREATININE, WHOLE BLOOD', 'CREATININE ISTAT', 'CREATININE VENOUS ISTAT'),
               RBC = c('RBC'),
               MCV = c('MCV'),
               MCH = c('MCH'),
               MCHC = c('MCHC'),
               WBC = c('WBC'),
               ANION_GAP = c('ANION GAP', 'ANION GAP ISTAT'),
               PLATELETS = c('PLATELETS'),
               MEAN_PLATELET_VOLUME = c('MEAN PLATELET VOLUME'),
               CALCIUM = c('CALCIUM(CA)'),
               EGFR = c('EGFR'),
               RDW = c('RDW', 'RDW-CV'),
               ABS_NEUTROPHILS = c('ABS NEUTROPHILS'),
               PCT_NEUTROPHILS = c('NEUTROPHILS', 'NEUTROPHILS %'),
               ABS_LYMPHOCYTES = c('ABS LYMPHOCYTES'),
               PCT_LYMPHOCYTES = c('LYMPHOCYTES'),
               ABS_MONOCYTES = c('ABS MONOCYTES'),
               PCT_MONOCYTES = c('MONOCYTES'),
               ABS_BASOPHILS = c('ABS BASOPHILS'),
               PCT_BASOPHILS = c('BASOPHILS'),
               ABS_EOSINOPHILS = c('ABS EOSINOPHILS'),
               PCT_EOSINOPHILS = c('EOSINOPHILS'),
               PHOSPHORUS = c('PHOSPHORUS'),
               ALBUMIN = 'ALBUMIN',
               TOTAL_BILIRUBIN = 'TOTAL BILIRUBIN',
               AST = c('ASPARTATE AMINOT.(AST)', 'ASPARTATE AMINOTRANSF. (AST)'),
               ALT = c('ALANINE AMINOTRANS(ALT)'),
               ALP = c('ALKALINE PHOSPHATASE'),
               TOTAL_PROTEIN = 'TOTAL PROTEIN',
               INR = c('INR', 'INR (POC)', 'INR ISTAT'),
               PROTHROMBIN_TIME = c('PROTHROMBIN TIME'),
               LACTATE = c('LACTATE BLOOD', 'LACTATE WHOLE BLOOD', 'LACTATE ARTERIAL (POC)', 'LACTATE ISTAT', 'LACTATE VENOUS (POC)'))
}

lab_names <- character(nrow(labs))
for (i in seq_along(cml)) {
   print(i)
   w <- which(labs$COMPONENT_NAME %in% cml[[i]])
   lab_names[w] <- names(cml)[i]
}
labs$LAB <- lab_names
rm(w, i, lab_names)

# remove less common lab results
labs <- labs %>% filter(LAB != '') # 8.9 million --> 6 million
labs <- labs %>% filter(!is.na(RESULT_VALUE)) # only 1,750 were missing

# handle non-numeric types
w <- grep('<|>', labs$RESULT_VALUE) # handle cutoffs - its mostly EGFR and creatinine, ~ 100K total
labs$RESULT_VALUE[w] <- gsub('<|>', '', labs$RESULT_VALUE[w])
w <- grep('%', labs$RESULT_VALUE) # handle percentages, ~ 3.5K
labs$RESULT_VALUE[w] <- gsub('%', '', labs$RESULT_VALUE[w])
w <- which(is.na(as.numeric(labs$RESULT_VALUE))) # REMAINING ~23K
# labs[w,] %>% count(RESULT_VALUE, sort=TRUE) %>% View() # all words (pending, result cancelled, unavailable, etc.)
labs <- labs[-w,]
rm(w)

# convert to numeric
labs$RESULT_VALUE <- as.numeric(labs$RESULT_VALUE)

# the same lab might be expressed in different units of measurement...LOL
units_3u <- c('x10e+09/l', '10e+9/l', '10*9/l',
              'x10e*3', 'x10e3/ul', 'x10e+03/ul', '10s3/ul', '10s3ul', 'x10e3/cumm', '10^3', '10[3]',
              'th/cumm', 'th/ul', 'th/mm3', 'thous/mcl', '/cmm',
              '1000/mcl',
              'k/ul')
units_6u <- c('x10e+12/l', '10e+12/l', '10*12/l',
              '/cmm',
              '10^3/ul',
              'x10e*6','x10e+06/ul', 'x10e6/ul', '10s6/ul', '10^6', '10[6]',
              'mil/ul', 'mill/mcl', 'm/ul',  'mil/mm3', 'mill/ul')
units_mmol <- c('mmol/l', 'meq/l', 'mm/l')
units_cells <- c('cells/mcl', 'cells/ul')
units_iuL <- c('iu/l', 'u/l')
units_rate <- c('ml/min/1.73**2', 'ml/min/1.73m**2', 'ml/min', 'ml/min/1.73m2', 'ml/min/1.73')
units_g_dl <- c('g/dl', 'gm/dl', 'g%')
units_secs <- c('seconds', 'sec.', 'sec', 'secs.')
units_fl <- c('fl', 'um3')

# prep units - make lowercase
labs <- labs %>% mutate(CLEAN_REFERENCE_UNIT = tolower(REFERENCE_UNIT))
labs <- labs %>%
   mutate(CLEAN_REFERENCE_UNIT = case_when(
      LAB == 'INR' & !is.na(CLEAN_REFERENCE_UNIT) ~ 'no units',
      LAB == 'RBC' & CLEAN_REFERENCE_UNIT %in% units_6u ~ '10^6/uL',
      !LAB == 'RBC' & CLEAN_REFERENCE_UNIT %in% units_3u ~ '10^3/uL',
      CLEAN_REFERENCE_UNIT %in% units_mmol ~ 'mM',
      CLEAN_REFERENCE_UNIT %in% units_cells ~ 'cells/uL',
      CLEAN_REFERENCE_UNIT %in% units_iuL ~ 'u/L',
      CLEAN_REFERENCE_UNIT %in% units_rate ~ 'ml/min/1.73m^2',
      CLEAN_REFERENCE_UNIT %in% units_g_dl ~ 'g/dL',
      CLEAN_REFERENCE_UNIT %in% units_secs ~ 'seconds',
      CLEAN_REFERENCE_UNIT %in% units_fl ~ 'fL',
      .default = CLEAN_REFERENCE_UNIT
   ))
labs <- labs %>% # THIS IS VERY ROWS ~ 31 in total
   mutate(CLEAN_REFERENCE_UNIT = case_when(
      LAB %in% c('EGFR', 'UREA_NITROGEN') & CLEAN_REFERENCE_UNIT == 'u' ~ 'marked for removal',
      LAB == 'LACTATE' & CLEAN_REFERENCE_UNIT == 'mg/dl' ~ 'marked for removal',
      LAB == 'PCT_LYMPHOCYTES' & CLEAN_REFERENCE_UNIT == 'absolute' ~ 'marked for removal',
      LAB == 'PCT_NEUTROPHILS' & CLEAN_REFERENCE_UNIT %in% c('absolute', 'fL') ~ 'marked for removal',
      LAB == 'POTASSIUM' & CLEAN_REFERENCE_UNIT == '.' ~ 'marked for removal',
      LAB == 'RDW' & CLEAN_REFERENCE_UNIT == 'fL' ~ 'marked for removal',
      LAB == 'TOTAL_BILIRUBIN' & CLEAN_REFERENCE_UNIT == 'mg/ml' ~ 'marked for removal',
      .default = CLEAN_REFERENCE_UNIT
   )) 
labs <- labs %>% filter(CLEAN_REFERENCE_UNIT != 'marked for removal' | is.na(CLEAN_REFERENCE_UNIT))
labs <- labs %>%
   mutate(RESULT_VALUE = case_when(
      CLEAN_REFERENCE_UNIT == 'cells/uL' ~ RESULT_VALUE / 1000,
      LAB == 'ABS_BASOPHILS' & CLEAN_REFERENCE_UNIT == '%' ~ RESULT_VALUE / 100,
      LAB == 'LACTATE' & CLEAN_REFERENCE_UNIT == 'mg/dl' ~ RESULT_VALUE * 0.111,
      .default = RESULT_VALUE
   ))
labs <- labs %>%
   mutate(CLEAN_REFERENCE_UNIT = case_when(
      grepl('^ABS_', LAB) & CLEAN_REFERENCE_UNIT %in% c('%', 'cells/uL', '10^3/uL') ~ '1000 cells/uL',
      .default = CLEAN_REFERENCE_UNIT
   ))
labs <- labs[-which(labs$LAB == 'LACTATE' & is.na(labs$CLEAN_REFERENCE_UNIT)),] # 1 instance of this

rm(units_3u, units_6u, units_cells, units_fl, units_g_dl, units_iuL, units_mmol, units_rate, units_secs)

# look at the value distributions when measurement unit is missing
info <- labs %>%
   summarise(n = n(),
             q25 = quantile(RESULT_VALUE, 0.25),
             median = median(RESULT_VALUE),
             q75 = quantile(RESULT_VALUE, 0.75),
             max = max(RESULT_VALUE),
             IQR = q75 - q25,
             xmax = q75 + IQR * 6,
             xmin = q25 - IQR * 6,
             .by = c(LAB, CLEAN_REFERENCE_UNIT)) %>%
   arrange(LAB, desc(n))
info$xmin <- sapply(info$xmin, function(x) max(c(0, x)))
units <- info %>%
   filter(!is.na(CLEAN_REFERENCE_UNIT)) %>%
   select(LAB, CLEAN_REFERENCE_UNIT)
units <- setNames(units$CLEAN_REFERENCE_UNIT, units$LAB)

# impute missing measurement unit
w <- which(is.na(labs$CLEAN_REFERENCE_UNIT)) # ~162K
labs$CLEAN_REFERENCE_UNIT[w] <- unname(units[labs$LAB[w]])

rm(info, w, units)

labs %>%
   count(LAB, CLEAN_REFERENCE_UNIT) %>%
   arrange(desc(n)) %>%
   View()

# calculate delay between order and result
labs <- labs %>% mutate(RESULT_DELAY = as.numeric(lubridate::as.duration(RESULT_DATE - ORDER_DATE)) / 86400)
x <- labs$RESULT_DELAY
sum(x >= 0) / length(x) # 98.6%
sum(x < 0 & x > -1) / length(x) # 0.9%
sum(x < 0) / length(x) # 1.4%
hist(x[x <= 0], breaks=diff(range(x[x <= 0])), ylim=c(0, 1e2))

labs %>%
   summarise(q01 = quantile(RESULT_DELAY, 0.01),
             q50 = quantile(RESULT_DELAY, 0.50),
             q99 = quantile(RESULT_DELAY, 0.99),
             .by = c(LAB)) %>%
   filter(q01 < -1)

labs <- labs %>%
   rename(UNITS = CLEAN_REFERENCE_UNIT) %>%
   select(PERSON_ID, ORDER_DATE, LAB, RESULT_VALUE, UNITS) %>%
   distinct() # 6,046,298 --> 3,438,025

labs <- labs %>%
   arrange(PERSON_ID, LAB, ORDER_DATE)

save(labs, file = paste0(data_path_name, 'LABS_cleaned_blood_2017subset.Rdata'))















