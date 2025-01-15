library(dplyr)
load(file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/MRSA_bacteremia/LabResults_raw.Rdata')

### CLEAN UP LAB RESULTS
# after done, make a general function that does this

labs <- labs %>%
   mutate(across(contains('DATE'), ~ strptime(., format='%m/%d/%Y %T'))) %>%
   arrange(PERSON_ID, ORDER_DATE, RESULT_DATE)



{
   cml <- list(
      GLUCOSE = c('GLUCOSE(BEDSIDE TEST)', 'GLUCOSE IN SER/PLAS (MG/DL)', 'GLUCOSE, POC', 'GLUCOSE WHOLE BLOOD', 
                  'WHOLE BLOOD GLUCOSE (POCT)', 'GLUCOSE ISTAT', 'GLUCOSE, WHOLE BLOOD ISTAT'),
      HEMOGLOBIN = c('HGB', 'HEMOGLOBIN-ARTERIAL', 'HEMOGLOBIN ISTAT', 'TOTAL HGB(VENOUS)', 'HEMOGLOBIN (POCT)'),
      HEMOGLOBIN_A1C = c('HEMOGLOBIN A1C'),
      HEMATOCRIT = c('HEMATOCRIT(HCT)', 'HEMATOCRIT ISTAT', 'HEMATOCRIT (POCT)', 'HEMATOCRIT DERIVED', 'HEMATOCRIT(HCT) MANUAL PCV &&',
                     'HEMATOCRIT-BODY FLUID (HCT)', 'HEMATOCRIT VENOUS POC', 'HEMATOCRIT DERIVED - MIXED VEN'),
      POTASSIUM = c('POTASSIUM(K)', 'POTASSIUM ISTAT', 'POTASSIUM(K) WHOLE BLOOD', 'POTASSIUM(K) WHOLE BLOOD &&', 'POTASSIUM (K) (POCT)', 'POTASSIUM K WHOLE BLOOD'),
      CARBON_DIOXIDE = c('CARBON DIOXIDE(CO2)'),
      CHLORIDE = c('CHLORIDE(CL)', 'CHLORIDE, WHOLE BLOOD', 'CHLORIDE ISTAT'),
      SODIUM = c('SODIUM(NA)', 'SODIUM(NA) WHOLE BLOOD', 'SODIUM ISTAT', 'SODIUM (NA) (POCT)', 'SODIUM NA WHOLE BLOOD', 'SODIUM VENOUS POC'),
      UREA_NITROGEN = c('UREA NITROGEN', 'UREA NITROGEN ISTAT', 'UREA NITROGEN, WHOLE BLOOD IST'), # 'BLOOD UREA NITROGEN', 'UREA NITROGEN URINE'
      CREATININE = c('CREATININE', 'CREATININE, WHOLE BLOOD', 'CREATININE ISTAT', 'CREATININE VENOUS ISTAT'),
      RBC = c('RBC'),
      MCV = c('MCV'),
      MCH = c('MCH'),
      MCHC = c('MCHC'),
      WBC = c('WBC'),
      ANION_GAP = c('ANION GAP', 'ANION GAP ISTAT'),
      PLATELETS = c('PLATELETS'),
      MEAN_PLATELET_VOLUME = c('MEAN PLATELET VOLUME'),
      CALCIUM = c('CALCIUM(CA)', 'CALCIUM I-STAT ARTERIAL'),
      EGFR = c('EGFR', 'EGFR (CKD-EPI 2021), SER/PLAS/WHOLE BLD', 'EGFR AFRICAN AMERICAN'),
      RDW = c('RDW', 'RDW-CV'),
      ABS_NEUTROPHILS = c('ABS NEUTROPHILS'),
      PCT_NEUTROPHILS = c('NEUTROPHILS', 'NEUTROPHILS %'),
      ABS_LYMPHOCYTES = c('ABS LYMPHOCYTES'),
      PCT_LYMPHOCYTES = c('LYMPHOCYTES'),
      ATYP_PCT_LYMPHOCYTES = c('ATYPICAL LYMPHOCYTES', 'ATYPICAL LYMPHOCYTES @'),
      ATYP_ABS_LYMPHOCYTES = c('ABS ATYPICAL LYMPHOCYTES'),
      ABS_MONOCYTES = c('ABS MONOCYTES'),
      PCT_MONOCYTES = c('MONOCYTES'),
      ABS_BASOPHILS = c('ABS BASOPHILS', 'BASOPHILS, ABS'),
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
      LACTATE = c('LACTATE BLOOD', 'LACTATE WHOLE BLOOD', 'LACTATE ARTERIAL (POC)', 'LACTATE ISTAT', 'LACTATE VENOUS (POC)')
   )
}

cns <- labs %>% 
   count(COMPONENT_NAME, sort=TRUE) %>%
   mutate(LAB = NA)
for (i in seq_along(cml)) {
   w <- which(cns$COMPONENT_NAME %in% cml[[i]])
   cns$LAB[w] <- names(cml)[i]
}
rm(i, w, cml)

# which component names have not been captured yet?
cns %>% filter(is.na(LAB))

# remove uncaptured component names?
cns <- cns %>% filter(!is.na(LAB))
cns <- setNames(cns$LAB, cns$COMPONENT_NAME)

# modify original dataframe to include common lab names
labs <- labs %>%
   mutate(LAB = unname(cns[COMPONENT_NAME])) %>%
   relocate(LAB, .before=COMPONENT_NAME) %>%
   filter(!is.na(LAB)) %>%
   filter(!is.na(RESULT_VALUE))
rm(cns)

labs %>% count(LAB, sort=TRUE)# %>% View()

labs <- labs %>%
   mutate(REFERENCE_UNIT = tolower(REFERENCE_UNIT))


# handle cutoffs
labs$CUTOFF <- 0
w <- grep('<', labs$RESULT_VALUE)
labs$RESULT_VALUE[w] <- gsub('<', '', labs$RESULT_VALUE[w])
labs$CUTOFF[w] <- -1

w <- grep('>', labs$RESULT_VALUE)
labs$RESULT_VALUE[w] <- gsub('>', '', labs$RESULT_VALUE[w])
labs$CUTOFF[w] <- 1

rm(w)


# what are each labs list of units?
ls <- unique(labs$LAB)
for (l in ls) {
   u <- labs %>% filter(LAB == l) %>% count(REFERENCE_UNIT) %>% pull(REFERENCE_UNIT)
   if (any(grepl('%', u))) cat(l, u, '\n')
}
rm(ls, l, u)

# ABS BASOPHILS SHOULDNT HAVE PERCENTAGE, remove those
w <- which(labs$LAB == 'ABS_BASOPHILS' & labs$REFERENCE_UNIT == '%')
labs <- labs[-w,]

# PCT_NEUTROPHILS shouldnt have absolute counts, remove them
w <- which(labs$LAB == 'PCT_NEUTROPHILS' & labs$REFERENCE_UNIT == 'absolute')
labs <- labs[-w,]

# PCT_LYMPHOCYTES shouldnt have absolute counts, remove them
w <- which(labs$LAB == 'PCT_LYMPHOCYTES' & labs$REFERENCE_UNIT == 'absolute')
labs <- labs[-w,]


# CLEAN UP REFERENCE UNITS
units <- list(
   'pcnt' = c('%', '% of total hgb'),
   'thou' = c('x10e+09/l', '10e+9/l', '10*9/l',
                  'x10e*3', 'x10e3/ul', 'x10e+03/ul', '10s3/ul', '10s3ul', 'x10e3/cumm', '10^3', '10[3]',
                  'thousand/ul', 'th/cumm', 'th/ul', 'th/mm3', 'thous/mcl', '/cmm', 'thou',
                  '1000/mcl',
                  'k/ul'),
   'mill' = c('x10e+12/l', '10e+12/l', '10*12/l',
                  '/cmm',
                  '10^3/ul',
                  'x10e*6','x10e+06/ul', 'x10e6/ul', '10s6/ul', '10^6', '10[6]',
                  'mil/ul', 'mill/mcl', 'm/ul',  'mil/mm3', 'mill/ul', 'million/ul', 'mil'),
   'mmol' = c('mmol/l', 'meq/l', 'mm/l', 'mmoll'),
   'cells' = c('cells/mcl', 'cells/ul'),
   'iuL' = c('iu/l', 'u/l'),
   'rate' = c('ml/min/1.73**2', 'ml/min/1.73m**2', 'ml/min', 'ml/min/1.73m2', 'ml/min/1.73'),
   'g_dl' = c('g/dl', 'gm/dl', 'g%'),
   'mg_dl' = c('mg/dl'),
   'secs' = c('seconds', 'sec.', 'sec', 'secs.'),
   'fl' = c('fl', 'um3', 'ume3', 'um'),
   'pg' = c('pg', 'uug'),
   'none' = c('therapeutic')
)
setdiff(unique(labs$REFERENCE_UNIT), unlist(units))

x <- labs %>% filter(REFERENCE_UNIT == '%')
x %>% count(LAB)
labs %>% 
   filter(LAB == 'ABS_EOSINOPHILS') %>% #count(CLEAN_REF_UNIT)
   mutate(val = as.numeric(RESULT_VALUE)) %>% 
   summarise(n=n(), 
             mean=mean(val, na.rm=T), 
             median=median(val, na.rm=T), 
             .by=REFERENCE_UNIT)
rm(x)


# remove x10
# INR: no ref unit (therapeutic)
w <- which(labs$REFERENCE_UNIT == 'x10') # 19
labs <- labs[-w,]
rm(w)


# CLEAN UP UNITS
labs$CLEAN_REF_UNIT <- NA_character_
for (i in seq_along(units)) {
   labs$CLEAN_REF_UNIT[labs$REFERENCE_UNIT %in% units[[i]]] <- names(units)[i]
}
rm(i)
labs <- labs %>%
   mutate(CLEAN_REF_UNIT = case_when(
      grepl('^PCT_', LAB) & is.na(REFERENCE_UNIT) ~ 'pcnt',
      LAB == 'INR' ~ 'none',
      LAB == 'EGFR' & is.na(REFERENCE_UNIT) ~ 'rate',
      .default = CLEAN_REF_UNIT
   ))
labs <- labs %>% filter(!is.na(CLEAN_REF_UNIT))

labs <- labs %>% relocate(CLEAN_REF_UNIT, .before=REFERENCE_UNIT)

# handle percentages
labs %>% count(grepl('%', RESULT_VALUE), CLEAN_REF_UNIT == 'pcnt')
w <- grep('%', labs$RESULT_VALUE) # 1243
labs$RESULT_VALUE[w] <- gsub('%', '', labs$RESULT_VALUE[w])

# handle some words??
labs <- labs %>% mutate(RESULT_VALUE = tolower(RESULT_VALUE))
w <- which(labs$RESULT_VALUE == 'less than 1')
labs$RESULT_VALUE[w] <- 1
rm(w)
# remove those rows
labs <- labs %>% filter(grepl('^[0-9.]+$', RESULT_VALUE))

# Make result values into actual numbers!!!
labs <- labs %>% mutate(RESULT_VALUE = as.numeric(RESULT_VALUE))

# Scale some of the result values according to their units
labs <- labs %>%
   mutate(RESULT_VALUE = case_when(
      CLEAN_REF_UNIT == 'cells' ~ RESULT_VALUE / 1000,
      .default = RESULT_VALUE
   ))
rm(units)




# EGFR stuff
if (FALSE) {
   e <- labs %>%
      filter(grepl('EGFR', COMPONENT_NAME)) %>% 
      select(PERSON_ID, COMPONENT_NAME, ORDER_DATE, RESULT_VALUE, CUTOFF) %>% 
      distinct() %>%
      mutate(ORDER_DAY = lubridate::as_date(ORDER_DATE))
   
   e %>%
      summarise(n=n(), 
                min=min(RESULT_VALUE, na.rm=T),
                mean=mean(RESULT_VALUE, na.rm=T), 
                median=median(RESULT_VALUE, na.rm=T),
                sd=sd(RESULT_VALUE, na.rm=T),
                max=max(RESULT_VALUE, na.rm=T),
                frac_cutoff=sum(CUTOFF == 1) / n,
                frac_high=sum(RESULT_VALUE > 60) / n,
                .by=COMPONENT_NAME)
   
   # most patients have either EGFR_CK OR EGFR AND EGFR_AA!!
   # EGFR_CK does not ever cap at 59 or 60 
   ew <- e %>%
      tidyr::pivot_wider(
         id_cols = c(PERSON_ID, ORDER_DATE),
         names_from = COMPONENT_NAME,
         values_from = RESULT_VALUE,
         values_fn = mean
      ) %>%
      rename(EGFR_AA = `EGFR AFRICAN AMERICAN`,
             EGFR_CK = `EGFR (CKD-EPI 2021), SER/PLAS/WHOLE BLD`) %>%
      mutate(ratioA = EGFR_AA / EGFR)
   ew %>% pull(ratioA) %>% hist(., breaks=1000)
   
   ew %>% count(hasNorm = !is.na(EGFR), hasAA = !is.na(EGFR_AA))
   ew %>% count(hasNorm = !is.na(EGFR), hasCK = !is.na(EGFR_CK))
   
   e %>% filter(COMPONENT_NAME == 'EGFR') %>% pull(RESULT_VALUE) %>% hist(ylim=c(0,10000), xlim=c(0,175))
   e %>% filter(COMPONENT_NAME == 'EGFR (CKD-EPI 2021), SER/PLAS/WHOLE BLD') %>% pull(RESULT_VALUE) %>% hist(xlim=c(0,175))
   
   rm(e, ew)
}
labs <- labs %>% filter(COMPONENT_NAME != 'EGFR AFRICAN AMERICAN')

# MORPHOLOGY - "Includes:" is most common, then "Slight", "Observations:"
labs %>% filter(COMPONENT_NAME == 'RBC MORPHOLOGY') %>% count(RESULT_VALUE, sort=TRUE)
# GLUCOSE (PRESENCE) - mostly "Neg" or "Negative", some are numbers
labs %>% filter(COMPONENT_NAME == 'GLUCOSE (PRESENCE) IN URINE BY TEST STRIP') %>% count(RESULT_VALUE, sort=TRUE)
# HEMATOCRIT CAPILLARY (POC) - mmol/L - everything else was %, only use this if HEMATOCRIT % missing??
# HEMOGLOBIN (others)
# ABS ATYPICAL LYMPHOCYTES - in X10E*3, X10E+09/L, x10e3/uL
# CALCIUM - split
# 



################################################
labs <- labs %>% mutate(ORDER_DAY = lubridate::as_date(ORDER_DATE))
#labsOG <- labsOG %>% mutate(as.numeric(lubridate::as.duration(RESULT_DATE - ORDER_DATE)) / 86400)
save(labs, file='~/Desktop/EHR/EHR work/RdataFiles/causal_prep/MRSA_bacteremia/LabResults_cleaned.Rdata')


library(dplyr)
load(file='~/Desktop/EHR/EHR work/RdataFiles/causal_prep/MRSA_bacteremia/LabResults_cleaned.Rdata')


# which labs have multiple measurements per day?
labsS <- labs %>% count(PERSON_ID, LAB, RESULT_VALUE, ORDER_DAY)
lab2 <- labsS %>% filter(n > 1L)
lab1 <- labsS %>% filter(n == 1L)

df <- full_join(
   x = lab2 %>% count(LAB, sort=TRUE),
   y = lab1 %>% count(LAB, sort=TRUE),
   by = join_by(LAB)
) %>%
   rename(n1 = n.y, n2 = n.x) %>%
   mutate(d = n1 / n2) %>%
   arrange(desc(d))

plot(x = df$d, y = seq_len(nrow(df)), log='x')
abline(v = 1, lty = 2)

df %>% arrange(d)
df %>% arrange(desc(d))
rm(df)


# how different are same day measurements?
lab2 %>%
   filter(LAB == 'GLUCOSE') %>%
   reframe(n = n(),
           range = abs(diff(range(RESULT_VALUE))),
           .by = c(PERSON_ID, LAB, ORDER_DAY)) %>% 
   filter(range > 0) %>%
   pull(range) %>% hist(breaks=500, xlim=c(0, 500))

x <- labs %>%
   filter(LAB == 'GLUCOSE', PERSON_ID == '1000002142', lubridate::year(ORDER_DAY) == 2022, lubridate::month(ORDER_DAY) == 5)

t <- lubridate::decimal_date(x$ORDER_DATE)
t <- t - min(t)
t <- t * 365
g <- x$RESULT_VALUE
plot(x=t, y=g, ylim=c(50, 400), pch=16)
#points(x=t, y=g^1.02)
lines(x=t, y=g, lwd=1.2)
lines(x=t, y=stats::filter(g, rep(1/4, 4)), col='blue', lwd=2)
lines(x=t, y=stats::filter(g^1.02, rep(1/4, 4)), col='blue', lwd=2)
lines(x=t, y=stats::filter(g^1.05, rep(1/4, 4)), col='blue', lwd=2)
abline(h = 200, lty=3)
rm(t, g, x)
rm(lab1, lab2)



load(file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/MRSA_bacteremia/MRSA_BSI_DAP_VAN.Rdata')
dfx <- dfx %>%
   select(PERSON_ID, ORDER_DAY, ORDER_DATE, TRT, time, FIRST_ABX_DAY_REL_ORDER, FIRST_ABX_TIME_REL_ORDER) %>%
   mutate(ABX_START_TIME = ORDER_DATE + 86400 * FIRST_ABX_TIME_REL_ORDER) %>%
   mutate(JOIN_START = ORDER_DAY + FIRST_ABX_DAY_REL_ORDER - 7,
          JOIN_END = ORDER_DAY + FIRST_ABX_DAY_REL_ORDER + 5)

# prep labs
labsS <- labs %>%
   count(PERSON_ID, ORDER_DAY, ORDER_DATE, RESULT_DATE, RESULT_VALUE, LAB) %>%
   rename(LAB_ORDER_DAY = ORDER_DAY)


# for each lab, what is the rate of missingness?
uniq_labs <- unique(labs$LAB) # 41
uniq_labs <- grep('ATYP', uniq_labs, value=TRUE, invert=TRUE)
for (l in seq_along(uniq_labs)) {
   cat(l, uniq_labs[l], '\n')
   df <- dfx %>%
      left_join(
         x = .,
         y = labsS %>% filter(LAB == uniq_labs[l]) %>% select(!c(ORDER_DATE, LAB_ORDER_DAY)),
         by = join_by(
            PERSON_ID,
            between(y$RESULT_DATE, x$JOIN_START, x$JOIN_END)
         )
      ) %>%
      mutate(X = as.numeric(lubridate::as.duration(RESULT_DATE - ABX_START_TIME)) / 86400) %>%
      select(!c(JOIN_START, FIRST_ABX_DAY, RESULT_VALUE, RESULT_DATE)) %>%
      mutate(n = case_when(
         is.na(n) ~ 1L,
         .default = as.integer(n > 0L),
      )) %>%
      mutate(Xr = round(X)) %>%
      arrange(Xr) %>%
      mutate(Xr = case_when(
         Xr < 0L ~ gsub('-', 'n', Xr),
         Xr > 0L ~ paste0('p', Xr),
         Xr == 0L ~ 'z',
         .default = 'missing'
      )) %>% 
      summarise(hasPre = any(X <= 0, na.rm=T),
                hasH4 = any(X <= 4/24, na.rm=T), 
                hasH8 = any(X <= 8/24, na.rm=T),
                hasH16 = any(X <= 16/24, na.rm=T),
                hasD1 = any(X <= 1, na.rm=T), 
                hasD3 = any(X <= 3, na.rm=T),
                .by=c(PERSON_ID, ORDER_DAY, TRT))
   
   names(df)[names(df) == 'hasPre'] <- paste0('pre_', uniq_labs[l])
   names(df)[names(df) == 'hasH4'] <- paste0('h4_', uniq_labs[l])
   names(df)[names(df) == 'hasH8'] <- paste0('h8_', uniq_labs[l])
   names(df)[names(df) == 'hasH16'] <- paste0('h16_', uniq_labs[l])
   names(df)[names(df) == 'hasD1'] <- paste0('d1_', uniq_labs[l])
   names(df)[names(df) == 'hasD3'] <- paste0('d3_', uniq_labs[l])
   dfx <- full_join(x = dfx, y = df, by=join_by(PERSON_ID, ORDER_DAY, TRT))
}
rm(df, l)

# creat missDF from dfx
flags <- c('pre', 'h4', 'h8', 'h16', 'd1', 'd3')
missDF <- array(NA, dim=c(nrow(dfx), length(uniq_labs), length(flags)), dimnames=list(dfx$TRT, uniq_labs, flags))
for (flag in flags) missDF[, , flag] <- as.matrix(dfx %>% select(contains(flag)))   
rm(flag)



lab_miss_rates <- apply(missDF, 3, function(mat) apply(mat, 2, function(x) 1 - (sum(x) / length(x))))

trt_diff_miss_rates <- array(NA, dim=c(2, length(uniq_labs), length(flags)), dimnames=list(c('iDAP', 'sDAP'), uniq_labs, flags))
for (flag in seq_along(flags)) {
   for (lab in seq_along(uniq_labs)) {
      t <- missDF[, lab, flag]
      t <- table(t, names(t))
      t <- t['FALSE',] / colSums(t)
      t <- c(t['iDAP'], t['sDAP']) - t['VAN']
      trt_diff_miss_rates[names(t), lab, flag] <- t
   }
}
rm(flag, lab, t)

o <- order(lab_miss_rates[,'pre'])
lab_miss_rates <- lab_miss_rates[o, ]
trt_diff_miss_rates <- trt_diff_miss_rates[, o, ]
rm(o)

col_vec <- colorRampPalette(c('darkblue', 'blue', 'red', 'pink3'))(n=length(flags))

{
   layout(matrix(1:3, nrow=1), widths=c(2,1,1))
   par(mar = c(3, 9, 3, 1), mgp=c(1.5, 0.5, 0), tck=-0.015, bg='white', cex=1, oma=c(0,0,0,1))
   plot(NA, xlim = c(0, 1), ylim=c(0.5, nrow(lab_miss_rates)+0.5), yaxt = 'n', ylab = '', xlab='Missing rate', xaxs='i', yaxs='i', main='Presence of lab results near abx treatment start')
   axis(side=2, at=1:nrow(lab_miss_rates), labels=stringr::str_to_sentence(gsub('_', ' ', rownames(lab_miss_rates))), las=1)
   abline(v = seq(0.2, 0.8, 0.2), lty=3, lwd=1.2)
   abline(v = seq(0.1, 0.9, 0.2), lty=3, lwd=0.6)
   for (flag in seq_along(flags)) {
      points(y = seq_len(nrow(lab_miss_rates)), x = lab_miss_rates[, flags[flag]], pch=16, col=col_vec[flag])      
   }
   legend('bottomright', legend=flags, pch=16, col=col_vec, title='Timing')
   
   par(mar=c(3, 2, 3, 1))
   for (t in 1:2) {
      plot(NA, xlim=c(-0.48, 0.48), ylim=c(0.5, ncol(trt_diff_miss_rates)+0.5), yaxt='n', ylab='', xlab='missing rate difference', yaxs='i', main=rownames(trt_diff_miss_rates)[t])
      abline(v = 0, lty=3, lwd=1.2)
      abline(v = seq(-0.4, 0.4, 0.1), lty=3)
      for (flag in seq_along(flags)) {
         points(y = seq_len(ncol(trt_diff_miss_rates)), x = trt_diff_miss_rates[t, , flags[flag]], pch=15, col=col_vec[flag])
      }
      arrows(x0=-0.15, x1=0.15, y0=ncol(trt_diff_miss_rates)+1, code=3, xpd=NA, lwd=1.5, length=0.1)
      text(x=0.4 * c(-1, 1), y=ncol(trt_diff_miss_rates)+1, labels=c('> data in DAP', '> data in VAN'), xpd=NA)
   }
}
rm(col_vec, t, flag, trt_diff_miss_rates, lab_miss_rates)

# lessons learned:
# VAN and DAP have the same rate of missingness at all times relative to their blood culture order time
# but since iDAP start on abx later than VAN group, so the pre-treatment missingness is higher
# the question is, 

# tons of lab results roll in between 8 and 16 hours after first antibiotic is administered
# iDAP group has lower rate of missingness than VAN for many common labs
# the delta between iDAP and VAN shrinks over time; therefore, iDAP group has more data EARLY than VAN group
# they end up with nearly identical amounts, but DAP has and gets the data earlier

# this means that if I keep only pre-treatment data, many more VAN individuals would have missing values

# let me visualize this with cumulative fraction of patients with data
# for a given lab, have 3 lines: 1 for each treatment group

### BELOW PLOTS CONFIRM THE SUSPICION

load(file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/MRSA_bacteremia/MRSA_BSI_DAP_VAN.Rdata')
dfx <- dfx %>%
   select(PERSON_ID, ORDER_DAY, ORDER_DATE, TRT, time, FIRST_ABX_DAY, FIRST_VAN_TIME, FIRST_DAP_TIME) %>%
   mutate(TRT = case_when(
      TRT == 'VANCOMYCIN' ~ 'VAN', 
      TRT == 'DAPTOMYCIN' ~ 'iDAP', 
      TRT == 'VANCOMYCIN,DAPTOMYCIN' ~ 'sDAP'
   )) %>%
   mutate(ABX_START_TIME = case_when(
      TRT == 'iDAP' ~ ORDER_DATE + 86400 * FIRST_DAP_TIME,
      .default = ORDER_DATE + 86400 * FIRST_VAN_TIME
   )) %>%
   mutate(JOIN_START = ORDER_DAY + FIRST_ABX_DAY - 7,
          JOIN_END = ORDER_DAY + FIRST_ABX_DAY + 5)

trt_cumsum <- setNames(vector('list', length(uniq_labs)), uniq_labs)
for (l in seq_along(uniq_labs)) {
   cat(l, uniq_labs[l], '\n')
   df <- dfx %>%
      left_join(
         x = .,
         y = labsS %>% filter(LAB == uniq_labs[l]) %>% select(!c(ORDER_DATE, LAB_ORDER_DAY)),
         by = join_by(
            PERSON_ID,
            between(y$RESULT_DATE, x$JOIN_START, x$JOIN_END)
         )
      ) %>%
      mutate(X = as.numeric(lubridate::as.duration(RESULT_DATE - ORDER_DATE)) / 86400) %>%
      group_by(PERSON_ID, ORDER_DAY) %>%
      slice_min(X, with_ties = FALSE) %>%
      ungroup() %>%
      arrange(X)
   t <- table(round(df$X * 24) / 24, df$TRT)
   t <- apply(t, 2, function(x) cumsum(x) / sum(x))
   trt_cumsum[[l]] <- t
}

diffs <- sapply(trt_cumsum, function(mat) max(mat[,'iDAP'] - mat[,'VAN']))
o <- order(diffs)
trt_cumsum <- trt_cumsum[o]


col_vec <- c('VAN' = "#0000FF", 'iDAP' = "#ef5675", 'sDAP' = "#ffa600")

{
   par(mfrow=c(6,7), mar = c(0, 0, 0, 0), mgp=c(1.5, 0.5, 0), tck=-0.02, oma=c(3, 3, 2, 1))
   for (l in 1:41) {
      plot(NA, xlim=range(as.numeric(rownames(t))), ylim=c(0,1), xlab='', ylab='', yaxt='n', xaxt='n')
      text(x=-8, y=1, adj=c(0,1), labels=names(trt_cumsum)[l])
      abline(v=c(-1,1), lty=3)
      abline(v=0, lty=3, lwd=2)
      
      # y - axis
      if ((l - 1) %% 7 == 0) {
         axis(side=2, las=1)
         title(ylab='Fraction with data', line=2, xpd=NA)
      }
      
      # x - axis
      if (l >= 36L) {
         axis(side=1)
         title(xlab='Days', xpd=NA)
      }
      
      # data 
      t <- trt_cumsum[[l]]
      for (i in 1:3) lines(x=as.numeric(rownames(t)), y=t[,i], col=col_vec[colnames(t)[i]], lwd=3)
   }
   legend('right', inset=c(-0.7, 0), legend=names(col_vec), col=col_vec, lwd=3, xpd=NA)
}
rm(t, i, l, diffs, o, trt_cumsum, col_vec, flags, uniq_labs, missDF, df, dfx)







library(dplyr)
source(file = '~/Desktop/WrightLab/UsefulRfuncs/newBarplotFxn.R')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/MRSA_bacteremia/LabResults_cleaned.Rdata')
labsS <- labs %>%
   count(PERSON_ID, ORDER_DAY, ORDER_DATE, RESULT_DATE, RESULT_VALUE, LAB) %>%
   rename(LAB_ORDER_DAY = ORDER_DAY)
load(file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/MRSA_bacteremia/MRSA_BSI_DAP_VAN.Rdata')
dfx <- dfx %>% select(!contains('POST'))
dfx <- dfx %>% select(!contains('PRE'))

# join all labs with rest of data
df <- dfx %>%
   select(-time) %>% rename(time = time_censored) %>%
   select(PERSON_ID, ORDER_DAY, ORDER_DATE, TRT, time, status, FIRST_ABX_DAY_REL_ORDER, FIRST_ABX_TIME_REL_ORDER, CATEGORY, AGE, FEMALE, year) %>%
   mutate(ABX_START_TIME = ORDER_DATE + 86400 * FIRST_ABX_TIME_REL_ORDER) %>%
   mutate(JOIN_START = ABX_START_TIME - 86400 * 3,
          JOIN_END = ABX_START_TIME + 86400 * 3) %>%
   left_join(
      x = .,
      y = labsS %>% select(!c(ORDER_DATE, LAB_ORDER_DAY)) %>% filter(!grepl('^ATYP_', LAB)),
      by = join_by(
         PERSON_ID,
         between(y$RESULT_DATE, x$JOIN_START, x$JOIN_END)
      )
   ) %>%
   mutate(X = as.numeric(lubridate::as.duration(RESULT_DATE - ABX_START_TIME)) / 86400) %>%
   select(-JOIN_START, -JOIN_END) %>%
   arrange(PERSON_ID, ORDER_DAY, LAB, RESULT_DATE)

df %>% pull(X) %>% hist(breaks=300)


df2 <- df %>%
   mutate(FLAG = case_when(
      !is.na(X) & X <= 0 ~ 'PRE',
      !is.na(X) & X > 0 ~ 'POST',
      .default = NA
   )) %>%
   mutate(LAB = case_when(
      FLAG == 'PRE' ~ paste0('PRE_', LAB),
      FLAG == 'POST' ~ paste0('POST_', LAB),
      .default = NA
   )) %>%
   select(PERSON_ID, ORDER_DAY, TRT, time, status, FIRST_ABX_DAY_REL_ORDER, LAB, RESULT_VALUE, X, FLAG) %>%
   distinct()

df2 %>% group_by(PERSON_ID, ORDER_DAY)

df2 <- df2 %>%
   mutate(Xn = ifelse(FLAG == 'PRE', X + 3, X)) %>%
   mutate(weight = case_when(
      FLAG == 'PRE' ~ 0.5 + (0.5 / 3) * Xn,
      FLAG == 'POST' ~ 1,
      .default = NA
   )) %>%
   summarise(
      MEAN = mean(RESULT_VALUE),
      WMEAN = weighted.mean(RESULT_VALUE, weight),
      .by = c(PERSON_ID, ORDER_DAY, TRT, LAB, time, status, FIRST_ABX_DAY_REL_ORDER, FLAG)
   )

df2 %>% filter(FLAG == 'POST') %>% summarise(abs(diff(range(WMEAN - MEAN))))
df2 %>% filter(FLAG == 'PRE') %>% summarise(abs(diff(range(WMEAN - MEAN))))

dfw <- df2 %>%
   select(-MEAN) %>%
   rename(VALUE = WMEAN) %>%
   tidyr::pivot_wider(
      id_cols = c(PERSON_ID, ORDER_DAY, TRT, time, status, FIRST_ABX_DAY_REL_ORDER),
      values_from = VALUE,
      names_from = LAB
   ) %>%
   select(-`NA`)



keep_labs <- character()
for (flag in c('PRE', 'POST')) {
   crl <- sapply(dfw %>% select(contains(flag)), function(x) sum(!is.na(x)) / length(x))
   crp <- apply(dfw %>% select(names(crl)), 1, function(x) sum(!is.na(x)) / length(x))
   crl <- sort(crl)
   
   par(mfrow=c(1,2), mar=c(4, 2, 3, 1), mgp=c(2, 0.5, 0), tck=-0.015, oma=c(0,11,0,0), cex=0.9)
   b <- barplot(crl, horiz=TRUE, names.arg=NA, xlim=c(0,1), main='By lab', xlab='Completion rate', xpd=NA)
   abline(v = c(0.2, 0.5, 0.8), lty=3)
   text(x=-0.02, y=b, adj=1, labels=gsub(paste0(flag, '_'), '', names(crl)), xpd=NA)
   barplot(sort(crp), horiz=TRUE, names.arg=NA, xlim=c(0,1), main='By patient', xlab='Completion rate')
   abline(v = c(0.2, 0.5, 0.8), lty=3)
   mtext(text=paste0(flag, ' antibiotic treatment'), side=3, outer=TRUE, line=-2, font=2)
   
   if (flag == 'PRE') keep_labs <- c(keep_labs, names(crl[crl > 0.4]))
   if (flag == 'POST') keep_labs <- c(keep_labs, names(crl[crl > 0.8]))
   
   # sicker patients have more labs present
   dfw %>%
      mutate(r = crp) %>%
      coxph(formula = Surv(time, status) ~ r, data=.) %>%
      summary()
}
lab_vals <- dfw %>% select(PERSON_ID, ORDER_DAY, !!grep('^PRE_', names(dfw), value=TRUE))
save(lab_vals, file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/MRSA_bacteremia/LabValues_DF_MRSA.Rdata')




library(survival)
dfw <- dfw %>% mutate(VAN = TRT == 'VAN') %>% relocate(VAN, .after=TRT)
keep_labs <- grep('^PRE', names(dfw[-1:-7]), value=TRUE)
x <- numeric()
presDF <- data.frame(
   row.names = keep_labs,
   trt_i = x,
   trt_s = x,
   sur_i = x,
   sur_s = x
)
#labDF <- presDF
rm(x)
for (l in seq_along(keep_labs)) {
   cat(l, keep_labs[l], '\n')
   # treatment model
   i <- glm(formula = paste0('VAN ~ ', keep_labs[l]), data=dfw, family='binomial', subset=dfw$TRT != 'sDAP') %>% summary()
   s <- glm(formula = paste0('VAN ~ ', keep_labs[l]), data=dfw, family='binomial', subset=dfw$TRT != 'iDAP') %>% summary()
   labDF$trt_i[l] <- i$coefficients[keep_labs[l], 'Pr(>|z|)']
   labDF$trt_s[l] <- s$coefficients[keep_labs[l], 'Pr(>|z|)']
   
   # outcome model
   i <- coxph(formula = as.formula(paste0('Surv(time, status) ~ ', keep_labs[l])), data=dfw, subset=which(dfw$TRT != 'sDAP')) %>% summary()
   s <- coxph(formula = as.formula(paste0('Surv(time, status) ~ ', keep_labs[l])), data=dfw, subset=which(dfw$TRT != 'iDAP')) %>% summary()
   labDF$sur_i[l] <- i$coefficients[keep_labs[l], 'Pr(>|z|)']
   labDF$sur_s[l] <- s$coefficients[keep_labs[l], 'Pr(>|z|)']
   
   # presence associated with treatment
   presDF$trt_i[l] <- fisher.test(table(is.na(dfw[[keep_labs[l]]][dfw$TRT != 'sDAP']), dfw$TRT[dfw$TRT != 'sDAP']))$p.value
   presDF$trt_s[l] <- fisher.test(table(is.na(dfw[[keep_labs[l]]][dfw$TRT != 'iDAP']), dfw$TRT[dfw$TRT != 'iDAP']))$p.value

   
   # presence associated with outcome
   presDF$sur_i[l] <- fisher.test(table(is.na(dfw[[keep_labs[l]]][dfw$TRT != 'sDAP']), dfw$time[dfw$TRT != 'sDAP'] < 30))$p.value
   presDF$sur_s[l] <- fisher.test(table(is.na(dfw[[keep_labs[l]]][dfw$TRT != 'iDAP']), dfw$time[dfw$TRT != 'iDAP'] < 30))$p.value
}
dfw <- dfw %>% select(-VAN)

rownames(labDF)[labDF$trt_i < 0.1]
rownames(labDF)[labDF$trt_s < 0.1]

im <- dfw %>%
   filter(TRT != 'sDAP') %>%
   select(PERSON_ID, TRT, !!keep_labs) %>%
   mutate(across(!!keep_labs, ~ tidyr::replace_na(., mean(., na.rm=T)))) %>%
   mutate(trt = TRT == 'VAN') %>% relocate(trt, .after=TRT) %>%
   glm(formula=paste0('trt ~ ', paste(keep_labs, collapse=' + ')), 
       data=., 
       family=binomial(link='logit')) %>%
   summary()
im <- im$coefficients[-1, 'Pr(>|z|)']
labDF$trt_im <- im[match(rownames(labDF), names(im))]

sm <- dfw %>%
   filter(TRT != 'iDAP') %>%
   select(PERSON_ID, TRT, !!keep_labs) %>%
   mutate(across(!!keep_labs, ~ tidyr::replace_na(., mean(., na.rm=T)))) %>%
   mutate(trt = TRT == 'VAN') %>% relocate(trt, .after=TRT) %>%
   glm(formula=paste0('trt ~ ', paste(keep_labs, collapse=' + ')), 
       data=., 
       family=binomial(link='logit')) %>%
   summary()
sm <- sm$coefficients[-1, 'Pr(>|z|)']
labDF$trt_sm <- sm[match(rownames(labDF), names(sm))]

im <- dfw %>%
   filter(TRT != 'sDAP') %>%
   select(PERSON_ID, time, status, !!keep_labs) %>%
   mutate(across(!!keep_labs, ~ tidyr::replace_na(., mean(., na.rm=T)))) %>%
   coxph(formula=as.formula(paste0('Surv(time, status) ~ ', paste(keep_labs, collapse=' + '))), data=.) %>%
   summary()
im <- im$coefficients[, 'Pr(>|z|)']
labDF$sur_im <- im[match(rownames(labDF), names(im))]

sm <- dfw %>%
   filter(TRT != 'iDAP') %>%
   select(PERSON_ID, time, status, !!keep_labs) %>%
   mutate(across(!!keep_labs, ~ tidyr::replace_na(., mean(., na.rm=T)))) %>%
   coxph(formula=as.formula(paste0('Surv(time, status) ~ ', paste(keep_labs, collapse=' + '))), data=.) %>%
   summary()
sm <- sm$coefficients[, 'Pr(>|z|)']
labDF$sur_sm <- sm[match(rownames(labDF), names(sm))]


im <- dfw %>%
   filter(TRT != 'sDAP') %>%
   select(PERSON_ID, TRT, !!keep_labs) %>%
   mutate(across(!!keep_labs, ~ !is.na(.))) %>%
   mutate(trt = TRT == 'VAN') %>%
   glm(formula=paste0('trt ~ ', paste(keep_labs, collapse=' + ')), 
       data=., 
       family=binomial(link='logit')) %>%
   summary()
im <- im$coefficients[-1, 'Pr(>|z|)']
presDF$trt_im <- im[match(rownames(presDF), gsub('TRUE', '', names(im)))]


sm <- dfw %>%
   filter(TRT != 'iDAP') %>%
   select(PERSON_ID, TRT, !!keep_labs) %>%
   mutate(across(!!keep_labs, ~ !is.na(.))) %>%
   mutate(trt = TRT == 'VAN') %>%
   glm(formula=paste0('trt ~ ', paste(keep_labs, collapse=' + ')), 
       data=., 
       family=binomial(link='logit')) %>%
   summary()
sm <- sm$coefficients[-1, 'Pr(>|z|)']
presDF$trt_sm <- sm[match(rownames(presDF), gsub('TRUE', '', names(sm)))]

im <- dfw %>%
   filter(TRT != 'sDAP') %>%
   select(PERSON_ID, time, status, !!keep_labs) %>%
   mutate(across(!!keep_labs, ~ !is.na(.))) %>%
   coxph(formula=as.formula(paste0('Surv(time, status) ~ ', paste(keep_labs, collapse=' + '))),
         data=.) %>%
   summary()
im <- im$coefficients[, 'Pr(>|z|)']
presDF$sur_im <- im[match(rownames(presDF), gsub('TRUE', '', names(im)))]

sm <- dfw %>%
   filter(TRT != 'iDAP') %>%
   select(PERSON_ID, time, status, !!keep_labs) %>%
   mutate(across(!!keep_labs, ~ !is.na(.))) %>%
   coxph(formula=as.formula(paste0('Surv(time, status) ~ ', paste(keep_labs, collapse=' + '))),
       data=.) %>%
   summary()
sm <- sm$coefficients[, 'Pr(>|z|)']
presDF$sur_sm <- sm[match(rownames(presDF), gsub('TRUE', '', names(sm)))]

p_order <- c('trt_i', 'trt_im', 'sur_i', 'sur_im', 'trt_s', 'trt_sm', 'sur_s', 'sur_sm')
labDF <- labDF[p_order]
presDF <- presDF[p_order]

w <- grep('^PRE_', rownames(labDF))
labDF <- labDF[w,]
presDF <- presDF[w,]

presDF[c('PRE_MCHC', 'PRE_MCV'), c('trt_im', 'trt_sm', 'sur_sm', 'sur_im')] <- presDF['PRE_MCH', c('trt_im', 'trt_sm', 'sur_sm', 'sur_im')]
presDF[c('PRE_RDW', 'PRE_WBC'), c('trt_im', 'trt_sm', 'sur_sm', 'sur_im')] <- presDF['PRE_PLATELETS', c('trt_im', 'trt_sm', 'sur_sm', 'sur_im')]
#presDF[c('POST_MCHC', 'POST_MCV', 'POST_PLATELETS', 'POST_RDW', 'POST_WBC'), c('trt_im', 'trt_sm', 'sur_sm', 'sur_im')] <- presDF['POST_MCH', c('trt_im', 'trt_sm', 'sur_sm', 'sur_im')]


x <- as.matrix(cbind(labDF, presDF))
x <- round(x, 3)
width <- 0.5
xpos <- rep(c(1:8, 9.5:16.5), each=nrow(x))
ypos <- rep(nrow(x):1, times=16)

# which labs show up in which places?
{
   par(mfrow=c(1,1), mar=c(10, 15, 4, 2), oma=c(0,0,0,0), cex=0.9)
   plot(NA, xlim=c(0.5, 16.5), ylim=c(0.5, nrow(x)+0.5), xaxs='i', yaxs='i', axes=F, ann=F)
   rect(xleft=xpos-width, xright=xpos+width,
        ybottom=ypos-width, ytop=ypos+width,
        xpd=NA, col=ifelse(x < 0.01, 'darkblue', ifelse(x < 0.05, 'blue', ifelse(x < 0.1, 'lightblue', 'white'))))
   abline(v=c(4.5, 13), lwd=4)
   text(x=xpos, y=ypos, labels=ifelse(x > 0.1, '', ifelse(x < 0.01, '<0.01', x)), col=ifelse(x < 0.05, 'white', 'black'), xpd=NA)
   
   text(x=0.4, y=nrow(x):1, adj=1, xpd=NA, labels=gsub('^PRE_', '', rownames(labDF)))
   text(x=c(mean(1:8), mean(9.5:16.5)), y=(nrow(x)+0.5)*1.075, xpd=NA, labels=c('Lab values', 'Lab presence'))
   text(x=c(mean(1:4), mean(5:8), mean(9.5:12.5), mean(13.5:16.5)), y=(nrow(x)+0.5)*1.05, xpd=NA, labels=rep(c('iDAP', 'sDAP'), times=2))
   text(x=c(seq(1.5,7.5,2), seq(10, 16, 2)), y=(nrow(x)+0.5)*1.025, xpd=NA, labels=rep(c('treatment', 'survival'), times=2))
}

# which variables are associated with treatment AND survival??
lab_names <- rownames(labDF)
conf_labs <- list()

# iDAP
trt <- lab_names[(labDF[, 'trt_i'] < 0.1) | (labDF[, 'trt_im'] < 0.1)]
sur <- lab_names[(labDF[, 'sur_i'] < 0.1) | (labDF[, 'sur_im'] < 0.1)]
conf_labs$iDAP <- union(trt, sur)

# sDAP
trt <- lab_names[(labDF[, 'trt_s'] < 0.1) | (labDF[, 'trt_sm'] < 0.1)]
sur <- lab_names[(labDF[, 'sur_s'] < 0.1) | (labDF[, 'sur_sm'] < 0.1)]
conf_labs$sDAP <- union(trt, sur)

lab_names <- rownames(labDF)
conf_labs_pres <- list()

# iDAP
trt <- lab_names[(presDF[, 'trt_i'] < 0.1) | (presDF[, 'trt_im'] < 0.1)]
sur <- lab_names[(presDF[, 'sur_i'] < 0.1) | (presDF[, 'sur_im'] < 0.1)]
conf_labs_pres$iDAP <- union(trt, sur)

# sDAP
trt <- lab_names[(presDF[, 'trt_s'] < 0.1) | (presDF[, 'trt_sm'] < 0.1)]
sur <- lab_names[(presDF[, 'sur_s'] < 0.1) | (presDF[, 'sur_sm'] < 0.1)]
conf_labs_pres$sDAP <- union(trt, sur)


# what is the missingness of these lab values?
x <- dfw %>%
   filter(TRT != 'sDAP') %>%
   select(!!conf_labs$iDAP)
crl <- sapply(x, function(x) sum(!is.na(x)) / length(x))
crl <- sort(crl)
crp <- sort(apply(x, 1, function(x) sum(!is.na(x)) / length(x)))
par(mfrow=c(1,2))
b <- barplot(crl, horiz=TRUE, xlim=c(0,1), names.arg=NA, xlab='Completion rate')
text(x=-0.05, y=b, adj=1, xpd=NA, labels=names(crl))
barplot(crp, horiz=TRUE, xlim=c(0,1), xlab='Completion rate')
mtext(text='iDAP', outer=TRUE, line=-2)

x <- dfw %>%
   filter(TRT != 'iDAP') %>%
   select(!!conf_labs$sDAP)
crl <- sapply(x, function(x) sum(!is.na(x)) / length(x))
crl <- sort(crl)
crp <- sort(apply(x, 1, function(x) sum(!is.na(x)) / length(x)))
par(mfrow=c(1,2))
b <- barplot(crl, horiz=TRUE, xlim=c(0,1), names.arg=NA, xlab='Completion rate')
text(x=-0.05, y=b, adj=1, xpd=NA, labels=names(crl))
barplot(crp, horiz=TRUE, xlim=c(0,1), xlab='Completion rate')
mtext(text='sDAP', outer=TRUE, line=-2)


dfw <- dfw %>% select(PERSON_ID, ORDER_DAY, TRT, time, status, !!unique(unlist(conf_labs)))
crp <- apply(dfw %>% select(!!unique(unlist(conf_labs))), 1, function(x) sum(!is.na(x)) / length(x))
dfw$LAB_COMP_RATE <- crp


dfw %>%
   filter(TRT != 'sDAP') %>%
   coxph(formula = Surv(time, status) ~ TRT + LAB_COMP_RATE, data=.) %>%
   summary()



# tasks:
#     turn each day into 1 number, turn multiple days into 1 number?
#     take lab results prior to 16 hours after first abx administration
#     relationship between abx admin delay (1 or 2 days?) and value of labs?










