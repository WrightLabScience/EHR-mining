# Parse filtering counts during analytic dataframe building:


library(dplyr)

path <- '~/Desktop/EHR-mining/Scripts/AnalysisScripts/BuildCohortDatasets/info/treatment/cat_output_df_building/'
lf <- list.files(path)
lf <- lf[c(1,4,6,5,8,7,2,3)]

labels <- c(
   "Index blood culture positive for requested pathogen(s)",
   "Removed units resistant to either treatment", 
   "Removed units with adjacent infections that contraindicate either treatment", 
   "Removed units >7 day AST result delay", 
   "Removed units that died <4 days after blood culture", 
   "Removed units <18 years of age", 
   "Removed units without any inpatient encounters data", 
   "Removed units discharged <4 days after blood culture", 
   "Removed units with no antibiotics <48 hours after blood culture", 
   "Removed units not treated with either treatment",
   "Removed patients that did not match specific treatment definition"
)

filterCounts <- vector('list', length(lf))
names(filterCounts)[1:5] <- gsub('([A-Z])_([a-z]+)_([A-Z]+)_([A-Z]+)', '<<i>\\1. \\2</i>  - \\3 vs. \\4>', gsub('_treatment.txt', '', lf[1:5], fixed=TRUE))
names(filterCounts)[6] <- '<<i>S. aureus</i>  (MRSA) - VAN vs. switch to DAP>'
names(filterCounts)[7] <- '<<i>E. faecalis</i>  - VAN vs. switch to AMP>'
names(filterCounts)[8] <- '<<i>E. faecium</i>  (VRE) - DAP vs. LZD>'

for (i in seq_along(lf)) {
   r <- readLines(con=paste0(path, lf[i]))
   
   pat <- '^(.+): ([0-9]+) ?$'
   steps_df <- tibble(
      step = 1:(length(r)-2L),
      # label = gsub(pat, '\\1', r[1:(length(r)-2L)]),
      label = labels,
      n = as.integer(gsub(pat, '\\2', r[1:(length(r)-2L)]))
   )
   
   filterCounts[[i]] <- list('df' = steps_df, 'trt' = setNames(as.integer(gsub(pat, '\\2', r[(length(r)-1):length(r)])),
                                                               paste0('Received ', gsub(pat, '\\1', r[(length(r)-1):length(r)]))))
}

save(filterCounts, file='~/Desktop/EHR-mining/Scripts/AnalysisScripts/BuildCohortDatasets/info/filterCounts.Rdata')







