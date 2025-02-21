library(dplyr)
Rdata_file_path <- '~/Desktop/EHR/EHR work/RdataFiles/R01/OutcomePredictionModel/'
load(file = paste0(Rdata_file_path, 'data_feat_outcomes_processed_all.Rdata'))
source('~/Desktop/EHR-mining/Scripts/CleaningScripts/featurizeXdayOutcome.R')


# what about a composite outcome that considers that some patients won't get readmitted
# because they died, and vice versa
times <- data %>% 
   select(contains('_time'), matches('^d[2369]')) %>%
   mutate(across(.cols = matches('^d[2369]'),
                 .fns = ~ as.integer(.) - 1L))

minEventTime <- function(x, y, z) {
   if (all(is.na(c(x, y, z))))
      return(NA_integer_)
   return(min(c(x, y, z), na.rm=T))  
}

data <- data %>%
   rowwise() %>%
   mutate(event_time = minEventTime(mortality_time, readmit_time, BSIrecur_time)) %>%
   relocate(event_time, .after=BSIrecur_time)

data <- XdayOutcome(data, col_name='event_time', includes_pid = FALSE)





# get names of target variables
targets <- grep('^d[2369]', names(times), value=TRUE)
rates <- sapply(targets, function(x) sum(times[[x]]) / nrow(times) * 100)


par(mar=c(3, 10, 2, 1), mgp=c(2, 0.5, 0), tck=-0.015)
b <- barplot(rates, horiz=TRUE, names.arg=rep('', length(targets)), xlim = c(0, 100))
text(x=-20, y=sort(tapply(b, gsub('^d[0-9]{2}', '', targets), mean)), labels=unique(gsub('^d[0-9]{2}', '', targets)), xpd=NA, adj=1, font=2)
text(x=-2.5, y=b, labels=unique(gsub('^d([0-9]{2}).+', '\\1 days', targets)), xpd=NA, adj=1)
abline(h = b[seq(length(b) / 4, length(b)-1, length(b) / 4)] + unique(round(as.vector(diff(b)), 1) / 2), lty = 2, lwd = 2)
text(x=rates+1, y=b, adj=0, labels=paste0(round(rates, 1), '%'))



























