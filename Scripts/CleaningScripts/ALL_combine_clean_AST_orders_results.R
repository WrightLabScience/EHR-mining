###############################################################
################# JOIN AST ORDERS AND RESULTS #################
###############################################################
library(dplyr)
library(tidyr)

load(file = '~/Desktop/EHR/EHR work/RdataFiles/AST_orders_clean.Rdata')  # 3,074,667
load(file = '~/Desktop/EHR/EHR work/RdataFiles/AST_results_clean.Rdata') # 1,860,785

astrDF <- astrDF %>% rename(RESULT_DAY = RESULT_DATE)
astoDF <- astoDF %>% mutate(RESULT_DAY = as.Date(substr(RESULT_DATE, 1, 10)))


# first join on BUG for those orders where bug is present 
astDF <- left_join(x = astrDF,
                   y = astoDF %>% filter(!is.na(BUG)), # order HAS bug - 1,523,058
                   by = join_by(PERSON_ID, ORDER_PROC_ID, BUG, RESULT_DAY)) %>%
   relocate(ORDER_DATE, RESULT_DATE, BLOOD, .after=ORDER_PROC_ID) %>%
   relocate(RESULT_DAY, .after=RESULT_DATE)

astDF %>% count(is.na(ORDER_DATE)) # 366,253 still need an order date

# then join without bug on the remaining records in ast that do not have an order_date
# with those orders in asto that do not have a bug
astrDF_unM <- astDF %>% # 366,253
   filter(is.na(ORDER_DATE)) %>%
   select(!c(ORDER_DATE, RESULT_DATE, BLOOD))
astoDF_unM <- astoDF %>% # 1,545,343 <-- 1,551,609 (out of 3,074,667)
   filter(is.na(BUG)) %>% 
   select(-BUG) %>%
   distinct()

astDF_unM <- left_join(x = astrDF_unM,
                       y = astoDF_unM,
                       by = join_by(PERSON_ID, ORDER_PROC_ID, RESULT_DAY)) %>%
   relocate(ORDER_DATE, RESULT_DATE, BLOOD, .after=ORDER_PROC_ID) %>%
   relocate(RESULT_DAY, .after=RESULT_DATE)

astDF_unM <- astDF_unM %>% filter(!is.na(ORDER_DATE)) # 356,225 (from 366,253)
astDF_M   <- astDF     %>% filter(!is.na(ORDER_DATE)) # 1,443,587

# combine initially matched and unmatched rows = 1,846,711
astDF <- rbind(astDF_M, astDF_unM) # 1,850,757
astDF <- astDF %>% arrange(PERSON_ID, ORDER_DATE, RESULT_DATE)

rm(astDF_M, astDF_unM, astoDF_unM, astrDF_unM)
gc()

length(unique(astDF$PERSON_ID))  # 613,829
length(unique(astoDF$PERSON_ID)) # 613,865
length(unique(astrDF$PERSON_ID)) # 615,872
length(unique(astDF$ORDER_PROC_ID))  # 1,548,048
length(unique(astoDF$ORDER_PROC_ID)) # 1,548,185
length(unique(astrDF$ORDER_PROC_ID)) # 1,556,083


rm(astoDF, astrDF)
print(object.size(astDF), units='Mb') # 2,028,7 Mb


###################################################################################
save(astDF, file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_clean_ASTs.Rdata')
###################################################################################





load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/plots_path_name.Rdata')

# RESULT_DELAY
astDF <- astDF %>% mutate(RESULT_DELAY = as.numeric(lubridate::as.duration(RESULT_DATE - ORDER_DATE)) / 86400)
x <- astDF$RESULT_DELAY
length(x[x < 1]) / length(x) * 100  # 0.03%
length(x[x > 60]) / length(x) * 100 # 0.51%
length(x[x > 5]) / length(x) * 100  # 11%
length(x[x > 14]) / length(x) * 100 # 1.8%
length(x[x > 30]) / length(x) * 100 # 0.94%

median(x) # 2.7
median(x[astDF$BLOOD]) # 3.6

length(x[!astDF$BLOOD][x[!astDF$BLOOD] > 5]) / length(x[!astDF$BLOOD]) * 100  # 8.9%
length(x[x > 5]) / length(x) * 100  # 11.1%
length(x[astDF$BLOOD][x[astDF$BLOOD] > 5]) / length(x[astDF$BLOOD]) * 100  # 34.8%

# which bugs are responsible for the increased delay of blood cultures?
ecb <- astDF$RESULT_DELAY[astDF$BLOOD & astDF$BUG == 'Escherichia coli']
ec <- astDF$RESULT_DELAY[!astDF$BLOOD & astDF$BUG == 'Escherichia coli']
median(ec)  # 2.20
median(ecb) # 2.78 (~14 hours)
rm(ec, ecb)

sab <- astDF$RESULT_DELAY[astDF$BLOOD & astDF$BUG == 'Staphylococcus aureus']
sa <- astDF$RESULT_DELAY[!astDF$BLOOD & astDF$BUG == 'Staphylococcus aureus']
median(sa)  # 2.75
median(sab) # 3.16 (~9.5 hours)
rm(sa, sab)

{
   pdf(file = paste0(plots_path_name, 'AST_result_delay.pdf'))
   par(mfrow = c(2, 1), mar=c(4,3,2,1), mgp=c(2, 0.5, 0), tck=-0.015)
   hist(x, xlim=c(0, 12), breaks=diff(range(x))*24, xlab='Days', main='Time between AST order and result - overall')
   hist(x[astDF$BLOOD],  xlim=c(0, 12), breaks=diff(range(x))*24, xlab='Days', main='Time between AST order and result - blood cultures')
   dev.off()
}
rm(x)
astDF <- astDF %>% select(-RESULT_DELAY)



# WHICH BUGS?
t <- astDF %>% select(BUG, BLOOD) %>% table()
n <- 18
t <- t(head(t[order(t[,1], decreasing=TRUE), 2:1], n))
colnames(t) <- gsub('^([A-Z])[a-z]+ ([a-z]+)$', '\\1\\. \\2', colnames(t))
colnames(t)[colnames(t) == 'Coagulase Negative Staph'] <- 'Coag Neg Staph'
x <- barplot(rep(1,n), horiz=TRUE, plot=FALSE)
{
   pdf(file = paste0(plots_path_name, 'BugCounts.pdf'), height=6)
   par(mar=c(4, 8, 1, 2), mgp=c(2, 0.5, 0), tck=-0.015)
   barplot(t, horiz=TRUE, names.arg=rep('', n), xlim=c(0,6e5), xlab='Number of times isolated',
           legend.text = c('blood', 'non-blood'))
   text(x=-10000, y=x, adj=c(1, 0.5), labels=colnames(t), xpd=NA)
   dev.off()
}
rm(t, x, n)



num_year <- astDF %>%
   mutate(year = substr(ORDER_DATE,1,4)) %>%
   select(year) %>%
   table()
num_year <- num_year[-length(num_year)]

num_year_b <- astDF %>%
   filter(BLOOD) %>%
   mutate(year = substr(ORDER_DATE,1,4)) %>%
   select(year) %>%
   table()
num_year_b <- num_year_b[-length(num_year_b)]

{
   pdf(file = paste0(plots_path_name, 'NumberCulturesPerYear.pdf'), height=5)
   par(mar = c(4, 5.5, 1, 1), mgp=c(2.2, 0.6, 0), tck=-0.01)
   plot(x = as.integer(names(num_year)), y = num_year, type = 'b', pch=16,
        xlab = 'Year', ylab='', yaxt='n')
   title(ylab = 'Number of cultures', line=4)
   points(x = as.integer(names(num_year_b)), y = num_year_b, type='b', pch=16, col='red')
   axis(side = 2, at=c(11000, seq(0, 150000, 25000)), las=1, gap.axis=1e-11)
   abline(h = 11000, lty=2)
   legend('topleft', legend = c('overall', 'blood'), pch=16, lty=1, col=c('black', 'red'))
   dev.off()
}
rm(num_year, num_year_b)




order_hour <- table(substr(astDF$ORDER_DATE, 12,13))
order_hour <- order_hour[names(order_hour) != '']
result_hour <- table(substr(astDF$RESULT_DATE, 12,13))
result_hour <- result_hour[names(result_hour) != '']
names(order_hour) <- gsub('^0', '', names(order_hour))
names(result_hour) <- gsub('^0', '', names(result_hour))
{
   pdf(file = paste0(plots_path_name, 'TimeOfDayASTordersResults.pdf'), height=4, width=9)
   par(mfrow=c(1,2), mgp=c(2, 0.6, 0), mar=c(4, 3, 2, 1))
   barplot(order_hour, main='AST orders', xlab = 'Hour of day')
   barplot(result_hour, main='AST results', xlab = 'Hour of day')
   dev.off()
}
rm(order_hour, result_hour)







