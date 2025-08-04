library(dplyr)
source('~/Desktop/WrightLab/UsefulRfuncs/newBarplotFxn.R')
source('~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/CreateNamedAbxAbbreviations.R')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/R01/bact_outcomes_variables.Rdata')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/R01/bact_features_variables.Rdata')
featuresDF <- full_join(
   x = outcomesDF %>% select(PERSON_ID, ORDER_DAY),
   y = featuresDF,
   by = join_by(PERSON_ID, ORDER_DAY)
)

# how many patients have data available for each set of features
df <- data.frame(
   var = featuresDF %>% select(-PERSON_ID, -ORDER_DAY) %>% names,
   pcnt_complete = sapply(X = featuresDF %>% select(-PERSON_ID, -ORDER_DAY), 
                          FUN = function(x) sum(!is.na(x)) / length(x)),
   row.names = NULL
)
df$group <- gsub('D_([A-z0-9]+_[A-z]+(1|C|X))_.+', '\\1', df$var)
df <- df[order(df$pcnt_complete), ]
x <- tapply(df$pcnt_complete, df$group, FUN = function(x) unique(round(x, 2)))
par(mfrow=c(1,1), mar=c(4,4,2,1), tck=-0.015, mgp=c(1.75, 0.5, 0))
b <- barplot(df$pcnt_complete, horiz=TRUE, names.arg=NA, xlim=c(0,1))
y <- tapply(b, df$group, mean)
text(x = x[order(x)], y = y[order(x)], labels = unique(df$group), adj=0)
rm(x, b, y)


# impute missing data as 0 and add flag per group
w <- sapply(unique(df$group), function(g) grep(g, names(featuresDF))[1])
missDF <- tibble(data.frame(sapply(featuresDF[w], function(x) as.integer(is.na(x)))))
names(missDF) <- gsub('^(D_[A-z0-9]+_[A-z]+(1|C|X))_.+', '\\1_missing', names(missDF))
featuresDF <- tibble(cbind(featuresDF, missDF))
rm(df, missDF, w)
featuresDF <- featuresDF %>%
   mutate(
      across(.cols = !c(PERSON_ID, ORDER_DAY),
             .fns = ~ ifelse(is.na(.), 0L, .)),
   )


# how sparse is the matrix? how many 0s/1s? ~ 10% non-zeros
x <- apply(X = featuresDF %>% select(-PERSON_ID, -ORDER_DAY, -contains('missing')),
               MARGIN = 1,
               FUN = sum)
plotBarplot(x)
sum(x > 50) # 14,332 (out of ~40K)
median(x) # 27 is median
mean(x) # 54
rm(x)


# which features are most common?
x <- sapply(featuresDF %>% select(-PERSON_ID, -ORDER_DAY, -contains('missing')), sum)
tail(sort(x), n=24)
hist(x)
sum(x > 10000) # 69
median(x) # 2398
mean(x) # 4206
rm(x)


# how sparse (imbalanced) are the outcomes?
x <- sort(sapply(outcomesDF %>% select(-PERSON_ID, -ORDER_DAY), sum))
data.frame(tail(x, n=24))
barplot(x, horiz=TRUE, names.arg=NA, xlim=c(0, nrow(featuresDF)))
rm(x)


# let's do a large-scale 
featuresDF <- featuresDF %>% arrange(PERSON_ID, ORDER_DAY)
outcomesDF <- outcomesDF %>% arrange(PERSON_ID, ORDER_DAY)

res <- array(
   NA,
   dim = c(
      length(featuresDF) - 2,
      length(outcomesDF) - 2,
      2
   ),
   dimnames = list(
      featuresDF %>% select(-PERSON_ID, -ORDER_DAY) %>% names,
      outcomesDF %>% select(-PERSON_ID, -ORDER_DAY) %>% names,
      c('est', 'pval')
   )
)
for (f in 3:length(featuresDF %>% select(-PERSON_ID, -ORDER_DAY))) {
   cat(f, '')
   feature <- names(featuresDF)[f]
   fvec <- featuresDF[[f]]
   for (o in 3:length(outcomesDF %>% select(-PERSON_ID, -ORDER_DAY))) {
      cat(o, '')
      outcome <- names(outcomesDF)[o]
      ovec <- outcomesDF[[o]]
      
      tab <- table(fvec, ovec)
      if (!(nrow(tab) == 2L & ncol(tab) == 2L))
         next
      ftest <- fisher.test(tab)
      res[feature, outcome, 'est'] <- unname(ftest$est)
      res[feature, outcome, 'pval'] <- ftest$p
   }
   cat('\n')
}
rm(f, o, fvec, ovec, feature, outcome, ftest, tab)
save(res, file='~/Desktop/EHR/EHR work/RdataFiles/R01/correlations.Rdata')
load(file='~/Desktop/EHR/EHR work/RdataFiles/R01/correlations.Rdata')



resDF <- tibble(data.frame(
   feature = rep(rownames(res), times=ncol(res)),
   outcome = rep(colnames(res), each=nrow(res)),
   est = as.vector(res[,,'est']),
   pval = as.vector(res[,,'pval'])
))

resDF %>% count(is.infinite(est), is.na(est), est == 0) # 135 inf, 3596 NA, 2515 0
resDF <- resDF %>% filter(!is.infinite(est), !is.na(est), est != 0)

w <- which(resDF$pval < 2.2e-16) # 315 are 0, 19,657 are < 2.2E-16
resDF$pval[w] <- 2.2e-16
rm(w)


par(mfrow=c(2,2))
hist(log(resDF$est), breaks=100, xlim=c(-5, 5))
abline(v = 0, lty=3, lwd=3)
plot(NA, xlim=c(1, prod(dim(res)[-3])), ylim=c(0.002, 2000), log='y')
abline(h = 1, lty=2)
points(sort(resDF$est), pch=16, cex=0.2)

hist(log(resDF$pval), breaks=100, xlim=c(-36, 0))
plot(NA, xlim=c(1, prod(dim(res)[-3])), ylim=range(resDF$pval), log='y')
abline(h = 0.01 / nrow(resDF), lty=2)
points(sort(resDF$pval), pch=16, cex=0.2)

par(mfrow=c(1,1))
plot(x = log2(resDF$est), y = -log10(resDF$pval),
     xaxs='i', yaxs='i', 
     xlab='log2 odds ratio', ylab='-log10 p-value',
     xlim=c(-3, 5),# xlim=range(log2(resDF$est)),
     ylim=range(-log10(resDF$pval)),
     pch = 16, cex=0.25)
abline(h = -log10(0.001), lty=1, lwd=2)

rm(res, resDF)




























