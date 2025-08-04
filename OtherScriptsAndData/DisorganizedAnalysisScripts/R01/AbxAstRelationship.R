library(dplyr)
source('~/Desktop/WrightLab/UsefulRfuncs/newBarplotFxn.R')
load(file = '~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/abbr_named.Rdata')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_clean_ASTs.Rdata')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/ALL_CLEANED_2017_2023_AbxAdmin.Rdata')
abx <- abxDF %>%
   select(PERSON_ID, START_DAY, ABX) %>%
   distinct() %>% # 4,938,854
   mutate(ABX_ABBR = unname(abbr[ABX])) %>%
   mutate(ABX_ABBR = ifelse(is.na(ABX_ABBR), ABX, ABX_ABBR)) %>%
   select(-ABX) %>%
   rename(ABX = ABX_ABBR)
rm(abxDF); gc()



# Remove any Coagulase-negative Staph species if they appear in blood cultures
astDF <- astDF %>% filter(BLOOD)
astDF <- astDF %>% filter(!(grepl('Staph', BUG) & !grepl('Staphylococcus aureus|Staphylococcus lugdunensis', BUG))) # 94,808
astDF <- astDF %>% filter(!grepl('Cryptococcus|Aspergillus|Candida', BUG)) # 92,090
astDF <- astDF %>% filter(lubridate::year(ORDER_DAY) %in% 2017:2023) # 52,487

# how many patients have abx admin events prior to their AST?
# get unique ASTs
ast <- astDF %>% 
   select(PERSON_ID, ORDER_DAY) %>% 
   distinct() %>% # 48,314
   arrange(PERSON_ID, ORDER_DAY)

ast <- ast %>%
   group_by(PERSON_ID) %>%
   mutate(DAYS_SINCE = as.integer(ORDER_DAY - lag(ORDER_DAY))) %>%
   ungroup()

hist(ast$DAYS_SINCE, xlim=c(0, 730), breaks=2000, ylim=c(0,80))
hist(ast$DAYS_SINCE, xlim=c(0, 30), breaks=2000, ylim=c(0,500))

ast <- ast %>% filter(is.na(DAYS_SINCE) | DAYS_SINCE > 30L) # 40,367


# join with Abx Admin data
df <- ast %>%
   select(-DAYS_SINCE) %>%
   left_join(
      x = .,
      y = abx,
      by = join_by(
         PERSON_ID, y$START_DAY < x$ORDER_DAY
      )
   ) %>%
   mutate(X = as.integer(ORDER_DAY - START_DAY)) %>%
   mutate(pastABX = !is.na(ABX)) %>%
   select(-ABX, -START_DAY) %>%
   distinct()

plotBarplot(df$X, xmin=0, xmax=365)


# how many have before the X days leading up to their AST order?
getPcnt <- function(d, w) {
   x <- df[w,] %>% filter(is.na(X) | X >= d) %>% select(-X) %>% distinct()
   sum(x$pastABX) / nrow(x)
}
d <- 0:365
years <- 2017:2023
p <- matrix(NA_real_, nrow=length(d)+1, ncol=length(years))
Ns <- sapply(years, function(y) df %>% filter(lubridate::year(ORDER_DAY) >= y) %>% select(PERSON_ID, ORDER_DAY) %>% distinct() %>% nrow())
for (i in seq_along(years)) {
   cat(i, '-\n')
   w <- which(lubridate::year(df$ORDER_DAY) >= years[i])
   for (j in c(0, seq_along(d))) {
      cat(j, '')
      p[j, i] <- getPcnt(d=j, w)
   }
}

par(mar=c(5,4,3,3))
plot(NA,
     ylim = c(0.5, 0.8),
     xlim = c(1, 450),
     ylab = 'Fraction',
     xlab = 'Days since AST order',
     main = 'How many AST orders are preceded by abx admin by X days?')
for (i in seq_len(ncol(p))) {
   lines(x = seq_len(nrow(p)),
         y = p[,i])
}
text(x = 368,
     y = p[nrow(p)-1,],
     labels = paste0(years, ' - ', Ns),
     adj = 0)



# A lot of the abx admin are "recent"
# meaning they were within 30d, 90d, 1 year of AST order
# therefore, restricting years of AST orders to ensure there
# are ample data from years prior for abx admin is unecessary

# let's start by keeping every AST >= 2018 
# and then look at abx admin more than 30 days before that
df <- ast %>%
   filter(lubridate::year(ORDER_DAY) >= 2018) %>%
   select(-DAYS_SINCE) %>%
   mutate(ORDER_DAY_m30 = ORDER_DAY - 30) %>%
   left_join(
      x = .,
      y = abx,
      by = join_by(
         PERSON_ID, ORDER_DAY_m30 > START_DAY
      )
   ) %>%
   mutate(X = as.integer(ORDER_DAY - START_DAY))


plotBarplot(df$X, xmin=0, xmax=365)

df %>% filter(is.na(ABX)) %>% group_by(PERSON_ID, ORDER_DAY) # 12,826 / 35,067 have no antibiotics


# how many days ago were all of these antibiotics prescribed and how many different?
df %>% group_by(PERSON_ID, ORDER_DAY)

windows <- 42 * 2^(1:6)

t <- tapply(X = df$X,
            INDEX = df$PERSON_ID,
            FUN = function(x) {
               sapply(windows, function(r) sum(x <= r))
            })
m <- matrix(unlist(t), nrow=length(windows))
m[is.na(m)] <- 0
rownames(m) <- windows

plot(NA,
     xlim = range(windows), log = 'x',
     ylim = c(0, max(m)),
     ylab = 'Number of Antibiotic-Days <= window',
     xlab = 'Window')
apply(m, 2, FUN = function(x) {
   points(x=windows, y=x, pch=16, cex=0.2, col='#00000055')
   lines(x=windows, y=x, lwd=0.2, col='#00000055')
})









# add ESBL flag to predict
ast <- ast %>%
   left_join(
      x = .,
      y = astDF %>% select(PERSON_ID, ORDER_DAY, ESBL) %>% distinct(),
      by = join_by(
         PERSON_ID,
         ORDER_DAY
      )
   )

df <- ast %>%
   filter(lubridate::year(ORDER_DAY) >= 2018) %>%
   select(-DAYS_SINCE) %>%
   mutate(ORDER_DAY_m30 = ORDER_DAY - 30) %>%
   left_join(
      x = .,
      y = abx,
      by = join_by(
         PERSON_ID, ORDER_DAY_m30 > START_DAY
      )
   ) %>%
   mutate(X = as.integer(ORDER_DAY - START_DAY)) %>%
   select(-ORDER_DAY_m30, -START_DAY)


cutoffs <- quantile(df$X, seq(0, 1, 0.2), na.rm=T, names=F)
bins <- mapply(cutoffs[-length(cutoffs)],
               cutoffs[-1],
               FUN = function(c1, c2) {
                  which(df$X >= c1 & df$X < c2)
               })
df$Y <- NA
for (b in seq_along(bins)) {
   df$Y[bins[[b]]] <- b
}
df %>% summarise(n(), min(X), max(X), .by=Y) %>% arrange(Y)
df <- df %>%
   mutate(X = 1) %>% 
   distinct()

dfw <- df %>%
   tidyr::pivot_wider(
      id_cols = c(PERSON_ID, ORDER_DAY, ESBL),
      names_from = c(ABX, Y),
      values_from = X,
      values_fill = 0
   ) %>%
   select(-`NA_NA`)

f <- sapply(dfw %>% select(-PERSON_ID, -ORDER_DAY, -ESBL), 
            function(x) sum(x) / length(x))
f <- sort(f)
plot(f, log='y')
# dfw <- dfw %>% select(PERSON_ID, ORDER_DAY, ESBL, !!names(f[f > 0.01]))


pvals <- sapply(dfw %>% select(-PERSON_ID, -ORDER_DAY, -ESBL), 
                function(x) fisher.test(table(x, dfw$ESBL))$p.value)
plot(sort(pvals), log='y')
abline(h = 0.01)

names(dfw) <- gsub('/|,| |-', '_', names(dfw))
# large model, takes a while to run
fit <- glm(formula = as.formula(paste0('ESBL ~ ',
                                       paste(dfw %>% 
                                                select(-PERSON_ID, -ORDER_DAY, -ESBL) %>% 
                                                names,
                                             collapse = ' + '))),
           family = binomial(link = 'logit'),
           data = dfw)
summary(fit)


# assess model predictions on training set
preds <- predict(fit, 
                 newdata=dfw %>% 
                    select(-PERSON_ID, -ORDER_DAY, -ESBL),
                 type='response')
labels <- dfw$ESBL
prevalence <- sum(labels == 1L) / length(labels)

o <- order(preds, decreasing=TRUE)
preds <- preds[o]
labels <- labels[o]
rm(o)

wil <- wilcox.test(preds[labels == 1L], preds[labels == 0L])
W <- wil$statistic
m <- sum(labels == 1L)
n <- sum(labels == 0L)

# plot(NA, xlim=c(0, 1), ylim=c(0, 100))
# lines(density(preds[dfw$ESBL == 0L]))
# lines(density(preds[dfw$ESBL == 1L]))

thresholds <- c(0, unique(preds), 1) # 16,037

tp_cumsum <- cumsum(labels == 1)
fp_cumsum <- cumsum(labels == 0)
total_pos <- sum(labels == 1)
total_neg <- sum(labels == 0)

sens <- tp_cumsum / total_pos  # Sensitivity
spec <- (total_neg - fp_cumsum) / total_neg  # Specificity
prec <- tp_cumsum / (tp_cumsum + fp_cumsum)  # Precision

# Handle edge cases for precision (avoid division by zero)
prec[is.nan(prec)] <- 0

# plot
par(mfrow=c(2,1), mgp=c(2, 0.5, 0), tck=-0.015)
plot(x=1 - spec, y=sens, xlim=c(0,1), ylim=c(0,1), type='l', xaxs='i', yaxs='i')
abline(a=0, b=1, lty=3)
plot(x=sens, y=prec, xlim=c(0,1), ylim=c(0,1), type='l', xaxs='i', yaxs='i')
abline(h = prevalence, lty=3)

(W - m * (m + 1) / 2) / (m * n)
# 0.6376054 
#



thresholds <- c(0, unique(preds), 1)
sens <- spec <- prec <- rep(NA, length(thresholds))
for (i in seq_along(thresholds)) {
   if (i %% 100 == 0) cat(i, '')
   pred_labels <- as.integer(preds > thresholds[i])
   conf_mat <- table(labels, pred_labels)[2:1,,drop=FALSE]
   if (ncol(conf_mat) == 1L) {
      if (colnames(conf_mat) == '1') {
         sens[i] <- 1
         spec[i] <- 0
         prec[i] <- conf_mat[1,1] / sum(conf_mat[,1])
         next
      } else {
         sens[i] <- 0
         spec[i] <- 1
         prec[i] <- 0
         next
      }
   } else {
      conf_mat <- conf_mat[, 2:1]
   }
   sens[i] <- conf_mat[1,1] / sum(conf_mat[1,])
   spec[i] <- conf_mat[2,2] / sum(conf_mat[2,])
   prec[i] <- conf_mat[1,1] / sum(conf_mat[,1])
}

plot(x = 1 - spec, y = sens, xlim=c(0,1), ylim=c(0,1), type='l', xaxs='i', yaxs='i')
abline(a=0, b=1, lty=3)
plot(x = sens, y = prec, xlim=c(0,1), ylim=c(0,1), type='l', xaxs='i', yaxs='i')



tp_cumsum <- cumsum(labels == 1)
fp_cumsum <- cumsum(labels == 0)
total_pos <- sum(labels == 1)
total_neg <- sum(labels == 0)









