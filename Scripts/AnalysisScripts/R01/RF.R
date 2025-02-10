library(dplyr)
library(randomForest)
source('~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/CreateNamedAbxAbbreviations.R')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/R01/bact_outcomes_variables.Rdata')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/R01/bact_features_variables.Rdata')
featuresDF <- full_join(
   x = outcomesDF %>% select(PERSON_ID, ORDER_DAY),
   y = featuresDF,
   by = join_by(PERSON_ID, ORDER_DAY)
) %>%
   mutate(
      across(.cols = !c(PERSON_ID, ORDER_DAY),
             .fns = ~ ifelse(is.na(.), 0L, .)),
   )
target_var <- 'ESBL'


df <- featuresDF
df[target_var] <- outcomesDF[target_var]
X <- df %>% select(-PERSON_ID, -ORDER_DAY, -!!target_var) %>% as.matrix()
labels <- df[[target_var]]

# 500+ features is too many
# how correlated are different features??
cor_mat <- matrix(NA, ncol=ncol(X), nrow=ncol(X), dimnames=list(colnames(X), colnames(X)))
for (i in seq_len(ncol(X))) {
   cat(i, '\n')
   for (j in i:ncol(X)) {
      if (i == j)
         next
      cor_mat[i,j] <- cor(X[,i], X[,j])
   }
}
hist(as.vector(cor_mat), breaks=10)

# PCA ??
pca <- prcomp(X, center=TRUE, scale.=TRUE)
s <- summary(pca)
plot(cumsum(s$sdev^2 / sum(s$sdev^2)), xlim=c(0,100), ylim=c(0,1), pch=16)



prevalence <- sum(labels == 1L) / length(labels)

folds <- sample(rep(1:5, length.out = length(labels)), size=length(labels), replace=FALSE)

probs <- rep(NA_real_, length(labels))
for (i in 1:5) {
   cat(i, '\n')
   rf <- randomForest(
      x = X[folds != i,],
      y = as.factor(labels[folds != i]),
      classwt = c(1, 1 / prevalence),
      ntree = 100,
      mtry = 10,
      nodesize = 1,
      do.trace = TRUE,
      importance = TRUE,
      keep.forest = TRUE
   )
   probs[folds == i] <- predict(rf, newdata=X[folds == i,], type='prob')[,'1']
}


# sort the prediction probabilities and labels so decreasing
o <- order(probs, decreasing=TRUE)
probs <- probs[o]
labels <- labels[o]

# calculate metrics across every possible threshold
tp_cumsum <- cumsum(labels == 1)
fp_cumsum <- cumsum(labels == 0)
total_pos <- sum(labels == 1)
total_neg <- sum(labels == 0)

sens <- tp_cumsum / total_pos  # Sensitivity
spec <- (total_neg - fp_cumsum) / total_neg  # Specificity
prec <- tp_cumsum / (tp_cumsum + fp_cumsum)  # Precision
prec[is.nan(prec)] <- 0

fpr <- 1 - spec
auroc <- sum(( (fpr[-1] - fpr[-length(fpr)]) * (sens[-1] + sens[-length(sens)]) ) / 2)
auprc <- sum(( (sens[-1] - sens[-length(sens)]) * (prec[-1] + prec[-length(prec)]) ) / 2)

plot(x=fpr, y=sens, xlim=c(0,1), ylim=c(0,1), type='l', xaxs='i', yaxs='i', xlab='FPR', ylab='TPR', yaxt='n', xaxt='n')
text(x=0.01, y=0.95, adj=0, labels=round(auroc, 2))
abline(a=0, b=1, lty=3)






















