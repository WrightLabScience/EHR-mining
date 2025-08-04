library(dplyr)

df <- read.table(file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/Sepsis/SepsisAllAbx.csv', sep='\t', header = TRUE)
df <- tibble(df)
y_true <- as.matrix(df %>% select(TZP, CRO, FEP, MEM) %>% mutate(across(everything(), ~ as.integer(.))))
rm(df)

y_pred <- read.table(file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/Sepsis/NNpreds.csv', sep='\t', header = TRUE)
y_pred <- as.matrix(y_pred)

nsamp <- nrow(y_pred)
nclass <- 4


y_rand <- matrix(sample(0:1, size=nsamp * nclass, replace=TRUE), nrow=nsamp, ncol=nclass)

mat <- matrix(c(table(apply(y_true == y_pred, 1, sum)), 
                table(apply(y_true == y_rand, 1, sum))),
              nrow = 2, byrow=TRUE,
              dimnames = list(c('actual', 'random'), 0:4))

barplot(mat,
        beside=TRUE,
        legend=TRUE,
        xlab='Number of agreed classes per sample')

