processTableOne <- function(tOg) {
   pvals <- c(
      setNames(attributes(tOg$ContTable)$pValues[,'pNormal'], 
               rownames(attributes(tOg$ContTable)$pValues)),
      setNames(attributes(tOg$CatTable)$pValues[,'pApprox'], 
               rownames(attributes(tOg$CatTable)$pValues))  
   )
   pvals[is.na(pvals)] <- 1
   smds <- tableone::ExtractSmd(tOg)
   smds <- smds[,1]
   smds[is.na(smds)] <- 0
   pvals <- pvals[match(names(smds), names(pvals))]
   
   signs_cat <- sapply(attributes(tOg$CatTable)$xtabs,
                       function(t) {
                          if (dim(t)[1] == 1L) return(1L)
                          x <- t[2,] / colSums(t)
                          return(sign(diff(x)))
                       })
   signs_cont <- sapply(df[names(attributes(tOg$ContTable)$percentMissing)], 
                        function(t) {
                           x <- tapply(t, df$TRT, mean, na.rm=T)
                           return(sign(diff(x)))
                        })
   signs <- unlist(c(signs_cat, signs_cont))
   names(signs) <- gsub('\\.1', '', names(signs))
   signs <- signs[match(names(pvals), names(signs))]
   if (!all(names(signs) == names(pvals) & names(pvals) == names(smds))) {
      print('table one vectors are out of order!')
   }
   
   table1 <- data.frame(
      row.names = names(smds),
      diffs = smds * signs, 
      pvals
   )
   table1 <- table1[order(table1$diffs), ]
   print(range(table1$diffs))
   return(table1)
}
