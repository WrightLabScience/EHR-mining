plotKP <- function(df, cohort='', trt=c('', ''), ymin=0.5, censor_time=30, cex=1, col_vec=c('VAN'='#0000FF', 'iDAP'='#ef5675', 'sDAP'='#ffa600', 'DAP'='#FF0000'), conf_int=TRUE) {
   getVertices <- function(x, y) {
      x <- c(0, x)
      y <- c(1, y)
      if (x[length(x)] < 730) {
         x <- c(x, 730)
         y <- c(y, y[length(y)])
      }
      return(list(x = c(x[1], rep(x[2:length(x)], each=2)),
                  y = c(rep(y[1:(length(y)-1)], each=2), y[length(y)])))
   }
   plotSurv <- function(x, y, s, d, col_vec, conf_int=TRUE) {
      zval <- qnorm(1- (1-0.95)/2, 0,1)
      for (i in seq_along(d)) {
         xvals <- x[d[[i]]]
         yvals <- y[d[[i]]]
         svals <- s[d[[i]]]
         v <- getVertices(xvals, yvals)
         lines(x = v$x,
               y = v$y,
               col = col_vec[i], lwd=2)
         vu <- getVertices(xvals, yvals + zval * svals)
         vl <- getVertices(xvals, yvals - zval * svals)
         if (conf_int) polygon(x=c(vu$x, rev(vu$x)), y=c(vu$y, rev(vl$y)), col=paste0(col_vec[i], '20'), border=NA)
      }
   }
   sample_sizes <- table(df$TRT)
   fit <- survfit(Surv(time, status) ~ TRT, data = df)
   d <- split(seq_along(summary(fit)$strata), summary(fit)$strata)
   names(d) <- gsub('TRT=', '', names(d))
   s <- summary(fit)
   xvals <- s$time
   yvals <- s$surv
   ste <- s$std.err
   col_vec <- col_vec[match(names(d), names(col_vec))]
   sample_sizes <- sample_sizes[match(names(d), names(sample_sizes))]
   
   ypos <- seq(0.575, 0.52, length.out=3)
   
   plot(NA, ylim=c(ymin, 1), xlim=c(0, censor_time), xaxt='n', yaxt='n', xaxs='i', yaxs='i', xpd=NA, xlab='', ylab='')
   title(ylab='Surviving fraction', line=2.6, cex.lab=cex)
   title(xlab='Days', line=2, cex.lab=cex)
   title(main = cohort, line=1, cex.main=cex)
   plotSurv(xvals, yvals, ste, d, col_vec, conf_int=conf_int)
   axis(side=1, at=seq(0, censor_time, 15), cex.axis=cex)
   axis(side=2, las=1, cex.axis=cex)
   #abline(v=30, lty=2)
   #text(x=2, y=ypos[1], adj=0, labels=paste0('p = ', prettyNum(pval, digits=3, scientific=-1)), font=ifelse(pval<0.05, 2, 1))
   legend('bottomleft', text.col=col_vec, pch=NA, x.intersp=0, bty='n', cex=cex,
          #legend=paste0(names(col_vec), ' (n = ', prettyNum(sample_sizes, big.mark=','), ')'),
          legend=paste0(names(col_vec), ' (n = ', sample_sizes[names(col_vec)], ')'))
}
