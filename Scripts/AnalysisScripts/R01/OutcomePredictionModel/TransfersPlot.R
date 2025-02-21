
sites <- names(sort(table(data$ENC_FACILITY1)))

trans <- data %>%
   count(ENC_FACILITY1, ENC_FACILITY2) %>% 
   arrange(ENC_FACILITY1)

trans$y1 <- match(trans$ENC_FACILITY1, sites)
trans$y2 <- match(trans$ENC_FACILITY2, sites)

par(mar=c(3, 6, 1, 6), mgp=c(1.5, 0.4, 0), tck=-0.015)
plot(NA, xlim=c(1, 2), ylim=c(1, length(sites)), yaxt='n', xaxt='n', ylab='', xlab='')
axis(side=1, at=1:2, labels=c('Admit', 'Transfer'))
axis(side=2, at=seq_along(sites), labels=sites, las=1)
axis(side=4, at=seq_along(sites), labels=sites, las=1)
segments(x0 = 1, x1 = 2,
         y0 = trans$y1,
         y1 = trans$y2,
         lwd = log(trans$n + 1) / 1.2)
