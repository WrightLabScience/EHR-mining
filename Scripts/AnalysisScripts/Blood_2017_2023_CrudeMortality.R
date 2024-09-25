#### JOIN DEATH DATA WITH CONCORDANCE ANALYSIS DATA ####
library(dplyr)
load(file = '~/Desktop/EHR/EHR-mining/UsefulDataForCleaning/data_path_name.Rdata')
load(file = paste0(data_path_name, 'ASTs_AbxAdmin_blood_2017_2023_imputed_flagged.Rdata'))
load(file = paste0(data_path_name, 'CleanedDemographics.Rdata'))

empDF <- empDF %>% filter(substr(ORDER_DAY,1,4) == '2017')

length(unique(empDF$PERSON_ID)) # 4,929
length(unique(dth$PERSON_ID))   # 4,546
length(intersect(empDF$PERSON_ID, dth$PERSON_ID)) # 3,247



empDF <- empDF %>% # 4,024
   left_join(x = .,
             y = dth,
             by = join_by(PERSON_ID),
             relationship = 'many-to-one') %>%
   mutate(AGE = as.integer(ORDER_DAY - DOB) / 365) %>%
   mutate(SURV_TIME = as.integer(DEATH_DATE - ORDER_DAY)) %>%
   filter(!is.na(GENDER))

empDF %>% count(PATIENT_STATUS) # 2300 are alive, 1700 are deceased


x <- empDF$AGE
hist(x, breaks=diff(range(x)), xlab='Years', main = 'Age')

x <- empDF$SURV_TIME

d <- sum(empDF$PATIENT_STATUS == 'DECEASED', na.rm=T) 
sum(x <= 30, na.rm=T) / d  # ~26%
sum(x <= 90, na.rm=T) / d  # ~37%
sum(x <= 365, na.rm=T) / d # ~53%
d <- nrow(empDF)
sum(x <= 30, na.rm=T) / d  # ~11%
sum(x <= 90, na.rm=T) / d  # ~16%
sum(x <= 365, na.rm=T) / d # ~23

hist(x, breaks=diff(range(x, na.rm=T)), xlab='Days', main='Survival time after order')
hist(x/30, breaks=diff(range(x, na.rm=T))/30, xlab='Months', main='Survival time after order')


empDF %>% count(BUG, sort=TRUE) # looks normal

empDF <- empDF %>% mutate(DAYS_ALIVE = as.integer(DEATH_DATE - ORDER_DAY))

empDF %>% count(DEATH_DATE <= RESULT_DAY) # 148 didn't make it to their AST culture results day

t <- table(empDF$DAYS_ALIVE[empDF$DAYS_ALIVE >= 0 & empDF$DAYS_ALIVE <= 60])
any(diff(as.integer(names(t))) > 1)
b <- barplot(t)
x <- sapply(0:max(as.integer(names(t))), function(x) sum(empDF$DAYS_ALIVE[!is.na(empDF$DAYS_ALIVE)] <= x)) / max(t)
lines(x = b, y = x)
rm(b, t, x)



{
   t <- empDF %>%
      #filter(AGE < 60)
      filter(n() >= 30L, .by=BUG) %>% # 3,254 / 4,024
      mutate(M = DAYS_ALIVE <= 30) %>%
      mutate(M = ifelse(is.na(M), FALSE, M)) %>%
      mutate(M = case_when(
         DAYS_ALIVE <= 30 ~ 'm1',
         DAYS_ALIVE > 30 & DAYS_ALIVE <= 180 ~ 'm2',
         .default = 'alive'
      )) %>%
      summarise(n = n(),
                .by = c(BUG, M)) %>%
      tidyr::pivot_wider(names_from = M,
                         values_from = n) %>%
      select(BUG, m1, m2, alive)
   
   bugs <- t$BUG
   t <- as.matrix(t %>% select(m1, m2, alive))
   rownames(t) <- bugs; rm(bugs)
   t <- t[order(t[,'m1'] / rowSums(t)), ]
   t <- t(t)
   colnames(t) <- gsub('^([A-Z])[a-z]+ ([a-z]+)$', '\\1. \\2', colnames(t))
   colnames(t)[colnames(t) == 'Coagulase Negative Staph'] <- 'Coag. neg. Staph'
   colnames(t)[colnames(t) == 'Viridans Streptococci'] <- 'Viridans group Strep'
   n <- colSums(t)
   t <- apply(t, 2, function(x) x / sum(x) * 100)
   xvals <- apply(t, 2, function(x) {
      cs <- cumsum(c(0,x))
      cs[1:3] + diff(cs) / 2
   })
   xvals[nrow(xvals), ] <- xvals[nrow(xvals), ncol(xvals)]
   pcnt <- round(t)
   pcnt[, ncol(pcnt)] <- paste0(pcnt[, ncol(pcnt)], '%')
   bugs <- colnames(t)
    
   {
      par(mar=c(3, 9, 1.5, 4), mgp=c(1.5, 0.5, 0))
      b <- barplot(t, horiz=TRUE, names.arg=rep('',ncol(t)), xlab='%')
      title(main='Mortality of 2017 bloodstream infections', line=0.5)
      text(x=-3, adj=1, y=b, labels=bugs, xpd=NA)
      text(x=xvals, y=rep(b, each=3), labels=pcnt, col=c('white', 'black', 'black'), xpd=NA)
      text(x=xvals[, ncol(xvals)], y=b[length(b)]+diff(b)[1]/1.3, labels=c('<30d', '<6m', '>6m'), xpd=NA)
      text(x=102, y=b, labels=n, adj=0, xpd=NA)
      text(x=102, y=b[length(b)]+diff(b)[1]/1.1, labels='N', font=2, adj=0, xpd=NA)
   }
}
rm(t, tn, b, xvals, bugs, n, pcnt)





######## SURVIVAL ANALYSIS ########
DF <- empDF %>%
   select(PERSON_ID, GENDER, DOB, DEATH_DATE, FLAG0, ORDER_DAY, BUG, AGE) %>%
   group_by(PERSON_ID) %>%  slice_max(ORDER_DAY) %>% ungroup() %>%
   reframe(BUG = unique(BUG),
           FLAG0 = paste(sort(unique(FLAG0)), collapse=' + '),
           # DISC = case_when(
           #    any(FLAG0) == 'DISCORDANT' ~ TRUE,
           #    all(FLAG0 != 'DISCORDANT') & any(FLAG0 == 'CONCORDANT') ~ FALSE
           # ),
           .by = c(PERSON_ID, GENDER, DOB, DEATH_DATE, ORDER_DAY, AGE)) %>%
   mutate(DISC = case_when(
      grepl('DISCORDANT', FLAG0) ~ 'DISCORDANT',
      !grepl('DISCORDANT', FLAG0) & grepl('CONCORDANT', FLAG0) ~ 'CONCORDANT',
      .default = FLAG0
   )) %>%
   mutate(NoEmp = DISC == 'No empiric therapy given') %>%
   mutate(status = ifelse(is.na(DEATH_DATE), 1, 2))

DF$DEATH_DATE[is.na(DF$DEATH_DATE)] <- as.Date('2024-08-31', format='%Y-%m-%d')
DF <- DF %>% mutate(time = as.integer(DEATH_DATE - ORDER_DAY))

makeKPplots <- function(bug='') {
   w <- 1:nrow(DF)
   if (bug != '') w <- which(DF$BUG == bug)
   df <- DF[w,]
   
   if (length(unique(df$DISC)) == 1L) return()
   fit <- survival::survfit(survival::Surv(time, status) ~ DISC, data=df %>% filter(DISC %in% c('CONCORDANT', 'DISCORDANT')))
   pval <- survival::survdiff(survival::Surv(time, status) ~ DISC, data=df)$pvalue
   par(mfrow=c(1,1), oma=c(0,1,0,0.25), mar=c(3.5, 2, 1.5, 0.25), mgp=c(2,0.5,0), tck=-0.01)
   plot(fit, col=c('blue', 'red'), lwd=2, xaxs='i', yaxs='i', ylim=c(0,1), xlim=c(0, 730), xaxt='n', yaxt='n',
        xlab='Months', ylab='Surviving fraction', main='', xpd=NA)
   axis(side=1, at=seq(0, 730, 365/3), labels=seq(0,24,4))
   axis(side=2, las=1)
   abline(v = 365/3, lty=2)
   text(x=365*2-5, y=0.99, adj=c(1,1), labels=paste0('p = ', prettyNum(pval, digits=2, scientific=-1)))
   text(x=365*2-5, y=c(0.95, 0.91), adj=c(1,1), labels=paste0(c('concordant', 'discordant'), ' (n = ', ')'), col=c('blue', 'red'))
   
   par(new=TRUE, mar=c(6, 23, 11, 1), mgp=c(1,0.2,0), tck=-0.025)
   plot(fit, col=c('blue', 'red'), lwd=2, xaxs='i', yaxt='n', yaxs='i', ylim=c(0.75,1), xlim=c(0,30), xlab='Days', ylab='', main='')
   axis(side=2, las=1, at=seq(0.75,1,0.05))
   
   mtext(text=paste0('Empiric therapy concordance (Kaplan-meier curves) - ', ifelse(bug=='', 'all pathogens', bug)), 
         side=3, line=-1.125, outer=TRUE, font=2)
}






df <- empDF %>%
   select(PERSON_ID, GENDER, DOB, DEATH_DATE, FLAG0, ORDER_DAY) %>%
   group_by(PERSON_ID) %>% 
   slice_max(ORDER_DAY) %>%
   ungroup() %>%
   summarise(ORDER_DAY = last(ORDER_DAY),
             DISC = any(FLAG0 == 'DISCORDANT'),
             .by = c(PERSON_ID, GENDER, DOB, DEATH_DATE)) %>%
   mutate(status = ifelse(is.na(DEATH_DATE), 1, 2)) # 1 = censored (still alive = NA death date), 2 = dead

df$DEATH_DATE[is.na(df$DEATH_DATE)] <- as.Date('2024-08-31', format='%Y-%m-%d')
df <- df %>% mutate(time = as.integer(DEATH_DATE - ORDER_DAY))


fit <- survfit(Surv(time, status) ~ DISC, data=df)
lgr <- survdiff(Surv(time, status) ~ DISC, data=df)


{
   par(mfrow=c(1,1), oma=c(0,1,0,0.25), mar=c(3.5, 2, 1.5, 0.25), mgp=c(2,0.5,0), tck=-0.015)
   plot(fit, col=c('blue', 'red'), lwd=2, xaxs='i', yaxs='i', ylim=c(0,1), xlim=c(0, 730), xaxt='n', yaxt='n',
        xlab='Months', ylab='Surviving fraction', main='', xpd=NA)
   axis(side=1, at=seq(0, 730, 365/3), labels=seq(0,24,4))
   axis(side=2, las=1)
   text(x = c(400, 380), y=c(0.8, 0.67), labels=c('concordant', 'discordant'), srt=-8, col=c('blue', 'red'))
   text(x=max(summary(fit)$time)+50, y=0.985, adj=c(1,1), labels='Empiric therapy', font=2)
   abline(v = 365/3, lty=2)
   
   par(new=TRUE, mar=c(6, 25, 13, 1), mgp=c(1,0.2,0), tck=-0.025)
   plot(fit, col=c('blue', 'red'), lwd=2, xaxs='i', yaxt='n', yaxs='i', ylim=c(0.5,1), xlim=c(0,365/3), xlab='Days', ylab='', main='')
   axis(side=2, las=1, at=seq(0.5,1,0.1))
   
   mtext(text='Empiric therapy concordance (Kaplan-meier curves)', side=3, line=-1.125, outer=TRUE, font=2)
}










plot(NA, )



x <- 0:5
y <- seq(0.95, 0.8, length.out=length(x))
s <- 0.02

plot(x,y,pch=16, ylim=c(0.7, 1))
lines(x,y)
v <- getVertices(x, y)
lines(x = v$x,
      y = v$y,
      col = col)
vu <- getVertices(x, y + 1.96 * s)
vl <- getVertices(x, y - 1.96 * s)
lines(x = vu$x, y = vu$y, col = col, lty=2)
lines(x = vl$x, y = vl$y, col = col, lty=2)
polygon(x=c(vu$x, rev(vu$x)), y=c(vu$y, rev(vl$y)), col=paste0(col, '55'))





cox <- coxph(Surv(time, status) ~ DISC, data=df)
plot(survfit(cox))


# Lower median survival time with empirically discordant therapy
# 123 days vs. 162 days
empDF %>%
   filter(AGE > 18, AGE < 65) %>%
   filter(FLAG0 %in% c('CONCORDANT', 'DISCORDANT', 'No empiric therapy given')) %>%
   summarise(n=n(),
             median = median(SURV_TIME, na.rm=T),
             mean = mean(SURV_TIME, na.rm=T),
             sd = sd(SURV_TIME, na.rm=T),
             .by = FLAG0)

empDF %>%
   filter(AGE > 18, AGE < 65) %>%
   filter(FLAGT %in% c('CONCORDANT', 'DISCORDANT', 'No empiric therapy given')) %>%
   summarise(n = n(),
             median = median(SURV_TIME, na.rm=T),
             mean = mean(SURV_TIME, na.rm=T),
             sd = sd(SURV_TIME, na.rm=T),
             .by = FLAGT)


empDF %>%
   filter(AGE > 18, AGE < 65, BUG == 'Staphylococcus aureus') %>%
   mutate(X = OXACILLIN == 1L) %>%
   summarise(n = n(),
             median = median(SURV_TIME, na.rm=T),
             mean = mean(SURV_TIME, na.rm=T),
             sd = sd(SURV_TIME, na.rm=T),
             .by = X)

empDF %>%
   filter(AGE > 18, AGE < 65, BUG == 'Enterococcus faecium') %>%
   mutate(X = VANCOMYCIN == 1L) %>%
   summarise(n = n(),
             median = median(SURV_TIME, na.rm=T),
             mean = mean(SURV_TIME, na.rm=T),
             sd = sd(SURV_TIME, na.rm=T),
             .by = X)

empDF %>%
   filter(AGE > 18, AGE < 65, BUG %in% c('Staphylococcus aureus', 'Enterococcus faecium')) %>%
   filter(FLAG0 %in% c('CONCORDANT', 'DISCORDANT', 'No empiric therapy given')) %>%
   summarise(n = n(),
             median = median(SURV_TIME, na.rm=T),
             mean = mean(SURV_TIME, na.rm=T),
             sd = sd(SURV_TIME, na.rm=T),
             .by = c(FLAG0, BUG))





# patients not given empiric treatments are generally slightly younger
empDF %>%
   mutate(not_given = FLAG0 == 'No empiric therapy given') %>%
   summarise(median = median(AGE, na.rm=T),
             sd = sd(AGE, na.rm=T),
             .by = not_given)
# patients not given targeted therapies are basically same age
empDF %>%
   mutate(not_given = FLAGT == 'No empiric therapy given') %>%
   summarise(median = median(AGE, na.rm=T),
             sd = sd(AGE, na.rm=T),
             .by = not_given)



