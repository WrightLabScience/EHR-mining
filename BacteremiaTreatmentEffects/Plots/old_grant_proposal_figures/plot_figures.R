library(dplyr)
col_vec <- c('#FFaa66', '#3399FF')
grant_figs_dir <- '~/Desktop/EHR-mining/BuildCohortDatasets/Plots/grant_proposal_figures/'
if (!dir.exists(grant_figs_dir)) dir.create(grant_figs_dir)
Rdata_file_path <- '~/Desktop/EHR/EHR work/RdataFiles/BuildCohortDatasets/CleanedData/'
load(file = '~/Desktop/EHR-mining/BuildCohortDatasets/CohortInfo/treatment/cohort_info.Rdata')
source(file='~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/CreateNamedAbxAbbreviations.R')
treatment_var <- 'treatment'


#### Bacteremia Cohort Counts ####
{
   x <- list(
      'E. coli' = c('Ceftriaxone' = 3226, 'Pip-tazo' = 1870),
      'K. pneumoniae' = c('Pip-tazo' = 774, 'Ceftriaxone' = 735),
      'P. mirabilis' = c('Ceftriaxone' = 346, 'Pip-tazo' = 298),
      'P. aeruginosa' = c('Cefepime' = 500, 'Pip-tazo' = 364),
      #'EKPSK' = c('Ceftriaxone' = 5632, 'Meropenem' = 866),
      #'EKPSK' = c('Ceftriaxone' = 4503, 'Pip-tazo' = 3662),
      # 'S. agalactiae' = c('Vancomycin' = 288, 'Ceftriaxone' = 179), # switch?
      # 'S. pneumoniae' = c('Vancomycin' = 166, 'Ceftriaxone' = 166),
      'S. aureus (MSSA)' = c('Cefazolin' = 2553, 'Oxacillin' = 1028),
      #'S. aureus (MSSA)' = c('Cefazolin or Oxacillin' = 1279, 'Vancomycin' = 1632),
      #'S. aureus (MSSA)' = c('Cefazolin or Oxacillin + Vancomycin' = 2576, 'Vancomycin' = 1632),
      'S. aureus (MRSA)' = c('Vancomycin' = 2275, 'Daptomycin' = 434),
      'E. faecalis' = c('Vancomycin' = 766, 'Ampicillin' = 536),
      'E. faecium (VRE)' = c('Daptomycin' = 209, 'Linezolid' = 166)
   )
   t <- matrix(unlist(x), nrow=2, dimnames=list(NULL, names(x)))
   
   grams <- paste0('Gram ', ifelse(names(x) %in% c('E. coli', 'K. pneumoniae', 'P. mirabilis', 'P. aeruginosa', 'EKPSK'), 'nega', 'posi'), 'tives')
   cvec2 <- c('Gram negatives' = '#999999', 'Gram positives' = '#000000')
   
   pdf(file = paste0(grant_figs_dir, 'bacteremia_cohort_counts.pdf'), width=5.25, height=4, bg='white')
   par(mar=c(3, 9, 0.2, 4), tck=-0.01, mgp=c(1.5, 0.25, 0))
   b <- barplot(t, horiz=TRUE, beside=TRUE, names.arg=rep(NA, ncol(t)), col=rev(col_vec), yaxt='n', xlab='No. of cases')
   text(x=-75, y=apply(b,2,mean)+0.2, adj=1, xpd=NA, labels=colnames(t), font=3, col=cvec2[grams])
   xpos <- split(apply(b,2,mean), grams)
   text(x=-0.6 * max(unlist(x)), y=sapply(xpos, mean), labels=unique(grams), xpd=NA, col=cvec2, srt=90)
   text(x=t+25, y=b, adj=0, xpd=NA, labels=as.vector(sapply(x, names)), cex=0.8)
   dev.off()
   
   rm(b, t, x, xpos, cvec2, grams)
}
#### END ####


#### Use over time and space ####
{
   source(file='~/Desktop/EHR-mining/UsefulDataForCleaning/handle_UPMC_site_names.R')
   rm(site_info, code2cat, code2cat_fxn, code2fac, code2fac_fxn, code2mixed, code2mixed_fxn, fac2cat, fac2mixed, fac2cat_fxn, fac2mixed_fxn)
   
   pdf(file = paste0(grant_figs_dir, 'UseByYearAndFacility.pdf'), height=4.5, width=6.3)
   layout(mat=matrix(1:4, nrow=2, byrow=T), widths=c(1.5,1))
   par(tck=-0.01, mgp=c(1.5,0.3,0), cex.main=1, mar=c(2.5,3,2,1))
   panel_letter_idx <- 1L
   
   for (cohort in 6:7) {
      # Cohort information
      cohort_bug <- cohort_info[[cohort]]$bug
      cohort_abx <- cohort_info[[cohort]]$abx
      cohort_name <- paste0(paste(cohort_bug, collapse='_'),
                            ifelse('switch' %in% names(cohort_info[[cohort]]), '_switch', ''), '_',
                            paste(unname(abbr[unlist(strsplit(cohort_abx, split='+', fixed=TRUE))]), collapse='_'), '_',
                            treatment_var)
      
      # Get clean version of bug name(s) for plotting
      clean_bug_name <- gsub('([A-Z])_([a-z]+)(_[A-Z]+)?', '\\1. \\2', cohort_bug)
      if (any(grepl('([A-Z])_([a-z]+)_([A-Z]+)$', cohort_bug))) 
         clean_bug_name <- paste0(clean_bug_name, gsub('([A-Z])_([a-z]+)_?([A-Z]+)?', ' (\\3)', cohort_bug))
      clean_bug_name <- paste(clean_bug_name, collapse=', ')
      if ('switch' %in% names(cohort_info[[cohort]])) clean_bug_name <- paste0(clean_bug_name, ' switch')
      
      # Load clean data.
      load(file = paste0(Rdata_file_path, treatment_var, '/', cohort_name, '.Rdata'))
      
      # Create named color vector.
      names(col_vec) <- names(sort(table(df$treatment)))
      
      t <- t(table(floor(df$DEM_year_decimal) + 2016, df$treatment))
      line <- log10(max(colSums(t))) / 1.3
      ymax <- max(colSums(t))*1.15
      b <- barplot(t, col=col_vec, yaxt='n', ylab='', names.arg=rep(NA, ncol(t)), ylim=c(0, ymax), xlab='')#, xlab='Year (20xx)')
      mtext(text=paste0(LETTERS[panel_letter_idx], ifelse(cohort == 6, ' - Year', '')), side=3, line=ifelse(cohort == 6, 0.25, 0), at=-1.4, adj=0, font=2, cex=0.9)
      panel_letter_idx <- panel_letter_idx + 1L
      title(ylab='# cases', line=line, xpd=NA)
      text(x=b, y=t[1,]+0.06*max(colSums(t)), font=2, labels=paste0(round(t[1,] / colSums(t) * 100, 0), rep(c('%', ''), times=c(1L, ncol(t)-1L))), xpd=NA, col='white')
      text(x=b+0.05, y=-0.04*max(colSums(t)), labels=2016:2024, xpd=NA, adj=1, srt=25)
      axis(side=2, las=1)
      rm(t, b, line, ymax)
      
      
      # normalized use rate of each treatment by facility
      cat_vec <- sapply(df %>% select(contains('ENC_next')), which)
      names(cat_vec) <- gsub('ENC_next_facility_(.+)', '\\1', names(cat_vec))
      cat_vec <- sort(unlist(cat_vec))
      cat_vec <- gsub('[0-9]+', '', names(cat_vec))
      cat_vec <- mixed2cat_fxn(cat_vec)
      t <- t(table(cat_vec, df$treatment)[c('Rural', 'Community', 'Regional', 'Academic'), names(col_vec)])
      line <- log10(max(colSums(t))) / 1.3
      ymax <- max(colSums(t))*1.15
      b <- barplot(t, col=col_vec, yaxt='n', ylab='', names.arg=rep(NA, ncol(t)), ylim=c(0, ymax), legend=TRUE, args.legend=list(x='topleft', bty='n'))
      mtext(text=paste0(LETTERS[panel_letter_idx], ifelse(cohort == 6, ' - UPMC facility type', '')), side=3, line=ifelse(cohort == 6, 0.25, 0), at=-1.4, adj=0, font=2, cex=0.9)
      panel_letter_idx <- panel_letter_idx + 1L
      title(ylab='# cases', line=line, xpd=NA)
      yshift <- 0.04
      text(x=b+0.1, y=-(1/2)*2*yshift*max(colSums(t)), labels=colnames(t), srt=25, xpd=NA, adj=1)
      text(x=b, adj=c(0.5,0), font=2, y=t[1,]+0.015*max(colSums(t)), labels=paste0(round(t[1,] / colSums(t) * 100, 0), '%'), xpd=NA, col='white')
      axis(side=2, las=1)
      rm(t, b, line, yshift)
   }
   
   mtext(text=c('MRSA', 'E. faecalis'), side=3, line=c(-1, -14.8), outer=TRUE, at=0.5, font=c(2,4))
   
   dev.off()
   
   rm(df, cat_vec, clean_bug_name, cohort, cohort_abx, cohort_bug, cohort_name, mixed2cat, panel_letter_idx, ymax, mixed2cat_fxn)
}
#### END ####


#### PS overlap + covariate imbalance ####
{
   # add legend
   pdf(file = paste0(grant_figs_dir, 'PSoverlap_CovariateImbalance.pdf'), height=5.5, width=5)
   layout(mat=matrix(1:6, nrow=3)) #heights=c(1,1.75))
   par(tck=-0.01, cex.main=1)
   panel_letter_idx <- 1L
   cex_mtext <- 0.7
   
   for (cohort in 6:7) {
      # Cohort information
      cohort_bug <- cohort_info[[cohort]]$bug
      cohort_abx <- cohort_info[[cohort]]$abx
      cohort_name <- paste0(paste(cohort_bug, collapse='_'),
                            ifelse('switch' %in% names(cohort_info[[cohort]]), '_switch', ''), '_',
                            paste(unname(abbr[unlist(strsplit(cohort_abx, split='+', fixed=TRUE))]), collapse='_'), '_',
                            treatment_var)
      
      # Get clean version of bug name(s) for plotting
      clean_bug_name <- gsub('([A-Z])_([a-z]+)(_[A-Z]+)?', '\\1. \\2', cohort_bug)
      if (any(grepl('([A-Z])_([a-z]+)_([A-Z]+)$', cohort_bug))) 
         clean_bug_name <- paste0(clean_bug_name, gsub('([A-Z])_([a-z]+)_?([A-Z]+)?', ' (\\3)', cohort_bug))
      clean_bug_name <- paste(clean_bug_name, collapse=', ')
      if ('switch' %in% names(cohort_info[[cohort]])) clean_bug_name <- paste0(clean_bug_name, ' switch')
      
      # Load clean data.
      load(file = paste0(Rdata_file_path, treatment_var, '/', cohort_name, '.Rdata'))
      
      # Load propensity score data.
      load(file = paste0('~/Desktop/EHR-mining/BuildCohortDatasets/CohortInfo/treatment/propensity_score_info/', cohort_name, '.Rdata'))
      ps_model <- 'gbm'
      estimand <- cohort_info[[cohort]]$estimand
      
      # Create named color vector.
      names(col_vec) <- names(sort(table(df$treatment)))
      
      #### PS distribution before weighting ####
      par(mar=c(3.5,3,2,1), mgp=c(1.5,0.15,0))
      brks <- seq(0,1,0.02)
      z <- df$treatment == names(col_vec)[1]
      e <- ps[[ps_model]][[estimand]]$prop_scores$e
      h0 <- hist(e[!z], breaks=brks, plot=FALSE)
      h1 <- hist(e[z],  breaks=brks, plot=FALSE)
      hist(e[!z], col=paste0(col_vec[2], 'aa'), breaks=brks, xlim=c(0,1), ylim=c(0, max(c(h1$counts, h0$counts))), ylab='',
           xlab='', xpd=NA, yaxt='n', main='')
      title(xlab='Propensity score', line=1.2)
      mtext(text=paste0(LETTERS[panel_letter_idx], ifelse(cohort == 6, ' - Pre-adjustment Overlap', '')), side=3, at=-0.245, adj=0, font=2, cex=cex_mtext, line=0.25)
      title(ylab='Frequency', line=2)
      axis(side=2, las=1)
      hist(e[z], col=paste0(col_vec[1], 'aa'), breaks=brks, add=TRUE)
      legend('topright', legend=rev(names(col_vec)), col=rev(col_vec), pch=15, pt.cex=1.4, bty='n')
      rm(brks)
      panel_letter_idx <- panel_letter_idx + 1L
      #### END ####
      
      #### Love plot ####
      balancePlot <- function(baltab, mode, panel_letter_idx, estimand, plot_type) {
         stopifnot('Rownames are screwed up across bal.tab' = rownames(baltab$unw) == rownames(baltab$adj))
         stopifnot('Balance metric not valid (std.eff.sz or var.ratio)' = mode %in% c('std.eff.sz', 'var.ratio'))
         
         xmin <- ifelse(mode == 'std.eff.sz', 0, 1)
         xmax <- ifelse(mode == 'std.eff.sz', 0.65, 3)#max(sapply(balDF, function(v) max(v[!is.infinite(v)]))))
         
         balDF <- data.frame(sapply(baltab, '[[', mode))
         rownames(balDF) <- rownames(baltab$unw)
         names(balDF) <- c('unw', 'adj')
         if (mode == 'var.ratio') {
            w <- which(balDF$unw < 1)
            if (length(w) > 0L) balDF$unw[w] <- 1 / balDF$unw[w]
            w <- which(balDF$adj < 1)
            if (length(w) > 0L) balDF$adj[w] <- 1 / balDF$adj[w]
            w <- which(balDF$unw > xmax)
            if (length(w) > 0L) balDF$unw[w] <- xmax
            w <- which(balDF$adj > xmax)
            if (length(w) > 0L) balDF$adj[w] <- xmax
         } else if (mode == 'std.eff.sz') {
            balDF$unw <- abs(balDF$unw)
            balDF$adj <- abs(balDF$adj)
         }
         balDF <- balDF[order(abs(balDF$unw)), ]
         wna <- which(is.na(balDF$unw))
         if (length(wna) > 0L) balDF[wna,] <- 1
         
         if (plot_type == 'individual') {
            par(mar=c(2.5,1.5,0.5,0.6), mgp=c(1.5,0.15,0))
            pch_vec <- setNames(c(1, 16), names(balDF))
            yshift <- 1
            
            xshift <- xmax / 100 #ifelse(mode == 'std.eff.sz', 0.01, 0.05)
            xlab <- ifelse(mode == 'std.eff.sz', '|Standardized mean differences|', 'Variance ratios')
            
            plot(NA, ylim=c(1-yshift, nrow(balDF)+yshift), yaxt='n', yaxs='i', xlim=c(xmin-xshift, xmax+xshift), xaxs='i',
                 xlab=xlab, ylab='', xaxt='n')#, main=paste0(LETTERS[panel_letter_idx], ') Covariate balance - ', estimand))
            mtext(text=paste0(LETTERS[panel_letter_idx], ifelse(cohort == 6, ' - Covariate Imbalance', '')), side=3, at=-0.07, font=2, adj=0, cex=cex_mtext)
            title(ylab='Covariates', line=0.5)
            legend('bottomright', legend=names(pch_vec), pch=pch_vec)
            
            sapply(1:2, function(i) {
               bal_vals <- balDF[[i]]
               points(x=bal_vals, y=seq_len(nrow(balDF)), pch=pch_vec[i], cex=0.5)
            })
            
            axis(side=1, at=seq(0,1,0.2))
            abline(v=0.2, lty=2)
            
         } else if (plot_type == 'distributions') {
            # brks <- seq(0,1,0.05)
            # hu <- hist(balDF$unw, breaks=brks, plot=FALSE)
            # ha <- hist(balDF$adj, breaks=brks, plot=FALSE)
            # hist(balDF$unw, breaks=brks, col='#87df22aa', xlim=c(0,1), ylim=c(-max(c(hu$counts, ha$counts)), max(c(hu$counts, ha$counts))))
            # hist(balDF$adj, breaks=brks, col='#dd65aaaa', add=TRUE)
            par(mar=c(2,2.5,1,0.6), mgp=c(1.5,0.25,0))
            names(balDF)[names(balDF) == 'unw'] <- 'unweighted'
            names(balDF)[names(balDF) == 'adj'] <- 'inverse weighted'
            vioplot::vioplot(balDF, ylim=c(xmin, xmax), yaxt='n')
            axis(side=2, las=1)
            title(main=ifelse(mode=='std.eff.sz', '|Standardized mean differences|', 'Variance ratio'), line=0.25, xpd=NA)
         }
      }
      
      balancePlot(ps[[ps_model]][[estimand]]$bal.tab, 'std.eff.sz', panel_letter_idx, estimand, 'distributions')
      balancePlot(ps[[ps_model]][[estimand]]$bal.tab, 'var.ratio', panel_letter_idx, estimand, 'distributions')
   }
   
   mtext(text=c('MRSA', 'E. faecalis'), side=3, line=-1, outer=TRUE, at=c(0.25,0.75), font=c(2,4), cex=0.8)
   
   dev.off()
   
   rm(ps, cohort, cohort_abx, cohort_bug, cohort_name, e, estimand, panel_letter_idx, w, wna, z, balancePlot, df, h0, h0, h1, ha, hu, cex_mtext, clean_bug_name, ps_model)
}
#### END ####


#### ATE estimates ####
{
   load(file = '~/Desktop/EHR-mining/BuildCohortDatasets/CohortInfo/treatment/ATE_estimates_list.Rdata')
   
   pdf(file = paste0(grant_figs_dir, 'ATE_estimates.pdf'), height=2.1, width=7.5)
   par(mfrow=c(1,2), tck=-0.01, mgp=c(1.75,0.3,0), mar=c(1.7,3,1.2,3.75), cex.main=1)
   
   for (i in 1:2) {
      col_vec <- ATE_estimates[[i+5]]$info$ABX
      coefs <- ATE_estimates[[i+5]]$COEFS$mort[c('unw', 'gbm'),]
      rownames(coefs) <- c('unweighted', 'inverse weighted')
      
      ylim <- 4
      xshift <- 0.3
      plot(NA, xlim=c(1-xshift, nrow(coefs)+xshift), xaxs='i', ylim=c(1/ylim, ylim), log='y', xaxt='n', ylab='Odds ratio', xlab='', font.main=4, yaxt='n')
      title(main = c('MRSA', 'E. faecalis')[i], font.main=c(2,4)[i], line=0.25)
      axis(side=2, at=c(0.25, 0.5, 1, 2, 4), labels=as.character(c(0.25, 0.5, 1, 2, 4)), las=1)
      abline(h=1, lty=2, lwd=0.5)
      axis(side=1, at=seq_len(nrow(coefs)), labels=rownames(coefs), gap.axis=-1)
      
      labs <- paste0('favors\n', names(col_vec))
      yshift <- ylim ^ 0.4
      text(x=nrow(coefs)+xshift+0.08, y=c(1/yshift, yshift), xpd=NA, labels=labs, col=col_vec, adj=0)
      sapply(1:2, function(i) {
         arrows(x0=nrow(coefs)+xshift, 
                y0=c(ylim^-0.5, 1)[i], 
                y1=c(1, ylim^0.5)[i],
                col=col_vec[i], code=i,
                length=0.15, lwd=2.5, xpd=NA)
      })
      
      yshift <- 1.2
      points(x=seq_len(nrow(coefs)), coefs[,'OR'], pch=16)
      arrows(x0=seq_len(nrow(coefs)), y0=coefs[,'lowerOR'], y1=coefs[,'upperOR'], angle=90, code=3, length=0.05)
      text(x=seq_len(nrow(coefs)), y=coefs[,'upperOR']*yshift, labels=paste0('p=', unname(formatC(coefs[,'pval'], format='g', digits=2, flag='#'))), 
           font=ifelse(coefs[,'pval'] < 0.05, 2, 1), cex=0.9)
      text(x=seq_len(nrow(coefs)), y=coefs[,'lowerOR']/yshift, xpd=NA, cex=0.9,
           labels=paste0(c('', 'a'), 'OR=', round(coefs[,'OR'], 2)))
   }
   
   dev.off()
   
   rm(ATE_estimates, coefs, i, labs, xshift, ylim, yshift)
}
#### END ####


#### Causal Decision Trees ####
load(file = paste0('~/Desktop/EHR-mining/BuildCohortDatasets/CohortInfo/treatment/CATE_estimates_list.Rdata'))
load(file = paste0('~/Desktop/EHR-mining/BuildCohortDatasets/CohortInfo/treatment/ATE_estimates_list.Rdata'))

for (i in 1:2) {
   col_vec <- ATE_estimates[[i+5]]$info$ABX
   coefs <- CATE_list[[i+5]]$ATEs
   coefs <- coefs[order(rownames(coefs)),]
   bush <- CATE_list[[i+5]]$bush
   cohort_name <- gsub('*', '', names(CATE_list)[i+5], fixed=TRUE)
   num_leaves <- nrow(coefs)
   
   pdf(file=paste0(grant_figs_dir, 'CausalTrees_', cohort_name, '.pdf'), height=6, width=7)
   layout(mat=matrix(1:2,nrow=2), height=c(0.6,0.5))
   par(oma=c(0,4,2,3))
   
   rpart.plot::rpart.plot(bush, roundint=FALSE)
   mtext(text=c('MRSA', 'E. faecalis')[i], side=3, line=-2, outer=TRUE, at=0.05, adj=0, cex=1, font=c(2,4)[i])
   # title(main = cohort_name, line=4, xpd=NA, font.main=4)
   
   ylim <- 60
   par(mar=c(4,0,1,0), tck=-0.015, mgp=c(2.5, 0.4, 0))
   plot(NA, xlim=c(1, 4), yaxs='i', ylim=c(1/ylim, ylim), log='y', xaxt='n', ylab='Odds ratio', xlab='', font.main=4, yaxt='n', xpd=NA)
   sq <- c(2, 4)
   # axis(side=2, at=c(1/rev(sq), 1, sq), labels=paste0(rep(c('1/', ''), times=c(length(sq), length(sq)+1L)), c(rev(sq), 1, sq)), las=1)
   axis(side=2, at=c(1/rev(sq), 1, sq), labels=c(0.25, 0.5, 1, 2, 4), las=1)
   abline(h=1, lty=2, lwd=0.5)
   
   yscale <- 15
   axis(side=4, at=1/(ylim/yscale), labels=paste0('favors\n', names(col_vec)[1]), las=1, col.axis=col_vec[1], tick=F)
   axis(side=4, at=ylim/yscale, labels=paste0('favors\n', names(col_vec)[2]), las=1, col.axis=col_vec[2], tick=F)
   
   
   if (num_leaves == 2) xpos <- c(1.96, 3.05)
   if (num_leaves == 3) xpos <- c(1.25, 2.51, 3.75)
   if (num_leaves == 4) xpos <- c(1.25, 2.09, 2.93, 3.77)
   
   if (i == 7) xpos[2:4] <- xpos[2:4] * 0.99
   if (i == 3 | i == 4) xpos[1] <- xpos[1] + 0.03
   if (i == 3) xpos[2:3] <- xpos[2:3] * 1.005
   
   axis(side=1, at=xpos, labels=paste0('leaf: ', seq_len(num_leaves)), xpd=NA, font=3)
   axis(side=3, at=xpos, labels=paste0('leaf: ', seq_len(num_leaves)), xpd=NA, font=3, tick=F)
   #axis(side=1, at=xpos, labels=paste0('ess=', round(as.vector(t(coefs[,'num.trt.ess'])), 2)), col.axis=col_vec[1], line=2, tick=F)
   #axis(side=1, at=xpos, labels=paste0('ess=', round(as.vector(t(coefs[,'num.ctrl.ess'])), 2)), col.axis=col_vec[2], line=1, tick=F)
   
   points(x=xpos, y=coefs[,'OR'], pch=16, cex=1.2)
   arrows(x0=xpos, y0=coefs[,'lowerOR'], y1=coefs[,'upperOR'], code=3, angle=90, length=0.075, lwd=1.25, xpd=NA)
   text(x=xpos, y=coefs[,'upperOR']*1.2, labels=paste0('p=', unname(formatC(coefs[,'pval'], format='g', digits=2, flag='#'))), xpd=NA, font=ifelse(coefs[,'pval'] < 0.05, 2, 1), adj=c(0.5,0))
   text(x=xpos, y=coefs[,'lowerOR']/1.2, labels=paste0('aOR=', round(coefs[,'OR'], 2)), adj=c(0.5,1), xpd=NA)
   
   dev.off()
}

#### END ####






















