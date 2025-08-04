library(dplyr)
source(file='~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/CreateNamedAbxAbbreviations.R')
source(file='~/Desktop/EHR-mining/UsefulDataForCleaning/UPMC_facilities/handle_UPMC_site_names.R')
rm(site_info, code2cat, code2cat_fxn, code2fac, code2fac_fxn, code2mixed, code2mixed_fxn, fac2cat, fac2mixed, fac2cat_fxn, fac2mixed_fxn)
Rdata_file_path <- '~/Desktop/EHR/EHR work/RdataFiles/BuildCohortDatasets/CleanedData/'


setwd('~/Desktop/EHR-mining/BuildCohortDatasets/')
plot_path <- paste0(base_path, 'Plots/')
if (!dir.exists(plot_path))
   dir.create(plot_path)

# Load cohort_info list.
load(file = paste0(base_path, 'CohortInfo/cohort_info.Rdata'))


# Loop over each cohort (pathogen(s)/antibiotic(s)) and generate balance diagnostic plots.
for (cohort in seq_along(cohort_info)) {
   
   # Cohort information
   cohort_bug <- cohort_info[[cohort]]$bug
   cohort_abx <- cohort_info[[cohort]]$abx
   cohort_name <- cohort_info[[cohort]]$cohort_name
   clean_bug_name <- cohort_info[[cohort]]$clean_bug_name
   cat(cohort_name, '\n')
   
   # Get clean version of bug name(s) for plotting
   clean_bug_name <- gsub('([A-Z])_([a-z]+)(_[A-Z]+)?', '\\1. \\2', cohort_bug)
   if (any(grepl('([A-Z])_([a-z]+)_([A-Z]+)$', cohort_bug))) 
      clean_bug_name <- paste0(clean_bug_name, gsub('([A-Z])_([a-z]+)_?([A-Z]+)?', ' (\\3)', cohort_bug))
   clean_bug_name <- paste(clean_bug_name, collapse=', ')
   if ('switch' %in% names(cohort_info[[cohort]])) clean_bug_name <- paste0(clean_bug_name, ' switch')
   
   # Load clean data.
   load(file = paste0(Rdata_file_path, treatment_var, '/', cohort_name, '.Rdata'))
   
   # Create named color vector.
   col_vec <- c('#FFaa66', '#3399FF')
   if (treatment_var == 'treatment') {
      names(col_vec) <- names(sort(table(df$treatment)))
   } else if (treatment_var == 'resistance') {
      names(col_vec) <- paste0(unname(abbr[cohort_abx]), '-', c('R','S'))
      df <- df %>%
         mutate(treatment = case_when(
            treatment == 1L ~ paste0(unname(abbr[cohort_abx]), '-R'),
            treatment == 0L ~ paste0(unname(abbr[cohort_abx]), '-S')
         ))
   }
   
   
   #### Plot Clinical Differences ####
   if (TRUE) {
      pdf(file = paste0(plot_path, 'ClinicalDifferences/', cohort_name, '_', treatment_var, '.pdf'), 
          width=5.5, height=6.5, bg='white')
      layout(matrix(c(rep(c(1,2), times=c(8,6)),rep(c(3,4,5), times=c(6,4,4)), rep(c(3,6,7), times=c(6,4,4))), nrow=3, byrow=T))
      par(tck=-0.015, oma=c(0,0,2,0), mar=c(3,4,2,1.5), mgp=c(1.5,0.3,0), cex.main=1)
      
      panel_letter_idx <- 1L
      # normalized use of each treatment by year and facility
      {
         textpos <- c('above', 'within')[2]
         par(mar=c(3,3,2,1))
         t <- t(table(floor(df$DEM_year_decimal) + 2016, df$treatment))
         line <- log10(max(colSums(t))) / 1.3
         ymax <- ifelse(textpos == 'above', max(colSums(t))*1.15, max(colSums(t)))
         b <- barplot(t, col=col_vec, yaxt='n', ylab='', names.arg=rep(NA, ncol(t)), ylim=c(0, ymax), main=paste0(LETTERS[panel_letter_idx], ') Use over time'), xlab='Year (20xx)')
         panel_letter_idx <- panel_letter_idx + 1L
         title(ylab='# cases', line=line, xpd=NA)
         if (textpos == 'above') {
            text(x=b, y=colSums(t)+0.05*max(colSums(t)), font=2, labels=paste0(round(t[1,] / colSums(t) * 100, 0), rep(c('%', ''), times=c(1L, ncol(t)-1L))), xpd=NA, col=col_vec[1])  
         } else if (textpos == 'within') {
            text(x=b, y=t[1,]+0.05*max(colSums(t)), font=2, labels=paste0(round(t[1,] / colSums(t) * 100, 0), rep(c('%', ''), times=c(1L, ncol(t)-1L))), xpd=NA, col='white')  
         }
         text(x=b, y=-0.03*max(colSums(t)), labels=16:24, xpd=NA, adj=c(0.5, 1))
         axis(side=2, las=1)
         rm(t, b, line, ymax)
         
         
         # normalized use rate of each treatment by facility
         par(mar=c(3,3,2,3))
         cat_vec <- sapply(df %>% select(contains('ENC_next')), which)
         names(cat_vec) <- gsub('ENC_next_facility_(.+)', '\\1', names(cat_vec))
         cat_vec <- sort(unlist(cat_vec))
         cat_vec <- gsub('[0-9]+', '', names(cat_vec))
         cat_vec <- mixed2cat_fxn(cat_vec)
         t <- t(table(cat_vec, df$treatment)[c('Rural', 'Regional', 'Community', 'Academic'), names(col_vec)])
         line <- log10(max(colSums(t))) / 1.3
         ymax <- ifelse(textpos == 'above', max(colSums(t))*1.15, max(colSums(t)))
         b <- barplot(t, col=col_vec, yaxt='n', ylab='', names.arg=rep(NA, ncol(t)), ylim=c(0, ymax), legend=TRUE, args.legend=list(x='topleft', bty='n'), #space=1,
                      main=paste0(LETTERS[panel_letter_idx], ') Use by UPMC facility category'))
         panel_letter_idx <- panel_letter_idx + 1L
         title(ylab='# cases', line=line, xpd=NA)
         yshift <- 0.04
         text(x=b, y=-(1/2)*2*yshift*max(colSums(t)), labels=colnames(t), srt=30, xpd=NA, adj=1)
         if (textpos == 'above') {
            text(x=b, adj=c(0.5,0), font=2, y=colSums(t)+0.03*max(colSums(t)), labels=paste0(round(t[1,] / colSums(t) * 100, 0), '%'), xpd=NA, col=col_vec[1])
         } else if (textpos == 'within') {
            text(x=b, adj=c(0.5,0), font=2, y=t[1,]+0.03*max(colSums(t)), labels=paste0(round(t[1,] / colSums(t) * 100, 0), '%'), xpd=NA, col='white')               
         }
         axis(side=2, las=1)
         rm(t, b, line, yshift)
      }
      
      # Binary flag prevalence
      {
         df$BUG_Hospital_acquired <- df$ENC_time_admit_order > 2
         bin_vars <- c(grep('^BUG_', names(df), value=TRUE),
                       grep('^ABX_early_', names(df), value=TRUE), 'ABX_last_2w', 'ICD_Sepsis_1w', 'ENC_vasopressors',
                       'ENC_bacteremia_1yr', 'ENC_transfer_2_quatern', 'ENC_transfer_during_trt', 'ENC_admit_source', 'ENC_admit_type', 'DEM_female')
         if (any(cohort_bug == 'S_aureus_MRSA')) bin_vars <- c('DEM_MRSA_1yr', bin_vars)
         bin_vars <- bin_vars[bin_vars %in% names(df)]
         t <- sapply(X = bin_vars,
                     FUN = function(var) {
                        t <- table(df[[var]], df$treatment)
                        if (var == 'ENC_admit_source') t <- table(df[[var]] == 'NursingHome', df$treatment)
                        if (var == 'ENC_admit_type') t <- table(df[[var]] == 'ED', df$treatment)
                        t['TRUE',] / colSums(t) * 100
                     })
         colnames(t)[colnames(t) == 'ENC_admit_source'] <- 'Nurs. home admit'
         colnames(t)[colnames(t) == 'ENC_admit_type'] <- 'ED admit'
         colnames(t) <- gsub('^(ENC|DEM)_', '', colnames(t))
         colnames(t) <- gsub('_', ' ', colnames(t))
         for (var in c('female', 'transfer 2 quatern', 'bacteremia 1yr')) {
            if (var %in% colnames(t)) {
               colnames(t)[colnames(t) == var] <- stringr::str_to_sentence(var)
            }
         }
         par(mar=c(3,7.5,2,0))
         plot(NA, xlim=c(-5,105), xaxs='i', ylim=c(0.5,length(bin_vars)+0.5), yaxs='i', xlab='%', yaxt='n', ylab='', 
              main=paste0(LETTERS[panel_letter_idx], ') Binary variable prevalence', paste(rep(' ', 28L), collapse='')))
         panel_letter_idx <- panel_letter_idx + 1L
         abline(v = c(0, 100), lty=2, lwd=0.5)
         sapply(seq_len(ncol(t)), function(i) points(x=t[,i], y=i+c(-1,1)*0.05, col=paste0(col_vec[rownames(t)], 'aa'), pch=16, cex=1.5))
         bug_names <- colnames(t)[grepl('BUG', colnames(t)) & colnames(t) != 'BUG Hospital acquired']
         text(x=-10, y=seq_len(ncol(t)), adj=1, xpd=NA, labels=gsub('^BUG ', '', colnames(t)), font=ifelse(colnames(t) %in% bug_names, 3, 1))
         rm(t, bin_vars)
      }
      
      # Continuous variable distributions
      {
         cont_vars <- c('DEM_age', 'ABX_firstAbxTime', 'ENC_days_inpt_1yr', 'ABX_numEmpAbxDoses')
         mains <- c('Age', 'First abx time', 'Inpt days last year', 'Empiric abx doses')
         xlabs <- c('Years', 'Hours', 'Days', 'Doses')
         xmins <- c(18, -48, 0, 0)
         xmaxs <- c(110, 48, 60, 15)
         bws <- c(5, 5, 5, 5)
         par(mar=c(3,0.75,2,0.75))
         for (i in seq_along(cont_vars)) {
            var <- cont_vars[i]
            plot(NA, xlim=c(ifelse(is.na(xmins[i]), min(df[[var]]), xmins[i]), ifelse(is.na(xmaxs[i]), max(df[[var]]), xmaxs[i])),
                 xaxs='i', yaxt='n', ylim=c(0, max(density(df[[var]], bw=bws[i])$y)*1.2), xlab=xlabs[i], ylab='', xaxt='n',
                 main=paste0(LETTERS[panel_letter_idx], ') ', mains[i]))
            axis(side=1, gap.axis = 1e-10)
            panel_letter_idx <- panel_letter_idx + 1L
            sapply(X = names(col_vec),
                   FUN = function(trt) {
                      v = df[[var]][df$treatment == trt]
                      lines(density(v, bw=bws[i]), col=col_vec[trt], lwd=2)
                      abline(v = median(v), col=col_vec[trt], lwd=2, lty=2)
                   })
         }
         rm(cont_vars, mains, xmaxs, xmins, xlabs, bws, var, i)
      }
      
      # Final plot stuff
      {
         xshift1 <- 0.02
         xshift2 <- (nchar(names(col_vec)[1])-3)/200
         mtext(text='vs.', font=2, outer=TRUE, at=0.5-xshift2, side=3, line=0, cex=0.7, adj=0.5)
         mtext(text=names(col_vec)[2], font=2, outer=TRUE, at=0.5-xshift1-xshift2, side=3, line=0, cex=0.7, adj=1, col=col_vec[2])
         mtext(text=names(col_vec)[1], font=2, outer=TRUE, at=0.5+xshift1-xshift2, side=3, line=0, cex=0.7, adj=0, col=col_vec[1])
         mtext(text=clean_bug_name, font=4, outer=TRUE, at=0.5, side=3, line=1, cex=0.7)
         rm(xshift1, xshift2)
      }
      
      dev.off()
   }
   #### END ####
   
   
   #### Propensity score balance diagnostics
   load(file = paste0('~/Desktop/EHR-mining/BuildCohortDatasets/CohortInfo/', treatment_var, '/propensity_score_info/', cohort_name, '.Rdata'))
   
   # if (treatment_var == 'treatment') df$treatment <- as.integer(df$treatment == names(sort(table(df$treatment)))[1])
   estimands <- cohort_info[[cohort]]$estimand
   for (estimand in estimands) {
      panel_letter_idx <- 1L
      pdf(file = paste0(plot_path, 'BalanceDiagnostics/', cohort_name, ifelse(only_sepsis_flags[s], '_sep', ''), '_', estimand, '.pdf'), width=6.5, height=8, bg='white')
      layout(mat=matrix(1:9, nrow=3), heights=c(1,1.5,1.5))
      par(tck=-0.015, oma=c(0,0,2,0), mgp=c(1.25,0.2,0), cex.main=1)
      
      for (ps_model in c('gbm', 'randomForest', 'logistic')) {
         #### PS distribution before weighting ####
         par(mar=c(3.5,3,3,1))
         brks <- seq(0,1,0.02)
         z <- df$treatment == names(col_vec)[1]
         e <- ps[[ps_model]][[estimand]]$prop_scores$e
         h0 <- hist(e[!z], breaks=brks, plot=FALSE)
         h1 <- hist(e[z],  breaks=brks, plot=FALSE)
         hist(e[!z], col=paste0(col_vec[2], 'aa'), breaks=brks, xlim=c(0,1), ylim=c(0, max(c(h1$counts, h0$counts))), ylab='',
              xlab='Propensity score', xpd=NA, main=paste0(LETTERS[panel_letter_idx], ') Pre-adjustment overlap\n', ps_model), yaxt='n')
         title(ylab='Frequency', line=2)
         axis(side=2, las=1)
         hist(e[z], col=paste0(col_vec[1], 'aa'), breaks=brks, add=TRUE)
         rm(brks)
         panel_letter_idx <- panel_letter_idx + 1L
         #### END ####
         
         #### Love plot ####
         balancePlot <- function(baltab, mode, panel_letter_idx, estimand) {
            stopifnot('Rownames are screwed up across bal.tab' = rownames(baltab$unw) == rownames(baltab$adj))
            stopifnot('Balance metric not valid (std.eff.sz or var.ratio)' = mode %in% c('std.eff.sz', 'var.ratio'))
            
            balDF <- data.frame(sapply(baltab, '[[', mode))
            rownames(balDF) <- rownames(baltab$unw)
            names(balDF) <- c('unw', 'adj')
            if (mode == 'var.ratio') {
               w <- which(balDF$unw < 1)
               if (length(w) > 0L) balDF$unw[w] <- 1 / balDF$unw[w]
               w <- which(balDF$adj < 1)
               if (length(w) > 0L) balDF$adj[w] <- 1 / balDF$adj[w]
            } else if (mode == 'std.eff.sz') {
               balDF$unw <- abs(balDF$unw)
               balDF$adj <- abs(balDF$adj)
            }
            balDF <- balDF[order(abs(balDF$unw)), ]
            wna <- which(is.na(balDF$unw))
            if (length(wna) > 0L) balDF[wna,] <- 1
            
            pch_vec <- setNames(c(1, 16), names(balDF))
            yshift <- 1
            
            xmin <- ifelse(mode == 'std.eff.sz', 0, 1)
            xmax <- ifelse(mode == 'std.eff.sz', 1, 4)#max(sapply(balDF, function(v) max(v[!is.infinite(v)]))))
            xshift <- xmax / 100 #ifelse(mode == 'std.eff.sz', 0.01, 0.05)
            xlab <- ifelse(mode == 'std.eff.sz', '|Standardized mean differences|', 'Variance ratios')
            
            plot(NA, ylim=c(1-yshift, nrow(balDF)+yshift), yaxt='n', yaxs='i', xlim=c(xmin-xshift, xmax+xshift), xaxs='i', 
                 xlab=xlab, ylab='', xaxt='n')#, main=paste0(LETTERS[panel_letter_idx], ') Covariate balance - ', estimand))
            title(ylab='Covariates', line=0.5)
            legend('bottomright', legend=names(pch_vec), pch=pch_vec)
            
            sapply(1:2, function(i) {
               bal_vals <- balDF[[i]]
               if (mode == 'var.ratio')
                  bal_vals[bal_vals > xmax] <- xmax
               points(x=bal_vals, y=seq_len(nrow(balDF)), pch=pch_vec[i], cex=0.5)
            })
            
            if (mode == 'std.eff.sz') {
               xpos <- xmin+0.26
               ypos <- seq(60, 24, length.out=6L)
               text(x=xpos, y=ypos[-1], adj=0, labels=colnames(ps[[ps_model]][[estimand]]$ess_mat))
               xshift <- c(0.33, 0.52)
               for (i in seq_along(xshift)) text(x=xpos+xshift[i], y=ypos, adj=c(0, 0.5), labels=c(names(balDF)[i], round(ps[[ps_model]][[estimand]]$ess_mat[i,1:2], 1), round(ps[[ps_model]][[estimand]]$ess_mat[i,3:5], 3)))
               vdash <- c(0.1, 0.2)
               xaxtcks <- seq(0,1,0.2)
               xaxlabs <- as.character(xaxtcks)
               
            } else if (mode == 'var.ratio') {
               vdash <- c(1.5, 2)
               xaxtcks <- seq(xmin, xmax, 1)
               xaxlabs <- paste0(rep(c('', '>'), times=c(length(xaxtcks)-1, 1)), xaxtcks)
            }
            axis(side=1, at=xaxtcks, labels=xaxlabs)
            abline(v=vdash, lty=2)
         }
         
         par(mar=c(2.5,2.4,0.5,0.6))
         balancePlot(ps[[ps_model]][[estimand]]$bal.tab, 'std.eff.sz', panel_letter_idx, estimand); panel_letter_idx <- panel_letter_idx + 1L
         balancePlot(ps[[ps_model]][[estimand]]$bal.tab, 'var.ratio', panel_letter_idx, estimand)
         #### END ####
      }
      
      # Final plot stuff
      xshift1 <- 0.02
      xshift2 <- (nchar(names(col_vec)[1])-3)/200
      mtext(text='vs.', font=2, outer=TRUE, at=0.5-xshift2, side=3, line=0, cex=0.7, adj=0.5)
      mtext(text=names(col_vec)[2], font=2, outer=TRUE, at=0.5-xshift1-xshift2, side=3, line=0, cex=0.7, adj=1, col=col_vec[2])
      mtext(text=names(col_vec)[1], font=2, outer=TRUE, at=0.5+xshift1-xshift2, side=3, line=0, cex=0.7, adj=0, col=col_vec[1])
      mtext(text=paste0(clean_bug_name, ' - ', ifelse(only_sepsis_flags[s], 'septic', ''), ' bacteremia - ', estimand), font=4, outer=TRUE, at=0.5, side=3, line=1, cex=0.7)
      rm(xshift1, xshift2)
      
      dev.off()
      
      # # Top 5 biggest differences before adjustment
      # if (ps_model == 'gbm') {
      #    top10 <- bal.table(ps$gbm[[estimand]])$unw[c('tx.mn', 'ct.mn', 'std.eff.sz')]
      #    w <- which(is.na(top10$std.eff.sz))
      #    if (length(w) > 0L)
      #       top10$std.eff.sz[w] <- 1  
      # } else if (ps_model == 'randomForest') {
      #    top10 <- ps$randomForest[[estimand]]$bal.tab$unw
      #    w <- which(is.na(top10$std.eff.sz))
      #    if (length(w) > 0L)
      #       top10$std.eff.sz[w] <- 1  
      # }
      # top10 <- head(top10[order(abs(top10$std.eff.sz), decreasing=TRUE), ], n=10)
      # 
      # ypos <- 64
      # text(x=xpos, y=ypos+27/9, adj=c(0,0.5), labels='Unadjusted top 10 diffs', font=2)
      # text(x=xpos, y=seq(ypos, ypos-27, length.out=10), adj=c(0,0.5), col=col_vec[as.integer(top10$std.eff.sz < 0) + 1],
      #      labels=gsub('_', ' ', gsub('^(ICD|LAB|ENC|ENC_next)+_', '', rownames(top10))))
   }
}



















