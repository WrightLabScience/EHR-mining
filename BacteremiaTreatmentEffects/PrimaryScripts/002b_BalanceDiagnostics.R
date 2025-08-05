# This script takes the output of script 002a and assesses the degree of 
# covariate balance across treatment groups before and after weighting.
#     statistical summaries of balance (smd, var ratios, KS statistic)
#     effective sample sizes
#     clinically-relevant variable balance visualization

setwd('~/Desktop/EHR-mining/BacteremiaTreatmentEffects/')

library(dplyr)
source(file='~/Desktop/EHR-mining/UsefulDataForCleaning/antibiotic_names/CreateNamedAbxAbbreviations.R')
# source(file='~/Desktop/EHR-mining/UsefulDataForCleaning/UPMC_facilities/handle_UPMC_site_names.R')
# rm(site_info, code2cat, code2cat_fxn, code2fac, code2fac_fxn, code2mixed, code2mixed_fxn, fac2cat, fac2mixed, fac2cat_fxn, fac2mixed_fxn, mixed2cat, mixed2cat_fxn)
Rdata_file_path <- '~/Desktop/EHR/EHR work/RdataFiles/BuildCohortDatasets/CleanedData/'

# Load cohort_info list.
load(file = 'CohortInfo/cohort_info.Rdata')

# Create plotting directories.
dir.create('Plots', showWarnings = FALSE)
sapply(X = c('/002b_BalanceDiagnostics', '/002b_ClinicalDifferences'),
       FUN = function(dir_name) dir.create(paste0('Plots/', dir_name), showWarnings = FALSE))


# Loop over cohorts.
for (cohort in seq_along(cohort_info)) {
   
   # Cohort name.
   clean_bug_name <- cohort_info[[cohort]]$clean_bug_name
   cohort_name <- cohort_info[[cohort]]$cohort_name
   cat(cohort_name, '\n')
   
   # Load previously calculated covariate balancing data.
   load(file = paste0('CohortInfo/balancing_weights_info/', cohort_name, '.Rdata'))
   bal_mat <- balance_list$summary_mat
   bal_mat <- bal_mat[,grep('(AT(T|C)|COMB)m', colnames(bal_mat), invert=TRUE)]
   bal_mat['Geom. Mean Var. Ratio',] <- bal_mat['Geom. Mean Var. Ratio',] - 1
   
   # Load raw data.
   load(file = paste0(Rdata_file_path, cohort_name, '.Rdata'))
   
   # Prep plotting data.
   col_vec <- setNames(balance_list$trt_df$col_vec, 
                       balance_list$trt_df$abx)
   
   
   # Plot balance diagnostics.
   tryCatch(
      {
         # Setup plot space.
         pdf(file = paste0('Plots/002b_BalanceDiagnostics/', cohort_name, '.pdf'), width=10, height=5.5)
         par(mfrow=c(2,3), tck=-0.015, mgp=c(1.5, 0.4, 0), mar=c(3, 6, 2, 2), cex.main=1)
         
         # Effective sample sizes.
         b <- barplot(bal_mat[c('ess.ctrl', 'ess.trt'),], beside=TRUE, horiz=TRUE, names.arg=rep(NA,ncol(bal_mat)), main='Effective sample sizes', 
                      xlab='n', col=col_vec, xlim=c(0, max(bal_mat[c('ess.ctrl', 'ess.trt'),]) * 1.02))
         axis(side=2, at=apply(b,2,mean), labels=colnames(bal_mat), las=1, tick=FALSE)
         abline(v=bal_mat[c('ess.ctrl', 'ess.trt'),'raw'], lty=2, xpd=FALSE)
         
         # Bug/cohort name.
         plot.new()
         title(main=clean_bug_name, line=-4, font.main=4, cex.main=1.5, xpd=NA)
         legend('center', legend=names(col_vec), pt.bg=col_vec, pch=22, pt.cex=3, bty='n', cex=1.5)
         
         # Propensity score distributions (I used to show this before after weighting, but now using non-PS weighting methods).
         plotPSdist <- function(e, z) {
            brks <- seq(0,1,0.02)
            h0 <- hist(e[z == 0L], breaks=brks, plot=FALSE)
            h1 <- hist(e[z == 1L], breaks=brks, plot=FALSE)
            hist(e[z == 0L], col=paste0(col_vec[1], 'aa'), breaks=brks, xlim=c(0,1), ylim=c(0, max(c(h1$counts, h0$counts))), 
                 ylab='', xlab='', xpd=NA, yaxt='n', main='Estimator: GBM')
            hist(e[z == 1L], col=paste0(col_vec[2], 'aa'), breaks=brks, add=TRUE)
            title(xlab='Propensity score', line=1.2)
            title(ylab='Frequency', line=2)
            axis(side=2, las=1)
         }
         plotPSdist(e=balance_list$prop_scores, z=balance_list$treat_vec)
         
         # Balance statistics.
         # xmaxs <- pmin(apply(bal_mat[1:3,], 1, max) * 1.1, c(0.5,0.5,2))
         xmaxs <- bal_mat[1:3,'raw'] * 1.1
         xmins <- c(0,0,0)
         for (i in seq_len(nrow(bal_mat[1:3,]))) {
            v <- bal_mat[i,]
            b <- barplot(v, horiz=TRUE, names.arg=NA, main=rownames(bal_mat)[i], xlim=c(xmins[i], xmaxs[i]), xaxt='n')
            if (i == 3) {
               ax <- axisTicks(usr=c(0,xmaxs[i]), log=F)
               axis(side=1, at=ax, labels=ax+1)
            } else {
               axis(side=1)
            }
            axis(side=2, at=b, labels=names(v), las=1, tick=FALSE)
            abline(v=v['raw'], lty=2, xpd=FALSE)
            if (any(bal_mat[i,] > xmaxs[i])) {
               text(x=bal_mat[i,'raw']*1.025, y=b[bal_mat[i,] > xmaxs[i]], adj=0, xpd=NA, labels=round(bal_mat[i, bal_mat[i,] > xmaxs[i]], 2)+1)
            }
         }
         
         dev.off()
      },
      error = function(e) {
         message(e)
         message('Error between opening pdf and closing. Closing now...')
         dev.off()
      }
   )
   
   
   # Plot clinically relevant variables balance.
   tryCatch(
      {
         pdf(file = paste0('Plots/002b_ClinicalDifferences/', cohort_name, '.pdf'),
             width=5.5, height=6.5, bg='white')
         layout(matrix(c(rep(c(1,2),   times=c(8,6)),
                         rep(c(3,4,5), times=c(6,4,4)), 
                         rep(c(3,6,7), times=c(6,4,4)),
                         rep(c(3,8,9), times=c(6,4,4))), 
                       nrow=4, byrow=T),
                heights=c(1,2/3,2/3,2/3))
         par(tck=-0.015, oma=c(0,0,2,0), mar=c(3,4,2,1.5), mgp=c(1.25,0.25,0), cex.main=1)
         
         panel_letter_idx <- 1L
         
         # Treatment by year.
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
               text(x=b, y=t[1,]+0.05*max(colSums(t)), font=2, labels=paste0(round(t[1,] / colSums(t) * 100, 0), rep(c('%', ''), times=c(1L, ncol(t)-1L))), xpd=NA, col='black')  
            }
            text(x=b, y=-0.03*max(colSums(t)), labels=16:24, xpd=NA, adj=c(0.5, 1))
            axis(side=2, las=1)
            rm(t, b, line, ymax)
         }
            
         # Treatment by facility.
         {
            par(mar=c(3,3,2,3))
            
            cat_vec <- sapply(df %>% dplyr::select(contains('ENC_admit_facility')), which)
            names(cat_vec) <- gsub('ENC_admit_facility_(.+)', '\\1', names(cat_vec))
            cat_vec <- sort(unlist(cat_vec))
            cat_vec <- gsub('[0-9]+', '', names(cat_vec))
            
            t <- t(table(cat_vec, df$treatment)[c('Specialty', 'Rural', 'Regional', 'Community', 'Academic'), names(col_vec)])
            line <- log10(max(colSums(t))) / 1.3
            ymax <- ifelse(textpos == 'above', max(colSums(t))*1.15, max(colSums(t)))
            b <- barplot(t, col=col_vec, yaxt='n', ylab='', names.arg=rep(NA, ncol(t)), ylim=c(0, ymax),#, legend=TRUE, args.legend=list(x='topleft', bty='n'), #space=1,
                         main=paste0(LETTERS[panel_letter_idx], ') Use by UPMC facility category'))
            panel_letter_idx <- panel_letter_idx + 1L
            title(ylab='# cases', line=line, xpd=NA)
            yshift <- 0.04
            text(x=b, y=-(1/2)*2*yshift*max(colSums(t)), labels=colnames(t), srt=30, xpd=NA, adj=1)
            if (textpos == 'above') {
               text(x=b, adj=c(0.5,0), font=2, y=colSums(t)+0.03*max(colSums(t)), labels=paste0(round(t[1,] / colSums(t) * 100, 0), rep(c('%', ''), times=c(1L, ncol(t)-1L))), xpd=NA, col=col_vec[1])
            } else if (textpos == 'within') {
               text(x=b, adj=c(0.5,0), font=2, y=t[1,]+0.03*max(colSums(t)), labels=paste0(round(t[1,] / colSums(t) * 100, 0), rep(c('%', ''), times=c(1L, ncol(t)-1L))), xpd=NA, col='black')               
            }
            axis(side=2, las=1)
            rm(t, b, line, yshift)
         }
         
         # Binary flag prevalence.
         {
            # df$BUG_Hospital_acquired <- df$ENC_time_admit_order > 2
            bin_vars <- c(grep('^BUG_', names(df), value=TRUE),
                          grep('^ABX_early_', names(df), value=TRUE),
                          'ABX_last_2w', 'ICD_Sepsis_1w', 'ENC_bacteremia_1yr', 'ENC_transfer_2_quatern', 
                          'ENC_transfer_during_trt', 'ENC_admit_source', 'ENC_admit_type', 'DEM_female')
            
            if (grepl('MRSA', cohort_name))
               bin_vars <- c('DEM_MRSA_1yr', bin_vars)
            
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
            
            par(mar=c(2.5,7.5,2,0.05))
            plot(NA, xlim=c(-5,105), xaxs='i', ylim=c(0.5,length(bin_vars)+0.5), yaxs='i', xlab='%', yaxt='n', ylab='',
                 main=paste0(LETTERS[panel_letter_idx], ') Binary variable prevalence', paste(rep(' ', 28L), collapse='')))
            panel_letter_idx <- panel_letter_idx + 1L
            abline(v = c(0, 100), lty=2, lwd=0.5)
            sapply(seq_len(ncol(t)), function(i) points(x=t[,i], y=i+c(-1,1)*0.05, col=paste0(col_vec[rownames(t)], 'aa'), pch=16, cex=1.5))
            bug_names <- colnames(t)[grepl('BUG', colnames(t)) & colnames(t) != 'BUG Hospital acquired']
            text(x=-10, y=seq_len(ncol(t)), adj=1, xpd=NA, labels=gsub('^BUG ', '', colnames(t)), font=ifelse(colnames(t) %in% bug_names, 3, 1))
            rm(t, bin_vars)
         }
         
         # Continuous variable distributions.
         {
            cont_vars <- c('DEM_age', 'ABX_firstAbxTime', 'ENC_days_inpt_1yr', 'ABX_numEmpDoses', 'ENC_vasopressors', 'ENC_time_admit_order')
            mains <- c('Age', 'First abx time', 'Inpt days last year', 'Empiric abx doses', 'Vasopressor doses', 'Days b/w admit culture')
            xlabs <- c('Years', 'Hours', 'Days', 'Doses', 'Doses', 'Days')
            xmins <- c(18, -48, 0, 0, -0.25, min(df$ENC_time_admit_order))
            xmaxs <- c(110, 48, 60, 15, 15, max(df$ENC_time_admit_order))
            bws <- c(5, 5, 5, 5, 5, 5)
            par(mar=c(2.5,0.75,2,0.75))
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
                         abline(v = mean(v), col=col_vec[trt], lwd=1.5, lty=2)
                      })
            }
            rm(cont_vars, mains, xmaxs, xmins, xlabs, bws, var, i)
         }
         
         # Final plot stuff.
         xshift1 <- 0.02
         xshift2 <- (nchar(names(col_vec)[1])-3)/200
         mtext(text='vs.', font=2, outer=TRUE, at=0.5-xshift2, side=3, line=0, cex=0.7, adj=0.5)
         mtext(text=names(col_vec)[2], font=2, outer=TRUE, at=0.5-xshift1-xshift2, side=3, line=0, cex=0.7, adj=1, col=col_vec[2])
         mtext(text=names(col_vec)[1], font=2, outer=TRUE, at=0.5+xshift1-xshift2, side=3, line=0, cex=0.7, adj=0, col=col_vec[1])
         mtext(text=clean_bug_name, font=4, outer=TRUE, at=0.5, side=3, line=1, cex=0.7)
         rm(xshift1, xshift2)
      
         dev.off()
      },
      error = function(e) {
         message(e)
         message('Error between opening pdf and closing. Closing now...')
         dev.off()
      }
   )
}




extra_weights <- list(
   'E. faecium (VRE)' = character(),
   'E. faecalis' = character(),
   'E. faecalis switch' = c('EB_ATT', 'EB_ATC'),
   'S. aureus (MRSA)' = c('EB_ATT', 'CBPS_ATT', 'SBW_ATT'),
   'S. aureus (MRSA) switch' = c('EB_ATT', 'CBPS_ATT', 'SBW_ATT'),
   'S. aureus' = c('EB_ATT', 'EB_ATC'),
   'ampC Producers' = c('EB_ATC', 'EB_ATT'),
   'P. aeruginosa' = c('EB_ATT', 'EB_ATC'),
   'P. mirabilis' = character(),
   'K. pneumoniae' = c('EB_ATT'),
   'E. coli' = c('EB_ATT', 'EB_ATC')
)
saveRDS(extra_weights, file='~/Desktop/EHR-mining/BuildCohortDatasets/CohortInfo/treatment/propensity_score_info/which_weighting_schemes_per_cohort.rds')


















