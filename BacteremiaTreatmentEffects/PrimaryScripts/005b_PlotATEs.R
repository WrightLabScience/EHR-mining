

# Load ATE estimates and plot
setwd('~/Desktop/EHR-mining/BuildCohortDatasets/')
load(file = paste0('CohortInfo/cohort_info.Rdata')) # Cohort_info list
load(file = paste0('CohortInfo/ATE_estimates_list.Rdata')) # ATE estimates list
plot_path <- 'Plots/ATEestimates/'
if (!dir.exists(plot_path))
   dir.create(plot_path)

for (cohort in seq_along(ATE_estimates)) {
   
   # Cohort info.
   cohort_name <- cohort_info[[cohort]]$cohort_name
   clean_bug_name <- cohort_info[[cohort]]$clean_bug_name
   cat(cohort_name, '\n')
   
   # Load propensity score info for extra stats (ess, col_vec).
   load(file = paste0('CohortInfo/propensity_score_info/', cohort_name, '.Rdata'))
   
   col_vec <- setNames(balance_list$trt_df$col_vec, 
                       balance_list$trt_df$abx)
   #ess_mat <- ATE_estimates[[cohort]]$info$ESS
   bug_name <- names(ATE_estimates)[cohort]
   if (bug_name == 'S. aureus') bug_name <- paste0(bug_name, ' (MSSA)')
   
   
   tryCatch(
      {
         pdf(file = paste0(plot_path, bug_name, '.pdf'), height=7.5, width=5)
         par(mfrow=c(3,1), tck=-0.015, mgp=c(2, 0.25, 0), mar=c(7.25, 3.5, 1.25, 3))
         cex <- 1
         
         # Loop over outcomes and plot
         for (o in seq_along(ATE_estimates[[cohort]])) {
            if (o == 2) {
               o <- 1
               # plot.new()
               # next
            }
            
            # Get treatment effect estimates, 95% CI, and p-values for this outcome.
            coefs <- ATE_estimates[[cohort]][[o]]
            mode <- c('uni', 'mlt', 'both')[3]
            if (mode != 'both')
               coefs <- coefs[c('unw uni', grepv(paste0('wtd ', mode), rownames(coefs))),]
            coefs <- coefs[c(1, order(gsub('^wtd (uni|mlt) ([A-z0-9]{3})([A-z0-9]+)?( - [A-z]+)?$', '\\2', rownames(coefs)[-1]))+1), ]
            coefs <- coefs[grep('.+ (AT(T|C)|COMB)m', rownames(coefs), invert=TRUE), ]
            
            # Setup plot area.
            ylim <- max(c(abs(coefs[,'lowerCI']), coefs[,'upperCI'])) * 1.15
            xshift <- nrow(coefs) * 0.02#9
            plot(NA, xlim=c(1-xshift, nrow(coefs)+xshift), xaxs='i', ylim=c(-ylim, ylim), xaxt='n', ylab='% difference', xlab='', font.main=4, yaxt='n')
            axis(side=2, las=1)
            abline(h=0, lty=2, lwd=0.5)
            # abline(v=c(1.5,3.5,5.5), lty=2, lwd=0.5)
            
            axis(side=1, at=seq_len(nrow(coefs)), labels=gsub('wtd ', '', rownames(coefs)), las=2)
            # axis(side=1, at=seq_len(nrow(coefs)), labels=c('raw', gsub('^wtd uni (.+) - [A-z]+$', '\\1', rownames(coefs)[-1])), gap.axis=-10)
            # axis(side=1, at=seq_len(nrow(coefs)), labels=rep('',nrow(coefs)))
            # text(x=seq_len(nrow(coefs)), y=-ylim*1.2, labels=rownames(coefs), adj=1, srt=30, xpd=NA)
            
            ypos <- ylim * 1.2
            text(x=1-xshift, y=ypos, font=2, adj=0, xpd=NA, labels=c('30-day mortality', 'Composite', 'Length of stay* > median')[o])
            text(x=nrow(coefs)+xshift, y=ypos, font=4, adj=1, xpd=NA, labels=bug_name)
            
            labs <- paste0('favors\n', names(col_vec))
            text(x=nrow(coefs)+xshift+nrow(coefs)*0.01, y=ylim * c(1,-1) * 0.5, xpd=NA, labels=labs, col=col_vec, adj=0)
            
            yshift <- 0.085
            points(x=seq_len(nrow(coefs)), coefs[,'estimate'], pch=16, col=ifelse(coefs[,'pvalue'] < 0.05, 'red', 'black'))
            arrows(x0=seq_len(nrow(coefs)), y0=coefs[,'lowerCI'], y1=coefs[,'upperCI'], angle=90, code=3, length=0.05, col=ifelse(coefs[,'pvalue'] < 0.05, 'red', 'black'))
            # text(x=seq_len(nrow(coefs)), y=coefs[,'upperCI']+yshift*ylim, labels=paste0(rep(c('p=',''),times=c(1L,nrow(coefs)-1L)), unname(formatC(coefs[,'pvalue'], format='g', digits=2, flag='#'))), 
            #      font=ifelse(coefs[,'pvalue'] < 0.05, 2, 1), cex=cex)
            # text(x=seq_len(nrow(coefs)), y=coefs[,'lowerCI']-yshift*ylim, xpd=NA, cex=cex, labels=paste0(rep(c('IRD=',''),times=c(1L,nrow(coefs)-1L)), round(coefs[,'estimate'], 1), '%'))
         }
         
         dev.off()
      },
      error = function(e) {
         message('Error between opening pdf and closing. Closing now...')
         dev.off()
      }
   )
}



























