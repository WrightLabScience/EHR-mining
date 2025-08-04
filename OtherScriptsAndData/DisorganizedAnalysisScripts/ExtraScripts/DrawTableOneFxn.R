drawTableOne <- function(p, trt, mar_left=13, xpos=max(abs(p$diffs))*1.1, col='black', cex=1, add=FALSE, addlegend=FALSE) {
   if (!add) {
      par(mar=c(3, mar_left, 1, 1.5), mgp=c(1.5, 0.25, 0), tck=-0.01)
      plot(NA, xlim=xpos*c(-1,1), cex.axis=cex, ylim=c(0, nrow(p)+1), yaxs='i', yaxt='n', ylab='', xlab='Standardized mean difference', cex.lab=cex)
      segments(x0=0, y0=0, y1=nrow(p)+3, lty=2, xpd=NA)
      abline(h=seq_len(nrow(p)), lty=3, lwd=0.5)
      abline(v=c(-0.1, 0.1), lty=2, lwd=0.75)
      if (addlegend)
         legend('bottomright', inset=c(-0.52,0), xpd=NA, cex=cex-0.2, legend=c('not sig.', 'p < 0.1', 'p < 0.05'), pch=c(1, 16, 18), col=col)
      rownames(p)[grep('^LAB_', rownames(p))] <- stringr::str_to_title(gsub('LAB', '', gsub('_', ' ', rownames(p)[grep('^LAB_', rownames(p))])))
      # rownames(p) <- stringr::str_to_title(gsub('_', ' ', rownames(p)))
      # rownames(p) <- gsub('^(Icd|Lab)', '', rownames(p))
      # make_cap_vars <- which(rownames(p) %in% c('Rbc', 'Mcv', 'Esbl', 'Mch', 'Mchc', 'Wbc') | grepl('Emp [A-z]+', rownames(p)))
      # rownames(p)[make_cap_vars] <- toupper(rownames(p)[make_cap_vars])
      # rownames(p)[rownames(p) == 'Aki 1w'] <- 'AKI 1w'
      text(x=-1.1*xpos, y=seq_len(nrow(p)), cex=cex, labels=rownames(p), xpd=NA, adj=1)
      #arrows(x0=-0.15*xpos, x1=0.15*xpos, y0=nrow(p)+1.5, length=0.1, lwd=2, code=3, xpd=NA)
   }
   pch_vec <- ifelse(p$pvals < 0.1, 16, 1)
   pch_vec <- ifelse(p$pvals < 0.05, 18, pch_vec)
   points(x=p$diffs, y=seq_len(nrow(p)), cex=cex, pch=pch_vec, col=col)
   text(x=0.7*xpos*c(-1, 1), y=nrow(p)+1.5, cex=cex, labels=paste0(' > in ', trt, ' group'), xpd=NA)
}
