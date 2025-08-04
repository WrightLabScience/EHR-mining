library(png)

png_list <- setNames(
   list(
      readPNG(source='~/Desktop/MSthesis/Final/MRSA_VAN_DAP.png'),
      readPNG(source='~/Desktop/MSthesis/Final/Efaec_VAN_AMP.png')
   ),
   c(
      'MRSA - VAN vs. switch to DAP',
      'E. faecalis - VAN vs. switch to AMP'
   )
)

{
   pdf(file = '~/Desktop/MSthesis/Final/flowcharts.pdf')
   par(oma=c(0,0,1.5,0), mar=c(0,0,0,0))
   plot(NA, xlim=c(1,1652+50), ylim=c(1,1560), ann=F, axes=F)
   
   xshift <- c(0, 826+50)
   for (i in seq_along(png_list)) {
      nm <- names(png_list)[i]
      rasterImage(png_list[[i]], xleft=1+xshift[i], ybottom=1, xright=826+xshift[i], ytop=1560)
   }
   text(x=c(mean(c(1,826)), mean(c(827+50, 826*2+50))), y=1650, labels=names(png_list), font=4, xpd=NA)
   
   dev.off()
}
