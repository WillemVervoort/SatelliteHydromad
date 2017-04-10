# Utility function to plot ET data
# Willem Vervoort/Joseph Guillaume
# September 2015
# ---------------------------------------------------
plot.ET <- function(caldata,ModelFit, main="ETa predicted and observed",
                    xlab="Date", ylab="Actual ET (mm/day)", col="red",
                    lwd=4, lty=2,ylim=c(0,max(caldata$aET)+1)) {
  plot(caldata$aET[caldata$aET>0,], 
       xlab=xlab, ylab=ylab, col=col,
       lwd=lwd, lty=lty,ylim=ylim, 
       main = main)
  plot.time <- unique(caldata$et.period)
  totals=aggregate(coredata(ModelFit$U$ET),
                     list(date=coredata(caldata$et.period)),sum)
  lines(zoo(totals[,2],order.by=totals[,1]))
  legend("topleft",c("MODIS ET", "Predicted aET"),
         lwd=c(lwd,1),col=c(col,1),lty=c(lty,1))
}
