# Utility function to plot ET data
# Willem Vervoort/Joseph Guillaume
# September 2015
# ---------------------------------------------------
plot.ET <- function(caldata,ModelFit) {
  plot(caldata$aET[caldata$aET>0,], 
       xlab="Date", ylab="Actual ET (mm/day)", col="red",
       lwd=4, lty=2,ylim=c(0,max(caldata$aET)+1), 
       main = "ETfun using Viney")
  plot.time <- unique(caldata$et.period)
  totals=aggregate(coredata(ModelFit$U$ET),
                     list(date=coredata(caldata$et.period)),sum)
  lines(zoo(totals[,2],order.by=totals[,1]))
  legend("topleft",c("MODIS ET", "Predicted aET"),
         lwd=c(3,1),col=c("red",1),lty=c(2,1))
}
