# --------------------------------------------------
# Calibration fo GR4J (or any hydromad model using MODISET)
# Satellite Hydromad
# Willem Vervoort/Joseph Guillaume
# September 2015
# -------------------------------

# This script is an example script that assumes that you have
# already downloaded MODIS data and developed a timeseries (zoo object possibly)
# to merge with the normal hydromad data set
# this means run: Extract.points() downloaded MODIS data and MODIS.transform()

# packages needed
library(hydromad)
library(lattice)
library(latticeExtra)
library(zoo)

Today <- format(Sys.Date(),"%Y%m%d")
basedir <- "C:/Users/rver4657/ownCloud/working/SatelliteHydromad"
setwd(basedir)

# source the ETa.merge function
source("ETa.merge.R")

#----------------------------
##   1.     READ IN DATA       ####
#-----------------------------
# the standard Cotter data file, 2000 - 2008
load("Cotter.rdata")

# load MODIS ET data 8 day cycle
#MODISET <- read.csv("20150913_AverageCorinCatchmentET.csv")
## convert to zoo
# Cot_MODISET <- zoo(MODISET[,2], order.by=as.Date(MODISET[,1]))
## generalise this data for later use
# save(Cot_MODISET,file="CotterMODISET.rdata")
load("CotterMODISET.rdata")

# ************************************************

# -----------------------------------------------
# Fitting hydromad GR4J using satellite data
# 2. Define new objective functions
# -------------------------------------------------

# use Viney's objective function(includes Bias), to fit 
# see http://hydromad.catchment.org/#hydromad.stats
hydromad.stats("viney" = function(Q, X, ...) {
  hmadstat("r.squared")(Q, X, ...) -
    5*(abs(log(1+hmadstat("rel.bias")(Q,X)))^2.5)})

# This is one way to do this, but uses the NSE
# Define an objective function that aggregates ET based on specified periods
# use buildTsObjective, but this uses the NSE
hydromad.stats("ETfun" = function(...,DATA,U) {
  .(buildTsObjective(DATA$aET, groups=DATA$et.period, 
                     FUN=sum))(DATA$aET,U$ET,...)})

# Also define an objective function that aggregates ET and fits
# # Use buildTsObjective to create a new objective function
# Define an objective function that aggregates ET based on specified periods
hydromad.stats("ETaggrViney" = function(..., DATA,U) {
  #This should be mean if the observed aET is repeated for each point in the period
  # using sum because inserted 0 values in data (ETa.merge)
  # inserted "coredata" statement after discussion with J Guillaume
  # relates to how bias is calculated
  aET.fin <- aggregate(DATA$aET,list(date=coredata(DATA$et.period)),sum)
  ET.fin <- aggregate(U$ET,list(date=coredata(DATA$et.period)),sum)
  #This can be any objective function
  obj <- hmadstat("viney")(coredata(aET.fin),coredata(ET.fin))
  return(obj)
})

# **************************************************************

# -----------------------------------------------------------
# 3.Do the standard fitting as test
# -----------------------------------------------------------
# Calibration data set
# Don't use bushfire set (so avoid 2003)
data.cal <- window(Cotter, start = "2005-01-01",
                   end = "2008-12-31")

# Define the model, important to define return_state=T
Cotter_mod <- hydromad(DATA=data.cal,
                   sma = "gr4j", routing = "gr4jrouting", 
                   x1 = c(100,1500), x2 = c(-30,20), x3 = c(5,500), 
                   x4 = c(0.5,10), etmult=c(0.01,0.5), 
                   return_state=TRUE)

# Fit without the MODIS data, traditional fit
# Using shuffled complex evolution algorithm for fitting
Cotter_fit<- fitBySCE(Cotter_mod,  
                     objective=~hmadstat("viney")(Q,X))
# Extract the coefficients and the summary
summary(Cotter_fit)
xyplot(Cotter_fit)
# ***************************************************

# ----------------------------------------------------------------
# 4. Including the ET data
# using ETa.merge()
Flow.Modis.zoo <- ETa.merge(Flowdata=Cotter,ETdata=Cot_MODISET)

# remake the calibration data
data.modis.cal <- window(Flow.Modis.zoo, start = "2005-01-01",end = "2008-12-31")

# Because we have rebuild data.cal, redefine the model
Cotter_mod_M <- hydromad(DATA=data.modis.cal,
                   sma = "gr4j", routing = "gr4jrouting", 
                   x1 = c(100,1500), x2 = c(-30,20), 
                   x3 = c(5,500), x4 = c(0.5,10), 
                   etmult=c(0.01,0.5), 
                   return_state=TRUE)

# Evaluate the model using the objective
# using equal weighting between ETa and Q
w=0.5
Cotter_Fit_B <- fitBySCE(Cotter_mod_M, 
                           objective=~w*hmadstat("viney")(X,Q) +
                             (1-w)*hmadstat("ETfun")(DATA=DATA,U=U))


summary(Cotter_Fit_B)
# Plot ET calibration
plot(data.modis.cal$aET[data.modis.cal$aET>0,], xlab="Date", 
     ylab="Actual ET (mm/day)", col="red",
     lwd=4, lty=2,ylim=c(0,max(data.modis.cal$aET)+1), 
     main = "ETfun using NSE")
plot.time <- time(window(Cot_MODISET,start="2005-01-01",end="2008-12-26"))
lines(zoo(aggregate(Cotter_Fit_B$U$ET,
                    list(date=data.modis.cal$et.period),sum),
          order.by=plot.time))
legend("topleft",c("MODIS ET", "Predicted aET"),
       lwd=c(3,1),col=c("red",1),lty=c(2,1))
# ******************************************************

# ------------------------------------------------------------
# 5. Now use the ETaggregate objective function (Viney)
# Fit the model again, using equal weighting
w = 0.5
Cotter_Fit_B_Viney <- fitBySCE(Cotter_mod_M, 
                         objective=~w*hmadstat("viney")(X,Q) +
                           (1-w)*hmadstat("ETaggrViney")(DATA=DATA,U=U))


summary(Cotter_Fit_B_Viney)
coef(Cotter_Fit_B_Viney)
hmadstat("viney")(Q=data.modis.cal$Q,X=Cotter_Fit_B_Viney$fitted.values)
hmadstat("ETaggrViney")(DATA=data.modis.cal,U=Cotter_Fit_B_Viney$U)
# the Q calibration
xyplot(Cotter_Fit_B_Viney)

# the ET calibration
plot(data.modis.cal$aET[data.modis.cal$aET>0,], 
     xlab="Date", ylab="Actual ET (mm/day)", col="red",
     lwd=4, lty=2,ylim=c(0,max(data.modis.cal$aET)+1), 
     main = "ETfun using Viney")
lines(zoo(aggregate(Cotter_Fit_B_Viney$U$ET,
                    list(date=data.modis.cal$et.period),sum),
          order.by=plot.time))
legend("topleft",c("MODIS ET", "Predicted aET"),
       lwd=c(3,1),col=c("red",1),lty=c(2,1))
# ***********************************************************
