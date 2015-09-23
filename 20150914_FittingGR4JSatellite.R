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
# Utility functions
source("ETfit.Utilities.R")

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


# -----------------------------------------------------------
# 2.Do the standard fitting as test
# ************************************************

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
# 3. Including the ET data
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
plot.ET(caldata=data.modis.cal,Cotter_Fit_B)

# ******************************************************

# ------------------------------------------------------------
# 4. Now use the ETaggregate objective function (Viney)
# Fit the model again, using equal weighting
w = 0.5
Cotter_Fit_B_Viney <- fitBySCE(Cotter_mod_M, 
                         objective=~w*hmadstat("viney")(X,Q) +
                           (1-w)*hmadstat("ETaggrViney")(DATA=DATA,U=U))


summary(Cotter_Fit_B_Viney)
coef(Cotter_Fit_B_Viney)
# Calculate the performance measures
hmadstat("viney")(Q=data.modis.cal$Q,X=Cotter_Fit_B_Viney$fitted.values)
hmadstat("ETaggrViney")(DATA=data.modis.cal,U=Cotter_Fit_B_Viney$U)
# Show the Q calibration (standard)
xyplot(Cotter_Fit_B_Viney)

# Show the ET calibration
plot.ET(caldata=data.modis.cal,Cotter_Fit_B_Viney)

# ***********************************************************
