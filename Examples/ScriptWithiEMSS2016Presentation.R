# Presentation iEMSS 2016
# some examples from Hydromad
library(hydromad)
# plot the basic data
data(Cotter)
xyplot(Cotter)
# define a sub set
Cotter.cal = window(Cotter, start = "1981-01-01",end = "1990-12-31")

# define a model
CMod <- hydromad(Cotter.cal, sma="gr4j", routing="gr4jrouting",
                 etmult=0.15,x1 = c(100,1500), x2 = c(-30,20), 
                 x3 =c(5,500), x4 = c(0.5,10))
print(CMod)

# fot the model
CotterFit <- fitByOptim(CMod,objective=~hmadstat("r.squared")(Q,X),
                          samples=500,method="PORT")
summary(CotterFit)
# make a plot
xyplot(CotterFit, with.P=TRUE, 
       xlim=as.Date(c("1981-01-01", "1990-01-01")))

## Time permitting
# # You can update the model with a different dataset or with different pars
# Cotter.val = window(Cotter, start = "1991-01-01",end = "1996-12-31")
# # update the model fit
# sim.val = update(CotterFit, newdata = Cotter.val)
# xyplot(sim.val, with.P=TRUE, 
#        xlim=as.Date(c("1991-01-01", "1993-01-01")))
# allMods = runlist("Calibration"=CotterFit, "Validation"= sim.val)
# round(summary(allMods),2)

# different objective function (logged)
CotterFit_log = fitByOptim(CMod,objective=hmadstat("r.sq.log"),
                            samples=500,method="PORT")
# use a runList
allMods = runlist(calib.rsq=CotterFit, calib.log=CotterFit_log)
# create summary
round(summary(allMods),2)

# Satellite ET extension
#1. Load packages, data, 2. define objective functions
#source("setup.R")
source("functions/leapfun.R")
source("functions/ETfit.objectives.R")
source("functions/ETa.merge.R")
#source("functions/ExtractPoints.R")
#source("functions/MODIS_to_Timeseries.R")
source("functions/plot.ET.R")
source("functions/runlistMODISfit.R")
load("Examples/CotterMODISET.rdata")
load("Examples/Cotter.rdata")
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

# Because we have rebuilt data.cal, redefine the model
Cotter_mod_M <- hydromad(DATA=data.modis.cal,
                         sma = "gr4j", routing = "gr4jrouting", 
                         x1 = c(100,1500), x2 = c(-30,20), 
                         x3 = c(5,500), x4 = c(0.5,10), 
                         etmult=c(0.01,0.5), 
                         return_state=TRUE)
hydromad.options(trace=TRUE)
options(warn=1)


# ******************************************************
# fit only on ET data
Cotter_Fit_ET <- fitBySCE(Cotter_mod_M,
                  objective=~hmadstat("ETfun")(Q=Q,X=X,DATA=DATA,model=Cotter_mod_M))
summary(Cotter_Fit_ET)
coef(Cotter_Fit_ET)
# Show the ET calibration
plot.ET(caldata=data.modis.cal,Cotter_Fit_ET, main="ETfun using NSE")


# ------------------------------------------------------------
# 5. Now use the ETaggregate objective function (Viney)
# Fit the model again, using equal weighting
w = 0.5
Cotter_Fit_B_Viney <- fitBySCE(Cotter_mod_M, 
                               objective=~w*hmadstat("viney")(Q,X) +
                                 (1-w)*hmadstat("ETaggr")(Q=Q,X=X,DATA=DATA,
                                                               model=model))


summary(Cotter_Fit_B_Viney)
coef(Cotter_Fit_B_Viney)

# Calculate the performance measures
hmadstat("viney")(Q=data.modis.cal$Q,X=Cotter_Fit_B_Viney$fitted.values)
# this does not work yet
# hmadstat("ETaggr")(Q=data.modis.cal$Q,X=Cotter_Fit_B_Viney$fitted.values,
#                    DATA=data.modis.cal,
#                    model=Cotter_Fit_B_Viney)
# Show the Q calibration (standard)
xyplot(Cotter_Fit_B_Viney)

# Show the ET calibration
plot.ET(caldata=data.modis.cal,Cotter_Fit_B_Viney, main="ETfun using Viney")

# **************************************************

# ----------------------------------------------------------
# 6. Show how this can be done all together and generate a runlist
test <- FitMODbySCE(mod=Cotter_mod_M, FIT_Q=T, 
                    FIT_Q_ET = T, FIT_ET = T)
round(summary(test, 
              items = c("rel.bias", "r.squared",
                        "r.sq.sqrt", "r.sq.log")),2)

# plot
xyplot(test)






