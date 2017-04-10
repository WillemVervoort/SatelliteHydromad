# --------------------------------------------------
# Calibration of GR4J (or any hydromad model using MODISET)
# Satellite Hydromad
# Willem Vervoort/Joseph Guillaume
# September 2015
# -------------------------------

# working directory should in "examples"
#1. Load packages, data, 2. define objective functions
source("setup.R")

# -----------------------------------------------------------
# 2.Do the standard fitting as test
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

# Because we have rebuilt data.cal, redefine the model
Cotter_mod_M <- hydromad(DATA=data.modis.cal,
                   sma = "gr4j", routing = "gr4jrouting", 
                   x1 = c(100,1500), x2 = c(-30,20), 
                   x3 = c(5,500), x4 = c(0.5,10), 
                   etmult=c(0.01,0.5), 
                   return_state=TRUE)
hydromad.options(trace=TRUE)
options(warn=1)

# Evaluate the model using the objective
# using equal weighting between ETa and Q
w=0.5
Cotter_Fit_B <- fitBySCE(Cotter_mod_M, 
                           objective=~w*hmadstat("viney") +
                             (1-w)*hmadstat("ETfun")(Q=Q,X=X,DATA=DATA,model=Cotter_mod_M))

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
                           (1-w)*hmadstat("ETaggr")(Q=Q,X=X,DATA=DATA,
                                                    model=model))

summary(Cotter_Fit_B_Viney)
coef(Cotter_Fit_B_Viney)

# Calculate the performance measures
hmadstat("viney")(Q=data.modis.cal$Q,X=Cotter_Fit_B_Viney$fitted.values)
# does not work at the moment
# hmadstat("ETaggr")(Q=data.modis.cal$Q,X=Cotter_Fit_B_Viney$fitted.values,
#                    DATA=data.modis.cal,
#                    model=Cotter_Fit_B_Viney)
# Show the Q calibration (standard)
xyplot(Cotter_Fit_B_Viney)

# Show the ET calibration
plot.ET(caldata=data.modis.cal,Cotter_Fit_B_Viney)

# **************************************************

# ----------------------------------------------------------
# 5. Show how this can be done all together and generate a runlist
test <- FitMODbySCE(mod=Cotter_mod_M, FIT_Q=T, 
                    FIT_Q_ET = T, FIT_ET = T)
summary(test, items = c("rel.bias", "r.squared","r.sq.sqrt", "r.sq.log"))

# plot
xyplot(test)

# **************************************************

