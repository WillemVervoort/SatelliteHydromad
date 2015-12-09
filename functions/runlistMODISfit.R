## MODIS ET fitting of models in one function
# # Calibration of GR4J (or any hydromad model using MODISET)
# Satellite Hydromad
# Willem Vervoort/Joseph Guillaume
# October 2015
# -------------------------------
#1. Load packages, data, 2. define objective functions
source("setup.R")

# Calibration data set
# Don't use bushfire set (so avoid 2003)
# Merge the ET MODIS data with the original Cotter data
# Maybe make this into a Cotter.M_ET dataset?
Flow.Modis.zoo <- ETa.merge(Flowdata=Cotter,ETdata=Cot_MODISET)

# the calibration data
data.modis.cal <- window(Flow.Modis.zoo, start = "2005-01-01",end = "2008-12-31")

# Create mother of a function that fits all three
# first define the hydromad object
# using GR4J but could be any Hydromad model
Cotter_mod_M <- hydromad(DATA=data.modis.cal,
                         sma = "gr4j", routing = "gr4jrouting", 
                         x1 = c(100,1500), x2 = c(-30,20), 
                         x3 = c(5,500), x4 = c(0.5,10), 
                         etmult=c(0.01,0.5), 
                         return_state=TRUE)

# Define a function that fits all three (FLOW_MODIS fit by SCE)
FM_fitBySCE <- function(mod, w=0.5, FIT_base=TRUE,
                        FIT_QET=TRUE,
                        FIT_Aggr=TRUE, Objfun=hmadstat("viney")) {
  #
  # mod is a hydromad model specified for fitting
  # this model might include ET data as well as flow data
  # w is the weighting between ET and flow data in the calibration
  # FIT_base is whether to fit just on flow data (TRUE or FALSE)
  # FIT_ETFun is whether to fit the ETFun function based on NSE
  # FIT_Aggr is whether to fit the ETAggr function with any obj fun
  # Objfun specifies the objective function to use in both the Q fit, 
  # JointQandET and the ETAggr 
  
  # 1. Fit Q only
  if (FIT_base==T) {
    base_fit<- fitBySCE(mod,  
                          objective=Objfun)
  }
  # 2. Fit Q and ET using JointQandET
  if (FIT_QET==T) {
    browser()
    # use master function
    ETQET_fit <- fitBySCE(mod, 
                           objective=~hmadstat("JointQandET")(Q,X,w,
                                                DATA=DATA,U=U,
                                                 objf = Objfun))
    
  }
  ## also include just fitting on ET and no fit on Q
  # 3. Fit ET alone using ETAggr and obj fun
  if (FIT_Aggr==T) {
    ETAggr_fit <- fitBySCE(mod,
                          objective=hmadstat("ETaggrViney")(DATA=DATA,U=U,
                                                objf = Objfun))
  }
  out <- runlist(
    "Base" = if(!exists("base_fit")) base_fit, 
    "ET Aggr fit" = if(!exists("ETAggr_fit")) ETAggr_fit,
    "ET QET fit" = if(!exists("ETQET_fit")) ETQET_fit)
return(out)
  
}

test <- FM_fitBySCE(mod=Cotter_mod_M, FIT_base=F, 
                    FIT_Aggr = T, FIT_QET = F)

