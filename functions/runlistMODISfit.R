## MODIS ET fitting of models in one function
# # Calibration of GR4J (or any hydromad model using MODISET)
# Satellite Hydromad
# Willem Vervoort/Joseph Guillaume
# October 2015
# -------------------------------
#1. Load packages, data, 2. define objective functions
# source("setup.R")

# Define a function that fits all three (FLOW_MODIS fit by SCE)
FitMODbySCE <- function(mod, w = 0.5, FIT_Q = TRUE,
                        FIT_Q_ET = TRUE,
                        FIT_ET = TRUE, Objfun=hmadstat("viney")) {
  #
  # mod is a hydromad model specified for fitting
  # this model might include ET data as well as flow data
  # w:  weighting between ET and Q if both are fitted in the calibration
  # FIT_Q is whether to fit just on flow data (TRUE or FALSE)
  # FIT_Q_ET is whether to fit both ET and flow data
  # FIT_ET is whether to fit just on ET data
  # Objfun specifies the objective function to use in both the Q fit, 
  # this is either JointQandET and ETAggr 
  
  # 1. Fit Q only
  if (FIT_Q==T) {
    base_fit<- fitBySCE(mod,  
                          objective=Objfun)
  }
  # 2. Fit Q and ET using JointQandET
  if (FIT_Q_ET==T) {
    #browser()
    # use master function
    ETQET_fit <- fitBySCE(mod, 
                           objective=~hmadstat("JointQandET")(Q,X,w=w,
                                                DATA=DATA, model=model,
                                                 objf = Objfun))
    
  }
  ## also include just fitting on ET and no fit on Q
  # 3. Fit ET alone using ETAggr and obj fun
  if (FIT_ET==T) {
    ETAggr_fit <- fitBySCE(mod,
                          objective=~hmadstat("ETaggr")(Q=Q,X=X,
                                                        DATA=DATA,
                                                        model=model,
                                                objf = Objfun))
  }
  out <- runlist(
    "Fit only Q" = if(FIT_Q==T) base_fit, 
    "Fit only ET" = if(FIT_ET==T) ETAggr_fit,
    "Fit ET and Q" = if(FIT_Q_ET==T) ETQET_fit)
  class(out) <- c("Mrunlist",class(out))
  out
}

# The summary function, adapted from the original runlist in hydromad
summary.Mrunlist <-
  function(object, ..., FUN = summary, items = NULL)
  {
    stopifnot(is.list(object))
    if (length(object) == 0)
      return(NULL)
    ## extract elements from summary which are single numbers

    cc <- lapply(object, function(x, ...) {
      tmp <- FUN(x, ...) 
      if (is.null(items)) {
        tmp <- tmp[unlist(lapply(tmp, function(z) {
          is.numeric(z) && !is.matrix(z) &&
            (length(z) == 1)
        }))]
      } else {
        tmp <- tmp[items]
      }
      unlist(tmp)
    }, ...)
    ## pad out missing entries with NAs
    ## find the set of all names
    allnms <- na.omit(unique(unlist(lapply(cc, names))))
    ans <- matrix(NA_real_,
                  nrow = length(object),
                  ncol = length(allnms),
                  dimnames = list(names(object), allnms))
    for (i in 1:NROW(ans))
      if(length(na.omit(names(cc[[i]])))!=0) {
        ans[i, names(cc[[i]])] <- cc[[i]]
      }
    ans <- as.data.frame(ans)
    class(ans) <- c("summary.runlist", class(ans))
    ans
  }



