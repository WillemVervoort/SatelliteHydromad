# --------------------------------------------------
# Utility Function to check for leap years
# Satellite Hydromad
# Willem Vervoort/Joseph Guillaume
# September 2015
# -------------------------------
leap.fun <- function(start.y,end.y) {
  # Utility function to deal with leap years
  # This function generates the "5" and "6" values
  # needed to find the dates of all the MODIS tile observations
  # find the first year
  year.b <- as.numeric(substr(start.y,1,4))
  year.e <- as.numeric(substr(end.y,1,4))
  # calculate when the next leap year will occur
  leap.count <-  (year.b/4 - floor(year.b/4))*4
  # Generate the series of "5" and "6"
  if (leap.count == 0) {
    extra.y <- (year.e-year.b)/4 - floor((year.e-year.b)/4) 
    leap1 <- rep(c(6,rep(5,3)),floor((year.e-year.b)/4))
    if (extra.y == 0) leap <- leap1 else {
      if(extra.y>1) leap <- c(leap1,6,rep(5,(extra.y-1))) else {
        leap<- c(leap1,6)
      }
    }
  } else {
    extra.y <- (year.e-year.b)/4 - floor((year.e-year.b)/4) 
    leap1 <- c(rep(5,leap.count),rep(c(6,rep(5,3)),floor((year.e-year.b)/4)))
    if (extra.y == 0) leap <- leap1 else {
      if(extra.y>1) leap <- c(leap1,6,rep(5,(extra.y-1))) else {
        leap<- c(leap1,6)
      }
    }
    
  }
  # for testing
  #return(list(leap.out=leap,leap.count.out=leap.count))
  return(leap)
}

# test
#test <- leap.fun(start.y="2000-01-01",end.y="2008")
