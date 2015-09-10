# set working dir
setwd("u:/my-workspace/rver4657/lwsc3007")

require(hydromad)

load("20150319_HydromadCalibration.RData")

# skip over all the MODISTOOLS stuff and packages
# Only read in stations
Stations <- read.csv("Station_locations.csv")

# start at Step 3 in the document

# The old flow file
head(Karuah)
# read in satellite data
# list of files
x1 <- list.files (pattern= (".asc"))
# find number of rows of one file
n <- nrow(read.csv(x1[1],header=F))

# create storage
Store <- data.frame(Year = numeric(length=n),
                    Jdate = numeric(length=n),
                    ET = numeric(length=n))
Store1 <- list()

# read in data into Store dataframe
# put for each location into Store1 list
for (i in 1:length(x1)) {
  Mdata <- read.csv(x1[i],header=F)
  # substring the Year from the eigth column
  Store[,1] <- as.numeric(substr(Mdata[,8],2,5))
  # substring the Julian Date from the eigth column
  Store[,2] <- as.numeric(substr(Mdata[,8],6,8))
  # get the ET values from column 11
  Store[,3] <- Mdata[,11]/10 ##  0.1 scaling factor (see MODIS read_me)
  # Put in each list element
  Store1[[i]] <- Store
}
# put names of stattions
names(Store1) <- Stations$Station

# However, we are just going to use Crabapple
# just use crabapple, so select this out of the list
ET.crab <- Store1["Crabapple"][[1]]

# this calculates the difference between the julian dates, 
# but this means the last day gets -360, for non-leap years should be 5
days.cor <- c(diff(ET.crab$Jdate),5)
# Now we need to replace the -360 values, 2000 is a leap year
# there are 14 years alltogether
leap <- c(rep(c(6,rep(5,3)),3),6,5)
days.cor[days.cor==-360] <- leap
# check
days.cor

## get dates for each 8-day tile
for (i in 1:(length(days.cor)-1)) {
  if (i ==1) {
    dates.tiles <- seq.Date(as.Date("2000-01-01"),
                            by=days.cor[i],length.out=2)
  } else {
    dates.tiles <- c(dates.tiles,
      (seq.Date(as.Date(dates.tiles[length(dates.tiles)]),
                       by=days.cor[i],length.out=2))[2])
  }
}
# you could make a plot to check
plot(as.Date(dates.tiles),ET.crab$ET,type="l")

# Step 4 fitting hydromad
data.cal <- window(Karuah, start = "2000-01-01",
                   end = "2005-12-31")

# Define the model, important to define return_state=T
# parameter guesses based on week 3
Crab.Q <- hydromad(DATA=data.cal,
                   sma = "gr4j", routing = "gr4jrouting", 
                   x1 = c(100,1500), x2 = c(-30,20), x3 = c(5,500), x4 = c(0.5,10), 
                   etmult=c(0.01,0.5), 
                   return_state=TRUE)
# use Viney's objective function(includes Bias), to fit 
# see http://hydromad.catchment.org/#hydromad.stats
hydromad.stats("viney" = function(Q, X, ...) {
  hmadstat("r.squared")(Q, X, ...) -
    5*(abs(log(1+hmadstat("rel.bias")(Q,X)))^2.5)})

# Using shuffled complex evolution algorithm for fitting
Crabfit.Q<- fitBySCE(Crab.Q,  
                     objective=~hmadstat("viney")(Q,X))
# you can extract the coefficients and the summary
# do this yourself
summary(Crabfit.Q)
xyplot(Crabfit.Q)

### STEP 5
# Make Modis ET data set into a zoo series
aET <- zoo(ET.crab$ET, order.by=dates.tiles)

# create a new vector of dates for the tiles.
# However we want to expand to make sure each date is repeated 8 times to later aggregate the predicted ET and make it into a zoo data
# we use dates.tiles and days.cor (the number of days between tiles)
date.sum <- zoo(rep(dates.tiles, days.cor),
                order.by= seq.Date(as.Date(dates.tiles[1]),
                                   by=1,length=sum(days.cor)))

# now merge ETa with Karuah
Flow.zoo <- merge(data.cal ,aET=aET, et.period=date.sum,all=T)

# remake the calibration data
data.cal <- window(Flow.zoo, start = "2000-01-01",end = "2005-12-31")
# insert 0 values for missing aET data
data.cal[is.na(data.cal$aET),"aET"] <- 0

# This is one way to do this, but uses the NSE
# Define an objective function that aggregates ET based on specified periods
# use buildTsObjective, but this uses the NSE
# This is ADVANCED R!
hydromad.stats("ETfun" = function(...,DATA,U) {
  .(buildTsObjective(DATA$aET, groups=DATA$et.period, 
                     FUN=sum))(DATA$aET,U$ET,...)})

# redefine the model
Crab.Q <- hydromad(DATA=data.cal,
                   sma = "gr4j", routing = "gr4jrouting", 
                   x1 = c(100,1500), x2 = c(-30,20), x3 = c(5,500), x4 = c(0.5,10), 
                   etmult=c(0.01,0.5), 
                   return_state=TRUE)

# Evaluate the model using the objective
w=0.5
Crabfit.both_a <- fitBySCE(Crab.Q, 
                 objective=~w*hmadstat("viney")(X,Q) +
                 (1-w)*hmadstat("ETfun")(DATA=DATA,U=U))


summary(Crabfit.both_a)
# Plot ET calibration
plot(data.cal$aET[data.cal$aET>0,], xlab="Date", 
     ylab="Actual ET (mm/day)", col="red",
     lwd=4, lty=2,ylim=c(0,max(data.cal$aET)+1), 
     main = "ETfun using NSE")
lines(zoo(aggregate(Crabfit.both_a$U$ET,
                    list(date=data.cal$et.period),sum),
          order.by=dates.tiles))
legend("topleft",c("MODIS ET", "Predicted aET"),
       lwd=c(3,1),col=c("red",1),lty=c(2,1))
