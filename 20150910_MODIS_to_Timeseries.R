
# Transform MODIS Satellite data to 8 day time series for modelling
# this should be a function
# Dir where MODIS data is stored (extracted with MODISTOOLS), so asc files
# This function can also work with tif files

# Function to check for leap years
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

# MODIS transform function
MODIS.transform <- function(MODISdir="MODIS",patt=".asc",
                            start="2000-01-01", end="2008-12-31",
                            xy.loc=NULL, merge.FUN=mean) 
{
  # this function transforms the MODIS data downloaded with MODIS tools
  # to a proper zoo timeseries rather than Julian dates
  # it also does catchment averaging
  # MODISdir is the directory where the MODIS downloads can be found
  # patt is the pattern to look for (Either tif or "asc")
  # start and end are the start and end dates
  # xy.loc are the xy locations of the downloaded sites 
  # (only needed if patt = ".tif")
  # This is two column data.frame with lon and lat
  # merge.FUN is the function to use to merge the data
  # if merge.FUN is FALSE, than all the series are returned
  #packages
  require(maptools,quietly=T)
  require(raster,quietly=T)
  require(rgdal,quietly=T)


  # This first looks at the list of files
  x1 <- list.files(path=MODISdir,
                 pattern=patt)

  #patt can also be "*tif" will get to this later
  if (patt == ".asc") {
    # for the asc files, each files stores all the values in time for 1 location
    n <- nrow(read.csv(paste(MODISdir,x1[1],sep="/"),header=F))
    # Create storage for the data
    Store <- data.frame(Year = numeric(length=n),
                        Jdate = numeric(length=n),
                        ET = numeric(length=n))
    # list to store the different locations
    Store1 <- list()
    
    for (i in 1:length(x1)) {
      Mdata <- read.csv(paste(MODISdir,x1[i],sep="/"),header=F)
      Store[,1] <- as.numeric(substr(Mdata[1:n,8],2,5))
      Store[,2] <- as.numeric(substr(Mdata[1:n,8],6,8))
      Store[,3] <- Mdata[1:n,11]/10 ##  0.1 scaling factor (see MODIS read_me)
      Store1[[i]] <- Store
    }
  } # end if patt = tif
  # in case the files are tif files and need to be first rastered
  if (patt==".tif") {
    # for the tif files, each individual file has all locations for one time
    # define locations as spatial points
    loc <- SpatialPoints(cbind(xy.loc$lon, xy.loc$lat), proj4string=CRS("+proj=longlat"))
    EandN <- spTransform(loc, CRS("+init=epsg:28355"))
    
    n <- nrow(xy.loc)
    # run through the files and extract values
    Store <- data.frame(Year = numeric(length=n),
                        Jdate = numeric(length=n),
                        ET = numeric(length=n))
    # list to store the different times
    Store1 <- list()

    for(i in 1:length(x1)){
      #i=1
#      test <- readGDAL(paste(MODISdir,x1[i],sep="/"))
      Mdata <- raster(paste(MODISdir,x1[i],sep="/"))
      
      Store[,1] <- as.numeric(substr(strsplit(x1[i],"[.]")[[1]][2],2,5))
      Store[,2] <- as.numeric(substr(strsplit(x1[i],"[.]")[[1]][2],6,8))
      Store[,3] <- extract(Mdata,EandN)
      Store1[[1]] <- Store
    }
  } # end if patt = tif
    

  # rbind to one data.frame  
  ETunlist <- do.call(rbind,Store1)
  # # you can check distribution
  # hist(test$ET, main="Histogram of raw ET values across all pixels", 
  #      xlab="8 daily sum ET (mm)")

  # Now aggregate across the catchment
  ET.mean <- aggregate(ETunlist[,3],list(Jdate=ETunlist$Jdate, 
                        Year=ETunlist$Year), merge.FUN,na.rm=T)
  # Julian Days
  Jdays <- ET.mean$Jdate
  # work out number of days between observations (not always 8!)
  days.cor <- c(diff(as.numeric(Jdays)),5)
  # Now we need to replace the -360 values at the end of the year
  # run leap.fun to account for leap years
  leap <- leap.fun(start,end)
  # correct the -360 values
  days.cor[days.cor==-360] <- leap
  # check
  #days.cor
    
  ## get dates for each 8-day tile
  for (i in 1:(length(days.cor)-1)) {
    if (i ==1) {
      dates.tiles <- seq.Date(as.Date(start),
                              by=days.cor[i],length.out=2)
    } else {
      dates.tiles <- c(dates.tiles,
                       (seq.Date(as.Date(dates.tiles[length(dates.tiles)]),
                                   by=days.cor[i],length.out=2))[2])
    }
  }
  # you could make a plot to check
  #plot(as.Date(dates.tiles),ET.mean[,3],type="l")
  return(data.frame(Dates=dates.tiles,ET=ET.mean$x))
}


# Example
Today <- format(Sys.Date(),"%Y%m%d")
basedir <- "C:/Users/rver4657/ownCloud/working/SatelliteHydromad"
setwd(basedir)
# read in file with xy locations (or create a data.frame)
xy.loc <- read.csv("CorinPoints.csv")

# run function with defaults
ET.out <- MODIS.transform()

# write file away
write.csv(ET.out,paste(Today,"AverageCorinCatchmentET.csv",sep="_"),row.names=F)
  
