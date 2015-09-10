
# Transform MODIS Satellite data to 8 day time series for modelling
# this should be a function
# Dir where MODIS data is stored (extracted with MODISTOOLS), so asc files
# This function might have to be adjusted if you are working with tif files

getwd()
xy.loc <- read.csv("CorinPoints.csv")

# Function to check for leap years
leap.fun <- function(start.y,end.y) {
  # find the first year
  year.b <- as.numeric(substr(start.y,1,4))
  year.e <- as.numeric(substr(end.y,1,4))
  # calculate when the next leap year will occur
  leap.count <-  (year.b/4 - floor(year.b/4))*4
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
  return(list(leap.out=leap,leap.count.out=leap.count))
         #return(leap)
}

# test
#test <- leap.fun(start.y="2000-01-01",end.y="2008")

MODIS.transform <- function(MODISdir="MODIS",patt=".asc",
                            start="2000-01-01", end="2008-12-31"
                            xy.loc)
  #patt can also be "*tif" will get to this later
  
  if (patt == ".asc") {
  # This first looks at the list of files
  x1 <- list.files(path=MODISdir,
                   pattern=patt)

  n <- nrow(read.csv(paste(MODISdir,x1[1],sep="/"),header=F))

  Store <- data.frame(Year = numeric(length=n),
                    Jdate = numeric(length=n),
                    ET = numeric(length=n))
  Store1 <- list()

  for (i in 1:length(x1)) {
    Mdata <- read.csv(paste(MODISdir,x1[i],sep="/"),header=F)
    Store[1:n,1] <- as.numeric(substr(Mdata[1:n,8],2,5))
    Store[1:n,2] <- as.numeric(substr(Mdata[1:n,8],6,8))
    Store[1:n,3] <- Mdata[1:n,11]/10 ##  0.1 scaling factor (see MODIS read_me)
    Store1[[i]] <- Store
  }


  ETunlist <- do.call(rbind,Store1)
  # # check distribution
  # hist(test$ET, main="Histogram of raw ET values across all pixels", 
  #      xlab="8 daily sum ET (mm)")

  ET.mean <- aggregate(ETunlist[,3],list(Jdate=ETunlist$Jdate, Year=ETunlist$Year), 
                     mean,na.rm=T)
  
  Jdays <- ET.mean$Jdate
  
  days.cor <- c(diff(as.numeric(Jdays)),5)
  # Now we need to replace the -360 values
  
  # run leap.fun
  leap <- leap.fun(start,end)
  # correct the -360 values
  days.cor[days.cor==-360] <- leap
  # check
  #days.cor
  
  # GOT TO HERE ON 20150910
  
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
  plot(as.Date(dates.tiles),Store1[[1]][,3],type="l")
  
  data.m <- do.call(cbind,Store1) 
  data.m2 <- data.m[,-grep("Jdate",names(data.m))]
  data.m3 <- data.m2[,-grep("Year",names(data.m2))]
  data.ET <- data.m3
  data.ET$MFRSGL <- data.ET$MFRAA.ET
  data.ET$Dates <- dates.tiles
  data.ET$Year <- data.m[,1]
  
  
  days.cor <- c(rep(8,45),6,rep(c(rep(8,45),5),3),rep(8,45),6,
              rep(c(rep(8,45),5),3),rep(8,45),6)
  # length(days.cor)
  # # should equal
  # nrow(z.mean)
  ET.mean.cor <- ET.mean$x/days.cor
  #z.mean.cor is ET/days.cor so it is the ET for one day, rather than the sum of 8 days - however, it is still at 8-day intervals

  plot(ET.mean$x, type="p", ylab="ET (mm)",ylim=c(0,16)) 
  # this is the 8 day sum of ET
points(ET.mean.cor,lty=2,col="red") #this is the 8 day sum divided by 8 (or 5/6)
legend("topleft",c("8 day sum of ET", "daily values"),
       col=c(1,"red"), pch=1)

# reset z.mean for rest of script
ET.mean.raw <- ET.mean
ET.mean <- ET.mean.cor


## get dates for each 8-day tile
for (i in 1:length(days.cor)) {
  if (i == 1) {
    dates.tiles <- seq.Date(as.Date("2000-01-01"),
                            by=days.cor[i]-1,length.out=2)
  } else {
    dates.tiles <- c(dates.tiles,
                     (seq.Date(as.Date(dates.tiles[length(dates.tiles)]),
                               by=days.cor[i],length.out=2))[2])
  }
}

modis.days <- ET.mean.raw ## Julian days the modis image was taken
dates <- as.character(seq.Date(as.Date("2000-01-01"),
                               to=as.Date("2008-12-31"),by=1))
Jdates <- c(1:length(dates[grep("2000",dates)]),
            1:length(dates[grep("2001",dates)]),
            1:length(dates[grep("2002",dates)]), 
            1:length(dates[grep("2003",dates)]),
            1:length(dates[grep("2004",dates)]),
            1:length(dates[grep("2005",dates)]),
            1:length(dates[grep("2006",dates)]),
            1:length(dates[grep("2007",dates)]),
            1:length(dates[grep("2008",dates)]))
full.days <- 1:length(Jdates)
# create a dataframe to combine
JDates2 <- data.frame(count=full.days,Julian=Jdates,dates=dates,
                      year=as.numeric(substr(dates,1,4)))

# find matching Modis days
modis.days2 <- JDates2[(JDates2$Julian %in% modis.days$Jdate & 
                          JDates2$year %in% modis.days$Year),] 

## Loess smooth
smooth.interpET <- loess.smooth(modis.days2$count, ET.mean, span = 0.04, 
                                degree = 2, evaluation = length(Jdates))$y
# smooth.interp is the daily, lumped, aET value, with length of 4383.



## Plot smooth  and linear spline
plot(full.days,smooth.interpET, type = "l", xaxt="n", lwd=2,
     ylim = c(0,3), ylab = "MODIS ET (mm)", xlab = "Time")       
m.days <- as.Date(paste(modis.days$Jdate,modis.days$Year,sep="-"),"%j-%Y")
points(modis.days2$dates,ET.mean,col="red")
ticks <- cumsum(rep(365, times = 10))-365
axis(1, at = ticks, labels = 2000:2009) 
legend("topleft",c("smooth MODIS ET", "MODIS 8 day values"),lwd=c(1,NA),
       col=c(1,"red"), pch=c(NA,1))


# THIS BELOW SHOULD DO 8 DAY EXTRACTION


Jdays <- Store1[[1]]$Jdate

days.cor <- c(diff(as.numeric(Jdays)),5)
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
plot(as.Date(dates.tiles),Store1[[1]][,3],type="l")

data.m <- do.call(cbind,Store1) 
data.m2 <- data.m[,-grep("Jdate",names(data.m))]
data.m3 <- data.m2[,-grep("Year",names(data.m2))]
data.ET <- data.m3
data.ET$MFRSGL <- data.ET$MFRAA.ET
data.ET$Dates <- dates.tiles
data.ET$Year <- data.m[,1]

## write to file
write.csv(data.ET, file=paste("x:/vervoort/research/cotter/MODISET/",today,"_modisETa.csv",sep=""),
          row.names=F)
