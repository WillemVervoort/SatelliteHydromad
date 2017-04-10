# packages needed
if(!require(hydromad)) 
  install.packages("hydromad", repos="http://hydromad.catchment.org",
                   dependencies = "Depends")
library(hydromad)
if(!require(lattice)) install.packages("lattice")
library(lattice)

# Load functions
source("../functions/leapfun.R")
source("../functions/ETfit.objectives.R")
source("../functions/ETa.merge.R")
#source("functions/ExtractPoints.R")
#source("functions/MODIS_to_Timeseries.R")
source("../functions/plot.ET.R")
source("../functions/runlistMODISfit.R")

#----------------------------
##   1.     READ IN DATA       ####
#-----------------------------
# the standard Cotter data file, 2000 - 2008
load("Cotter.rdata")

# If not already done, download MODIS data and 
# develop a timeseries (zoo object possibly)
# to merge with the normal hydromad data set

if(!file.exists("Corinpoints.csv")){
  # Example
  # Use data in folder Shapefiles_testdata
  Points <- Extract.points(wdir="Shapefiles_Testdata",
                           catchment_file="corin_pro",modis_eg_file="test1")
  head(Points)
  # write away for later use
  write.csv(Points,"Corinpoints.csv")
}

if(!dir.exists("MODIS")){
  source("DownloadMODISdata.R")
}

if(!file.exists("AverageCorinCatchmentET.csv")){

  # read in file with xy locations (or create a data.frame)
  xy.loc <- read.csv("CorinPoints.csv")
  
  # run function with defaults
  ET.out <- MODIS.transform()
  
  # write file away
  write.csv(ET.out,"AverageCorinCatchmentET.csv",row.names=F)
}


# load MODIS ET data 8 day cycle
if(!file.exists("CotterMODISET.rdata")){
  MODISET <- read.csv("AverageCorinCatchmentET.csv")
  # convert to zoo
  Cot_MODISET <- zoo(MODISET[,2], order.by=as.Date(MODISET[,1]))
  # generalise this data for later use
  save(Cot_MODISET,file="CotterMODISET.rdata")
}
load("CotterMODISET.rdata")
