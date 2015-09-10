# Demonstrate downloading satellite data using MODIS tools
# Example catchment  Corin

require(MODISTools)
# you migth need to set your proxy
Sys.setenv(http_proxy="web-cache.usyd.edu.au:8080")

basedir <- "C:/Users/rver4657/ownCloud/working/SatelliteHydromad"
setwd(basedir)
# read in file with xy locations (or create a data.frame)
xy.loc <- read.csv("CorinPoints.csv")
# Following the MODIS tools manual
coords <- data.frame(lat=xy.loc$lat, 
                     long=xy.loc$lon,
                     start.date=rep(2000,nrow(xy.loc)), 
                     end.date=rep(2008,nrow(xy.loc)), 
                     transect=1:nrow(xy.loc))

# We need to figure out the name of the product, you can use GetProducts()
GetProducts()
# and check out the data bands that are in the product
GetBands(Product="MOD16A2")

# # Now use MODISSubsets (This can take very long)
 MODISSubsets(LoadDat=coords, Product = "MOD16A2",
                Bands=c("ET_1km","ET_QC_1km"), StartDate=T, Size=c(0,0),
                SaveDir=paste(getwd(),"MODIS",sep="/"))
