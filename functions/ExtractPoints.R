# ------------------------------
# Satellite Hydromad project
# Willem Vervoort/Joseph Guillaume
# Sept 2015
# Utility function to 
# Extract the latitude and longitude of points in a shape file (catchment)
# This requires:
# a shape file of the area
# a test image from MODIS or a different satellite

#Install packages if missing
for(pkg in c("raster","maptools","rgeos","rgdal"))
  if(!require(pkg,character.only=TRUE)) install.packages(pkg)
library(raster)
library(maptools)
library(rgeos)
library(rgdal)

# You don't need this function if you have the lat and long of your locations
Extract.points <- function(wdir=getwd(), catchment_file, modis_eg_file,show.plot=TRUE) {
  # wdir is the working directory, required
  # catchment_file is the base name for the catchment shapefile
  # modis_eg_file is the name of a test MODIS satellite image of the area
  
  # Reading the shape file of the catchment
  shape <- readOGR(wdir,catchment_file)
  
  # reading in the test MODIS file
  modis <- raster(paste(wdir,modis_eg_file,sep="/")) 
  if(show.plot){
    image(modis)
    lines(shape)
  }
  
  # Now clipping to the catchment boundary
  modis.clip <- crop(modis, extent(shape)) 
  modis.clip <- mask(modis.clip, shape)

  # Now extracting the points
  p <- rasterToPoints(modis.clip) # it contains the cell values as well
    
  if(show.plot){
    image(modis.clip)
    lines(shape)
    points(p)
  }
  
  # only taking the x(lon) and y(lat) values
  return(data.frame(lon = p[,1], lat = p[,2]))
}

# # Example
# # Use data in folder Shapefiles_testdata
#  Points <- Extract.points(wdir="Shapefiles_Testdata",
#                 catchment_file="corin_pro",modis_eg_file="test1")
# head(Points)
# # write away for later use
# # write.csv(Points,"Corinpoints.csv")