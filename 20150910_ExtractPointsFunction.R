# Extract the latitude and longitude of points in a shape file (catchment)
# This requires:
# a shape file of the area
# a test image from MODIS or a different satellite


# You don't need this function if you have the lat and long of your locations
Extract.points <- function(wdir=getwd(), shape_name="corin_pro", test_image,
                           proj = "+proj=longlat +ellps=WGS84") {
    # wdir is the working directory, required
    # shape_name is the base name for the shapefile
    # test_image is the name of a test MODIS satellite image of the area
    # proj is the projection of the shape file
    # Following packages are required
    # make sure all packages are updated
    # this needs an error message? using try?
    try({require(raster)
    require(maptools)
    require(rgeos)
    require(rgdal)})

    # Reading the shape file of the catchment
    shape <- readShapePoly(paste(wdir,paste(shape_name,"shp",sep="."),sep="/"))

    # setting up the projection of the shapefile
    crs(shape) <- proj

    # reading in the test MODIS file
    modis <- raster(paste(wdir,test_image,sep="/")) 
#     image(modis)
#     plot(shape, add = T)

    # Now clipping to the catchment boundary
    modis.clip <- crop(modis, extent(shape)) 
    modis.clip <- mask(modis.clip, shape)

    image(modis.clip)
    plot(shape, add = T)
    # Now extracting the points
    p <- rasterToPoints(modis.clip) # it contains the cell values as well

    #plotting 
    points(p)
  
    # only taking the x(lon) and y(lat) values
    return(data.frame(lon = p[,1], lat = p[,2]))
}

# Example
# Use data in folder Shapefiles_testdata
basedir <- "C:/Users/rver4657/ownCloud/working/SatelliteHydromad"


Points <- Extract.points(wdir=paste(basedir,"/Shapefiles_Testdata",sep=""),
               shape_name="corin_pro",test_image="test1",
               proj = "+proj=longlat +ellps=WGS84")
head(Points)
# write away for later use
# write.csv(Points,paste(basedir,"Corinpoints.csv",sep="/"))
