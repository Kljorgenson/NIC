## Compile and clip ERA5 climate data

library(ecmwfr)
library(ncdf4)
library(raster) # package for raster manipulation
library(sf) # package for geospatial analysis
library(terra)
library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)


#### Compile ERA5 data
# Followed code from: https://rpubs.com/Ajeep05/era5

# Make a blank stack (using stack()) to store rasterlayer data and use a brick() command to read it.
listfile <- list.files(pattern="Data/*.nc", full.names = TRUE)
listfile


## Covert to raster stack
d0 <- stack()

for (i in 1:length(listfile)) {
  
  #read data
  filename <- brick(listfile[i])
  
  d0 <- stack(d0,  filename)
  
  
}

nlayers(d0)

# Plot raster
plot(d0)

# Make raster stack into a list or rasters
r.list <- list(setNames(unstack(d0), names(d0)))
plot(d0[[1]]$X2024.01.01)


#### Clip ERA5 data to watersheds and average

## Read in watershed boundaries from shapefile
shape.watershed <- st_read("Data/Watershed/AK_discharge.shp")

plot(d0[[1]]$X2024.01.01)
plot(shape.watershed, add = T)


## Crop raster to watershed and create mask
# first clip to boundary of all watersheds to save processing time
r2 <- crop(d0, extent(shape.watershed) )
r.shed <- mask(r2, shape.watershed)

plot(d0[[1]][[1]])
plot(r.shed[[1]])
plot(shape.watershed, add = T, lwd = T, col = "transparent")


## Extract mean of raster values within each watershed to list object
temp = raster::extract(r.shed, shape.watershed, method="simple", fun=mean, sp=T, df=TRUE 
                       #exact = T #If TRUE the fraction of a cell that is covered is returned or used by fun
)

# Convert SpatialPolygonDataframe into data frame
temp_df<- as.data.frame(temp)
head(temp_df)

# Make data long
temp_df <- gather(temp_df, date, temp.K, 5:ncol(temp_df), factor_key = T) 
temp_df$date <- substr(temp_df$date, 2, 11)
temp_df$date <- as.Date(temp_df$date, format = "%Y.%m.%d", tz = "UTC") # Data are origionally in UTC
temp_df$date <- format(temp_df$date, tz = "America/Anchorage")
head(temp_df)

# Units are in Kelvin - conver to Celcius
temp_df$temp <- temp_df$temp.K-273.15

temp_df %>% ggplot(aes(date, temp, color = as.factor(gridcode))) + geom_line()


## Export data
write.csv(temp_df, "ERA5-daily.csv")
