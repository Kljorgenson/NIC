## Compile and clip ERA5 climate data

library(ecmwfr)
library(ncdf4)
library(raster) # package for raster manipulation
library(sf) # package for geospatial analysis
library(terra)
library(tidyr)
library(dplyr)
library(data.table)
library(stats)
library(stars)
library(foreach)
library(ggplot2)

#### Import ERA5-Land monthly data from CDS
# AK bounds:  71.4 ,-167 ,54.5,-126

wf_set_key(user ="287249", key = "eb85f0b0-0226-4f58-83e3-7e150b869bb4", service = "cds")

ph <- "C:/Users/kljorgenson/Documents/Repos/AK_discharge_data/ERA5-Land/" 
vr <- c("2m_temperature")
yr <- as.character(2004:2022)
mn <- c("01",  "02",  "03",  "04",  "05",  "06",  "07",  "08",  "09",  "10", "11", "12")
dy <- c(sprintf("0%d",1:9), 10:31)
tm <- c("00:00")

# Request monthly data
for (i in 1:length(yr)) {
  request <- list(
    "dataset_short_name" = 'reanalysis-era5-land-monthly-means',
    "product_type" = "monthly_averaged_reanalysis",
    "variable" = vr[1],
    "year" = yr[i],
    "month" = mn,
    "time" = '00:00',
    "area" = "71.4/-167/54.5/-126",
    "format" = "netcdf",
    "target" = paste0("era5_temp2m_",yr[i],".nc")
  )
  
  file <- wf_request(
    user     = "287249",   # user ID (for authentification)
    request  = request,  # the request
    transfer = TRUE,     # download the file
    path     = ph,       # store data in current working directory
    time_out = 200000,
  )
  
}

# Request hourly data
for (i in 1:length(yr)) {
  for (j in 1:length(mn)) {
    request <- list(
      "dataset_short_name" = "reanalysis-era5-land",
      "product_type" = "reanalysis",
      "variable" = vr[1],
      "year" = yr[i],
      "month" = mn[j],
      "day" = dy,
      "time" = tm,
      "area" = "71.4/-167/54.5/-126",
      "format" = "netcdf",
      "target" = paste0("era5_",yr[i],"_",mn[j],".nc")
    )
    
    file <- wf_request(
      user     = "287249",   # user ID (for authentification)
      request  = request,  # the request
      transfer = TRUE,     # download the file
      path     = ph,       # store data in current working directory
      time_out = 200000
    )
    
  }
}

#### Compile temperature (2m) data
# Followed code from: https://rpubs.com/Ajeep05/era5

# Make a blank stack (using stack()) to store rasterlayer data and use a brick() command to read it.
listfile <- list.files(path="C:/Users/kljorgenson/Documents/Repos/AK_discharge_data/ERA5", pattern="ERA-hourly-temp2m*", full.names = TRUE)
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
#plot(d0)

# Make raster stack into a list or rasters
r.list <- list(setNames(unstack(d0), names(d0)))
plot(d0[[1]]$X2002.01.01)


#### Clip ERA5 data to watersheds and average

## Read in watershed boundaries from shapefile
shape.watershed <- st_read("C:/Users/kljorgenson/Documents/AK_Discharge/Data/Watershed_export/AK_watersheds_58.shp")

plot(d0[[1]]$X2002.01.01)
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
write.csv(temp_df, "ERA5-Land-monthly.csv")



#### Compile daily total precipitation data

# Make a blank stack (using stack()) to store rasterlayer data and use a brick() command to read it.
listfile <- list.files(path="C:/Users/kljorgenson/Documents/Repos/AK_discharge_data/ERA5", pattern="ERA-hourly-ppt*", full.names = TRUE)
listfile


## Covert to raster stack
d0 <- stack()

for (i in 1:length(listfile)) {
  #read data
  filename <- brick(listfile[i])
  
  d0 <- stack(d0,  filename)
}

nlayers(d0)

# Make raster stack into a list or rasters
r.list <- list(setNames(unstack(d0), names(d0)))

## Clip ERA5 data to watersheds and average
## Read in watershed boundaries from shapefile
shape.watershed <- st_read("C:/Users/kljorgenson/Documents/Repos/AK_Discharge/Data/Watersheds_export/AK_watersheds.shp")

## Crop raster to watershed and create mask
# first clip to boundary of all watersheds to save processing time
r2 <- crop(d0, extent(shape.watershed) )
r.shed <- mask(r2, shape.watershed)

## Extract mean of raster values within each watershed to list object
ppt = raster::extract(r.shed, shape.watershed, method="simple", fun=mean, sp=T, df=TRUE 
                       #exact = T #If TRUE the fraction of a cell that is covered is returned or used by fun
)

# Convert SpatialPolygonDataframe into data frame
ppt_df<- as.data.frame(ppt)

# Make data long
ppt_df <- gather(ppt_df, datetime, ppt, 5:ncol(ppt_df)) 
ppt_df$date <- substr(ppt_df$datetime, 2, 11)
ppt_df$date <- as.Date(ppt_df$date, format = "%Y.%m.%d")
ppt_df$date <- format(ppt_df$date, tz = "America/Anchorage")
ppt_df$time <- substr(ppt_df$datetime, 13, 20)

# Take daily sum
ppt.dat <- ppt_df %>% group_by(Id, gridcode, Shape_Leng, Shape_Area, date) %>% summarise(ppt = sum(ppt))

## Export data
write.csv(ppt.dat, "Data/ERA5-daily-ppt-1998-2022-watersheds.csv")




#### Compile daily evaporation data

# Make a blank stack (using stack()) to store rasterlayer data and use a brick() command to read it.
listfile <- list.files(path="C:/Users/kljorgenson/Documents/Repos/AK_discharge_data/ERA5", pattern="ERA-hourly-evap*", full.names = TRUE)
listfile


## Covert to raster stack
d0.e <- stack()

for (i in 1:length(listfile)) {
  #read data
  filename <- brick(listfile[i])
  
  d0.e <- stack(d0.e,  filename)
}

nlayers(d0.e)

# Make raster stack into a list or rasters
r.list <- list(setNames(unstack(d0.e), names(d0.e)))

## Clip ERA5 data to watersheds and average
## Read in watershed boundaries from shapefile
shape.watershed <- st_read("C:/Users/kljorgenson/Documents/Repos/AK_Discharge/Data/Watersheds_export/AK_watersheds.shp")

## Crop raster to watershed and create mask
# first clip to boundary of all watersheds to save processing time
r2.e <- crop(d0.e, extent(shape.watershed) )
r.shed.e <- mask(r2.e, shape.watershed)

## Extract mean of raster values within each watershed to list object
evap = raster::extract(r.shed.e, shape.watershed, method="simple", fun=mean, sp=T, df=TRUE 
                      #exact = T #If TRUE the fraction of a cell that is covered is returned or used by fun
)

# Convert SpatialPolygonDataframe into data frame
evap_df<- as.data.frame(evap)

# Make data long
evap_df <- gather(evap_df, datetime, evaporation, 5:ncol(evap_df)) 
evap_df$date <- substr(evap_df$datetime, 2, 11)
evap_df$date <- as.Date(evap_df$date, format = "%Y.%m.%d")
evap_df$date <- format(evap_df$date, tz = "America/Anchorage")
evap_df$time <- substr(evap_df$datetime, 13, 20)

# Take daily sum
# Change time zone to AKST

evap.dat <- evap_df %>% group_by(Id, gridcode, Shape_Leng, Shape_Area, date) %>% summarise(evaporation = sum(evaporation))

## Export data
write.csv(evap.dat, "Data/ERA5-daily-evap-1998-2022-watersheds.csv")









#### Compile daily snowmelt data

# Make a blank stack (using stack()) to store rasterlayer data and use a brick() command to read it.
listfile <- list.files(path="C:/Users/kljorgenson/Documents/Repos/AK_discharge_data/ERA5", pattern="ERA-hourly-snowmelt*", full.names = TRUE)
listfile


## Covert to raster stack
d0 <- stack()

for (i in 1:length(listfile)) {
  #read data
  filename <- brick(listfile[i])
  
  d0 <- stack(d0,  filename)
}

nlayers(d0)

# Make raster stack into a list or rasters
r.list <- list(setNames(unstack(d0), names(d0)))

## Clip ERA5 data to watersheds and average
## Read in watershed boundaries from shapefile
shape.watershed <- st_read("C:/Users/kljorgenson/Documents/Repos/AK_Discharge/Data/Watersheds_export/AK_watersheds.shp")

## Crop raster to watershed and create mask
# first clip to boundary of all watersheds to save processing time
r2.melt <- crop(d0, extent(shape.watershed) )
r.shed.melt <- mask(r2.melt, shape.watershed)

## Extract mean of raster values within each watershed to list object
snwmlt = raster::extract(r.shed.melt, shape.watershed, method="simple", fun=mean, sp=T, df=TRUE 
                      #exact = T #If TRUE the fraction of a cell that is covered is returned or used by fun
)

# Convert SpatialPolygonDataframe into data frame
snwmlt_df<- as.data.frame(snwmlt)

# Make data long
snwmlt_df <- gather(snwmlt_df, datetime, snowmelt, 5:ncol(snwmlt_df)) 
snwmlt_df$date <- substr(snwmlt_df$datetime, 2, 11)
snwmlt_df$date <- as.Date(snwmlt_df$date, format = "%Y.%m.%d")
snwmlt_df$date <- format(snwmlt_df$date, tz = "America/Anchorage")
snwmlt_df$time <- substr(snwmlt_df$datetime, 13, 20)

# Take daily sum
snwmlt.dat <- snwmlt_df %>% group_by(Id, gridcode, Shape_Leng, Shape_Area, date) %>% summarise(snowmelt = sum(snowmelt))

## Export data
write.csv(snwmlt.dat, "Data/ERA5-daily-snowmelt-1998-2022-watersheds.csv")






#### Compile daily total precipitation data

# Make a blank stack (using stack()) to store rasterlayer data and use a brick() command to read it.
listfile <- list.files(path="C:/Users/kljorgenson/Documents/Repos/AK_discharge_data/ERA5", pattern="ERA-hourly-runoff*", full.names = TRUE)
listfile


## Covert to raster stack
d0 <- stack()

for (i in 1:length(listfile)) {
  #read data
  filename <- brick(listfile[i])
  
  d0 <- stack(d0,  filename)
}

nlayers(d0)

# Make raster stack into a list or rasters
r.list <- list(setNames(unstack(d0), names(d0)))

## Clip ERA5 data to watersheds and average
## Read in watershed boundaries from shapefile
shape.watershed <- st_read("C:/Users/kljorgenson/Documents/Repos/AK_Discharge/Data/Watersheds_export/AK_watersheds.shp")

## Crop raster to watershed and create mask
# first clip to boundary of all watersheds to save processing time
r2.run <- crop(d0, extent(shape.watershed) )
r.shed.run <- mask(r2.run, shape.watershed)

## Extract mean of raster values within each watershed to list object
runoff = raster::extract(r.shed.run, shape.watershed, method="simple", fun=mean, sp=T, df=TRUE 
                      #exact = T #If TRUE the fraction of a cell that is covered is returned or used by fun
)

# Convert SpatialPolygonDataframe into data frame
runoff_df<- as.data.frame(runoff)

# Make data long
runoff_df <- gather(runoff_df, datetime, runoff, 5:ncol(runoff_df)) 
runoff_df$date <- substr(runoff_df$datetime, 2, 11)
runoff_df$date <- as.Date(runoff_df$date, format = "%Y.%m.%d")
runoff_df$date <- format(runoff_df$date, tz = "America/Anchorage")
runoff_df$time <- substr(runoff_df$datetime, 13, 20)

# Take daily sum
runoff.dat <- runoff_df %>% group_by(Id, gridcode, Shape_Leng, Shape_Area, date) %>% summarise(runoff = sum(runoff))

## Export data
write.csv(runoff.dat, "Data/ERA5-daily-runoff-1998-2022-watersheds.csv")

save.image("ERA5_compilation.RData")
