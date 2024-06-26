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
# temp = average over area, monthly average
# runoff = average over area, monthly sum
# SWE = average over area, monthly average
# ice_depth = average over area, monthly average



# Make a blank stack (using stack()) to store rasterlayer data and use a brick() command to read it.
listfile <- list.files(path = "Data/", pattern="runoff_2023.nc", full.names = TRUE)
listfile


## Covert to raster stack
d0 <- stack()

for (i in 1:length(listfile)) {
  
  #read data
  filename <- brick(listfile[i])
  
  d0 <- stack(d0,  filename)
  
}

nlayers(d0)


# Make raster stack into a list of rasters
r.list <- list(setNames(unstack(d0), names(d0)))
plot(d0[[1]]$X2023.08.01)


#### Clip ERA5 data to watersheds and average

## Read in watershed boundaries from shapefile
shape.watershed <- st_read("Data/Watershed/AK_discharge.shp")

plot(d0[[1]]$X2023.08.01)
plot(shape.watershed)


## Extract mean of raster values within each watershed to list object
temp = raster::extract(d0, shape.watershed, method="simple", fun=mean, sp=T, df=TRUE)

# Convert SpatialPolygonDataframe into dataframe
temp_df<- as.data.frame(temp)
head(temp_df)

# Make data long
temp_df <- gather(temp_df, date, runoff, 5:ncol(temp_df), factor_key = T) 


# Select daily rows and take monthly average
temp_df$year = substr(temp_df$date, start = 2, stop = 5)
temp_df$month = substr(temp_df$date, start = 7, stop = 8)
temp_df <- temp_df %>% filter(month == "03") %>% group_by(year) %>% summarise(runoff = sum(runoff))


temp_df %>% ggplot(aes(year, temp)) + geom_line()


## Export data
write.csv(temp_df, "ERA5-Mar-runoff_23.csv", row.names = F)





#### Join all data together
temp1 <- read.csv("ERA5-Mar-temp_early.csv")
temp2 <-read.csv("ERA5-Mar-temp_2002-2022.csv")

temp <- rbind(temp1, temp2) %>% select(year, temp)


runoff <- read.csv("ERA5-Mar-runoff.csv")


# lake ice depth
ice <- read.csv("ERA5-Mar-icedepth.csv") %>% select(year, ice_depth)

# SWE
swe <- read.csv("ERA5-Mar-SWE.csv")

# Join all vars
vars <- full_join(temp, runoff)
vars <- full_join(vars, ice)
vars <- full_join(vars, swe)
head(vars)

vars %>% ggplot(aes(year, temp)) + geom_line() + geom_line(aes(year, runoff), color = "green") + geom_line(aes(year, ice_depth), color = "blue") +
  geom_line(aes(year, SWE), color = "yellow")


write.csv(vars, "covars.csv", row.names = F)
