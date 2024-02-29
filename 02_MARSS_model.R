#### MARSS models with ERA5 covariates
library(dataRetrieval)
library(ggplot2)
library(dplyr)
library(lubridate)
library(waterData)
library(MARSS)
library(discharge)
library(tidyverse)
library(data.table)

## Read in data
# Read in discharge data
flow.dat <- read.csv("C:/Users/kljorgenson/Documents/Repos/AK_discharge_data/Discharge_data/DFFT_2011-2022.csv")
# Add region to data
region <- read.csv("Data/USGS_stations.csv") %>% select(name, region, station) %>% unique() %>% rename("site" = "name")
dat <- left_join(flow.dat, region)
dat %>% ungroup() %>% select(site, region)  %>% unique() #%>% view()
unique(dat$site)

# Site metadata
names <- read.csv("Data/USGS_stations.csv") %>% filter(name %in% dat$site)

metadata <- read.csv("Data/USGS_station_metadata.csv") %>% filter(site_no %in% names$station)
head(metadata)
length(metadata$station_nm)

write.csv(metadata, "Selected_station_metadata.csv")

# Read in covariate data
ppt_df <- read.csv("C:/Users/kljorgenson/Documents/Repos/AK_discharge_data/ERA5/ERA5-daily-ppt-1998-2022-watersheds.csv") %>% select(gridcode, ppt, date)
names(ppt_df) <- c("station", "ppt", "date")

# Join discharge and covariate data
data <- left_join(dat, ppt_df) %>% filter(date <= "2022-04-25") %>% unique()
head(data)

# I don't have all these stations currently
data <- data %>% filter(is.na(ppt) == F)


data$date <- as.Date(data$date)


### MARSS models
# Daily
# test
data <- data %>% filter(site %in% c("CHENA_1", "SALCHA"))
data$season <- ifelse(data$jday >= 91 & data$jday <= 304, 1, 0)

vars<-matrix(data$resid.sig, nrow = length(unique(data$site)),byrow=TRUE)
dim(vars)

# seasonal covariate

vals <- unique(data$season)
TT <- length(data[data$site == "CHENA_1",17])
p <- length(vals)
c <- matrix(0, p, TT)
for(i in 1:p) c[i,] <- data[data$site == "CHENA_1",17] == vals[i]
rownames(c) <- vals
dim(c)

c <- data[data$site == "CHENA_1",17]
covars <- rbind(matrix(data$ppt, nrow = length(unique(data$site)),byrow=TRUE), c)
dim(covars)


any(is.na(vars)) 
any(is.na(covars)) # OK

# C matrix
l <- list(0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0)
ppt.l <- as.list(c("ppt1", l, "ppt2", l,"ppt3", l, "ppt4", l,"ppt5", l, "ppt6", l,"ppt7", l, "ppt8", l,"ppt9", l, "ppt10", l,"ppt11", l, "ppt12", l,"ppt13", l, "ppt14", l,"ppt15", l, "ppt16", l,"ppt17", l, "ppt18", l,"ppt19", l, "ppt20", l,"ppt21", l, "ppt22", l,"ppt23", l, "ppt24", l,"ppt25", l, "ppt26", l,"ppt27", l, "ppt28"))
ppt.l <- list("ppt1", 0,0,"ppt2")
ppt.l

C0 <- matrix(ppt.l, ncol = 2)
dim(C0)

C <- cbind(C0, as.list(rep("S", 2)))
dim(C)



## Model with separate processes for each river
mod3 = list(
  R = "diagonal and equal",
  B = "identity", # Mean regression
  A = "zero",
  U = matrix(rep("u", 2)),
  Z = "identity",
  c = covars,
  C = C,
  Q = "diagonal and equal"
)


## FIT MODEL
mod3.fit = MARSS(vars, model=mod3)
MARSSparamCIs(mod3.fit)

acf(t(as.data.frame(mod3.fit$ytT))[,12])
 plot(mod3.fit)

 
 
 # Monthly
 
 ## PDO data
 PDO <- read.table("Data/ersst.v5.pdo.dat.txt", skip = 1, header = T) %>% pivot_longer(cols = Jan:Dec, names_to = "Month", values_to = "PDO")
 PDO$date <- paste(PDO$Year, PDO$Month, "1") # Put day as 1st of month
 PDO$date <- strptime(PDO$date, format = "%Y %b %d", tz = "America/Anchorage")
 PDO$month <- format(as.Date(PDO$date), "%m")
 PDO$month <- as.numeric(PDO$month)
 names(PDO) <- c("year", "Month", "PDO", "Date", "month")
 head(PDO)
 

 dat <- data %>% group_by(site, year, month) %>% summarise(ppt = sum(ppt),resid.sig = mean(resid.sig))
 dat$season <- ifelse(dat$month %in% 5:10, 1, 0)
 cos.t <- cos(2 * pi * seq(dim(vars)[2])/12)
 sin.t <- sin(2 * pi * seq(dim(vars)[2])/12)
 c.Four <- rbind(cos.t, sin.t)
 
dat <- left_join(dat, PDO)
  
 vars<-matrix(dat$resid.sig, nrow = length(unique(dat$site)),byrow=TRUE)
 dim(vars)
 
 c <- dat[dat$site == "CHENA_1",6]
 
 covars <- rbind(matrix(dat$ppt, nrow = length(unique(dat$site)),byrow=TRUE), t(as.matrix(dat[dat$site == "CHENA_1",8])), c.Four)

 dim(covars)
 any(is.na(vars)) 
 any(is.na(covars)) # OK
 

 # C matrix
 n = 28 # number of sites
 l <- list(0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0)
 ppt.l <- as.list(c("ppt1", l, "ppt2", l,"ppt3", l, "ppt4", l,"ppt5", l, "ppt6", l,"ppt7", l, "ppt8", l,"ppt9", l, "ppt10", l,"ppt11", l, "ppt12", l,"ppt13", l, "ppt14", l,"ppt15", l, "ppt16", l,"ppt17", l, "ppt18", l,"ppt19", l, "ppt20", l,"ppt21", l, "ppt22", l,"ppt23", l, "ppt24", l,"ppt25", l, "ppt26", l,"ppt27", l, "ppt28"))
 #ppt.l <- list("ppt1", 0,0,"ppt2")
 ppt.l
 
 C0 <- matrix(ppt.l, ncol = n)
 dim(C0)
 
 C <- cbind(C0, as.list(rep("PDO", n)), as.list(paste0("seas1", 1:28)), as.list(paste0("seas2", 1:28)))
 dim(C)
 
 
 ## Model with separate processes for each river
 mod3 = list(
   R = "diagonal and equal",
   B = "identity", # Mean regression
   A = "zero",
   U = matrix(as.list(paste0("u", 1:n))),
   Z = "identity",
   c = covars,
   C = C,
   Q = "diagonal and equal"
 )
 
 
 ## FIT MODEL
 mod3.fit = MARSS(vars, model=mod3)
 MARSSparamCIs(mod3.fit)
 
 acf(t(as.data.frame(mod3.fit$ytT))[,12])
 plot(mod3.fit)
 
 
 ### Display model output from monthly model
 # Plot effect sizes
 coefs <- mod3.fit$coef
 
 
 # Plot origional observations and model fit
 
 
