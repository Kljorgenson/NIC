#### MARSS models with ERA5 covariates
library(ggplot2)
library(dplyr)
library(lubridate)
library(discharge)
library(tidyverse)
library(data.table)
library(R2jags)
library(rjags)

## Read in tower dates
dat <- read.csv("Raw data/NenanaIceClassic_1917-2021.csv") %>% rename(year = Year, doy = Decimal_doy)
head(dat)

## PDO data
PDO <- read.table("Raw data/PDO.dat.txt", skip = 1, header = T) %>% pivot_longer(cols = Jan:Dec, names_to = "Month", values_to = "PDO")

# Take yearly average and crop to length of variables
PDO <- PDO %>% filter(Month %in% c("Jan", "Feb", "Mar")) %>% group_by(Year) %>% summarise(PDO = mean(PDO)) %>% filter(Year >= 1917 & Year <= 2024) %>% rename(year = Year)


### Covariates
covars <- read.csv("covars.csv")
covars <- left_join(covars, PDO)

data <- left_join(covars, dat) %>% filter(year < 2024) # I am redoing 2021-2023 temp data
head(data)

## Plot doy
data %>% ggplot(aes(year, doy)) + geom_point()

## Simple linear model
lm <- lm(data$doy ~ data$temp + data$runoff + data$ice_depth + data$SWE + data$PDO)
summary(lm)

pred <- predict(lm, data.frame(temp = 13, runoff = 0.007, ice_depth = 11, SWE = 31.7, PDO = 0))
mean(pred)
sd(pred)



## I am not sure how to predict from this model

### Model with rjags
writeLines("
  model {
  
  ## Likelihood
  for (i in 1:N) {
    mu[i] <- alpha + B1*temp[i] + B2*runoff[i] + B3*ice_depth[i] + B4*SWE[i] + B5*PDO[i]
    doy[i]~dnorm(mu[i], tau)
  }
  
  ## Priors
  alpha ~ dnorm(130, 100)
  B1 ~ dnorm(0, 100)
  B2 ~ dnorm(0, 2)
  B3 ~ dnorm(0, 100)
  B4 ~ dnorm(0, 100)
  B5 ~ dnorm(0, 10)
  
  
  sigma~dunif(0, 100) # Residual standard deviation
  tau <- 1/(sigma*sigma) # Residual precision
  
}
", con = "mod.jags")


jagsdata <- with(data[,2:7], list(doy = doy , temp = temp, runoff = runoff, ice_depth = ice_depth, SWE = SWE, PDO = PDO, N = length(doy)))

# create a jags model object
jags <- jags.model("mod.jags",
                   data = jagsdata,
                   n.chains = 4,
                   n.adapt = 100)

# burn-in
update(jags, 1000)

# draw samples
samples <- jags.samples(jags, c("alpha", "B1", "B2", "B3", "B4", "B5"), 1000)
str(samples)

# extract posterior means from the mcarray object by marginalizing over chains and iterations (alternative: posterior modes)
posterior_means <- lapply(samples, apply, 1, "mean")


## Hmm I don't know how to predict the correct way
# Just mean, no error, predict for 2024
new <- data.frame(temp = 13, runoff = 0.007, ice_depth = 11, SWE = 31.7, PDO = 0)
posterior_means$alpha + posterior_means$B1*new[1,1] + posterior_means$B2*new[1,2] + posterior_means$B3*new[1,3] + posterior_means$B4*new[1,4]






### Better way?
# take our posterior means 
B <- as.matrix(unlist(posterior_means[c("alpha", "B1", "B2", "B3", "B4", "B5")]))

# create a model matrix from x
X <- cbind(1, data.frame(temp = 13, runoff = 0.007, ice_depth = 11, SWE = 31.7, PDO = 0))

# predicted outcomes are the product of our model matrix and estimates
y_hat <- X %*% B ## Hmm idk how to do this
head(y_hat)


### Model checks

traceplot(fit)













