#################################################
# Author: Robin Elahi
# Date: 151210
# Modifed IPM parameterized by modern data
# But using the slope of the growth function
# From historical growth data
#################################################

# rm(list=ls(all=TRUE)) # removes all previous material from R's memory

##### LOAD PACKAGES, DATA #####
# source("./R/bael_params.R")
# source("./R/ipmFunctions.R")

# Create params DF with modified growth function
# Get growth data
source("bael_growth.R")

# Summary of modeled parameters
summary(pastMod)
summary(pastModAll)
summary(presMod)

# Fitted coefficients of general linear model
fixef(pastMod) 
fixef(pastModAll) 

# Create new paramsDF  
paramsWA_hist <- paramsWA
paramsWA_histAll <- paramsWA

# Substitute growth slope from historical data
paramsWA_hist$growth.slope <- fixef(pastMod)[2] 
paramsWA_histAll$growth.slope <- fixef(pastModAll)[2] 

paramsWA_hist$growth.sd <- fixef(pastMod)[1] 
paramsWA_histAll$growth.sd <- fixef(pastModAll)[1] 

paramsWA_hist
paramsWA_histAll

##### RUN IPMS #####

min.size <- .9*min(hist0710$area, na.rm=T)
max.size <- 1.1*max(hist0710$area, na.rm=T)

# Will use larger size range for IPM because
# it is otherwise artificially truncated
min.size <- 0.02
max.size <- 2

binSize <- 0.02
binN <- (max.size - min.size)/binSize
binN

# Basic IPM, paramsWA
ipm0 <- bigmatrix(n = binN, params = paramsWA)
res0 <- popF(ipm0, binSize)
res0$max99

# Basic IPM, paramsWA_hist (historical growth slope, truncated)
ipm1 <- bigmatrix(n = binN, params = paramsWA_hist)
res1 <- popF(ipm1, binSize)
res1$max99

# Basic IPM, paramsWA_histAll (historical growth slope, all data)
ipm2 <- bigmatrix(n = binN, params = paramsWA_histAll)
res2 <- popF(ipm2, binSize)
res2$max99

# Check for eviction
plot(ipm2$meshpts, s.x(ipm1$meshpts, params), xlab = "Size", type = "l", ylab = "Survival probability", lwd = 12) # plot the survival model
points(ipm2$meshpts, apply(ipm1$P, 2, sum), col = "red", lwd = 3, cex = 0.5, pch = 19) # plot the column sums of the survival/growth matrix

