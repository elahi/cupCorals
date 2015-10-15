#################################################
# Author: Robin Elahi
# Date: 151015

# Modifed IPM parameterized by modern data
# But using the slope of the growth function
# From historical growth data
#################################################

# rm(list=ls(all=TRUE)) # removes all previous material from R's memory
source("./R/baelParamsWA.R")
source("./R/baelParamsCA.R")
source("./R/ipmFunctions.R")
paramsWA
paramsCA

### Create params DF with modified growth function
# Get growth data
source("bael_growth.R")

# Summary of modeled parameters
summary(pastMod)
summary(presMod)

# Fitted coefficients of general linear model
fixef(pastMod) 

# Create new paramsDF  
paramsWA_hist <- paramsWA

# Substitute growth slope from historical data
paramsWA_hist$growth.slope <- fixef(pastMod)[2] 

# Create one more new DF that also has historical intercept
paramsWA_hist2 <- paramsWA_hist
paramsWA_hist2$growth.int <- fixef(pastMod)[1] 

##########################################
##########################################
# MAX AND MIN SIZES
names(hist0710)
range(hist0710$area)

min.size <- .9*min(hist0710$area, na.rm=T)
max.size <- 1.1*max(hist0710$area, na.rm=T)

min.size; max.size

# Will use slightly larger size range for IPM because
# it is otherwise artificially truncated
min.size <- 0.02
max.size <- 2

binSize <- 0.02
binN <- (max.size - min.size)/binSize
binN

##########################################
##########################################
### RUN IPMS
# Basic IPM, paramsWA_hist
ipm0 <- bigmatrix(n = binN, params = paramsWA)
res0 <- popF(ipm0, binSize)

# Basic IPM, paramsWA_hist (historical growth slope)
ipm1 <- bigmatrix(n = binN, params = paramsWA_hist)
res1 <- popF(ipm1, binSize)

# Check for eviction
plot(ipm1$meshpts, s.x(ipm1$meshpts, params), xlab = "Size", type = "l", ylab = "Survival probability", lwd = 12) # plot the survival model
points(ipm1$meshpts, apply(ipm1$P, 2, sum), col = "red", lwd = 3, cex = 0.5, pch = 19) # plot the column sums of the survival/growth matrix

# Basic IPM, paramsWA_hist (historical growth slope and intercept)
ipm2 <- bigmatrix(n = binN, params = paramsWA_hist2)
res2 <- popF(ipm2, binSize)

### Plot 
range(hist0710$area)
hist(hist0710$area, breaks = seq(0.0, 2, 0.05), freq = FALSE, 
	xlab = xlab2, ylab="Probability density", col = "gray87", main = "", 
	ylim = c(-0.025, 2), xlim = c(0, 2), border = "gray90", las = 1) 
box()

# Original IPM with paramsWA
points(ipm0$meshpts , res0$ssd2, type="l", lty=1, lwd = 2, 
       col = "black")

# Modified growth slope
points(ipm1$meshpts , res1$ssd2, type="l", lty=1, lwd = 2, 
		col = "darkgray")

# Modified growth slope and intercept (doesn't change much)
points(ipm2$meshpts , res2$ssd2, type="l", lty=1, lwd = 2, 
		col = "red")

arrows(res0$max99, 0.6, res0$max99, 0.3, col = "black", 
       length = 0.1, lwd = 1.5, angle = 20)
arrows(res1$max99, 0.6, res1$max99, 0.3, col = "darkgray", 
 		length = 0.1, lwd = 1.5, angle = 20)
arrows(res2$max99, 0.6, res2$max99, 0.3, col = "red", 
 		length = 0.1, lwd = 1.5, angle = 20)		

# legend("topright", leg.txt, lwd = 2, bty = "n", 
# 	col = c("black", "darkgray"), lty = 1, 
# 	cex = 1, text.col = c("black", "darkgray"))


