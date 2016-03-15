#################################################
# Author: Robin Elahi
# Date: 151209
# Coral embryos
#################################################

# rm(list=ls(all=TRUE)) 

##### LOAD PACKAGES, DATA #####
library(ggplot2)
theme_set(theme_classic(base_size = 16))
library(dplyr)

source("./R/baelFunctions.R")
source("./R/metabolicTheory.R")
source("bael_temperature.R")

# Need to estimate the relationship between coral size and embryo #
# use data in Fadlallah Figure 5 and Table 6
# Figure 5 gives # of embryos per female, the raw data
# Table 6 assumes a 1:1 ratio of males:females, and thus divides the
# number of embryos by 2 to get mx per coral

##### AREA VOLUME RELATIONSHIPS #####
# upload empirical data on area volume relationships
# because I have calyx area, but Fadlallahs data are volume

avDat <- read.csv("./data/bael area volume.csv", 
                  header=TRUE, na.strings = "NA")

# Regression for ImageJ area vs calculated area
aaReg <- lm(areaIJmm2 ~ areaCALCmm2, data = avDat)
plot(areaIJmm2 ~ areaCALCmm2, data = avDat)
abline(a = 0, b = 1, col = "red", lty = 2)
abline(aaReg, col = "blue") 

# Regression for calculated volume vs ImageJ area (log both x and y)
avReg <- lm(log(areaIJmm2) ~ log(volCALCmm3), data = avDat)
summary(avReg)

plot(log(areaIJmm2) ~ log(volCALCmm3), data = avDat)
abline(a = 0, b = 0.66, col = "red", lty = 2) # expect a slope of 2/3
abline(avReg, col = "blue") 

##### LIFE TABLE DATA #####
# Upload bael life table
ltDat <- read.csv("./data/bael life table.csv", 
                  header=TRUE, na.strings="NA")
head(ltDat)
summary(ltDat)

# Calculate area for corals in life table based on regression
ltDat$areaR <- (exp(log(ltDat$Vol.mid)*avReg$coefficients[2] + 
                      avReg$coefficients[1]))/100
# biggest coral is 700 mm3  
areaCalcF(700)

# Compare to calculations in excel
head(ltDat)
plot(areaR ~ areaREGR, ltDat)

# Plot calculated area vs regressed area
plot(areaCALC ~ areaREGR, ltDat)
abline(a = 0, b = 1, col = "red") # calculated volume underestimates 

# Estimate the relationship between areaR and mx.coral
head(ltDat)
plot(mx.coral ~ Vol.mid, data = ltDat)
plot(mx.coral ~ areaR, data = ltDat)
ltDat

##### EMBRYOS VS SIZE #####
# Original California relationship
mxReg <- lm(mx.coral ~ areaR, data = ltDat[ltDat$Vol.mid > 75, ])
summary(mxReg)

plot(mx.coral ~ areaR, data = ltDat, xlim = c(0, 1.5), ylim = c(0, 50),
     xlab = "Size (cm2)", ylab = "No. of embryos", las = 1)
abline(mxReg, col = "darkgray")

xx <- seq(0, 1.2, by = 0.1)
mxReg$coefficients

# Adjust for temperature to obtain Washington relationship
yIntCA <- mxReg$coefficients[1]
slopeCA <- mxReg$coefficients[2]
xIntCA <- -mxReg$coefficients[1]/mxReg$coefficients[2]

# Get a_coef for the intercept of the embryo function
kelvin_CA <- as.numeric(embryoTemp) + 273.15
kelvin_WA <- as.numeric(sc_meanTemp_pres) + 273.15

y_intercept_CA <- mxReg$coefficients[1]
y_intercept_CA

embryo_A <- aCoef(y_intercept_CA, E, k, temp = kelvin_CA, dir = "neg")
y_intercept_WA <- ArrF(embryo_A, E, k, temp = kelvin_WA, dir = "neg")

xIntWA <- -y_intercept_WA/mxReg$coefficients[2]
xIntWA

sizeDiff <- xIntWA - xIntCA; sizeDiff

ltDat$areaRcorr <- ltDat$areaR + sizeDiff
mxRegWA <- lm(mx.coral ~ areaRcorr, data = ltDat[ltDat$Vol.mid > 75, ])
mxRegCA <- lm(mx.coral ~ areaR, data = ltDat[ltDat$Vol.mid > 75, ])

mxRegWA
mxRegCA

plot(mx.coral ~ areaR, data = ltDat, xlim = c(0, 1.5), ylim = c(0, 50),
     xlab = "Size (cm2)", ylab = "No. of embryos", las = 1, pch = 21, bg = 1)
abline(mxRegCA)

points(mx.coral ~ areaRcorr, data = ltDat, xlim = c(0, 1.5), ylim = c(0, 50),
       xlab = "Size (cm2)", ylab = "No. of embryos", las = 1)
abline(mxRegWA, lty = 2)
