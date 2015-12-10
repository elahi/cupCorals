#################################################
# Author: Robin Elahi
# Date: 151209
# Coral establishment probability (recruitment over 3 years)
#################################################

# rm(list=ls(all=TRUE)) 

##### LOAD PACKAGES, DATA #####
library(reshape2)

dat <- read.csv("./data/bael_ipmData.csv", header=TRUE, na.strings="NA")

source("./bael_embryos.R") # this also loads bael_functions
source("./R/graphicalParams.R")

### Three steps to calculate establishment probability 
# 1. Calculate the # of embryos produced in quadrat x over 3 years
# 2. Calculate the # of recruits that survived to time 3
# 3. Establishment probability = # of living recruits/ # of total embryos

##### CALCULATE NUMBER OF EMBRYOS PRODUCED IN EACH QUADRAT #####
# Use embryo function for washington
dat$embryo.no <- embryoFwa(dat$area)

# Calculate density for each quadrat and year
density.q <- tapply(dat$area, list(dat$quad, dat$date.no), FUN=densityF)
density.q

### Calculate number of embryos for each quadrat and year
embryo.q <- tapply(dat$embryo.no, list(dat$quad, dat$date.no), FUN=sum, na.rm=TRUE)
colnames(embryo.q) <- c("2007", "2008", "2009", "2010")
# remove the last column, because the embryo no. in 2010 is irrelevant
eq <- embryo.q[,-4] 
eqSum <- rowSums(eq)
eqMult3 <- 3*eq[, 1] # 3 years worth of embryos
eqMult2 <- 2*eq[, 1] # 2 years worth of embryos
eqMult1 <- eq[, 1] # 1 years worth of embryos 2007
plot(eqSum ~ eqMult3)
abline(a = 0, b = 1, col = "red")

##### CALCULATE NUMBER OF RECRUITS THAT SURVIVED TO TIME 3 (2010) #####
names(dat)
# Select the relevant columns
d <- dat[, which(names(dat) %in% c("quad", "date", "date.no", 
                                   "coral.id", "area", "feret", "code", "sizeOK", "surv", "growth", "recruit"))]
# Get data at time 3
unique(d$date.no)
dT03 <- droplevels(d[d$date.no == 39426 | d$date.no == 40515, ])
summary(dT03)

dT0 <- droplevels(dT03[dT03$date.no == 39426, ])
dT3 <- droplevels(dT03[dT03$date.no == 40515, ])

recruitDat <- d[d$recruit == 1,]

# create dataframe with code coral.id for recruitDat
recruitDat2 <- droplevels(recruitDat[, which(names(recruitDat) %in% 
                                               c("coral.id", "code"))])
dim(recruitDat2)
names(recruitDat2) <- c("coral.id", "recruitCode")
head(recruitDat2)
summary(recruitDat2)

dT3b <- merge(dT3, recruitDat2, all.x = TRUE)

# Calculate recruit density in T3
recruit.qT3 <- tapply(dT3b$recruitCode, list(dT3$quad, dT3$date.no), FUN=recruitF)
recruit.qT3

##### CALCULATE ESTABLISHMENT PROBABILITY AFTER 3 YEARS IN 2010 #####
# Sum the total number of embryos produced over three years
eqSum

recProbSum <- (recruit.qT3/eqSum) * 100; 
recProb2 <- recProbSum
recProb2[is.infinite(recProb2)] <- NaN
recProb2


