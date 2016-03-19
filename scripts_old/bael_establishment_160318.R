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

# Calculate density for each quadrat and year
density.q <- tapply(dat$area, list(dat$quad, dat$date.no), FUN=densityF)
density.q

str(dat)
unique(dat$date)

# Use embryo function for Washington and California
dat$embryosWA <- embryoFwa(dat$area)
dat$embryosCA <- embryoFca(dat$area)

# Calculate the number of embryos produced by each coral
dat$embryos <- embryoF(dat$area, xIntercept = xIntWA, mxRegression = mxRegWA)

# Remove 2010 data, and then get total number of embryos per quad

quad_embryosDF <- dat %>% filter(date.no != 40515) %>% 
  group_by(quad) %>%
  summarise(quadEmbryos = sum(embryos, na.rm = TRUE))
  

### Calculate number of embryos for each quadrat and year - WA
embryoWA.q <- tapply(dat$embryosWA, list(dat$quad, dat$date.no), FUN=sum, na.rm=TRUE)
colnames(embryoWA.q) <- c("2007", "2008", "2009", "2010")
# remove the last column, because the embryo no. in 2010 is irrelevant
eqWA <- embryoWA.q[,-4] 
eqWASum <- rowSums(eqWA)

### Calculate number of embryos for each quadrat and year - CA
embryoCA.q <- tapply(dat$embryosCA, list(dat$quad, dat$date.no), FUN=sum, na.rm=TRUE)
colnames(embryoCA.q) <- c("2007", "2008", "2009", "2010")
# remove the last column, because the embryo no. in 2010 is irrelevant
eqCA <- embryoCA.q[,-4] 
eqCASum <- rowSums(eqCA)

##### CALCULATE NUMBER OF RECRUITS THAT SURVIVED TO TIME 3 (2010) #####

### NEW WAY
source("R/get_estab_prob.R")

# These are the observed recruits every year
recruitDat <- dat %>% filter(recruit ==1) %>% select(coral.id, code) %>%
  rename(recruitCode = code)
head(recruitDat)

quad_recruitsDF <- dat %>% filter(date.no == 40515) %>% 
  left_join(., recruitDat, by = "coral.id") %>%
  group_by(quad) %>%
  summarise(quadRecruits = recruitF(recruitCode)) 

# Join embryos with recruits, then calculate establishment probability
quad_DF <- inner_join(quad_recruitsDF, quad_embryosDF) %>% 
  mutate(recProb = quadRecruits/quadEmbryos * 100)

quad_DF$recProb[is.infinite(quad_DF$recProb)] <- NaN
quad_DF

### OLD WAY
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
eqCASum; eqWASum

### For California (original)
recProbSum <- (recruit.qT3/eqCASum) * 100; 
recProb2 <- recProbSum
recProb2[is.infinite(recProb2)] <- NaN
recProbCA <- recProb2

### For Washington (temperature-adjusted)
recProbSum <- (recruit.qT3/eqWASum) * 100; 
recProb2 <- recProbSum
recProb2[is.infinite(recProb2)] <- NaN
recProb2
recProbWA <- recProb2

plot(recProbWA ~ recProbCA)

