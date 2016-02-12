#################################################
# Author: Robin Elahi
# Date: 151208
# Temperature manipulation of IPM 
# Figure 5
#################################################

rm(list=ls(all=TRUE)) # removes all previous material from R's memory

##### LOAD PACKAGES ETC #####
library(fields) # need for image.plot
library(ggplot2)
theme_set(theme_classic(base_size = 12))
library(dplyr)
library(RColorBrewer)

source("./R/bael_params.R") 
source("./R/ipmFunctions.R")
source("./R/metabolicTheory.R")
source("./R/multiplotF.R")
source("./R/modify_vital_rates.R")

# Temperatures
modTemp <- 9.25 + 273.15
hisTemp <- 8.6 + 273.15

#### SET GLOBAL PARAMETERS FOR THE SIMULATIONS #####

### Set IPM parameters for Washington population
params <- paramsWA

### Range of activation energies; 0.2 - 1.2
lowerEa <- 0.2; upperEa <- 1.2; Ea_increment <- 0.2

### Range of temperatures; 8-9.5 degrees C
lowerTemp <- 281.5; upperTemp <- 282.65; temp_increment <- 0.05;
originalTemp <- modTemp

##### MODIFYING VITAL RATES WITH THE ARRHENIUS EQUATION: GROWTH #####

growthGrid <- modifyVitalRates(slope = params$growth.slope, slopeTempEffect = "neg", 
                               intercept = params$growth.int, interceptTempEffect = "pos", 
                               lowerEa, upperEa, Ea_increment, 
                               lowerTemp, upperTemp, temp_increment,
                               originalTemp)
head(growthGrid)

##### MODIFYING VITAL RATES WITH THE ARRHENIUS EQUATION: SURVIVAL #####

survivalGrid <- modifyVitalRates(slope = params$surv.slope, slopeTempEffect = "neg", 
                               intercept = params$surv.int, interceptTempEffect = "neg", 
                               lowerEa, upperEa, Ea_increment, 
                               lowerTemp, upperTemp, temp_increment,
                               originalTemp)
head(survivalGrid)

# sample plot of the slope as a function of temperature and Ea ####
ggplot(survivalGrid, aes((Kelvin-273.15),
                  vec_slope, color = as.factor(Ea))) +
  geom_point(alpha = 0.2, size = 2) + theme_bw() + 
  xlab("Temperature (C)") + ylab("") + 
  geom_smooth(se = FALSE) + 
  theme(legend.justification = c(1,1), legend.position = c(1, 1)) +
  guides(color = guide_legend(reverse=TRUE)) + 
  scale_color_discrete(name = "Activation\nEnergy") +
  theme(text = element_text(size = 18)) +
  theme(legend.title = element_text(size = 12)) + 
  theme(legend.text = element_text(size = 12))

##### COLLATE VITAL RATES INTO FINAL DATAFRAME FOR IPMS #####
simulationDF <- growthGrid %>% select(Kelvin:row) %>% 
  rename(simulation = row)

# Growth: constant intercept, but variable slope
growthDF <- cbind(simulationDF, growthGrid$const_int, growthGrid$vec_slope, 
                  survivalGrid$const_int, survivalGrid$const_slope) 

names(growthDF)[4:7] <- c("growth.int", "growth.slope", "surv.int", "surv.slope")
growthDF$parameter <- "Growth"

# Survival: constant intercept, but variable slope
survDF <- cbind(simulationDF, growthGrid$const_int, growthGrid$const_slope, 
                  survivalGrid$const_int, survivalGrid$vec_slope) 

names(survDF)[4:7] <- c("growth.int", "growth.slope", "surv.int", "surv.slope")
survDF$parameter <- "Survival"

# Survival and Growth: variable slopes for both
survGrowthDF <- cbind(simulationDF, growthGrid$const_int, growthGrid$vec_slope, 
                survivalGrid$const_int, survivalGrid$vec_slope) 

names(survGrowthDF)[4:7] <- c("growth.int", "growth.slope", "surv.int", "surv.slope")
survGrowthDF$parameter <- "Survival and growth"

# Combine all scenarios into one dataframe
masterDF <- rbind(growthDF, survDF, survGrowthDF)
masterDF$row <- seq(1:nrow(masterDF))

# Will only use growth slope scenario, 
# because sensitivity analysis suggests only growth matters
masterDF2 <- masterDF %>% filter(parameter == "Growth")

##### RUN IPMS AND EXTRACT POPULATION-LEVEL TRAITS #####
### ACCORDING TO VITAL RATES DEFINED ABOVE

###1.4.4 Make a kernel###
#min.size <- .9*min(c(growthDat$size,growthDat$sizeNext), na.rm=T)
#max.size <- 1.1*max(c(growthDat$size,growthDat$sizeNext), na.rm=T)

min.size <- 0.02
max.size <- 2.02

n <- 100 # number of cells in the matrix

b <- min.size+c(0:n)*(max.size-min.size)/n # boundary points (edges of cells defining the matrix)
b
y <- 0.5*(b[1:n]+b[2:(n+1)]) # mesh points (i.e., midpoints to be used in numerical integration)
y
h <- y[2]-y[1] # step size (i.e., width of cells)
h

### Set up a for loop to run the IPM simulations
loopDat <- masterDF2
AllReps <- unique(loopDat$row)
N <- length(AllReps); N

# Create a matrix to store the results
mat1 <- matrix(nrow = N, ncol = 4)
colnames(mat1) <- c("lambda", "maxSize95", "maxSize99", "meanSize")

# Loop over each row in masterDF
for(g in 1:N){
  row.g <- AllReps[g]
  
  param.g <- data.frame(
    growth.int = loopDat$growth.int[row.g], 
    growth.slope = loopDat$growth.slope[row.g], 
    growth.sd = params$growth.sd, 
    embryo.int = params$embryo.int, 
    embryo.slope = params$embryo.slope, 
    mature.size = - params$embryo.int/params$embryo.slope, 
    recruit.size.mean = params$recruit.size.mean, 
    recruit.size.sd = params$recruit.size.sd, 
    estab.prob.mean = params$estab.prob.mean, 
    surv.int = loopDat$surv.int[row.g], 
    surv.slope = loopDat$surv.slope[row.g])
  
  S <- s.x(y, params = param.g) 						# survival 	
  F <- h*outer(y, y, f.yx, params = param.g)  	# reproduction
  G <- h*outer(y, y, g.yx, params = param.g) 	# growth
  P <- G # placeholder; redefine P on the next line
  # fix eviction of offspring
  for(i in 1:(n/2)) {
    G[1,i] <- G[1,i] + 1-sum(G[,i])
    P[,i] <- G[,i]*S[i]
  }
  
  # fix eviction of large adults
  for(i in (n/2+1):n) {
    G[n,i] <- G[n,i] + 1 - sum(G[,i])
    P[,i]<-G[,i]*S[i]
  }
  K <- P + F  			# full matrix
  # calculate lambda
  lam <- Re(eigen(K)$values[1]) # returns dominant eigenvalue, which gives the asymptotic population growth rate
  # calculate max size (99% of stable size distribution)
  w.eigen <- Re(eigen(K)$vectors[,1])
  stable.dist <- w.eigen/sum(w.eigen)
  maxSize95 <- y[min(which(cumsum(stable.dist) > 0.95))] 
  maxSize99 <- y[min(which(cumsum(stable.dist) > 0.99))] 
  meanSize <- sum(y*stable.dist)
  # populate matrix with continuous variables
  mat1[g,] <- c(lam, maxSize95, maxSize99, meanSize)
}

head(mat1)
mat2 <- cbind(masterDF2, mat1)
mat2$Ea <- as.factor(mat2$Ea)

# Rename simulated dataframe
simDat <- mat2

##### PLOT SIMULATED DATA #####

# Plotting details
label1 <- expression(paste("Maximum size (", cm^2, ")"))
tempLab <- expression(paste("Temperature (", degree, "C)"))
regression <- geom_smooth(method = "lm", se = FALSE, alpha = 0.5, 
                          size = 0.4)
ULClabel <- theme(plot.title = element_text(hjust = -0.07, vjust = 1, 
                                            size = rel(1.5)))
theme_set(theme_classic(base_size = 12))

# Max size plot
maxSizePlot <- ggplot(data = simDat, aes((Kelvin-273.15), maxSize99, linetype = Ea)) +
  xlab(tempLab) + ylab(label1) +
  geom_point(alpha = 0.5, size = 0, color = "white") + 
  geom_smooth(se = FALSE, size = 0.7, color = "black") + 
  theme(legend.justification = "center", legend.position = c(0.9, 0.7)) +
  theme(legend.text = element_text(size = 10)) + 
  theme(legend.title = element_text(size = 10)) + 
  scale_linetype_discrete(name = "Activation\nenergy") + 
  # scale_colour_grey(start = 0.8, end = 0.2) + 
  guides(linetype = guide_legend(reverse=TRUE)) + 
  # guides(color = guide_legend(reverse=TRUE)) + 
  coord_cartesian(ylim = c(0.9, 1.7)) 
maxSizePlot

##### OBSERVED MAXIMUM SIZE: PRESENT #####
# Present: max of corals in 2007 or 2010
# Data from ipmData.csv:
dat <- read.csv("./data/bael_ipmData.csv", header=TRUE, na.strings="NA")
# Select the relevant columns
d <- dat[, which(names(dat) %in% 
                   c("quad", "date", "date.no", "coral.id", "area", 
                     "feret", "code", "sizeOK", "surv", "growth", "recruit"))]

dHIST <- droplevels(d[d$code != "angle" & d$code != "algae" &
                        d$code != "nv" & d$code != "dead", ])
dHIST <- droplevels(dHIST[complete.cases(dHIST$area), ]) # drop NAs

# include only 2007 and 2010: years for IPM
hist07 <- droplevels(dHIST[dHIST$date.no == 39426, ])
hist10 <- droplevels(dHIST[dHIST$date.no == 40515, ])
hist0710 <- rbind(hist07, hist10)
range(hist0710$area)
quantile(hist0710$area, c(0.95, 0.99))

##### OBSERVED MAXIMUM SIZE: PAST #####
# Past: max of corals in 1969 or 1972
# Past size data in histoData.csv:
dat <- read.csv("./data/bael_histoData.csv", header=TRUE, na.strings="NA")
head(dat)

# Initial data
unique(dat$ini.notes)
ini.dat <- droplevels(dat[dat$ini.notes != "angle" & dat$ini.notes != "fuzzy" &
                            dat$ini.notes != "gone" & dat$ini.notes != "nv" &
                            dat$ini.notes != "tentacles", ])

ini.dat <- droplevels(ini.dat[complete.cases(ini.dat$ini.area), ]) # drop NAs

# Final data
unique(dat$fin.notes)
fin.dat <- droplevels(dat[dat$fin.notes != "angle" & dat$fin.notes != "fuzzy" &
                            dat$fin.notes != "gone" & dat$fin.notes != "nv" &
                            dat$fin.notes != "tentacles" & dat$fin.notes != "algae" &
                            dat$fin.notes != "dead" & dat$fin.notes != "overgrown", ])

fin.dat <- droplevels(fin.dat[complete.cases(fin.dat$fin.area), ]) # drop NAs
unique(fin.dat$fin.notes)

ini.sc<- subset(ini.dat, site=="SC")
fin.sc <- subset(fin.dat, site == "SC")

ini.sc$tempC <- with(ini.sc, ifelse(time == "past", hisTemp-273.15, modTemp-273.15))
fin.sc$tempC <- with(fin.sc, ifelse(time == "past", hisTemp-273.15, modTemp-273.15))

# Initial data
ggplot(ini.sc, aes(tempC, ini.area)) +
  geom_violin(alpha = I(0.5), aes(color = time), 
              position = position_jitter(width = 0.01)) 

# Final data
ggplot(fin.sc, aes(tempC, fin.area)) +
  geom_violin(alpha = I(0.5), aes(color = time), 
              position = position_jitter(width = 0.01)) 
summary(ini.sc)

initialSize <- ini.sc %>% group_by(time) %>%
  summarise(maxSize = max(ini.area), 
            max99 = quantile(ini.area, 0.99), 
            max95 = quantile(ini.area, 0.95))
initialSize

finalSize <- fin.sc %>% group_by(time) %>%
  summarise(maxSize = max(fin.area), 
            max99 = quantile(fin.area, 0.99), 
            max95 = quantile(fin.area, 0.95))

initialSize
finalSize

sizeObs <- initialSize
sizeObs$tempC <- c(hisTemp-273.15, modTemp-273.15)

# substitute max size from largest observed size (to match past data)
range(hist0710$area)
quantile(hist0710$area, c(0.95, 0.99))

present_hist0710 <- data.frame(time = "present0710", 
                                  maxSize = max(hist0710$area), 
                                  max99 = quantile(hist0710$area, 0.99), 
                                  max95 = quantile(hist0710$area, 0.95), 
                                  tempC = 9.25)
present_hist0710
sizeObs

sizeObs <- rbind(sizeObs, present_hist0710)
sizeObs
##### GET PREDICTED MAX SIZES  #####
# Empirical historical and modern growth and survival functions 

# Predicted modern size (empirical growth and survival functions)

# what is the predicted size at the modern temp?
modTemp
maxSizeModPred <- simDat[simDat$Kelvin == modTemp &
                             simDat$Ea == 0.6, ]$maxSize99
maxSizeModPred

# what is the predicted size, using the IPM with the historical
# growth slope?
source("bael_IPM_historicGrowth.R")

res0$max99 # should be the same as maxSizeModPred
maxSizeModPred <- res0$max99

# truncated growth curve (historic data)
res1$max99
# full growth curve (historic data)
res2$max99

# use full growth curve
maxSizeHistPred <- res1$max99

# Create data frame with predicted size at modern temp (above)
# and predicted size at historic temperature, using the empirical
# growth function for past data

maxSizePred <- data.frame(
  era = c("Historic", "Modern"),
  tempC = c(hisTemp-273.15, modTemp-273.15),
  size = c(maxSizeHistPred, maxSizeModPred)
)

maxSizePred

##### FINAL PLOTS  #####
# observed points
sizeObs <- sizeObs %>% filter(time != "present")
maxObs <- geom_point(aes(tempC, max99, linetype = NULL), 
                     data = sizeObs,
                     size = 3, shape = 15, color = c("darkgray", 1)) 

maxObs2 <- geom_point(aes(tempC, max99, linetype = NULL), 
                     data = sizeObs,
                     size = 1, shape = 15, color = "white") 

# predicted points 
pPred <- geom_point(aes(tempC, size, linetype = NULL), 
                    data = maxSizePred,
                    size = 3, shape = 17, color = c("darkgray", 1))

pPred2 <- geom_point(aes(tempC, size, linetype = NULL), 
                    data = maxSizePred,
                    size = 1, shape = 17, color = "white")

# final plot
maxSizePlot + maxObs + maxObs2 + pPred + pPred2

ggsave("./figs/ipm_temp.pdf", width = 3.5, height = 3.5)


##### PLOT USING 99%ILE OBSERVED SIZE  #####

sizeObs
library(tidyr)

sizeObsLong <- gather(sizeObs, key = sizeType, value = size_cm, maxSize:max95)

sizeObsLong <- sizeObsLong %>% filter(sizeType == "max99") %>%
  filter(time != "present")
sizeObsLong

obsPoints <- geom_point(aes(tempC, size_cm, by = sizeType, linetype = NULL), 
                    data = sizeObsLong)
maxSizePred

maxSizePlot + pPred + pPred2 + obsPoints

### How different are these observed and predicted values?
sizeObs
maxSizeHistObs <- sizeObs$max99[1]
maxSizeModObs <- sizeObs$max99[2]

maxSizeHistObs  
maxSizeModObs
maxSizeHistPred  
maxSizeModPred

(maxSizeHistPred - maxSizeModPred)/maxSizeModPred

observedChange <- maxSizeHistObs - maxSizeModObs
observedChange/maxSizeHistObs
