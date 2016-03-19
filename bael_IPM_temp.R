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
source("./R/modify_vital_rates.R")

# Temperatures
sc_meanTemp_pres; sc_meanTemp_past; hms_meanTemp; embryoTemp
modTemp <- as.numeric(sc_meanTemp_pres) + 273.15
hisTemp <- as.numeric(sc_meanTemp_past) + 273.15

#### SET GLOBAL PARAMETERS FOR THE SIMULATIONS #####

### Set IPM parameters for Washington population
params <- paramsOptim

### Range of activation energies; 0.2 - 1.2
lowerEa <- 0.2; upperEa <- 1.2; Ea_increment <- 0.2

### Range of temperatures; 8-9.5 degrees C
lowerTemp <- 281.5; upperTemp <- 282.65; temp_increment <- 0.05;
originalTemp <- modTemp

##### MODIFYING VITAL RATES WITH THE ARRHENIUS EQUATION: GROWTH #####

growthGrid <- modifyVitalRates(slope = params$growth.slope, 
                               slopeTempEffect = "neg", 
                               intercept = params$growth.int, 
                               interceptTempEffect = "pos", 
                               lowerEa, upperEa, Ea_increment, 
                               lowerTemp, upperTemp, temp_increment,
                               originalTemp)
head(growthGrid)

##### MODIFYING VITAL RATES WITH THE ARRHENIUS EQUATION: SURVIVAL #####

survivalGrid <- modifyVitalRates(slope = params$surv.slope, 
                                 slopeTempEffect = "neg", 
                               intercept = params$surv.int, 
                               interceptTempEffect = "neg", 
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
