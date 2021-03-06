#################################################
# Author: Robin Elahi
# Date: 151013

# Temperature manipulation of IPM growth parameter
# Figure 5
#################################################

library(fields) # need for image.plot
library(ggplot2)
library(dplyr)

# rm(list=ls(all=TRUE)) # removes all previous material from R's memory

source("./R/baelParamsWA.R") 
source("./R/metabolicTheory.R")
source("./R/multiplotF.R")
params

# Temperatures
modTemp <- 9.25 + 273.15
hisTemp <- 8.6 + 273.15

################################################
#### MODIFYING THE IPM
#### GROWTH: INTERCEPT CONSTANT, SLOPE VARIABLE
#### ACCORDING TO ARRHENIUS EQUATION
################################################
params
growthSlope <- params$growth.slope
growthInt <- params$growth.int

# SET UP VECTOR OF ACTIVATION ENERGIES
vec_Ea <- seq(0.2, 1.2, by = 0.2) 

# SET UP VECTOR OF TEMPERATURE
# 0C = 273.15K
# vec_temp <- seq(280.15, 288.15, by = 0.5) # 7 - 15 degrees C

modTemp; hisTemp

vec_temp <- seq(281.5, 282.65, by = 0.05) # 8-9.5 degrees C

# EXPAND THE EaDF BY TEMP GRID
grid1 <- expand.grid(vec_temp, vec_Ea)
head(grid1)
names(grid1) <- c("Kelvin", "Ea")
dim(grid1)
grid1$row <- seq(from = 1, to = dim(grid1)[1], by = 1)
dim(grid1)
summary(grid1)

# Mean temperature for modern corals

# Get pre-expontial coefficient, a, for each Ea FOR INTERCEPT
grid1$intA <- aCoef(growthInt, E = grid1$Ea, k, temp = modTemp, dir = "pos")

# Get pre-expontial coefficient, a, for each Ea FOR SLOPE
grid1$slopeA <- aCoef(growthSlope, E = grid1$Ea, k, temp = modTemp, dir = "neg")

head(grid1)

# Calculate intercept for each temperature, assuming E = 0.65
# Here we assume that the intercept will increase with temperature
# ie, corals grow faster at warmer temps

# Calculate intercepts as a function of temp and Ea
grid1$vec_int <- ArrF(coef_a = grid1$intA, E = grid1$Ea, k, temp = grid1$Kelvin, dir = "pos")

# Calculate slopes as a function of temp and Ea
grid1$vec_slope <- ArrF(coef_a = grid1$slopeA, E = grid1$Ea, k, temp = grid1$Kelvin, dir = "neg")

head(grid1)

grid1$const_int <- rep(growthInt, length = nrow(grid1))
grid1$const_slope <- rep(growthSlope, length = nrow(grid1))

# PLOT PARAMETER VS TEMPERATURE FOR EACH EA
pSlope <- ggplot(grid1, aes((Kelvin-273.15),
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

# PLOT PARAMETER VS TEMPERATURE FOR EACH EA
pInt <- ggplot(grid1, aes((Kelvin-273.15),
                          vec_int, color = as.factor(Ea))) +
  geom_point(alpha = 0.2, size = 2) + theme_bw() + 
  xlab("Temperature (C)") + ylab("") + 
  geom_smooth(se = FALSE) + 
  theme(legend.justification = c(1,1), legend.position = c(1, 1)) +
  guides(color = guide_legend(reverse=TRUE)) + 
  scale_color_discrete(name = "Activation\nEnergy") +
  theme(text = element_text(size = 18)) +
  theme(legend.title = element_text(size = 12)) + 
  theme(legend.text = element_text(size = 12))

pSlope + ylab("Growth Slope")
pInt + ylab("Growth Intercept")

masterGrid <- grid1[, 1:3]

# for growth, I want a constant intercept, but variable slope
masterGrid$intGrowth <- grid1$const_int
masterGrid$slopeGrowth <- grid1$vec_slope
head(masterGrid)

################################################
#### MODIFYING THE IPM
#### CALCULATE LAMBDA AND MAX SIZE
#### ACCORDING TO GROWTH, SURV, EMBRYOS above
################################################

###1.4.4 Make a kernel###
min.size <- .9*min(c(growthDat$size,growthDat$sizeNext), na.rm=T)
max.size <- 1.1*max(c(growthDat$size,growthDat$sizeNext), na.rm=T)

min.size <- 0.02
max.size <- 2.02

n <- 100 # number of cells in the matrix

b <- min.size+c(0:n)*(max.size-min.size)/n # boundary points (edges of cells defining the matrix)
b
y <- 0.5*(b[1:n]+b[2:(n+1)]) # mesh points (i.e., midpoints to be used in numerical integration)
y
h <- y[2]-y[1] # step size (i.e., width of cells)
h

### Now calculate lambda and max size of the stable size distribution after running the ipm
# as a function of different temps
head(grid1)
# for loop
loopDat <- masterGrid
dim(loopDat)
AllReps <- unique(loopDat$row)
AllReps
N <- length(AllReps); N

######################################
######################################
### CHANGE GROWTH SLOPE
######################################
######################################

# matrix
mat1 <- matrix(nrow = N, ncol = 4)
colnames(mat1) <- c("lambda", "maxSize95", "maxSize99", "meanSize")
mat1
head(loopDat)
params

loopDat$intEmbryo[1]/loopDat$slopeEmbryo[1]

# Loop over each row in grid1; get lambda and maxSize
for(g in 1:N){
  row.g <- AllReps[g]
  
  param.g <- data.frame(
    growth.int = loopDat$intGrowth[row.g], 
    growth.slope = loopDat$slopeGrowth[row.g], 
    growth.sd = params$growth.sd, 
    embryo.int = params$embryo.int, 
    embryo.slope = params$embryo.slope, 
    mature.size = - params$embryo.int/params$embryo.slope, 
    recruit.size.mean = params$recruit.size.mean, 
    recruit.size.sd = params$recruit.size.sd, 
    estab.prob.mean = params$estab.prob.mean, 
    surv.int = params$surv.int, 
    surv.slope = params$surv.slope)
  
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
  # population matrix with continuous variables
  mat1[g,] <- c(lam, maxSize95, maxSize99, meanSize)
}

head(mat1)
mat2 <- cbind(masterGrid, mat1)
mat2$Ea <- as.factor(mat2$Ea)

###RENAME DATAFRAME
oneParam <- mat2
head(oneParam)


##########################################################
# Plotting model results
##########################################################
label1 <- expression(paste("Maximum size (", cm^2, ")"))
tempLab <- expression(paste("Temperature (", degree, "C)"))

regression <- geom_smooth(method = "lm", se = FALSE, alpha = 0.5, 
                          size = 0.4)
ULClabel <- theme(plot.title = element_text(hjust = -0.07, vjust = 1, 
                                            size = rel(1.5)))

range(oneParam$maxSize99)

theme_set(theme_classic(base_size = 12))

# Max size 
maxSizePlot <- ggplot(data = oneParam,
                      aes((Kelvin-273.15), maxSize99, color = Ea)) +
  xlab(tempLab) + ylab(label1) +
  geom_point(alpha = 0.5, size = 0) + 
  geom_smooth(se = FALSE, size = 0.7) + 
  theme(legend.justification = c(1,1), legend.position = c(1.1, 1.1)) +
  scale_color_discrete(name = "Activation\nenergy") + 
  guides(color = guide_legend(reverse=TRUE)) + 
  coord_cartesian(ylim = c(0.9, 1.75))
maxSizePlot


##########################################################
# Get the observed max sizes
##########################################################
# Past: max of corals in 1969 or 1972
# Present: max of corals in 2007 or 2010

# Present size data already stored here, from ipmData.csv:
names(hist0710)
range(hist0710$area)

# Past size data in histoData.csv:
dat <- read.csv("./data/bael_histoData.csv", header=TRUE, na.strings="NA")
head(dat)

### Initial data
unique(dat$ini.notes)
ini.dat <- droplevels(dat[dat$ini.notes != "angle" & dat$ini.notes != "fuzzy" &
                            dat$ini.notes != "gone" & dat$ini.notes != "nv" &
                            dat$ini.notes != "tentacles", ])

ini.dat <- droplevels(ini.dat[complete.cases(ini.dat$ini.area), ]) # drop NAs
dim(ini.dat)

### Final data
unique(dat$fin.notes)

fin.dat <- droplevels(dat[dat$fin.notes != "angle" & dat$fin.notes != "fuzzy" &
                            dat$fin.notes != "gone" & dat$fin.notes != "nv" &
                            dat$fin.notes != "tentacles" & dat$fin.notes != "algae" &
                            dat$fin.notes != "dead" & dat$fin.notes != "overgrown", ])

fin.dat <- droplevels(fin.dat[complete.cases(fin.dat$fin.area), ]) # drop NAs
unique(fin.dat$fin.notes)

dim(fin.dat)

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

initialSize <- ini.sc %>% group_by(time) %>%
  summarise(maxSize = max(ini.area))

finalSize <- fin.sc %>% group_by(time) %>%
  summarise(maxSize = max(fin.area))

initialSize
finalSize

sizeObs <- initialSize
sizeObs$tempC <- c(hisTemp-273.15, modTemp-273.15)

# substitute max size from largest observed size (to match past data)
sizeObs$maxSize[2] <- max(hist0710$area)
sizeObs

##########################################################
# Get the predicted mean and max sizes
# Using the historical and modern growth slopes 
##########################################################

# what is the predicted size at the modern temp?
modTemp
maxSizeModPred <- oneParam[oneParam$Kelvin == modTemp &
                             oneParam$Ea == 0.6, ]$maxSize99
maxSizeModPred

# what is the predicted size, using the IPM with the historical
# growth slope?
source("bael_IPM_historicGrowth.R")
maxSizeHistPred <- res1$max99

res0$max99 # should be the same as maxSizeModPred
maxSizeModPred <- res0$max99

# Create data frame with predicted size at modern temp (above)
# and predicted size at historic temperature, using the empirical
# growth function for past data

maxSizePred <- data.frame(
  era = c("Historic", "Modern"),
  tempC = c(hisTemp-273.15, modTemp-273.15),
  size = c(maxSizeHistPred, maxSizeModPred)
)

##########################################################
# Create ggplot objects
##########################################################
sizeObs
# observed points
maxObs <- geom_point(aes(tempC, maxSize), 
                     data = sizeObs,
                     size = 3, shape = 16, color = c("darkgray", 1))

# predicted points 
pPred <- geom_point(aes(tempC, size), 
                    data = maxSizePred,
                    size = 3, shape = 17, color = c("darkgray", 1))

##########################################################
# Create ggplot objects
##########################################################

maxSizePlot + maxObs + pPred

ggsave("./figs/ipm_temp.pdf", width = 3.5, height = 3.5)

##########################################
# How different are these observed and predicted values?
sizeObs
maxSizeHistObs <- sizeObs$maxSize[1]
maxSizeModObs <- sizeObs$maxSize[2]

maxSizeHistObs  
maxSizeModObs
maxSizeHistPred  
maxSizeModPred

observedChange <- maxSizeHistObs - maxSizeModObs
  
##########################################


