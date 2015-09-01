# RE 150521
# This script will calculate model components for an IPM
# FOR WASHINGTON POPULATION
# USES CORRECTED EMBRYO FUNCTION FOR
# ESTIMATATION OF FECUNDITY FUNCTION, 
# SIZE AT MATURITY, 
# ESTABLISHMENT PROBABILITY

library(ggplot2)
library(reshape2)
library(lme4)

#rm(list=ls(all=TRUE)) # removes all previous material from R's memory

source("./R/baelFunctions.R")
source("./R/graphicalParams.R")

dat <- read.csv("./data/bael_ipmData.csv", header=TRUE, na.strings="NA")

names(dat)
# Select the relevant columns
d <- dat[, which(names(dat) %in% c("quad", "date", "date.no", 
        "coral.id", "area", "feret", "code", "sizeOK", "surv", "growth", "recruit"))]

############################################
# get appropriate data for histogram in 2007-2010
dim(d)

dHIST <- droplevels(d[d$code != "angle" & d$code != "algae" &
                        d$code != "nv" & d$code != "dead", ])
dHIST <- droplevels(dHIST[complete.cases(dHIST$area), ]) # drop NAs

# examine size frequency distribution for each year
set_graph_pars(ptype = "panel4")
hist(dHIST[dHIST$date.no == 39426, ]$area, xlab = "Size", main = "", breaks = seq(0, 1.4, 0.05), ylim = c(0, 3), freq = FALSE, col =  "darkgray")
hist(dHIST[dHIST$date.no == 39797, ]$area, xlab = "Size", main = "", breaks = seq(0, 1.4, 0.05), ylim = c(0, 3), freq = FALSE, col =  "darkgray")
hist(dHIST[dHIST$date.no == 40160, ]$area, xlab = "Size", main = "", breaks = seq(0, 1.4, 0.05), ylim = c(0, 3), freq = FALSE, col =  "darkgray")
hist(dHIST[dHIST$date.no == 40515, ]$area, xlab = "Size", main = "", breaks = seq(0, 1.4, 0.05), ylim = c(0, 3), freq = FALSE, col =  "darkgray")

# include only 2007 and 2010: years for IPM
hist07 <- droplevels(dHIST[dHIST$date.no == 39426, ])
dim(hist07)
hist10 <- droplevels(dHIST[dHIST$date.no == 40515, ])
dim(hist10)
hist0710 <- rbind(hist07, hist10)

# plot
set_graph_pars(ptype = "panel1")
xlab2 <- expression(paste("Size (", cm^2, ")"))
hist(hist0710$area, breaks = seq(0.0, 1.3, 0.02), freq = FALSE, 
	xlab = xlab2, ylab="Probability density", col = "gray90", main = "", 
	ylim = c(0, 2.8), xlim = c(0, 1.3), border = "gray90") 

hist(hist0710$area, breaks = seq(0, 1.1, 0.05), freq = FALSE, 
	xlab = xlab2, ylab="Probability density", col = "gray", main = "", 
	cex.lab = 1.4, cex.axis = 1.2, ylim = c(0, 2.2), xlim = c(0, 1.1)) 
	
############################################
### Get appropriate data for size and survival 
#at t+3 as a function of size at time t 

dim(d)
unique(d$date.no)
dT03 <- droplevels(d[d$date.no == 39426 | d$date.no == 40515, ])
summary(dT03)

dT0 <- droplevels(dT03[dT03$date.no == 39426, ])
dT3 <- droplevels(dT03[dT03$date.no == 40515, ])
summary(dT3)

dim(dT03)
unique(dT03$coral.id) 
dim(dT0)
unique(dT0$coral.id) 
dim(dT3)
names(dT3)
unique(dT3$coral.id) 

# reshape the dataset to get survival as a function of area for each time point in columns
dW <- dcast(dT03, quad + coral.id ~ date.no, value.var = "area")
head(dW)

names(dW)[3:4] <- c("T0", "T3")
unique(dW$coral.id) 
dim(dW) # length matches the unique set of coral.IDs, that is good
summary(dW)

# remove columns with NA for T0 (because I only care about corals that I measured in 2007, not ones that showed up later either through recruitment or because the framer shifted)

dW2 <- droplevels(dW[complete.cases(dW$T0), ])
dim(dW2)
summary(dW2)

# now I would like to add code column for each time point
names(dW2)
names(dT0) # want to add sizeOK for T0
dW3 <- merge(dW2, dT0[, which(names(dT0) %in% c("coral.id", "code"))], by = "coral.id")
dim(dW3)
names(dW3)[5] <- "T0code"
head(dW3)

dW4 <- merge(dW3, dT3[, which(names(dT3) %in% c("coral.id", "code"))], by.x = "coral.id")
dim(dW4)
head(dW4)
names(dW4)[6] <- "T3code"
head(dW4)
summary(dW4)

# Remove rows that have angle, nv, and algae in the T0code column
unique(dW4$T0code)
dW4
summary(dW4) # should get a dataframe that is 119 rows long
dim(dW4)

dW5 <- droplevels(dW4[dW4$T0code != "angle" & 
             dW4$T0code != "nv" & 
             dW4$T0code != "algae" , ])

dim(dW5)
unique(dW5$T0code)

### Create survival dataset
# remove rows that are nv and algae for T3code column
unique(dW5$T3code)
survDat <- droplevels(dW5[dW5$T3code != "nv" & 
            dW5$T3code != "algae", ])
dim(survDat)

## NEED SURVIVAL FUNCTION
survDat$T3surv <- survF(survDat$T3code)
survDat
names(survDat)
names(survDat) <- c("coral.id", "quad", "size", "sizeNext", "T0code",
                    "T3code", "surv")
sum(survDat$surv) # 85 out of 104 surivived

### Create growth dataset
# remove rows that are nv and algae for T3code column
unique(dW5$T3code)
growthDat <- droplevels(dW5[dW5$T3code != "nv" & 
                            dW5$T3code != "algae" &
                            dW5$T3code != "angle" &
                            dW5$T3code != "dead" , ])
dim(growthDat)
head(growthDat)
names(growthDat) <- c("coral.id", "quad", "size", 
	"sizeNext", "T0code", "T3code")
growthDat$delta <- with(growthDat, sizeNext - size)

# check deltas
plot(delta ~ size, data = growthDat)
growthDat[growthDat$delta > 0.4, ]
growthDat[growthDat$delta < -0.1, ]

#########################################
# Need to estimate the relationship between coral size and embryo #
# use data in Fadlallah Figure 5 and Table 6
# Figure 5 gives # of embryos per female, the raw data
# Table 6 assumes a 50:50 ratio of males:females, and thus divides the
# number of embryos by 2 to get mx per coral

# upload empirical data on area volume relationships
# because I have calyx area, but Fadlallahs data are volume
avDat <- read.csv("./data/bael area volume.csv", 
                  header=TRUE, na.strings="NA")

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

# Upload bael life table
ltDat <- read.csv("./data/bael life table.csv", 
                  header=TRUE, na.strings="NA")
head(ltDat)
summary(ltDat)
# Calculate area for corals in life table based on regression
ltDat$areaR <- (exp(log(ltDat$Vol.mid)*avReg$coefficients[2] + 
               avReg$coefficients[1]))/100
ltDat$areaR

# biggest coral is 700 mm3  
areaCalcF(700)   # 0.77 cm2
  
# Compare to calculations in excel
head(ltDat)
plot(areaR ~ areaREGR, ltDat)
# Plot calculated area vs regressed area
plot(areaCALC ~ areaREGR, ltDat)
abline(a = 0, b = 1, col = "red") # calculated volume underestimates 
# regressed volume

# Estimate the relationship between areaR and mx.coral
head(ltDat)
plot(mx.coral ~ Vol.mid, data = ltDat)
plot(mx.coral ~ areaR, data = ltDat)
ltDat

# THIS ASSUMES THAT WA POPULATION FUNCTION IS THE SAME AS THE CA POPULATION - BUT I WILL CORRECT THIS
mxReg <- lm(mx.coral ~ areaR, data = ltDat[ltDat$Vol.mid > 75, ])
summary(mxReg)

plot(mx.coral ~ areaR, data = ltDat, xlim = c(0, 1.5), ylim = c(0, 50),
     xlab = "Size (cm2)", ylab = "No. of embryos", las = 1)
abline(mxReg, col = "darkgray")
xx <- seq(0, 1.2, by = 0.1)
mxReg$coefficients

xIntCA <- -mxReg$coefficients[1]/mxReg$coefficients[2]
xIntWA <- 16.28972/mxReg$coefficients[2]
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

#########################################
# Need to calculate establishment probability 
# 1. Calculate the # of embryos produced in quadrat x in time 0
# 2. Calculate the # of recruits that survived to time 3
# 3. Establishment probability = # of living recruits/ # of embryos in T0

# USE WA EMBRYO FUNCTION
#dat$embryo.no <- embryoFca(dat$area)
dat$embryo.no <- embryoFwa(dat$area)

# Sum the embryos per quadrat for each time point
tapply(dat$embryo.no, list(dat$quad, dat$date.no), FUN=sum, na.rm=TRUE)
############################################
# Calculate density 
density.q <- tapply(dat$area, list(dat$quad, dat$date.no), FUN=densityF)
density.q
#write.csv(data.frame(density.q), "./output/baelDensity.csv")

# 36 x 25 cm quadrat
quadArea <- 36*25*0.0001
adultDens <- density.q/quadArea
adultDens
adultDens2 <- data.frame(adultDens)
adultDens2
str(adultDens2)

# Calculate recruit density
recruit.q <- tapply(dat$code, list(dat$quad, dat$date.no), FUN=recruitF)
colnames(recruit.q) <- c("2007", "2008", "2009", "2010")
recruit.q
dim(recruit.q)

recruitDens <- recruit.q/quadArea
recruitDens

# Calculate recruitment success as # of recruits in year t+1 divided by # of embryos in year t
embryo.q <- tapply(dat$embryo.no, list(dat$quad, dat$date.no), FUN=sum, na.rm=TRUE)
colnames(embryo.q) <- c("2007", "2008", "2009", "2010")
embryo.q
dim(embryo.q)

rq <- recruit.q[,-1] # remove the first column, because it does not have recruitment data (for year 2007)
eq <- embryo.q[,-4] # remove the last column, because the embryo no. in 2010 is irrelevant
dq <- density.q[,-4] # remove the last column, because the density in 2010 is irrelevant
dq
eq

eqSum <- rowSums(eq)
eqSum
eqMult3 <- 3*eq[, 1] # 3 years worth of embryos
eqMult2 <- 2*eq[, 1] # 2 years worth of embryos
eqMult1 <- eq[, 1] # 1 years worth of embryos 2007
plot(eqSum ~ eqMult3)
abline(a = 0, b = 1, col = "red")

# calculate recruitment success
recruit.prop <- (rq/eq)*100
recruit.prop
recruit.prop[is.infinite(recruit.prop)] <- NaN
mean(recruit.prop, na.rm = TRUE)
sd(recruit.prop, na.rm = TRUE)

annualRecProb <- colMeans(recruit.prop, na.rm = TRUE)
# another way
apply(as.matrix(recruit.prop), 2, FUN = function(x) mean(x, na.rm = TRUE))
apply(as.matrix(recruit.prop), 2, FUN = function(x) sd(x, na.rm = TRUE))
annualRecProb

# OVERALL MEAN AND SD
mean(annualRecProb)
sd(annualRecProb)

str(recruit.prop)
hist(recruit.prop)

eq

# plot
par(mfrow=c(2,1))
plot(jitter(rq) ~ eq, las=1, 
     xlab = "Estimated no. of embryos per quadrat at time t", 
     ylab = "No. of recruits at time t + 1")

plot(jitter(recruit.prop) ~ eq, las=1, 
     xlab = "Estimated no. of embryos per quadrat at time t", 
     ylab = "Recruitment success (%) at time t + 1")

plot(recruit.prop ~ eq)
plot(recruit.prop ~ dq)

estabProp <- as.data.frame(recruit.prop)
estabProp

establishment.prob <- as.data.frame(recruit.prop)
head(establishment.prob)

############################################
############################################
# Get distribution of recruit size
head(d)
summary(d)
recruitDat <- d[d$recruit == 1,]
dim(recruitDat)
summary(recruitDat)
# 52 total recruits over the 3 year period

recruit.size.mean <- mean(recruitDat$area)
recruit.size.sd <- sd(recruitDat$area)

### Now do this for recruits at the end of T3
head(recruitDat)
# Compare recruitDat with data for T3
head(dT3)
dim(dT3)

# create dataframe with code coral.id for recruitDat
recruitDat2 <- droplevels(recruitDat[, which(names(recruitDat) %in% 
                                c("coral.id", "code"))])
dim(recruitDat2)
names(recruitDat2) <- c("coral.id", "recruitCode")
head(recruitDat2)
summary(recruitDat2)

dT3b <- merge(dT3, recruitDat2, all.x = TRUE)
#dT3b <- droplevels(dT3b[complete.cases(dT3b$code), ])

dim(dT3b)
head(dT3b)
summary(dT3b) # 52 recruits, and 30 had died

# how many of the 52 recruits were still alive at T3?
dT3c <- droplevels(dT3b[dT3b$recruitCode == "recruit", ]) # doesn't work?

dT3c <- droplevels(dT3b[dT3b$code == "dead" &
                          dT3b$recruitCode == "recruit", ]) # also doesn't work?

dT3c <- droplevels(dT3b[dT3b$code == "dead", ]) # works

dim(dT3c)
summary(dT3c) # so of the 30 that died, 3 were recruits
head(dT3c)

# after 3 years 49/52 recruits survived
# Get the size distribution of these 49 recruits in time 3
names(dT3b)
dT3d <- droplevels(dT3b[complete.cases(dT3b$recruitCode), ])
dim(dT3d)
summary(dT3d)

# Have to drop those covered by algae, at an angle, or dead guys
dT3e <- droplevels(dT3d[dT3d$code != "algae" &
                          dT3d$code != "angle" & 
                   dT3d$code != "dead",  ])

summary(dT3e)
hist(dT3e$area)

#corals bigger than 0.3
dT3e[dT3e$area > 0.2, ]

recruit.size.mean2 <- mean(dT3e$area)
recruit.size.sd2 <- sd(dT3e$area)

#write.csv(dT3e, "./output/recruitDatPresent.csv")

eqSum # sum of embryos from 2007 - 2010
eqMult3 # the total # of embryos per quadrat in 2007, multiplied by 3
# Now I need the total number of recruits per quadrat at T3

# Calculate recruit density in T3
recruit.qT3 <- tapply(dT3b$recruitCode, list(dT3$quad, dT3$date.no), FUN=recruitF)
recruit.qT3
str(recruit.qT3)
sum(recruit.qT3)


names(adultDens2) <- c("dens07", "dens08", "dens09", "dens10")

adultDens3 <- adultDens2

recruitDensT3 <- as.numeric(recruit.qT3)/quadArea
str(recruitDensT3)
adultDens3$recruitDens <- recruitDensT3

adultDens3
ggplot(adultDens3, aes(dens07, recruitDens)) +
	geom_point() + geom_smooth(method = "lm")
(summary(lm(recruitDens ~ dens07, data = adultDens3)))

#write.csv(adultDens3, "./output/adultRecDensPresent_150411.csv")
# plot density over time for 2007-2010
densDF <- melt(adultDens2, measure.vars = c("dens07", "dens08", "dens09", "dens10"), value.name = "density")
densDF$year <- c(rep(2007, 24), rep(2008, 24), rep(2009, 24), rep(2010, 24))
densDF

eq; eqSum
sum(eqSum > 0)

# based on 1 years worth of embryos
recProb1 <- (recruit.qT3/(eq[, 1]))*100 
recProb1
recProb1[is.infinite(recProb1)] <- NaN
recProb1

# based on sum of 3 years worth of embryos
recProbSum <- (recruit.qT3/eqSum) * 100; 
recProbSum
recProb2 <- recProbSum
recProb2[is.infinite(recProb2)] <- NaN
recProb2

# based on year1 multipled by 3
eq[, 1] * 3
recProbMult <- (recruit.qT3/(eq[,1]*3))*100
recProbMult
recProb3 <- recProbMult
recProb3[is.infinite(recProb3)] <- NaN
recProb3

# compare to annual probabilities
mean(recProb1, na.rm = TRUE) # 
sd(recProb1, na.rm = TRUE) # 

mean(recProb2, na.rm = TRUE) # 
sd(recProb2, na.rm = TRUE) # 

mean(recProb3, na.rm = TRUE) # 
sd(recProb3, na.rm = TRUE) #

plot(recProbSum ~ recProbMult)

set_graph_pars(ptype = "panel2")
hist(recProb1, breaks = 15, las=1, xlab = "Recruitment success (%)", main = NULL)

stockRecDat <- cbind(eq[, 1], recruit.qT3, recProb2)
srDat <- as.data.frame(stockRecDat)
names(srDat) <- c("embryos", "recruits", "recProb2")
srDat

# Stock-recruitment regressions
sr1 <- lm(recruits ~ embryos, data = srDat)
sr2 <- lm(recProb2 ~ embryos, data = srDat)
# using # of embryos from 2007 - 2010

# plot
set_graph_pars(ptype = "panel2")

plot(jitter(recruits) ~ embryos, las=1, 
     xlab = "Estimated no. of embryos\nper quadrat 2007-2009", 
     ylab = "No. of recruits in 2010", data = srDat)
abline(sr1, col = "black", lwd=2)

plot(jitter(recProb2) ~ embryos, las=1, 
     xlab = "Estimated no. of embryos\nper quadrat (2007-2010)", 
     ylab = "Recruitment success (%) in 2010", data = srDat)
abline(sr2, col = "black", lwd=2)

recProb2
hist(recProb2)
sd(recProb2/100, na.rm = TRUE)
#######################################################
########################################################
########################################################

###1.4.2 Build regressions for vital rate functions###
# set up a data frame of the model parameters
params <- data.frame(
  surv.int=NA,			# Intercept from logistic regression of survival
  surv.slope=NA,	# Slope from logistic regression of survival
  mature.size=NA,  # Size at maturity
  embryo.int=NA,			# Intercept from linear regression of embryo number
  embryo.slope=NA,		# Slope from linear regression of embryo number	
  embryo.sd=NA,   # Residual sd from the linear regression of embryo number
  recruit.size.mean=NA, # Mean recruit size
  recruit.size.sd=NA,   # Standard deviation of recruit size
  estab.prob.mean=NA, # Mean of establishment probability
  estab.prob.sd=NA,      # SD of establishment probability
  growth.int=NA,  	# Intercept from linear regression of growth
  growth.slope=NA,	# Slope from linear regression of growth
  growth.sd=NA		# Residual sd from the linear regression of growth    
)

# Next we build each of the vital rate regressions and store the coefficients

# 1. survival: logistic regression
surv.reg <- glm(surv ~ size, data = survDat, family = binomial())
summary(surv.reg)
length(survDat$size)

curve(plogis(1.14 + 0.98*x), ylim = c(0,1), lwd = 2)

# summary of logistic regression meaning
# the interpretation of the slope term is as follows
# for every one unit change in size (cm2), the log(odds survival) differ by 0.98
# therefore the odds of survival is exp(0.98)
exp(0.98)
exp(1.2*params$surv.slope + params$surv.int)

# use glmer to test significance of size
head(survDat)
survDat$tran <- substr(survDat$quad, 1, 3)
survMod0 <- glmer(surv ~ size + (1 | quad) + (1 | tran), 
		data = survDat, family = binomial)
survMod <- glmer(surv ~ size + (1 | quad), 
		data = survDat, family = binomial)
AIC(survMod0, survMod) # don't need the transect term
summary(survMod)		
summary(survMod)[[10]][2]
summary(surv.reg)

params$surv.int <- summary(survMod)[[10]][1]
params$surv.slope <- summary(survMod)[[10]][2]
curve(plogis(1.71 + 0.45*x), ylim = c(0,1), lwd = 2)
params

# 2. growth: linear regression
growth.reg <- lm(sizeNext ~ size, data = growthDat)
summary(growth.reg)
plot(resid(growth.reg) ~ growthDat$size, las=1, 
     cex=1.2, lwd=1.3, cex.axis=1, cex.lab=1.2, 
     ylab="Residuals", main = "Growth rate regression")
abline(a = 0, b = 0, lty = 2, lwd = 2, col = "darkgray")
# the fact that the residuals do not increase with size obviate the need to log transform size
head(growthDat)

# compare linear growth to power law growth
gReg1 <- nls(sizeNext ~ a + b*size, start = list(a = 0.2, b = 0.75), data = growthDat)
gReg1
gReg2 <- nls(sizeNext ~ a * size^b, start = list(a = 0.5, b = 0.5), data = growthDat)
gReg2

AIC(gReg1, gReg2) # use linear regression

# use glmer to test signifiance of size
head(growthDat)
growthDat$tran <- substr(growthDat$quad, 1, 3)
growthMod0 <- lmer(sizeNext ~ size + (1 | quad) + (1|tran), 
		data = growthDat)
growthMod <- lmer(sizeNext ~ size + (1 | quad) , 
		data = growthDat)
anova(growthMod0, growthMod)		
summary(growthMod)
summary(growthMod0)
summary(growth.reg)

params$growth.int <- summary(growthMod)[[10]][1]
params$growth.slope <- summary(growthMod)[[10]][2]
params$growth.sd <- sd(resid(growthMod))

# 2.1 log growth plots
lnGrowthReg <- lm(log(sizeNext) ~ log(size), data = growthDat)
summary(lnGrowthReg)

# 3. embryos: linear regression; CORRECTED for WA
embryo.reg <- lm(mx.coral ~ areaRcorr, data = ltDat[ltDat$Vol.mid > 75, ])
summary(embryo.reg)
params$embryo.int=coefficients(embryo.reg)[1]
params$embryo.slope=coefficients(embryo.reg)[2]
params$embryo.sd=sd(resid(embryo.reg))
params

embryo.reg$coefficients
params$mature.size <- -embryo.reg$coefficients[1]/embryo.reg$coefficients[2]
params

# 4. size distribution of recruits
params$recruit.size.mean <- recruit.size.mean2
params$recruit.size.sd <- recruit.size.sd2
params

# 5. establishment probability - using the sum of three years of embryos
params$estab.prob.mean <- mean(recProb2/100, na.rm = TRUE) 
params$estab.prob.sd <- sd(recProb2/100, na.rm = TRUE) # 
params

paramsWA <- params

################################################
################################################
### IPM diagnostics for growth
################################################
################################################
growthMod
xlab1 <- expression(paste("Size at time t (", cm^2, ")"))

# Residuals v fitted 
set_graph_pars(ptype = "panel1")
par(cex = 1.4)
#plot(resid(growthMod) ~ fitted(growthMod), xlab = "Fitted", ylab = "Residuals", col = "black", las = 1, main = "Growth model for IPM")

# Residuals v observed 
plot(resid(growthMod) ~ growthDat$size, xlab = xlab1, ylab = "Residuals", col = "black", las = 1, main = "Growth model for IPM", lwd = 2)
abline(h=0, col = "darkgray", lty = 2, lwd = 2)

################################################
################################################
### GLM for survival
################################################
################################################
head(survDat)
survMod <- glmer(surv ~ size + (1 | quad) , 
		data = survDat, family = binomial)
survNull <- glmer(surv ~ 1 + (1 | quad) , 
		data = survDat, family = binomial)		
summary(survMod)[[10]][2]
summary(survMod)
summary(surv.reg)
AIC(survMod) - AIC(survNull)

#params$surv.int <- summary(survMod)[[10]][1]
#params$surv.slope <- summary(survMod)[[10]][2]
params
