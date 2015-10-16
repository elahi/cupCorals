#################################################
# Author: Robin Elahi
# Date: 151016

# Coral densities for 1969 and 2007
# Figure 1
#################################################

library(ggplot2)
theme_set(theme_classic(base_size = 12))

library(reshape2)
library(nlme)
library(dplyr)

# rm(list=ls(all=TRUE)) # removes all previous material from R's memory

source("./R/graphicalParams.R")
source("./R/summarizeData.R")

###########################################	
### CALCULATE

dat <- read.csv("./data/bael_densityModWide.csv", header=TRUE, na.strings="NA")
names(dat)
head(dat)

lambdaDF <- with(dat, dens10/dens07)
lambdaDF

meanLambda <- mean(lambdaDF, na.rm = TRUE)
sdLambda <- sd(lambdaDF, na.rm = TRUE)
meanLambda; sdLambda

densDF <- melt(dat, measure.vars = 
                 c("dens07", "dens08", "dens09", "dens10"), 
               value.name = "density", id.vars = c("quad", "tran"))
densDF$year <- c(rep(2007, 24), rep(2008, 24), 
                 rep(2009, 24), rep(2010, 24))
densDF
densDF$logDens <- log(densDF$density + 1)
hist(densDF$density)
hist(densDF$logDens)
# write.csv(densDF, "modernDensity.csv")

# does density increase over time?
# use nlme
# random effects
rand1 <- ~ 1 | quad
rand2 <- ~ 1 | tran

boxplot(density ~ year, data = densDF)
boxplot(logDens ~ year, data = densDF)

lme1 <- lme(fixed = density ~ year, 
	data = densDF, method = "ML", 
	random = list(rand1, rand2))

# now add autocorrelation term
lme2 <- update(lme1, correlation = corAR1())
lme3 <- update(lme1, correlation = corCAR1())

AIC(lme1, lme2, lme3)
# need to include autocorrelation
plot(lme2)
summary(lme2) # year NS

nullModel <- update(lme2, . ~ 1)
summary(nullModel)
nullModel

AIC(lme2, nullModel)
anova(lme2, nullModel) # marginally significant

quadArea <- 36*25*0.0001
densDF$densM <- densDF$density/quadArea

summDat <- summarySE(densDF, measurevar = "densM", groupvars = "year")
summDat

ggplot(summDat, aes(x = year, y = densM)) + 
	theme_classic(base_size = 22) + xlab("Year") +
	ylab("densM per m2") + 
	geom_point(size = 5) + 
	geom_errorbar(aes(ymin = densM - ci, ymax = densM + ci), 
		width = 0.02) +
	scale_y_continuous(limits = c(0, 200))

densDF	

ggplot(densDF, aes(x = year, y = densM, color = quad)) + 
	theme_classic(base_size = 22) + xlab("Year") +
	ylab("densM per m2") + 
	geom_point(size = 3) + 
	geom_line() +
	theme(legend.position = "none") 

###########################################	
# load new file with historic data also

dat <- read.csv("./data/bael_densityMaster.csv", header=TRUE, na.strings="NA")
names(dat)
head(dat2)

# to compare densities of old and new:
# remove D3
dat2 <- dat[dat$quad != "D3", ]
dim(dat2)
dim(dat)

boxplot(densM ~ year, data = dat2)

ggplot(data = dat2, aes(as.factor(year), densM)) + geom_violin()

summ2 <- summarySE(dat2, measurevar = "densM", 
                   groupvars = c("tran", "year"))
summ2

boxplot(densM ~ year, data = summ2)
ggplot(data = summ2, aes(as.factor(year), densM)) + geom_violin()

summ3 <- summarySE(summ2, measurevar = "densM", groupvars = c("year"))
summ3$era <- c(rep("1969-1972", 2), rep("2007-2010", 4))
summ3

ggplot(summ3, aes(x = year, y = densM, color = era)) + 
	theme_classic(base_size = 22) + xlab("Year") +
	ylab("densM per m2") + 
	geom_point(size = 3) + 
	theme(legend.position = "none") + 
	geom_errorbar(aes(ymin = densM - sd, ymax = densM + sd), 
		width = 0.02) 

#
summ2$era <- c(rep("1969-1972", 8), rep("2007-2010", 24))

ggplot(summ2, aes(x = year, y = densM, color = tran)) + 
	theme_classic(base_size = 22) + xlab("Year") +
	ylab("densM per m2") + 
	geom_point(size = 3) + 
	theme(legend.position = "none") +
	geom_line() 

facet2 <- facet_wrap(~ era, scales = "free_x")

library(grid)
ggplot(summ2, aes(x = year, y = densM, color = tran)) + 
	theme_classic(base_size = 22) + xlab("Year") +
	ylab("densM per m2") + 
	geom_point(size = 3, color = "black") + 
	theme(legend.position = "none") +
	geom_line() + facet2 +
	theme(panel.margin = unit(3, "lines")) 
	
	+ 
	geom_errorbar(aes(ymin = densM - se, ymax = densM + se), 
		width = 0.02) 
	
	
# model for summ2
summ2
	
# what is lambda?
# calculate as Nt+1/Nt
#
histPre <- summ3[summ3$year == 1969, ]$densM
histPost <- summ3[summ3$year == 1972, ]$densM
modPre <- summ3[summ3$year == 2007, ]$densM
modPost <- summ3[summ3$year == 2010, ]$densM

histLambda <- histPost/histPre
modLambda <- modPost/modPre

histLambda
modLambda

summ4 <- summ3[c(1:3, 6), ]
summ4$census <- c("t", "t+3", "t", "t+3")
summ4

ylab1 <- expression(paste("Density (no. ", m^-2, ")"))
labelHist <- expression(paste(lambda, " = 1.09"))
labelMod <- expression(paste(lambda, " = 1.25"))

labelHist <- ("lambda == 1.09")
labelMod <- ("lambda == 1.25")

plot(1, 1, xlab = labelMod)

ggplot(summ4, aes(x = census, y = densM, shape = era, color = era)) + 
  xlab("Year") + ylab(ylab1) + 
	geom_point(size = 6, 
		position = position_dodge(0.3)) + 
	geom_errorbar(aes(ymin = densM - sd, ymax = densM + sd), 
		width = 0.0, position = position_dodge(0.3)) +
	scale_colour_manual(breaks = c("1969-1972", "2007-2010"), 
		values = c("darkgray", "black"))	 +
	scale_shape_manual(breaks = c("1969-1972", "2007-2010"), 
		values = c(18, 20)) +
	theme(legend.justification = c(1, 0), legend.position = c(1.05, 0.8)) +
	theme(legend.title = element_blank()) + 
	theme(legend.key.size = unit(0.8, "cm")) +
	annotate("text", x = 2.35, y = 107, label = labelMod, 
		parse = T, fontface = "bold", size = 8) +
	annotate("text", x = 2.2, y = 170, label = labelHist, 
		color = "darkgray", parse = T, size = 8)
	
summ2	
	
	