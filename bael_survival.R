#################################################
# Author: Robin Elahi
# Date: 151209
# Coral survival
#################################################

# rm(list=ls(all=TRUE)) 

##### LOAD PACKAGES, DATA #####
library(lme4)
library(ggplot2)
theme_set(theme_classic(base_size = 16))
library(AICcmodavg)
library(dplyr)

dat <- read.csv("./data/bael_survivalData.csv", header=TRUE, na.strings="NA")

source("./R/graphicalParams.R")
source("./R/multiplotF.R")

dat$ini.areaLN <- log(dat$ini.area)
dat$fin.areaLN <- log(dat$fin.area)

##### USE GLMER TO TEST THE EFFECT OF ERA ON SURVIVAL #####
# test whether survival is different among eras
lmerDat <- dat[dat$ini.area <= 0.96, ] # for some reason does not converge for area < 1

unique(lmerDat$quad)

Cand.mod <- list()
Cand.mod[[1]] <- glmer(survival ~ ini.area * time + (1|quad),
							data = lmerDat, family = binomial)

Cand.mod[[2]] <- glmer(survival ~ ini.area + time + (1|quad),
							data = lmerDat, family = binomial)							
							
Cand.mod[[3]] <- glmer(survival ~ ini.area + (1|quad),
							data = lmerDat, family = binomial)										

Cand.mod[[4]] <- glmer(survival ~ time + (1|quad),
							data = lmerDat, family = binomial)			
							
Cand.mod[[5]] <- glmer(survival ~ 1 + (1|quad),
							data = lmerDat, family = binomial)										
													
#generate AICc table with names

mod_numbers <- paste("Cand.mod", 1:length(Cand.mod), sep=" ")	
mod_text <- c("Era x Size", "Era + Size", "Size", "Era", "Null model")     

mod.aicctab <- aictab(cand.set= Cand.mod, modnames= mod_text, 
                      sort=TRUE, second.ord=TRUE) # second.ord =TRUE means AICc is used 
print(mod.aicctab, digits=2, LL=TRUE)
write.csv(mod.aicctab, "./output/survivalAIC.csv")

summary(Cand.mod[[1]])

# check normality and homogeneity of variances
globalMod <- Cand.mod[[1]] # 
																	
# check normality and homogeneity of variances
par(mfrow = c(1,2))
qqnorm(resid(globalMod))
qqline(resid(globalMod))
plot(resid(globalMod) ~ fitted(globalMod)); abline(h=0)

summary(globalMod)
fitted(globalMod)

# tables of estimates with 95% CI
se <- sqrt(diag(vcov(globalMod)))
tab <- cbind(Est = fixef(globalMod), 
					LL = fixef(globalMod)  - 1.96 * se, 
					UL = fixef(globalMod) + 1.96 * se)
tab					

# odds ratios, by exponentiation
exp(tab)

# caterpillar plot
ranef(globalMod, which = "quad", condVar = TRUE)

##### GET MODEL PARAMETERS BY ERA #####					
datPres <- lmerDat[lmerDat$time == "present", ]
datPast <- lmerDat[lmerDat$time == "past", ]

pastMod <- glmer(survival ~ ini.area + (1|quad),
							data = datPast, family = binomial)		

presMod <- glmer(survival ~ ini.area + (1|quad),
							data = datPres, family = binomial)	

summary(pastMod)
summary(presMod)

curve(plogis(0.5186 + 3.4766*x), ylim = c(0,1), lwd = 2, col = "darkgray")
curve(plogis(1.7106 + 0.4489*x), lwd = 2, col = 1, add = TRUE)

curve(plogis(0.5186 + 0.4489*x), ylim = c(0,1), lwd = 2, col = "red", add = TRUE)
curve(plogis(1.7106 + 3.4766*x), ylim = c(0,1), lwd = 2, col = "blue", add = TRUE)

survTrendPast <- ggplot(data.frame(x = c(0, 0.95)), aes(x)) + 
	stat_function(fun = function(x) plogis(x*3.4776 + 0.5186), 
	geom = "line", size = 1.5, color = "darkgray") +
	scale_y_continuous(limits = c(0, 1))

survTrendPres <- ggplot(data.frame(x = c(0, 0.95)), aes(x)) + 
	stat_function(fun = function(x) plogis(x*0.4489 + 1.7106), 
	geom = "line", size = 1.5, color = "black") +
	scale_y_continuous(limits = c(0, 1))

survTrendPast
survTrendPres

survTrendPast2 <- stat_function(fun = function(x) 
								plogis(x*3.4776 + 0.5186), 
								geom = "line", size = 1, color = "darkgray") 

survTrendPres2 <- stat_function(fun = function(x) 
								plogis(x*0.4489 + 1.7106), 
								geom = "line", size = 1, color = "black")

survTrendPast + survTrendPres2

##### FINAL FIGURE - SURVIVAL SCALING BY ERA #####
ylab1 <- "Survival at time t+3"
xlab1 <- expression(paste("Size at time t (", cm^2, ")"))
ULClabel <- theme(plot.title = element_text(hjust = -0.07, vjust = 0, size = rel(1.2)))
# ULClabel <- theme(plot.title = element_text(hjust = -0.15, vjust = 0, size = rel(1.2)))

theme_set(theme_classic(base_size = 12))

surv1 <- ggplot(lmerDat, aes(ini.area, survival, 
	color = time, shape = time)) +
	ylab(ylab1) + xlab(xlab1) + 
	theme(legend.justification = c(1, 0), legend.position = c(1, 0.2)) +
	theme(legend.title = element_blank()) + 
	geom_point(size = 4, alpha = 0.8, 
		position = position_jitter(h = 0.05)) +
	scale_colour_manual(breaks = c("past", "present"), 
		values = c("darkgray", "black")) +
		scale_shape_manual(breaks = c("past", "present"), 
		values = c(18, 20))

survPlot <- surv1 + labs(title = "A") + ULClabel + 
	survTrendPast2 + survTrendPres2 + theme(legend.position = "none")
	
survPlot

# RENAME lmerDat
dat_survival <- lmerDat
