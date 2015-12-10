#################################################
# Author: Robin Elahi
# Date: 150828

# Coral survival
# Figure 4
#################################################

library(lme4)
library(ggplot2)
library(AICcmodavg)

# rm(list=ls(all=TRUE)) # removes all previous material from R's memory

dat <- read.csv("./data/bael_survivalData.csv", header=TRUE, na.strings="NA")

source("./R/graphicalParams.R")
source("./R/multiplotF.R")

dat
names(dat)
summary(dat)
dim(dat)

dat$ini.areaLN <- log(dat$ini.area)
dat$fin.areaLN <- log(dat$fin.area)

##########################################################
##########################################################
# LMER; TEST SC DATA ONLY (TRUNCATED)
##########################################################
##########################################################

# lme4
# test whether survival is differrent among eras
lmerDat <- dat[dat$ini.area < 0.96, ]
#lmerDat <- datSC
dim(lmerDat)

unique(lmerDat$quad)

Cand.mod <- list()
Cand.mod[[1]] <- glmer(survival ~ ini.area*time + (1|quad),
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

library(sjPlot)

sjp.lmer(globalMod)
sjp.lmer(globalMod, type = "fe.cor")
sjp.lmer(globalMod, type = "re.qq")

sjp.lmer(globalMod, fade.ns = TRUE, free.scale = TRUE, 
	geom.colors = c(1,1), showValueLabels = FALSE, 
	sort.coef = TRUE) 

sjp.lmer(globalMod, fade.ns = TRUE, free.scale = TRUE, 
	geom.colors = c(1,1), showValueLabels = FALSE, 
	type = "fe") 	
					
##########################################################
##########################################################
# COMPARE PARAMETERS
##########################################################
##########################################################
glm1 <- glm(survival ~ ini.area*time, data = lmerDat, family = binomial)
summary(glm1)

head(dat)
datPres <- lmerDat[lmerDat$time == "present", ]
datPast <- lmerDat[lmerDat$time == "past", ]

pastMod <- glmer(survival ~ ini.area + (1|quad),
							data = datPast, family = binomial)		

presMod <- glmer(survival ~ ini.area + (1|quad),
							data = datPres, family = binomial)	

pastMod2 <- glm(survival ~ ini.area,
							data = datPast, family = binomial)		

presMod2 <- glm(survival ~ ini.area,
							data = datPres, family = binomial)	

summary(pastMod)
summary(pastMod2)

summary(presMod)
summary(presMod2)

curve(plogis(0.5186 + 3.4766*x), ylim = c(0,1), lwd = 2, col = "darkgray")
curve(plogis(1.7106 + 0.4489*x), lwd = 2, col = 1, add = TRUE)

curve(plogis(0.5186 + 0.4489*x), ylim = c(0,1), lwd = 2, col = "red", add = TRUE)
curve(plogis(1.7106 + 3.4766*x), ylim = c(0,1), lwd = 2, col = "blue", add = TRUE)


ggplot(data.frame(x = c(0, 0.96)), aes(x)) + 
	stat_function(fun = function(x)x^2, geom = "line")

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

##########################################################
##########################################################
# FIGURE - SURVIVAL SCALING BY ERA (TRUNCATED)
##########################################################
##########################################################
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


