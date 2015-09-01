#################################################
# Author: Robin Elahi
# Date: 150828

# Coral growth
# Figure 4
#################################################

library(lme4)
library(ggplot2)
library(AICcmodavg)

#rm(list=ls(all=TRUE)) # removes all previous material from R's memory

dat <- read.csv("./data/bael_growthData.csv", header=TRUE, na.strings="NA")
source("./R/graphicalParams.R")

dat
names(dat)
summary(dat)
dim(dat)

dat$ini.areaLN <- log(dat$ini.area)
dat$fin.areaLN <- log(dat$fin.area)

datSC <- dat
# Create dataset for SC present
datSCpresent <- datSC[datSC$time == "present", ]
# Create dataset for SC past
datSCpast <- datSC[datSC$time == "past", ]

dim(datSC)
dim(datSCpresent)
dim(datSCpast)
range(datSCpast$ini.area) # 1.669
range(datSCpresent$ini.area) # 0.828

# Truncate initial size range so that era's match at upper limit
presMaxSC <- max(datSCpresent[datSCpresent$time == "present", ]$ini.area)

set_graph_pars(ptype = "panel1")
plot(fin.area ~ ini.area, data = datSCpast)
#abline(v = presMaxSC)
points(fin.area ~ ini.area, data = datSCpresent, col = "red")
abline(v = 0.96, col = "black", lty = 3)
abline(a = 0, b = 1, col = "darkgray", lty = 2, lwd = 2)

datSCTrunc <- datSC[datSC$ini.area <= 0.96, ]

##########################################################
##########################################################
# LMER; TEST SC DATA ONLY (TRUNCATED)
##########################################################
##########################################################

lmerDat <- datSCTrunc
dim(lmerDat)
names(lmerDat)

# SIZE: LINEAR
# use varying slopes by quadrat? no - doesn't add anything
mod1 <- lmer(fin.area ~ ini.area*time +  (ini.area|quad), 
                      REML = FALSE, data=lmerDat)
summary(mod1)
mod2 <- lmer(fin.area ~ ini.area*time + (1|quad), REML = FALSE, data=lmerDat)    
summary(mod2)
AIC(mod1, mod2) # mod 1 is better; including 
lm1 <- lm(fin.area ~ ini.area*time, data=lmerDat)
summary(lm1)

Cand.mod <- list()                                        
Cand.mod[[1]] <- lmer(fin.area ~ ini.area*time + (1|quad), 
                      REML = FALSE, data=lmerDat)                     
Cand.mod[[2]] <- lmer(fin.area ~ ini.area + time + (1|quad), 
                      REML = FALSE, data=lmerDat)
Cand.mod[[3]] <- lmer(fin.area ~ ini.area  + (1|quad), 
                      REML = FALSE, data=lmerDat)
Cand.mod[[4]] <- lmer(fin.area ~ time + (1|quad), 
                      REML = FALSE, data = lmerDat)
Cand.mod[[5]] <- lmer(fin.area ~ 1 + (1|quad), 
                      REML = FALSE, data=lmerDat)  

#create a vector of names to trace back models in set
mod_numbers <- paste("Cand.mod", 1:length(Cand.mod), sep=" ")	
mod_text <- c("Era x Size", "Era + Size", "Size", "Era", "Null model")

#generate AICc table with names
mod.aicctab <- aictab(cand.set= Cand.mod, modnames= mod_text, 
                      sort=TRUE, second.ord=TRUE) # second.ord =TRUE means AICc is used 
print(mod.aicctab, digits=2, LL=TRUE)

write.csv(mod.aicctab, "./output/growthAIC.csv")

globalMod <- Cand.mod[[1]]
summary(globalMod)

##########################################################
##########################################################
# MODEL PARAMETERS BY ERA
##########################################################
##########################################################
datSCpastTrunc <- datSCpast[datSCpast$ini.area <= 0.96, ]
datSCpresent
datSCpastTrunc

presMod <- lmer(fin.area ~ ini.area + (1|quad), 
                      REML = FALSE, data = datSCpresent)

pastMod <- lmer(fin.area ~ ini.area + (1|quad), 
                      REML = FALSE, data = datSCpastTrunc)

summary(presMod)
summary(pastMod)
sd(resid(pastMod))

##########################################################
##########################################################
# FIGURE - GROWTH SCALING BY ERA (TRUNCATED)
##########################################################
##########################################################
ylab1 <- expression(paste("Size at time t+3 (", cm^2, ")"))
xlab1 <- expression(paste("Size at time t (", cm^2, ")"))
ULClabel <- theme(plot.title = element_text(hjust = -0.07, vjust = 0.05, size = rel(1.2)))
theme_set(theme_classic(base_size = 12))

dim(lmerDat)
lmerDat$fitted <- fitted(globalMod)
plot(fitted ~ fin.area, data = lmerDat)

size1 <- ggplot(lmerDat, aes(ini.area, fin.area, color = time, shape = time)) +
	ylab(ylab1) + xlab(xlab1) + 
	theme(legend.justification = c(1, 0), legend.position = c(1, 0.01)) +
	theme(legend.title = element_blank()) + 
	geom_point(size = 4, alpha = 0.8, 
		position = position_jitter(h = 0.05)) +
	scale_colour_manual(breaks = c("past", "present"), 
		values = c("darkgray", "black"), 
		labels = c("1969-1972", "2007-2010")) +
	scale_shape_manual(breaks = c("past", "present"), 
		values = c(18, 20), 
		labels = c("1969-1972", "2007-2010"))

size1
size1 + layer(stat = "smooth", stat_params = list(method = "lm", se = FALSE, size = 10))
size1 + stat_smooth(method = "lm", se = FALSE, size = 1)


# but i want to plot the lmer output, taking into account random effects
size2 <- ggplot(lmerDat, aes(ini.area, fitted, color = time)) +
	ylab("Fitted survival") + xlab(xlab1) + 
	theme(legend.justification = c(1, 0), legend.position = c(1, 0.2)) +
	theme(legend.title = element_blank()) + 
	geom_point(size = 4, alpha = 1, pch = 21, 
		position = position_jitter(h = 0.05)) +
	scale_colour_manual(breaks = c("past", "present"), 
		values = c("darkgray", "black"))

size2 + layer(stat = "smooth", stat_params = list(method = "lm", se = FALSE, size = 10))

size1 + geom_point(aes(y = fitted, color = time)) + 
	layer(stat = "smooth", stat_params = list(method = "lm", se = FALSE, size = 10))

summary(globalMod) # random effects contribute to 0 variance, 
#hence the perfect overlap between GLM and GLMER

size1 + stat_smooth(method = "lm", se = TRUE, size = 2, 
	aes(color = NULL), linetype = 1, color = 1) + 
	geom_line(aes(y = fitted, color = time), size = 1, linetype = 2) +
	labs(title = "B") + ULClabel +
	geom_abline(a = 0, b = 1, linetype = 3)

ULClabel <- theme(plot.title = element_text(hjust = -0.15, vjust = 0, size = rel(1.2)))

sizePlot <- size1 + 
	geom_smooth(method = "lm", se = FALSE, size = 1) + 
	labs(title = "B") + ULClabel + 
	geom_abline(a = 0, b = 1, linetype = 2, color = "gray") 

sizePlot

# Save lmerDat as different file
dat_growth <- lmerDat
