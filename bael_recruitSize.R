#################################################
# Author: Robin Elahi
# Date: 150828

# Coral recruit size
# Figure 4
#################################################

library(lme4)
library(ggplot2)
library(AICcmodavg)

#rm(list=ls(all=TRUE)) # removes all previous material from R's memory

source("./R/graphicalParams.R")

dat <- read.csv("./data/bael_recruitSizeData.csv", header=TRUE, na.strings="NA")

datPast <- dat[dat$time == "past", ]
datPres <- dat[dat$time == "present", ]

mean(datPast$area)
sd(datPast$area)

mean(datPres$area)
sd(datPres$area)

##########################################################
##########################################################
# LMER; TEST FOR SIZE DIFFERENCES BY ERA
##########################################################
##########################################################
names(dat)
lmerDat <- dat

# random effects
Cand.mod <- list()
Cand.mod[[1]] <- lmer(area ~ time + (1|quad),
								data = lmerDat)
								
Cand.mod[[2]] <- lmer(area ~ time + (1|quadOriginal),
								data = lmerDat)

Cand.mod[[3]] <- lmer(area ~ time + (1|quadOriginal) + (1|tran),
								data = lmerDat)

mod_numbers <- paste("Cand.mod", 1:length(Cand.mod), sep=" ")	
mod_text <- c("quad", "quadOrig", "quadOrig + tran")     

mod.aicctab <- aictab(cand.set= Cand.mod, modnames= mod_text, 
                      sort=TRUE, second.ord=TRUE) # second.ord =TRUE means AICc is used 
print(mod.aicctab, digits=2, LL=TRUE)

# fixed effects
Cand.mod <- list()
Cand.mod[[1]] <- lmer(area ~ time + (1|quad), 
								REML = FALSE, data = lmerDat)
								
Cand.mod[[2]] <- lmer(area ~ 1 + (1|quad),
								REML = FALSE, data = lmerDat)

mod_text <- c("Era", "Null model")     

mod.aicctab <- aictab(cand.set= Cand.mod, modnames= mod_text, 
                      sort=TRUE, second.ord=TRUE) # second.ord =TRUE means AICc is used 
print(mod.aicctab, digits=2, LL=TRUE)
write.csv(mod.aicctab, "./output/recruitSizeAIC.csv")

globalMod <- Cand.mod[[1]]
											
# check normality and homogeneity of variances
par(mfrow = c(1,2))
qqnorm(resid(globalMod))
qqline(resid(globalMod))
plot(resid(globalMod) ~ fitted(globalMod)); abline(h=0)

##########################################################
##########################################################
# FIGURE - RECRUIT SIZE BY ERA 
##########################################################
##########################################################
with(dat, table(time)) # sample sizes

boxplot(area ~ time, data = dat, notch = TRUE)

ylab1 <- expression(paste("Recruit size (", cm^2, ")"))

ULClabel <- theme(plot.title = element_text(hjust = -0.1, vjust = 0.05, size = rel(1.2)))
theme_set(theme_classic(base_size = 12))

#ULClabel <- theme(plot.title = element_text(hjust = -0.18, vjust = 0, size = rel(1.2)))

recruit1 <- ggplot(data = dat, aes(time, area, fill = time)) + 
	geom_boxplot() +
	ylab(ylab1) + xlab("Year") +
	scale_fill_manual(breaks = c("past", "present"), 
		values = c("darkgray", "white")) +
	theme(legend.position = "none") + 
	scale_x_discrete(labels = c("1972", "2010"))

recruitPlot <- recruit1 + labs(title = "C") + ULClabel 
recruitPlot

# Save lmerDat as different file
dat_recruitSize <- lmerDat
