#################################################
# Author: Robin Elahi
# Date: 151208
# Coral growth
# Figure 4
#################################################

# rm(list=ls(all=TRUE)) 

##### LOAD PACKAGES, DATA #####
library(lme4)
library(ggplot2)
theme_set(theme_classic(base_size = 16))
library(AICcmodavg)
library(dplyr)

dat <- read.csv("./data/bael_growthData.csv", header=TRUE, na.strings="NA")
source("./R/graphicalParams.R")
dat$ini.areaLN <- log(dat$ini.area)
dat$fin.areaLN <- log(dat$fin.area)

datSC <- dat
# Create dataset for SC present
datSCpresent <- datSC[datSC$time == "present", ]
# Create dataset for SC past
datSCpast <- datSC[datSC$time == "past", ]

# Truncate initial size range to reflect the largest observed size in 2007-2010 (1.0 cm2)
set_graph_pars(ptype = "panel1")
plot(fin.area ~ ini.area, data = datSCpast)
points(fin.area ~ ini.area, data = datSCpresent, col = "red")
abline(v = 1.0, col = "black", lty = 3)
abline(a = 0, b = 1, col = "darkgray", lty = 2, lwd = 2)

datSCTrunc <- datSC[datSC$ini.area <= 1.0, ]

##### USE LMER TO TEST THE EFFECT OF ERA ON CORAL GROWTH #####
### Use truncated historical dataset to match observed size range in modern dataset
lmerDat <- datSCTrunc

# Random effects are quadrat (original quadrats for historical study, not subquads)
# Don't need varying slopes by quadrat 

Cand.mod <- list()                                        
Cand.mod[[1]] <- lmer(fin.area ~ ini.area * time + (1|quad), 
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

# Save lmerDat as different file
dat_growth <- lmerDat

##### GET MODEL PARAMETERS BY ERA #####
datSCpastTrunc <- datSCpast[datSCpast$ini.area <= 1.0, ]

presMod <- lmer(fin.area ~ ini.area + (1|quad), 
                      REML = FALSE, data = datSCpresent)

pastMod <- lmer(fin.area ~ ini.area + (1|quad), 
                      REML = FALSE, data = datSCpastTrunc)

summary(presMod)
summary(pastMod)
sd(resid(pastMod))

##### FINAL FIGURE - GROWTH SCALING BY ERA #####
ylab_growth <- expression(paste("Size at time t+3 (", cm^2, ")"))

xlab_growth <- expression(paste("Size at time t (", cm^2, ")"))

ULClabel <- theme(plot.title = element_text(hjust = -0.15, vjust = 1, size = rel(1.2)))

size1 <- ggplot(dat_growth, aes(ini.area, fin.area, color = time, shape = time)) +
  ylab(ylab_growth) + xlab(xlab_growth) + 
  theme(legend.justification = c(0, 0), legend.position = c(0.5, 0)) +
  theme(legend.title = element_blank()) + 
  geom_point(size = 2.5, alpha = 0.6, 
             position = position_jitter(h = 0.05)) +
  scale_colour_manual(breaks = c("past", "present"), 
                      values = c("darkgray", "black"), 
                      labels = c("1969-1972", "2007-2010")) +
  scale_shape_manual(breaks = c("past", "present"), 
                     values = c(18, 20), 
                     labels = c("1969-1972", "2007-2010")) 

sizePlot <- size1 + 
  geom_smooth(method = "lm", se = FALSE, size = 0.75) + 
  # labs(title = "B") + ULClabel + 
  geom_abline(a = 0, b = 1, linetype = 2, color = "black", size = 0.5) 

sizePlot

ggsave("./figs/growthPlot.pdf", height = 3.5, width = 3.5)

