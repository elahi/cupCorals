#################################################
# Author: Robin Elahi
# Date: 150828

# Multipanel plot of survival, growth, and recruit size
# Figure 4
#################################################

rm(list=ls(all=TRUE)) # removes all previous material from R's memory

# get growth data
source("./bael_growth.R")
head(dat_growth)

# get survival data
source("./bael_survival.R")
head(dat_survival)

# get recruit size data
source("./bael_recruitSize.R")
head(dat_recruitSize)

# ggplot2 settings
theme_set(theme_classic(base_size = 12))

##########################################################
# FIGURE - SURVIVAL SCALING BY ERA
##########################################################
ylab_surv <- "Survival at time t+3"
xlab_surv <- expression(paste("Size at time t (", cm^2, ")"))
ULClabel <- theme(plot.title = element_text(hjust = -0.2, vjust = 1, size = rel(1.2)))

surv1 <- ggplot(dat_survival, aes(ini.area, survival, 
                             color = time, shape = time)) +
  ylab(ylab_surv) + xlab(xlab_surv) + 
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.2)) +
  theme(legend.title = element_blank()) + 
  geom_point(size = 2, alpha = 0.8, 
             position = position_jitter(h = 0.05)) +
  scale_colour_manual(breaks = c("past", "present"), 
                      values = c("darkgray", "black"), 
                      labels = c("1969-1972", "2007-2010")) +
  scale_shape_manual(breaks = c("past", "present"), 
                     values = c(18, 20), 
                     labels = c("1969-1972", "2007-2010"))

survPlot <- surv1 + labs(title = "A") + ULClabel + 
  survTrendPast2 + survTrendPres2 

# theme(legend.position = "none")
survPlot

##########################################################
# FIGURE - GROWTH SCALING BY ERA
##########################################################
ylab_growth <- expression(paste("Size at time t+3 (", cm^2, ")"))
xlab_growth <- expression(paste("Size at time t (", cm^2, ")"))
ULClabel <- theme(plot.title = element_text(hjust = -0.2, vjust = 1, size = rel(1.2)))

size1 <- ggplot(dat_growth, aes(ini.area, fin.area, color = time, shape = time)) +
  ylab(ylab_growth) + xlab(xlab_growth) + 
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.01)) +
  theme(legend.title = element_blank()) + 
  geom_point(size = 2, alpha = 0.8, 
             position = position_jitter(h = 0.05)) +
  scale_colour_manual(breaks = c("past", "present"), 
                      values = c("darkgray", "black"), 
                      labels = c("1969-1972", "2007-2010")) +
  scale_shape_manual(breaks = c("past", "present"), 
                     values = c(18, 20), 
                     labels = c("1969-1972", "2007-2010"))

sizePlot <- size1 + 
  geom_smooth(method = "lm", se = FALSE, size = 1) + 
  labs(title = "B") + ULClabel + 
  geom_abline(a = 0, b = 1, linetype = 2, color = "black", size = 0.5) + 
  theme(legend.position = "none")
  
sizePlot

##########################################################
# FIGURE - RECRUIT SIZE BY ERA
##########################################################
ylab_recruitSize <- expression(paste("Recruit size (", cm^2, ")"))
ULClabel <- theme(plot.title = element_text(hjust = -0.2, vjust = 1, size = rel(1.2)))

recruit1 <- ggplot(data = dat_recruitSize, aes(time, area, fill = time)) + 
  geom_boxplot() +
  ylab(ylab_recruitSize) + xlab("Year") +
  scale_fill_manual(breaks = c("past", "present"), 
                    values = c("darkgray", "white")) +
  theme(legend.position = "none") + 
  scale_x_discrete(labels = c("1972", "2010"))

recruitPlot <- recruit1 + labs(title = "C") + ULClabel 
recruitPlot


##########################################################
# Multi-panel plot
##########################################################
pdf("./figs/vitalRatesPlot.pdf", 7, 3.5)
multiplot(survPlot, sizePlot, recruitPlot, cols = 3)
dev.off()
