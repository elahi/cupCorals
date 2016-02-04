#################################################
# Author: Robin Elahi
# Date: 151203
# Multipanel plot
# Sensitivity analysis and growth rate comparison
#################################################

rm(list=ls(all=TRUE)) # removes all previous material from R's memory

##### LOAD PACKAGES ETC #####

library(grid)

# get sensitivity data
source("./bael_IPM_sensitivity.R")
head(plotDat)

# get growth data
source("./bael_growth.R")
head(dat_growth)

# ggplot2 settings
theme_set(theme_classic(base_size = 12))

##### FIGURE - SENSITIVITY ANALYSIS #####
str(plotDat)

# Relevel factors of variable
unique(plotDat$variable)
plotDat$variable <- as.factor(plotDat$variable)
print(levels(plotDat$variable))
plotDat$variable <- with(plotDat, factor(variable, levels(variable)[c(1, 3, 2)]))

ULClabel <- theme(plot.title = element_text(hjust = -0.1, vjust = 1, size = rel(1.2)))

sensitivityPlot <- ggplot(plotDat, aes(param, value)) +
  geom_bar(stat = "identity", color = "black", fill = "darkgray") +
  geom_hline(yintercept = c(0, 10), color = "gray", linetype = "dashed") + 
  # geom_hline(yintercept = 10, lty = 2, color = "darkgray") + 
  theme(legend.justification = c(1, 1), legend.position = c(0.95, 0.95)) +
  scale_fill_discrete(name = "Perturbation") +
  ylab("Percentage change") + xlab("Parameter") +
  scale_x_discrete(limits = rev(paramOrder), labels = rev(paramLabels)) +
  theme(strip.background = element_blank()) +
  # scale_y_continuous(limits = c(-5, 35)) + 
  theme(axis.text.y = element_text(size = 8)) + 
  facet_wrap(~ variable) + coord_flip() + 
  labs(title = "A") + ULClabel + 
  theme(panel.margin = unit(1.5, "lines"))

sensitivityPlot

##### FIGURE - GROWTH SCALING BY ERA #####
ylab_growth <- expression(paste("Size at time t+3 (", cm^2, ")"))

xlab_growth <- expression(paste("Size at time t (", cm^2, ")"))

ULClabel <- theme(plot.title = element_text(hjust = -0.15, vjust = 1, size = rel(1.2)))

size1 <- ggplot(dat_growth, aes(ini.area, fin.area, color = time, shape = time)) +
  ylab(ylab_growth) + xlab(xlab_growth) + 
  theme(legend.justification = c(0, 0), legend.position = c(0.65, 0.03)) +
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
  labs(title = "B") + ULClabel + 
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "black", size = 0.5) 
  
sizePlot

##### MULTI-PANEL PLOT #####
pdf("./figs/sensitivityGrowthPlot.pdf", width = 9, height = 3.5)
multiplot(sensitivityPlot, sizePlot, 
          layout = matrix(c(1, 1, 2), nrow = 1, byrow = TRUE))

dev.off()

##### MULTI-PANEL PLOT - 2 rows #####
pdf("./figs/sensitivityGrowthPlot_2rows.pdf", width = 5, height = 7)
multiplot(sensitivityPlot, sizePlot, 
          layout = matrix(c(1, 2), nrow = 2, byrow = TRUE))

dev.off()

