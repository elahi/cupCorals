#################################################
# Author: Robin Elahi
# Date: 151029

# Coral biomass scaling
# Data from Ken Sebens
#################################################

library(smatr)
library(ggplot2)
theme_set(theme_classic(base_size = 12))

source("./R/graphicalParams.R")

dat <- read.csv("./data/bael_CN_2015.csv")
dat

dat$massLN <- log(dat$biomass)
dat$areaLN <- log(dat$area)

###########
# Untransformed plot
# OLS regression for prediction
###########

head(dat)
dat$category <- rep("Empirical")

# change biomass from g to mg
biomass_mg <- dat$biomass * 1000
dat$biomass <- biomass_mg

# Regression
ols1 <- lm(biomass ~ area, data = dat)
summary(ols1)
ols1$coefficients

# predicted biomass for 1.00 cm2, and 1.67 cm2
area_given <- c(1, 1.67)

# predicted biomass for 99% max sizes (0.91, 1.44)
area_given <- c(0.91, 1.44)

predict_biomass <- function(area) {
  (area * ols1$coefficients[2] + 
                ols1$coefficients[1])
}

biomass_pred <- predict_biomass(area_given)
biomass_pred

# create dataframe for predicted biomass
predBiomass <- as.data.frame(cbind(area_given, biomass_pred))
predBiomass$era <- c("modern", "historic")

# Calculate percent change per degree
predBiomass
mass_Mod <- predBiomass$biomass_pred[1]
mass_Hist <- predBiomass$biomass_pred[2]

perChangeMass <- (mass_Hist - mass_Mod)/(mass_Hist) * 100
perChangeMass

# Temperature has increased by 0.65 C
# Normalize to one degree change to match Forster 2012
tempSizeResponse <- perChangeMass/0.65
part1 <- paste(round(tempSizeResponse, 1), 
      "% change in mass per", sep = "")
part1
text1 <- expression(paste("65% reduction in mass per", degree, "C"))
text1

label1 <- expression(paste("Surface area (", cm^2, ")"))
label2 <- "Biomass (mg)"

ggplot(data = dat, aes(area, biomass)) + 
  geom_point(shape = 1) + 
  geom_smooth(method = "lm", color = "black") + 
  xlab(label1) + ylab(label2) + 
  scale_x_continuous(limits = c(0, 1.5)) + 
  scale_y_continuous(limits = c(0, 140)) + 
  geom_point(data = predBiomass, 
             aes(area_given, biomass_pred, 
                 color = era, shape = era)) + 
  scale_color_manual(values = c("darkgray", "black")) +
  scale_shape_manual(values = c(18, 20)) + 
  theme(legend.position = "none")

# p + annotate("text", x = 0, y = 60, label = text1, 
#             size = 0.2)

ggsave("./figs/biomassPlot.pdf", height = 3.5, width = 3.5)  

###########
# Scaling plot
###########
mod0 <- ma(massLN ~ areaLN, data = dat, 
           slope.test = 1)
summary(mod0)

mod1 <- sma(massLN ~ areaLN, data = dat, 
            slope.test = 1)

summary(mod1)


label3 <- bquote(italic(R)^2 == 0.92)

set_graph_pars(ptype = "panel1")
par(cex = 1.5, las = 1)

plot(data = dat, massLN ~ areaLN, 
	xlim = c(-3, 1),  ylim = c(-8, -1), 
	xlab = label1, ylab = label2, lwd = 2)
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "gray")

text(0, -5, "y = 2.43 - 2.02", cex = 0.8)
text(0, -5.5, label3, cex = 0.8)

clip(min(areaLN), max(areaLN), -100, 100)
abline(mod1, lwd = 2)



# # Call: sma(formula = massLN ~ areaLN, slope.test = 1) 

# Fit using Standardized Major Axis 

# ------------------------------------------------------------
# Coefficients:
            # elevation    slope
# estimate    -2.016757 2.426102
# lower limit -2.412922 1.943084
# upper limit -1.620592 3.029189

# H0 : variables uncorrelated
# R-squared : 0.9246237 
# P-value : 9.1061e-06 

# ------------------------------------------------------------
# H0 : slope not different from 1 
# Test statistic : r= 0.9648 with 8 degrees of freedom under H0
# P-value : 6.4495e-06 
