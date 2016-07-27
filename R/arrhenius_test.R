#################################################
# Author: Robin Elahi
# Date: 160727
# Arrhenius tests
#################################################

library(ggplot2)
library(tidyr)

##### SET UP DATAFRAME TO STORE MODEL PARAMETERS #####
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

##### GROWTH FUNCTION #####
source("./bael_growth.R")

# Want the model fitted to modern data
growthMod <- presMod
summary(growthMod)
head(datSCpresent)

# Store vital rate parameters
params$growth.int <- summary(growthMod)[[10]][1]
params$growth.slope <- summary(growthMod)[[10]][2]
params$growth.sd <- sd(resid(growthMod))

datSCpresent %>%
  ggplot(aes(ini.area, fin.area)) +
  geom_point() + geom_smooth(method = "lm") + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")

datSCpresent %>%
  ggplot(aes(ini.area, delta)) +
  geom_point() + geom_smooth(method = "lm") + 
  geom_hline(yintercept = 0, linetype = "dashed")

# What is the average growth rate?
growth_mean = mean(datSCpresent$delta, na.rm = TRUE)
# What is the relevant temperature?
source("bael_temperature.R")
modTemp <- as.numeric(sc_meanTemp_pres) + 273.15

##### PLOT EXPONENTIAL FUNCTION #####
eq <- function(c1, a1, x1){
  c1 * exp(a1 * x1)
}

c1 = 2
a1 = -0.1
a2 = 0.1
x_seq <- seq(0, 25)
y_seq <- eq(c1, a1, x1 = x_seq)
y_seq2 <- eq(c1, a2, x1 = x_seq)

xydf <- data.frame(x_seq, y_seq, y_seq2)
xydfL <- gather(xydf, key = sequence, value = y, y_seq:y_seq2)

##' if a is negative, the function decreases exponentially
##' if a is positive, the function increases exponentially
xydfL %>% 
  ggplot(aes(x_seq, y, color = sequence)) + 
  geom_line()

##### CALCULATE NEW GROWTH RATES BASED ON ARRHENIUS #####

##' Arrhenius equation
##' rate = a * e ^ (-E / k * T) 
##' a = constant (specific to each rate)
##' E = activation energy
##' k = boltzmann constant
##' T = temp in Kelvin

### Constants ###
E <- 0.65 # activation energy
k <- 0.0000862 # boltzmann constant

##' Calculate a for growth rate
growth_mean # mean change in surface area over 3 years
modTemp # temperature in kelvin

##' Solve for a
##' a = R / (exp(-E / k * T))
##' 

# My solution
get_a <- function(rate, E, k, temp) {
    rate / (exp(-E/(k * temp)))
  } 
get_a(rate = growth_mean, E, k, temp = modTemp)
aCoef(rate = growth_mean, E, k, temp = modTemp, dir = "pos")

coef_growth <- aCoef(rate = growth_mean, E, k, temp = modTemp)
coef_growth


get_new_rate <- function(coef_a, E, k, temp) {
  coef_a * exp(-E/(k * temp))}

temp_seq <- seq(modTemp + 2, modTemp - 2, length.out = 20)

get_new_rate(coef_a = coef_growth, E, k, temp = temp_seq)
ArrF(coef_a = coef_growth, E, k, temp = temp_seq, dir = "pos")

growth_seq <- ArrF(coef_a = coef_growth, E, k, temp = temp_seq, dir = "pos")

growth_seq; temp_seq

theme_set(theme_bw())
qplot(temp_seq - 273.15, growth_seq, geom = c("line", "point")) + 
  geom_hline(yintercept = growth_mean, color = 'red') + 
  labs(y = "Growth (cm2 per 3 years)", x = "Temperature (C)") + 
  geom_point(aes(growth_mean))



