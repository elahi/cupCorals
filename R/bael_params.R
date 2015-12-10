#################################################
# Author: Robin Elahi
# Date: 151209
# Parameters for IPM
#################################################

rm(list=ls(all=TRUE)) 

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

##### SURVIVAL FUNCTION #####
source("./bael_survival.R")

# These are the relevant outputs
presMod
pastMod

# Rename relevant model
survMod <- presMod

# Store vital rate parameters
params$surv.int <- summary(survMod)[[10]][1]
params$surv.slope <- summary(survMod)[[10]][2]

params

##### GROWTH FUNCTION #####
source("./bael_growth.R")

# These are the relevant outputs
presMod
pastMod

# Rename relevant model
growthMod <- presMod

# Store vital rate parameters
params$growth.int <- summary(growthMod)[[10]][1]
params$growth.slope <- summary(growthMod)[[10]][2]
params$growth.sd <- sd(resid(growthMod))

params

##### SIZE DISTRIBUTION OF RECRUITS #####
dat <- read.csv("./data/bael_recruitSizeData.csv", header=TRUE, na.strings="NA")
datPres <- dat[dat$time == "present", ]

# Store vital rate parameters
params$recruit.size.mean <- mean(datPres$area)
params$recruit.size.sd <- sd(datPres$area)

params

##### CREATE FORK FOR CALIFORNIA AND WASHINGTON #####
# The embryo function is originally from California, 
# and was modified using the Arrhenius equation for Washington, 
# which has colder seawater

# The modified embryo-size function causes downstream changes to
# size at maturity and establishment probability

paramsWA <- params
paramsCA <- params

##### EMBRYO FUNCTION #####
# Note that the embryo function determines size at maturity and
# establishment probability

source("./bael_embryos.R")

# These are the relevant outputs
mxRegWA
mxRegCA

# Store vital rate parameters for California
paramsCA$embryo.int <- coefficients(mxRegCA)[1]
paramsCA$embryo.slope <- coefficients(mxRegCA)[2]
paramsCA$embryo.sd <- sd(resid(mxRegCA))

# Store vital rate parameters for Washington
paramsWA$embryo.int <- coefficients(mxRegWA)[1]
paramsWA$embryo.slope <- coefficients(mxRegWA)[2]
paramsWA$embryo.sd <- sd(resid(mxRegWA))

paramsCA
paramsWA

##### SIZE AT MATURITY #####
paramsCA$mature.size <- -paramsCA$embryo.int/paramsCA$embryo.slope
paramsWA$mature.size <- -paramsWA$embryo.int/paramsCA$embryo.slope

##### ESTABLISHMENT PROBABILITY #####
source("./bael_establishment.R")

paramsWA$estab.prob.mean <- mean(recProbWA/100, na.rm = TRUE) 
paramsWA$estab.prob.sd <- sd(recProbWA/100, na.rm = TRUE)

paramsCA$estab.prob.mean <- mean(recProbCA/100, na.rm = TRUE) 
paramsCA$estab.prob.sd <- sd(recProbCA/100, na.rm = TRUE)

##### COMPLETE DATAFRAME #####
paramsCA
paramsWA
