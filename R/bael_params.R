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

##### EMBRYO SLOPE AND EMBRYO SD #####

source("./bael_embryos_Ea.R")

# Embryo intercept will vary with temperature or Ea, but slope and sd will remain the same
params$embryo.slope <- coefficients(mxRegCA)[2]
params$embryo.sd <- sd(resid(mxRegCA))

##### ESTABLISHMENT PROBABILITY #####
fecCA <- get_fec(E = 0.65, k, kelvin = kelvin_CA, annual = FALSE)
fecWA <- get_fec(E = 0.65, k, kelvin = kelvin_WA, annual = FALSE)
fecOptim <- get_fec(E = 0.655, k, kelvin = kelvin_WA, annual = FALSE)

paramsCA <- params
paramsWA <- params
paramsOptim <- params

paramsCA$estab.prob.mean <- fecCA$estab.prob.mean
paramsCA$estab.prob.sd <- fecCA$estab.prob.sd
paramsCA$mature.size <- fecCA$mature.size
paramsCA$embryo.int <- fecCA$embryo.int

paramsWA$estab.prob.mean <- fecWA$estab.prob.mean
paramsWA$estab.prob.sd <- fecWA$estab.prob.sd
paramsWA$mature.size <- fecWA$mature.size
paramsWA$embryo.int <- fecWA$embryo.int

paramsOptim$estab.prob.mean <- fecOptim$estab.prob.mean
paramsOptim$estab.prob.sd <- fecOptim$estab.prob.sd
paramsOptim$mature.size <- fecOptim$mature.size
paramsOptim$embryo.int <- fecOptim$embryo.int


##### COMPLETE DATAFRAME #####
paramsCA
paramsWA
paramsOptim
