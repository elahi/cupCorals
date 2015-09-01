# RE 150312
# BAEL functions
# Modified the embryo function to reflect WA temperature
# (previously I used Fadlallah's value, but that was quantified at 13.1 C)

# Function to code survival
# Define survival as zero if it is has "dead" in notes column
survF <- function(code) {ifelse(code == "dead", 0, 1)}


####################
# Calculate area given a volume
areaCalcF <- function(x){
	exp(log(x)*avReg$coefficients[2] + 
               avReg$coefficients[1])/100
}
####################


###########################################
# Function to calculate the # of embryos per coral
# California embryo function, unadulterated
embryoFca <- function(area) {
  ifelse(area < xIntCA, 
         0, 
         area*mxReg$coefficients[2] - mxReg$coefficients[1])
}
matureFca <- function(area) {
  ifelse(area < xIntCA, 0, 1)
}


# NEW FUNCTION for modified WA regression
embryoFwa <- function(area) {
  ifelse(area < xIntWA, 
         0, 
         area*mxRegWA$coefficients[2] + mxRegWA$coefficients[1])
}
matureFwa <- function(area) {
  ifelse(area < xIntWA, 0, 1)
}

#########################################


############################################
# Function to calculate density
# do not count NAs or corals smaller than 0.01 (?)

densityF <- function(X1) {
  D1 <- X1 > 0.01 # do not include corals smaller than 0.01
  sum(D1, na.rm = TRUE) # sum up true values, do not include NAs
}
############################################

############################################
# Function to calculate recruit density
# define a recruit as one that is has "recruit" in notes column

recruitF <- function(X1) {
  D1 <- X1 == "recruit" 
  sum(D1, na.rm = TRUE) # sum up true values, do not include NAs
}
############################################
#IPM
###1.4.3 Define vital rate functions to describe life history###
# 1. survival probability function
s.x <-function(x, params) {
  u = exp(params$surv.int + params$surv.slope*x)
  return(u/(1+u))
}

# 2. growth function
g.yx <- function(xp, x, params) {
  dnorm(xp, mean = params$growth.int + params$growth.slope*x,
        sd  = params$growth.sd)
}

# 3. reproduction function
f.yx <- function(xp, x, params) {
  embryos = ifelse(x < params$mature.size, 0, 
        x*params$embryo.slope + params$embryo.int)
  embryos * params$estab.prob.mean *
  dnorm(xp, mean = params$recruit.size.mean, 
          sd = params$recruit.size.sd)
}


    
