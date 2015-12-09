#################################################
# Author: Robin Elahi
# Date: 151208
# Functions for cup coral calculations
#################################################

##### FUNCTION TO CODE SURVIVAL #####
# Define survival as zero if it is has "dead" in notes column
survF <- function(code) {ifelse(code == "dead", 0, 1)}

##### FUNCTION TO CALCULATE CORAL AREA GIVEN A VOLUME #####
# Calculate area given a volume
areaCalcF <- function(x){
	exp(log(x)*avReg$coefficients[2] + 
               avReg$coefficients[1])/100
}

##### FUNCTION TO CALCULATE NUMBER OF EMBRYOS PER CORAL #####

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

##### FUNCTION TO CALCULATE ADULT DENSITY #####
# do not count NAs or corals smaller than 0.01 (?)

densityF <- function(X1) {
  D1 <- X1 > 0.01 # do not include corals smaller than 0.01
  sum(D1, na.rm = TRUE) # sum up true values, do not include NAs
}

##### FUNCTION TO CALCULATE RECRUIT DENSITY #####
# define a recruit as one that is has "recruit" in notes column

recruitF <- function(X1) {
  D1 <- X1 == "recruit" 
  sum(D1, na.rm = TRUE) # sum up true values, do not include NAs
}
