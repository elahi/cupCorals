#################################################
# Author: Robin Elahi
# Date: 160314
# Coral embryos - testing a range of activation energies
#################################################

source("./embryos.R")

##### ESTIMATE WA X-INTERCEPT OVER A RANGE OF ACTIVATION ENERGIES #####

# Function to get new estimate of x-intercept for WA based on CA data
# a_coefficient is the same; based on california embryo-size function
xIntCA

embryo_A <- aCoef(y_intercept_CA, E, k, temp = kelvin_CA, dir = "neg")

get_x_int_WA <- function(rate, E, k, newKelvin) {
  # Get y-intercept for WA
  y_intercept_WA <- ArrF(embryo_A, E, k, temp = newKelvin, dir = "neg")
  # Calculate x-intercept for WA 
  xIntWA <- -y_intercept_WA/mxReg$coefficients[2]
  return(xIntWA)
}

get_x_int_WA(y_intercept_CA, E, k, newKelvin = kelvin_WA)

# We assume that size at maturity will be larger in the colder Washington population
# If we use an Ea = 0.65 (mean), then we get a size of ~ 0.4 cm2
# Can we explore a range of activation energies, because we don't know the sensitivity of this function to temperature?

vector_Ea <- seq(0.635, 0.67, by = 0.0005)
loopDat <- data.frame(Ea = vector_Ea, yIntWA = NA, xIntWA = NA)
loopDat

for (i in 1:length(vector_Ea)) {
  E.i <- loopDat[i, ]$Ea
  
  yIntWA.i <- ArrF(embryo_A, E = E.i, k, 
                   temp = kelvin_WA, dir = "neg")
  
  xIntWA.i <- get_x_int_WA(y_intercept_CA, E = E.i, k, 
                           newKelvin = kelvin_WA)
  
  loopDat[i, ]$yIntWA <- as.numeric(yIntWA.i)
  
  loopDat[i, ]$xIntWA <- as.numeric(xIntWA.i)
}

embryoInterceptsDF <- loopDat



sizeDiff <- xIntWA - xIntCA; sizeDiff

ltDat$areaRcorr <- ltDat$areaR + sizeDiff
mxRegWA <- lm(mx.coral ~ areaRcorr, data = ltDat[ltDat$Vol.mid > 75, ])
mxRegCA <- lm(mx.coral ~ areaR, data = ltDat[ltDat$Vol.mid > 75, ])

