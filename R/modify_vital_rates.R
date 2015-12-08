#################################################
# Author: Robin Elahi
# Date: 151207

# FUNCTION TO GET TEMPERATURE-MODIFIED PREDICTIONS
# OF VITAL RATES BASED ON THE ARRHENIUS EQUATION
#################################################

##### DESCRIPTION OF FUNCTION #####
### This function takes eleven arguments
# slope = the slope of the vital rate in question
# intercept = the intercept of the vital rate
# slopeTempEffect = the hypothesized effect of a temperature increase 
 # on the slope (takes "pos" or "neg" as arguments)
# interceptTempEffect = the hypothesized effect of a temperature increase on the slope 
 # on the intercept (takes "pos" or "neg" as arguments)
# lowerEa = minimum activation energy
# upperEa = maximum activation energy
# Ea_increment = incremental steps of Ea
# lowerTemp = minimum temperature
# upperTemp = maximum temperature
# temp_increment = incremental steps of temperature
# originalTemp = the original temperature at which the vital rate was measured

##### FUNCTION #####
modifyVitalRates <- function(slope, slopeTempEffect, 
                             intercept, interceptTempEffect, 
                             lowerEa, upperEa, Ea_increment, 
                             lowerTemp, upperTemp, temp_increment,
                             originalTemp) {
  
  # SET UP VECTOR OF ACTIVATION ENERGIES
  vec_Ea <- seq(lowerEa, upperEa, by = Ea_increment) 
  
  # SET UP VECTOR OF TEMPERATURE
  vec_temp <- seq(lowerTemp, upperTemp, by = temp_increment) 
  
  # EXPAND THE EaDF BY TEMP GRID
  grid1 <- expand.grid(vec_temp, vec_Ea)
  names(grid1) <- c("Kelvin", "Ea")
  grid1$row <- seq(from = 1, to = dim(grid1)[1], by = 1)
  
  # Get pre-expontial coefficient, a, for each Ea FOR INTERCEPT
  grid1$intA <- aCoef(intercept, E = grid1$Ea, k, temp = originalTemp, 
                      dir = interceptTempEffect)
  
  # Get pre-expontial coefficient, a, for each Ea FOR SLOPE
  grid1$slopeA <- aCoef(slope, E = grid1$Ea, k, temp = originalTemp, 
                        dir = slopeTempEffect)
  
  # Calculate intercepts as a function of temp and Ea
  grid1$vec_int <- ArrF(coef_a = grid1$intA, E = grid1$Ea, k, 
                        temp = grid1$Kelvin, dir = interceptTempEffect)
  
  # Calculate slopes as a function of temp and Ea
  grid1$vec_slope <- ArrF(coef_a = grid1$slopeA, E = grid1$Ea, k, 
                          temp = grid1$Kelvin, dir = slopeTempEffect)
  
  # Add in columns for constant intercepts and slopes
  grid1$const_int <- rep(intercept, length = nrow(grid1))
  grid1$const_slope <- rep(slope, length = nrow(grid1))
  
  return(grid1)
}

