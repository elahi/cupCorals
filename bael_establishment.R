#################################################
# Author: Robin Elahi
# Date: 151318
# Coral establishment probability 
#################################################

# rm(list=ls(all=TRUE)) 

##### LOAD PACKAGES, DATA #####
source("./bael_embryos_Ea.R") # this also loads bael_functions
source("R/get_fecundity_params.R")

### Test function

E = .655

kzTrue_annual <- get_fec(E = E, k, kelvin_WA, keepZeroes = T, annual = TRUE) %>%
  summarise(estab.mean = mean(estab.prob.mean), 
            estab.sd = sd(estab.prob.sd))

kzFalse_annual <- get_fec(E = E, k, kelvin_WA, keepZeroes = F, annual = TRUE) %>%
  summarise(estab.mean = mean(estab.prob.mean), 
            estab.sd = sd(estab.prob.sd))

kzTrue_all <- get_fec(E = E, k, kelvin_WA, keepZeroes = T, annual = FALSE) %>%
  select(., starts_with("estab"))

kzFalsee_all <- get_fec(E = E, k, kelvin_WA, keepZeroes = F, annual = FALSE) %>%
  select(., starts_with("estab"))


kzTrue_annual; kzFalse_annual; kzTrue_all; kzFalsee_all

get_fec(E = E, k, kelvin_WA, keepZeroes = T, annual = TRUE)
get_fec(E = E, k, kelvin_WA, keepZeroes = T, annual = FALSE)


