###############################################################
# Author: Robin Elahi
# Date: 160927
# Modeling mass transfer as a function of size
###############################################################

##' Following theory from Patterson 1992 
##' Using coefficients in Nakamura et al. 2001
##' Goal: duplicate Figure 1 in Nakamura 

library(dplyr)

c1 = 2 # constant
Dv = 2 * 10^-9 # diffusion coefficient of O2 in water
W = seq(0.01, 0.45, by = 0.01) # vector of organism's characteristic dimension
d = 0.5 # flow size exponent; here I just use 0.5 to represent mass transfer in laminar flow
rho = 998.4 # fluid density
U = 0.4 # water velocity; constant at 0.4 ms-1
mu = 0.001 # dynamic viscosity of water


dat <- data.frame(c1, Dv, W, d, rho, U, mu)
head(dat)

# Now calculate mass transfer coefficient

dat <- dat %>% 
  mutate(hm = (c1 * Dv * W^(d - 1) * rho^d * U^d) / mu^d)

plot(hm ~ W, data = dat, type = "b")

# Clearly, mass transfer for a given area is more rapid than to and from large organisms
