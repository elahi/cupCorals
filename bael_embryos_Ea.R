#################################################
# Author: Robin Elahi
# Date: 160314
# Coral embryos - testing a range of activation energies
#################################################

##### LOAD PACKAGES, DATA #####
source("./bael_embryos.R") # this also loads bael_functions
source("R/get_estab_prob.R")
library(tidyr)

# California derived parameters (from embryos script)
xIntCA; yIntCA; slopeCA
ltDat # life table

##' 1. Start with California embryo-size relationship
##' 2. Estimate a modified embryo-size relationship
##' 3. Extract size at new intercept, and get new regression coefficients
##' 4. Use the new embryo function to calculate establishment probability (mean and sd)

##### LOOP THROUGH A RANGE OF EA #####

vector_Ea <- seq(0.6, 0.7, by = 0.001)
loopDat <- data.frame(Ea = vector_Ea, 
                      estab.prob.mean = NA, 
                      estab.prob.sd = NA, 
                      embryo.int = NA, 
                      mature.size = NA, 
                      embryos.mean = NA, 
                      quad.n = NA)
loopDat

for (i in 1:length(vector_Ea)) {
  E.i <- loopDat[i, ]$Ea
  
  estabDF.i <- get_fecundity_params(E = E.i, 
                                    k = k, kelvin = kelvin_WA)

  loopDat[i, ]$Ea <- E.i
  loopDat[i, ]$estab.prob.mean <- estabDF.i$estab.prob.mean
  loopDat[i, ]$estab.prob.sd <- estabDF.i$estab.prob.sd
  loopDat[i, ]$embryo.int <- estabDF.i$embryo.int
  loopDat[i, ]$mature.size <- estabDF.i$mature.size
  loopDat[i, ]$embryos.mean <- estabDF.i$embryos.mean
  loopDat[i, ]$quad.n <- estabDF.i$quad.n
  
}

estabDF.i
loopDat

##### CLEAN UP RESULTS #####

# Remove data where mature.size exceeds the observed value
max(dat$area, na.rm = TRUE)
head(loopDat)
tail(loopDat)

# loopDat2 <- loopDat %>% filter(estab.prob.mean > 0)

loopLong <- loopDat %>% gather(key = parameter, value = value, 
                                estab.prob.mean:quad.n)

head(loopLong)
str(loopLong)
unique(loopLong$parameter)

# Reorder factor levels
loopLong$parameter <- factor(loopLong$parameter)
levels(loopLong$parameter)
loopLong$parameter <- factor(loopLong$parameter, 
                             levels(loopLong$parameter)[c(1,5,2,6,3,4)])

##### PLOT RESULTS #####

ggDat <- loopLong

ggplot(ggDat, aes(Ea, value)) +
  geom_point(size = 0.5) + geom_line(color = 'red') +
  facet_wrap(~ parameter, scales = "free") +
  ylab("Parameter value") + xlab("Activation energy (Ea)") +
  geom_point(data = loopLong[loopLong$Ea == 0.65, ],
             aes(Ea, value), color = "blue")

ggsave("figs/fecundityParams_Ea_sensitivity.pdf",
       height = 5, width = 8)




