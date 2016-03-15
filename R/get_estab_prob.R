#################################################
# Author: Robin Elahi
# Date: 160314
# Function to calculate fecundity related parameters that change with temperature or activation energy
#################################################


##### LOAD PACKAGES, DATA #####
source("./bael_embryos.R") # this also loads bael_functions
library(tidyr)
# California derived parameters from this embryos script
xIntCA; yIntCA; slopeCA
ltDat

# Store the coefficient A from arrhenius equation for the y intercept of the embryo-size regression
embryo_A <- aCoef(yIntCA, E, k, temp = kelvin_CA, dir = "neg")

# Load the ipm data - need this for estimating recruitment probability
dat <- read.csv("./data/bael_ipmData.csv", header=TRUE, na.strings="NA")

##' 1. Start with California embryo-size relationship
##' 2. Estimate a modified embryo-size relationship
##' 3. Extract size at new intercept, and get new regression coefficients
##' 4. Use the new embryo function to calculate establishment probability (mean and sd)

##' Function takes as input:
##' E = desired activation energy
##' k = Boltzmann's constant
##' kelvin = desired temperature

##### START FUNCTION #####

get_fecundity_params <- function(E, k, kelvin) {
  
  ### PART 1: ESTIMATE NEW EMBRYO-SIZE RELATIONSHIP
  # Calculate new y-intercept and x-intercept (size at maturity)
  yIntWA <- ArrF(embryo_A, E = E, k = k, temp = kelvin, dir = "neg")
  xIntWA <- -yIntWA/slopeCA

  ### PART 2: ESTIMATE NEW RECRUITMENT PROBABILITY
  # Calculate the number of embryos produced by each coral
  dat$embryos <- ifelse(dat$area < xIntWA, 0, 
                                  dat$area * slopeCA + yIntWA)
  
  # Remove 2010 data, and then get total number of embryos per quad
  quad_embryosDF <- dat %>% filter(date.no != 40515) %>% 
    group_by(quad) %>%
    summarise(quadEmbryos = sum(embryos, na.rm = TRUE))
  
  # These are the observed recruits every year
  recruitDat <- dat %>% filter(recruit ==1) %>% select(coral.id, code) %>%
    rename(recruitCode = code)

  quad_recruitsDF <- dat %>% filter(date.no == 40515) %>% 
    left_join(., recruitDat, by = "coral.id") %>%
    group_by(quad) %>%
    summarise(quadRecruits = recruitF(recruitCode)) 
  
  # Join embryos with recruits, then calculate establishment probability
  quad_DF <- inner_join(quad_recruitsDF, quad_embryosDF, by = "quad") %>% 
    mutate(recProb = quadRecruits/quadEmbryos)
  
  quad_DF$recProb[is.infinite(quad_DF$recProb)] <- NaN

  estabDF <- quad_DF %>% summarise(estab.prob.mean = mean(recProb, na.rm = TRUE), 
                                   estab.prob.sd = sd(recProb, na.rm = TRUE), 
                                   embryos.mean = mean(quadEmbryos, na.rm = TRUE))
  
  estabDF$embryo.int <- yIntWA
  estabDF$mature.size <- xIntWA

  return(estabDF)

}


get_fecundity_params(E = E, k = k, kelvin = kelvin_CA)


vector_Ea <- seq(0.5, 0.7, by = 0.005)
loopDat <- data.frame(Ea = vector_Ea, 
                      estab.prob.mean = NA, 
                      estab.prob.sd = NA, 
                      embryo.int = NA, 
                      mature.size = NA, 
                      embryos.mean = NA)
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
  
}

loopDat

# Remove data where mature.size exceeds the observed value
max(dat$area, na.rm = TRUE)
loopDat2 <- loopDat %>% filter(mature.size < 1.1 &
                                estab.prob.mean < 1)

loopLong <- loopDat2 %>% gather(key = parameter, value = value, 
                               estab.prob.mean:embryos.mean)

head(loopLong)
str(loopLong)
unique(loopLong$parameter)

# Reorder factor levels
loopLong$parameter <- factor(loopLong$parameter)
levels(loopLong$parameter)
loopLong$parameter <- factor(loopLong$parameter, 
                    levels(loopLong$parameter)[c(1,5,2,3,4)])


ggplot(loopLong, aes(Ea, value)) + 
  geom_point(size = 0.5) + geom_line(color = 'red') +
  facet_wrap(~ parameter, scales = "free") + 
  ylab("Parameter value") + xlab("Activation energy (Ea)") + 
  geom_point(data = loopLong[loopLong$Ea == 0.65, ], 
             aes(Ea, value), color = "blue")

ggsave("figs/fecundityParams_Ea_sensitivity.pdf", 
       height = 3.5, width = 5)



