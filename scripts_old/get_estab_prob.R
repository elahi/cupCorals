#################################################
# Author: Robin Elahi
# Date: 160314
# Function to calculate fecundity related parameters that change with temperature or activation energy
#################################################

##' Function takes as input:
##' E = desired activation energy
##' k = Boltzmann's constant
##' kelvin = desired temperature
##' keepZeroes = TRUE to keep quadrats where embryos = 0
##' keepZeroes = FALSE to remove quadrats where embryos = 0
##' (in both cases, when recruits > embryos, then the quadrat is removed)

# source("bael_embryos.R")

# California derived parameters (from embryos script)
xIntCA; yIntCA; slopeCA
ltDat # life table

### Constants ###
k <- 0.0000862 # boltzmann constant (metabolic theory)
E <- 0.65 #activation energy for development rate = 0.65. E for PLD (1/dev) = -0.65

# Store the coefficient A from arrhenius equation for the y intercept of the embryo-size regression
embryo_A <- aCoef(yIntCA, E, k, temp = kelvin_CA, dir = "neg")

# Load the histo ipm data - need this for estimating recruitment probability
# Select the relevant columns
dat <- read.csv("./data/bael_ipmData.csv", header=TRUE, na.strings="NA")

kelvin = kelvin_WA

get_fecundity_params <- function(E, k, kelvin, keepZeroes = TRUE) {
  
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
  
  # Join embryos with recruits,  
  quad_DF <- inner_join(quad_recruitsDF, quad_embryosDF, by = "quad")
  
  # Calculate establishment probability
  ifelse(keepZeroes == TRUE,
         quad_DF <- quad_DF %>% 
           mutate(recProb = ifelse(quadRecruits == 0 & quadEmbryos == 0, 0,
                                   ifelse(quadRecruits > quadEmbryos, NA, 
                                          quadRecruits/quadEmbryos))), 
         
         quad_DF <- quad_DF %>% 
           mutate(recProb = ifelse(quadRecruits == 0 & quadEmbryos == 0, NA,
                                   ifelse(quadRecruits > quadEmbryos, NA, 
                                          quadRecruits/quadEmbryos))))
  
  # Summarize data
  estabDF <- quad_DF %>% summarise(estab.prob.mean = mean(recProb, na.rm = TRUE), 
                                   estab.prob.sd = sd(recProb, na.rm = TRUE), 
                                   embryos.mean = mean(quadEmbryos, na.rm = TRUE), 
                                   quad.n = length((recProb)[!is.na(recProb)]))
  
  estabDF$embryo.int <- yIntWA
  estabDF$mature.size <- xIntWA
  
  estabDF
  
  return(estabDF)

}

get_fecundity_params(E = 0.65, k, kelvin_WA, keepZeroes = TRUE)
