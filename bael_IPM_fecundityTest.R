#################################################
# Author: Robin Elahi
# Date: 160315
# Testing the sensitivity of IPM results to variation in 
# Ea on x-intercept of fecundity function (which has
# downstream effects on size at maturity, establishment probability)
#################################################

rm(list=ls(all=TRUE)) 

##### LOAD PACKAGES ETC #####
source("./R/bael_params.R")
source("./R/ipmFunctions.R")
source("R/get_histo_ipm_data.R")

##### MAX AND MIN SIZES #####
# Will use slightly larger size range for IPM because
# it is otherwise artificially truncated
min.size <- 0.02
max.size <- 1.52

binSize <- 0.02
binN <- (max.size - min.size)/binSize
binN

##### IPMs WITH FOR LOOP THROUGH EAS #####
data_in <- loopDat2
N <- dim(data_in)[1]

# Lists to store the results of IPMs
ipm_out <- result <- vector("list", N) 
res_out <- result <- vector("list", N) 

# Dataframe to store Ea and the log-likelihood of the size distribution (ssd2) to the observed size-frequency distribution

# Use same bin sizes and range as IPM
histBreaks <- seq(min.size, max.size, binSize)

areaH_07 <- hist(hist10$area, breaks = histBreaks)
areaH_10 <- hist(hist10$area, breaks = histBreaks)
areaH_0710 <- hist(hist0710$area, breaks = histBreaks)

dns07 <- areaH_07$density
dns10 <- areaH_10$density
dns0710 <- areaH_0710$density

data_out <- data.frame(Ea = data_in$Ea, 
                       LL07 = NA, 
                       LL10 = NA, 
                       LL0710 = NA)

# Initialize parameter dataframe
paramDF <- paramsWA # the cells will be modified in the loop

for (i in 1:N){
  # Create new paramDF, with modified fecundity parameters
  params.i <- data_in[i, ]
  paramDF$estab.prob.mean <- params.i$estab.prob.mean
  paramDF$estab.prob.sd <- params.i$estab.prob.sd
  paramDF$embryo.int <- params.i$embryo.int
  paramDF$mature.size <- params.i$mature.size
  
  # Now run IPM
  ipm1 <- bigmatrix(n = binN, params = paramDF)
  res1 <- popF(ipm1, binSize)
  
  # Store log likelihood results in data_out
  Ea.i <- as.character(params.i$Ea)
  LL07.i <- - sum(dns07 * log(res1$ssd2))
  LL10.i <- - sum(dns10 * log(res1$ssd2))
  LL0710.i <- - sum(dns0710 * log(res1$ssd2))
  data_out[i, ] <- c(Ea.i, LL07.i, LL10.i, LL0710.i)

  # Store IPM results in lists  
  ipm_out[[i]] <- ipm1
  res_out[[i]] <- res1
  
}

head(data_out)
str(ipm_out[[1]])
str(res_out[[1]])

##### IPMs for WA and CA #####


### Modified IPM: WA temperatures, E = 0.65
ipm1 <- bigmatrix(n = binN, params = paramsWA)
res1 <- popF(ipm1, binSize)
res1[1:8]

### Base IPM: E = 0.65
ipm2 <- bigmatrix(n = binN, params = paramsCA)
res2 <- popF(ipm2, binSize)

### Modified IPM: Optimized E = 0.655
ipm3 <- bigmatrix(n = binN, params = paramsOptim)
res3 <- popF(ipm3, binSize)
res3[1:4]

##### CALCULATE LL #####
# Use 2010 data
area.h <- hist(hist0710$area, breaks = seq(min.size, max.size, binSize))
area.h
cnt <- area.h$counts 
dns <- area.h$density
sum(cnt) 
sum(dns); length(dns)
sum(dns * 0.05)

# using density
mod1LL <- - sum(dns * log(res1$ssd2)); mod1LL # LL for empirical data is -1.2

### Plot LL's
head(data_out)
names(data_out) <- c("Ea", "SizeDistribution_2007", 
                     "SizeDistribution_2010", "SizeDistribution_bothYears")

logLikLong <- data_out %>% gather(key = Year, value = LogLik, 
                                  SizeDistribution_2007:SizeDistribution_bothYears) %>% 
  mutate(LogLik = as.numeric(LogLik), Ea = as.numeric(Ea))

str(logLikLong)

##### PLOT LL RESULTS #####
bothYears <- logLikLong %>% filter(Year == "SizeDistribution_bothYears") %>%
  filter(LogLik < 50)

ggplot(bothYears, aes(Ea, LogLik)) +
  geom_point(size = 0.5) + geom_line(color = 'red') +
  ylab("Log likelihood") + xlab("Activation energy (Ea)") +
  geom_point(data = bothYears[bothYears$Ea == 0.65, ],
             aes(Ea, LogLik), color = "blue")

ggsave("figs/fecundityParams_Ea_LogLik_1.pdf",
       height = 3.5, width = 3.5)

##### MINIMIZE LOG-LIKELIHOOD, SELECT EA #####

# Get 5 lowest LL values for each Year, then take the median value
llDF <- logLikLong %>% group_by(Year) %>% 
  arrange(LogLik) %>% slice(1:5)

medianEa <- llDF %>% ungroup() %>% summarise(median = median(Ea))
medianEa

# This is the lowest LogLik (i.e. optimized)
llDF %>% ungroup() %>% arrange(LogLik) %>% slice(1:1)

medianEa_position <- which(data_out$Ea == as.numeric(medianEa))
