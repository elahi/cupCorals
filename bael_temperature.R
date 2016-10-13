###############################################################
# Author: Robin Elahi
# Date: 150611
# Calculating seawater temperatures for each site-year combo
###############################################################

##### NOTES #####
# Calculate the average temperature for the decade prior to 
# the last year of coral observations
# Shady Cove - 1962-1972; 2000-2010
# HMS - 1968-1978

### Why decade?
# Fadlallah estimates that to reach large sizes, it takes 5.6 - 10.6 years

##### LOAD PACKAGES, DATA #####
library(dplyr)
library(ggplot2)
#theme_set(theme_classic(base_size=12))

theme_set(theme_classic(base_size = 12) + 
            theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')))

source("./R/multiplotF.R")

##### WASHINGTON ANALYSIS #####
### At Shady Cove, I have temperature data every half hour at 15m depth
# from 2007-2011 - these were averaged to get daily means
#These data will be correlated with a longer-time series of daily temperature
# records at Race Rocks lighthouse, ~ 60km away

### Race Rocks daily SST ###
rr <- read.table("./data/racerockday.txt", 
                 header=FALSE, na.strings="99.9", skip = 3)
names(rr) <- c("year", "month", "day", "salinity", "tempC")
rr$dateR <- with(rr, as.Date(ISOdate(year, month, day)))

qplot(dateR, tempC, data = rr[rr$year > 1960, ])
rr$julien <- as.numeric(format(rr$dateR, "%j"))
qplot(julien, tempC, data = rr, color = year)

### Shady Cove daily mean SST ###
sc <- read.csv("./data/SC_temp daily means.csv", 
               header=TRUE)

sc$dateR <- with(sc, as.Date(sc$date, "%m/%d/%y"))

# extract year, month, day, julien day
sc$year <- as.numeric(format(sc$dateR, "%Y"))
sc$month <- as.numeric(format(sc$dateR, "%m"))
sc$day <- as.numeric(format(sc$dateR, "%d"))
sc$julien <- as.numeric(format(sc$dateR, "%j"))

# join two dataframes by dateR
corrDat <- inner_join(rr, sc, by = "dateR")
dim(corrDat)
head(corrDat)

### Examine correlation between two sites ###
ggplot(data = corrDat, aes(tempC, temp.mean)) +
  geom_point(alpha = 0.5) + geom_smooth(method = "lm") + 
  geom_abline(intercept = 0, slope = 1, color = "red", lwd = 1.5, lty = 2) + 
  xlab("Race Rocks temperature") + 
  ylab("Shady Cove temperature")

# compare variances
with(corrDat, sd(temp.mean, na.rm = TRUE))
with(corrDat, sd(tempC, na.rm = TRUE))

# Between 6 and 9 degrees there is approximately a 1:1 relationship
# However, at higher temperatures (ie summer), Shady Cove temps are underpredicted
# by the Race Rocks data. This is likely due to stratification in the water column, 
# because Race Rocks temperature was taken at the surface. 

# Regression
lm1 <- lm(temp.mean ~ tempC, data = corrDat)
summary(lm1)

# Use regression to predict SC temperatures for 10-year chunks
# 1962-1972, 2000-2010
dim(rr)
rr$tempSC <- 1.3976 + 0.8105 * rr$tempC
tail(rr)

### Use daily temperatures
summary(rr)

##' The last census dates:
##' past = January 1972 
##' present = December 2010
##' I want temperatures for a 10 year period prior to the last census date, so
##' past = 1971
##' present = 2010

# calculate average temps for estimated temperatures at Shady Cove
pastSC <- rr %>% filter(year > 1961 & year < 1972)
presSC <- rr %>% filter(year > 2000 & year < 2011)

pastSC %>% summarise(meanTemp = mean(tempSC, na.rm = TRUE), 
            sdTemp = sd(tempSC, na.rm = TRUE))

presSC %>% summarise(meanTemp = mean(tempSC, na.rm = TRUE), 
            sdTemp = sd(tempSC, na.rm = TRUE))

pastSC$era <- "1962-1971"
presSC$era <- "2001-2010"
longSC <- rbind(pastSC, presSC)

longSC %>% distinct(era, year)

plot_temp_boxplot <- ggplot(data = longSC, aes(era, tempSC)) +
  geom_violin(fill = "gray") + 
  geom_boxplot(width = 0.3, notch = TRUE, color = "black") + 
  ylab(expression(paste("Daily temperature (", degree, "C)"))) + 
  xlab("")

### Use monthly temperatures
# Calculate monthly averages for predicted SC temperature
sc_monthly <- rr %>% group_by(year, month) %>%
  summarise(monthC = mean(tempSC, na.rm = TRUE)) %>%
  ungroup() 
head(sc_monthly)
sc_monthly %>% filter(year > 1962 & year < 2010) %>%
  qplot(x = year, y = monthC, data = .) + geom_smooth()

pastSC_monthly <- sc_monthly %>% filter(year > 1962 & year < 1972) %>%
  summarise(meanTemp = mean(monthC, na.rm = TRUE), 
            sdTemp = sd(monthC, na.rm = TRUE))

presSC_monthly <- sc_monthly %>% filter(year > 2002 & year < 2010) %>%
  summarise(meanTemp = mean(monthC, na.rm = TRUE), 
            sdTemp = sd(monthC, na.rm = TRUE))

pastSC_monthly; presSC_monthly

### Use annual temperatures
# Calculate annual averages for predicted SC temperature
sc_yearly <- rr %>% group_by(year) %>%
  summarise(tempC = mean(tempSC, na.rm = TRUE)) %>%
  ungroup() 
head(sc_yearly)

# plot annual trend for study duration
sc_yearly %>% filter(year > 1962 & year < 2010) %>%
  qplot(x = year, y = tempC, data = .) + geom_smooth()

sc_yearly %>% filter(year > 1962 & year < 1972) %>%
  summarise(meanTemp = mean(tempC, na.rm = TRUE), 
            sdTemp = sd(tempC, na.rm = TRUE))

sc_yearly %>% filter(year > 2002 & year < 2010) %>%
  summarise(meanTemp = mean(tempC, na.rm = TRUE), 
            sdTemp = sd(tempC, na.rm = TRUE))

### Compare with original temperature data at Shady Cove
glimpse(sc)
qplot(dateR, temp.mean, data = sc)
scM <- sc %>% group_by(year, month) %>%
  summarise(monthC = mean(temp.mean, na.rm = TRUE)) %>%
  ungroup() 
scM
scM$dateR <- with(scM, ISOdate(year, month, day = 15))
qplot(dateR, monthC, data = scM, geom = "line")

# original data at SC
sc %>% filter(year > 2007 & year < 2011) %>%
  summarise(tempC = mean(temp.mean, na.rm = TRUE))
# estimated data at SC
rr %>% filter(year > 2007 & year < 2011) %>%
  summarise(tempC = mean(tempSC, na.rm = TRUE))
# Very close

##### WASHINGTON TEMPERATURE PLOTS #####
ULClabel <- theme(plot.title = element_text(hjust = -0.1, vjust = 1, 
                                            size = rel(1.2)))

### Annual temperatures - time-series
plot_temp_annual <- sc_yearly %>% filter(year > 1961 & year < 2011) %>%
  ggplot(data = ., aes(year, tempC)) + 
  geom_point() + geom_smooth(color = "black") + 
  ylab(expression(paste("Annual temperature (", degree, "C)"))) + 
  xlab("") + 
  ULClabel + labs(title = "A") 
  
### Daily temperatures - era comparison
plot_temp_daily <- ggplot(data = longSC, aes(era, tempSC)) +
  geom_violin(fill = "gray") + 
  geom_boxplot(width = 0.3, notch = TRUE, color = "black") + 
  ylab(expression(paste("Daily temperature (", degree, "C)"))) + 
  xlab("") + 
  ULClabel + labs(title = "B")

### save as pdf
# pdf("./figs/sc_temp_plot.pdf", width = 7, height = 3.5)
# multiplot(plot_temp_annual, plot_temp_daily,cols = 2)
# dev.off()

##### CALIFORNIA ANALYSIS #####
### Hopkins Marine Station handheld SST from caretaker's lodge
# Single measurement daily
# Get corrected temperature data (up to 2004)
hms_corr <- read.table("./data/HMStemp.corrected.txt", 
                       skip = 13, header = TRUE, na.strings = "NaN")
head(hms_corr)
hms_corr$dateR <- as.Date(with(hms_corr, ISOdate(year, month, day)))
str(hms_corr)
hms <- hms_corr %>% select(dateR, new) %>% rename(SST = new)

# get month and year to calculate averages
hms$month <- strftime(hms$dateR, "%m")
hms$year <- strftime(hms$dateR, "%Y")
hms$mo_yr <- with(hms, paste(month, year, sep = "_"))
hms$SST <- as.numeric(hms$SST)
str(hms)

# get data relevant to Fadlallah's data
fadTemp <- filter(hms, year > 1968 & year < 1978)
fadTemp$study <- "Fadlallah"
head(fadTemp)

##### GET TEMPERATURES FOR IPMS #####

sc_meanTemp_past <- pastSC %>% summarise(meanTemp = mean(tempSC, na.rm = TRUE))

sc_meanTemp_pres <- presSC %>% summarise(meanTemp = mean(tempSC, na.rm = TRUE))

hms_meanTemp <- fadTemp %>% summarise(meanTemp = mean(SST, na.rm = TRUE))

# get data relevant to Fadlallah's embryo-size function
# Use years 1976-1980 - these were the years during which Fadlallah did fieldwork, based on the 1983 paper - he does not say when the embryo study was done, only stating that it was conducted over an 18 month period

embryoTemp <- filter(hms, year > 1975 & year < 1980) %>%
  summarise(meanTemp = mean(SST, na.rm = TRUE)) # 4 year period, 

sc_meanTemp_pres; sc_meanTemp_past; hms_meanTemp; embryoTemp
