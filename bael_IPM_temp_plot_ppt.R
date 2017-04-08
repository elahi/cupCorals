#################################################
# Author: Robin Elahi
# Date: 151208
# Temperature manipulation of IPM - PLOTTING
# Figure 6
#################################################

source("bael_IPM_temp.R")
source("bael_temperature.R")
source("./R/multiplotF.R")
source("R/get_histo_ipm_data.R")


##### PLOT SIMULATED DATA #####
# Plotting details
label1 <- expression(paste("Maximum size (", cm^2, ")"))
tempLab <- expression(paste("Temperature (", degree, "C)"))
regression <- geom_smooth(method = "lm", se = FALSE, alpha = 0.5, 
                          size = 0.4)
ULClabel <- theme(plot.title = element_text(hjust = -0.07, vjust = 1, 
                                            size = rel(1.5)))
theme_set(theme_classic(base_size = 12))

# Max size plot
max99Plot <- ggplot(data = simDat, aes((Kelvin-273.15), maxSize99, linetype = Ea)) +
  xlab(tempLab) + ylab(label1) +
  geom_point(alpha = 0.5, size = 0, color = "white") + 
  geom_smooth(se = FALSE, size = 0.7, color = "black") + 
  theme(legend.justification = "center", legend.position = c(0.9, 0.7)) +
  theme(legend.text = element_text(size = 10)) + 
  theme(legend.title = element_text(size = 10)) + 
  scale_linetype_discrete(name = "Activation\nenergy") + 
  # scale_colour_grey(start = 0.8, end = 0.2) + 
  guides(linetype = guide_legend(reverse=TRUE)) + 
  # guides(color = guide_legend(reverse=TRUE)) + 
  coord_cartesian(ylim = c(0.9, 1.7)) 

max95Plot <- ggplot(data = simDat, aes((Kelvin-273.15), maxSize95, linetype = Ea)) +
  xlab(tempLab) + ylab(label1) +
  geom_point(alpha = 0.5, size = 0, color = "white") + 
  geom_smooth(se = FALSE, size = 0.7, color = "black") + 
  theme(legend.justification = "center", legend.position = c(0.9, 0.7)) +
  theme(legend.text = element_text(size = 10)) + 
  theme(legend.title = element_text(size = 10)) + 
  scale_linetype_discrete(name = "Activation\nenergy") + 
  # scale_colour_grey(start = 0.8, end = 0.2) + 
  guides(linetype = guide_legend(reverse=TRUE)) + 
  # guides(color = guide_legend(reverse=TRUE)) + 
  coord_cartesian(ylim = c(0.5, 1.7)) 
max95Plot

##### OBSERVED MAXIMUM SIZE: PRESENT #####
# Size frequencies in 2007 + 2010
sizeSummary_modern <- hist0710 %>%
  summarise(maxSize = max(area), 
            max99 = quantile(area, 0.99), 
            max95 = quantile(area, 0.95),
            median = median(area), 
            n = length(area))

##### OBSERVED MAXIMUM SIZE: PAST #####
# Size frequences in 1969 + 1972
# Past size data in histoData.csv:
dat <- read.csv("./data/bael_histoData.csv", header=TRUE, na.strings="NA")

# Initial data
ini.dat <- droplevels(dat[dat$ini.notes != "angle" & dat$ini.notes != "fuzzy" &
                            dat$ini.notes != "gone" & dat$ini.notes != "nv" &
                            dat$ini.notes != "tentacles", ])

ini.dat <- droplevels(ini.dat[complete.cases(ini.dat$ini.area), ]) # drop NAs

# Final data
fin.dat <- droplevels(dat[dat$fin.notes != "angle" & dat$fin.notes != "fuzzy" &
                            dat$fin.notes != "gone" & dat$fin.notes != "nv" &
                            dat$fin.notes != "tentacles" & 
                            dat$fin.notes != "algae" &
                            dat$fin.notes != "dead" & 
                            dat$fin.notes != "overgrown", ])

fin.dat <- droplevels(fin.dat[complete.cases(fin.dat$fin.area), ]) # drop NAs

ini.sc <- subset(ini.dat, site=="SC")
fin.sc <- subset(fin.dat, site == "SC")

# Combine past and present SC data 
iniDF <- ini.sc %>% select(time, ini.area) %>% dplyr::rename(area = ini.area)
finDF <- fin.sc %>% select(time, fin.area) %>% dplyr::rename(area = fin.area)
histDF <- rbind(iniDF, finDF)

histDF %>% group_by(time) %>%
  summarise(maxSize = max(area), 
            max99 = quantile(area, 0.99), 
            max95 = quantile(area, 0.95), 
            median = median(area), 
            n = length(area))

sizeSummary_historic <- histDF %>% filter(time == "past") %>%
  summarise(maxSize = max(area), 
            max99 = quantile(area, 0.99), 
            max95 = quantile(area, 0.95), 
            median = median(area), 
            n = length(area))

sizeObs <- rbind(sizeSummary_modern, sizeSummary_historic)
sizeObs$tempC <- c(modTemp-273.15, hisTemp-273.15)
sizeObs$era <- c("present", "past")
sizeObs

##### GET PREDICTED MAX SIZES  #####
# Predicted modern size (empirical growth functions)

# what is the predicted size at the modern temp?
modTemp
maxSizeModPred <- simDat[simDat$Kelvin == 282.40 &
                           simDat$Ea == 0.6, ]$maxSize99
maxSizeModPred

# what is the predicted size, using the IPM with the historical
# growth slope?
source("bael_IPM_historicGrowth.R")

res0$max99 # should be the same as maxSizeModPred
maxSizeModPred <- res0$max99

# truncated growth curve (historic data)
res1$max99
# full growth curve (historic data)
res2$max99

# use truncated growth curve
maxSizeHistPred <- res1$max99

# Create data frame with predicted size at modern temp (above)
# and predicted size at historic temperature, using the empirical
# growth function for past data

maxSizePred <- data.frame(
  era = c("Historic", "Modern"),
  tempC = c(hisTemp-273.15, modTemp-273.15),
  size = c(maxSizeHistPred, maxSizeModPred)
)

maxSizePred

##### FINAL PLOTS  #####
presentObs <- sizeObs %>% filter(era == "present")
pastObs <- sizeObs %>% filter(era == "past")

# Upper edge is max size, lower edge is 99%ile size
rectPresent <- annotate("rect", 
                        xmin = (as.numeric(sc_meanTemp_pres) - 0.02), 
                        xmax = (as.numeric(sc_meanTemp_pres) + 0.02), 
                        ymin = presentObs$max99, 
                        ymax = presentObs$maxSize,
                        alpha = .2, fill = "black") 

rectPast <- annotate("rect", 
                     xmin = (as.numeric(sc_meanTemp_past) - 0.02), 
                     xmax = (as.numeric(sc_meanTemp_past) + 0.02), 
                     ymin = pastObs$max99, 
                     ymax = pastObs$maxSize,
                     alpha = .2, fill = "black")


max99Plot + rectPresent + rectPast + 
  geom_point(aes(tempC, size, linetype = NULL), 
             data = maxSizePred,
             size = 3, shape = c(23, 21), fill = 1) + 
  geom_point(aes(tempC, size, linetype = NULL), 
             data = maxSizePred,
             size = 1.5, shape = c(23, 21), fill = "darkgray")

ggsave("./figs/ipm_temp.pdf", width = 3.5, height = 3.5)

##### FINAL PLOTS - PPT  #####

# Max size plot
max99Plot <- ggplot(data = simDat, aes((Kelvin-273.15), maxSize99, linetype = Ea)) +
  xlab(tempLab) + ylab(label1) +
  geom_point(alpha = 0.5, size = 0, color = "white") + 
  geom_smooth(se = FALSE, size = 0.7, color = "black") + 
  theme(legend.justification = "center", legend.position = c(0.9, 0.7)) +
  theme(legend.text = element_text(size = 10)) + 
  theme(legend.title = element_text(size = 10)) + 
  scale_linetype_discrete(name = "Activation\nenergy") + 
  # scale_colour_grey(start = 0.8, end = 0.2) + 
  guides(linetype = guide_legend(reverse=TRUE)) + 
  # guides(color = guide_legend(reverse=TRUE)) + 
  coord_cartesian(ylim = c(0.9, 1.7)) 

# Max size plot - for ppt
simDat2 <- simDat %>% filter((Kelvin-273.15) < 9.26)

theme_set(theme_bw(base_size = 18))

max99Plot2 <- ggplot(data = simDat2, aes((Kelvin-273.15), maxSize99, color = Ea)) +
  xlab(tempLab) + ylab(label1) +
  geom_point(alpha = 0.5, size = 0, color = "white") + 
  geom_smooth(se = FALSE, size = 0.7) + 
  theme(legend.text = element_text(size = 14)) + 
  theme(legend.title = element_text(size = 14)) + 
  scale_shape_discrete(name = "Era") + 
  scale_color_discrete(name = "Activation\nenergy") + 
  guides(color = guide_legend(reverse=TRUE)) + 
  coord_cartesian(ylim = c(0.9, 1.65)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

max99Plot2

max99Plot2 + rectPresent + rectPast + 
  geom_point(aes(tempC, size, color = NULL, shape = era), 
             data = maxSizePred, size = 4, fill = 1)

ggsave("./figs/ipm_temp_ppt1.pdf", width = 7, height = 5)


max99Plot2_blank <- ggplot(data = simDat2, aes((Kelvin-273.15), maxSize99, color = Ea)) +
  xlab(tempLab) + ylab(label1) +
  geom_point(alpha = 0.5, size = 0, color = "white") + 
  #geom_smooth(se = FALSE, size = 0.7) + 
  theme(legend.text = element_text(size = 14)) + 
  theme(legend.title = element_text(size = 14)) + 
  scale_shape_discrete(name = "Era") + 
  scale_color_discrete(name = "Activation\nenergy") + 
  guides(color = guide_legend(reverse=TRUE)) + 
  coord_cartesian(ylim = c(0.9, 1.65)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

max99Plot2_blank + rectPresent + #rectPast + 
  geom_point(aes(tempC, size, color = NULL, shape = era), 
             data = maxSizePred, size = 4, fill = 1)

ggsave("./figs/ipm_temp_ppt2.pdf", width = 7, height = 5)

max99Plot2_blank + rectPresent + rectPast + 
  geom_point(aes(tempC, size, color = NULL, shape = era), 
             data = maxSizePred, size = 4, fill = 1)

ggsave("./figs/ipm_temp_ppt3.pdf", width = 7, height = 5)


### How different are these observed and predicted values?
sizeObs
maxSizeHistObs <- sizeObs$max99[1]
maxSizeModObs <- sizeObs$max99[2]

maxSizeHistObs  
maxSizeModObs
maxSizeHistPred  
maxSizeModPred

(maxSizeHistPred - maxSizeModPred)/maxSizeModPred

observedChange <- maxSizeHistObs - maxSizeModObs
observedChange/maxSizeHistObs

