#################################################
# Author: Robin Elahi
# Date: 151016

# Coral densities for 1969 and 2007
# Figure x
#################################################

library(ggplot2)
theme_set(theme_classic(base_size = 10))
library(reshape2)
library(nlme)
library(lme4)
library(dplyr)


rm(list=ls(all=TRUE)) # removes all previous material from R's memory

source("./R/graphicalParams.R")
source("./R/summarizeData.R")
source("./R/multiplotF.R")

#################################################
# DENSITY
#################################################
dat <- read.csv("./data/bael_densityMaster.csv", header=TRUE, na.strings="NA")
names(dat)

### Remove D3 because this subquadrat was not sampled after 3 years
dat2 <- dat[dat$quad != "D3", ]

### Quick plots of data at the subquadrat and quadrat scale
# Note that this analysis treats subquads as replicates, not appropriate
boxplot(densM ~ year, data = dat2)
ggplot(data = dat2, aes(as.factor(year), densM)) + geom_violin()
tail(dat2)

facet1 <- facet_wrap(~ era, scales = "free_x")
my_label <- expression(paste("Density (no. ", m^-2, ")"))
unique(dat2$quad)

ULClabel <- theme(plot.title = element_text(hjust = -0.2, vjust = 1, 
                                            size = rel(1.2)))
head(dat2)

library("grid") # needed for unit call in ggplot2
densityPlot <- ggplot(dat2, aes(x = year, y = densM, by = quad)) + 
  xlab("Year") + ylab(my_label) + 
  geom_point(size = 2, alpha = 0.5, aes(color = era, shape = era)) + 
  theme(legend.position = "none") +
  geom_line(size = 0.2, alpha = 0.3, aes(color = era)) + facet1 + 
  theme(strip.background = element_blank(), 
        strip.text.x = element_blank()) + 
  scale_colour_manual(breaks = c("historic", "modern"), 
                      values = c("darkgray", "black")) +
  scale_shape_manual(breaks = c("historic", "modern"), 
                     values = c(18, 20)) + 
  theme(panel.margin = unit(1.5, "lines"))

p2 <- densityPlot + ULClabel + labs(title = "B                    Shady Cove") 
p2

###########################
# Models for density
# lme4
# test whether density is different among eras
# compare 1969 with 2007 for simplicity
head(dat)
tail(dat) # use transect as a random effect, 
# b/c responses are quad densities
lmerDat <- dat %>% filter(year == "2007" | year == "1969")
summary(lmerDat)
unique(lmerDat$quad)
with(lmerDat, table(era, quad))

# need to calculate average density for historical quads
meanDensData <- dat %>% filter(year == "2007" | year == "1969") %>%
  group_by(quad1, tran, era, year) %>% 
  summarise(meanDens = mean(densM, na.rm = TRUE))
meanDensData

# No matter whether I use raw subquad densities, or mean quad densities
# for historic study, there is no difference between eras
Cand.mod <- list()
Cand.mod[[1]] <- lmer(densM ~ era + (1|tran),
                      REML = FALSE, data = lmerDat)
Cand.mod[[2]] <- lmer(densM ~ 1 + (1|tran),
                      REML = FALSE, data = lmerDat)		

# OR
Cand.mod <- list()
Cand.mod[[1]] <- lmer(meanDens ~ era + (1|tran),
                      REML = FALSE, data = meanDensData)
Cand.mod[[2]] <- lmer(meanDens ~ 1 + (1|tran),
                      REML = FALSE, data = meanDensData)		

# Present the more conservative statistical results in manuscript 
# (meanDensData)

#generate AICc table with names

mod_text <- c("Era", "Null model")							
mod.aicctab <- aictab(cand.set= Cand.mod, modnames= mod_text, 
                      sort=TRUE, second.ord=TRUE) # second.ord =TRUE means AICc is used 
print(mod.aicctab, digits=2, LL=TRUE)
write.csv(mod.aicctab, "./output/densityAIC.csv")

summary(Cand.mod[[1]])
bestMod <- Cand.mod[[1]]

# check normality and homogeneity of variances
par(mfrow = c(1,2))
qqnorm(resid(bestMod))
qqline(resid(bestMod))
plot(resid(bestMod) ~ fitted(bestMod)); abline(h=0)

# Plot this nicely
meanDensPlot <- ggplot(data = meanDensData, aes(as.factor(year), meanDens)) + 
  xlab("") + ylab(my_label) + 
  geom_jitter(size = 3, alpha = 1, 
              position = position_jitter(width = 0.05), 
              aes(color = era, shape = era)) + 
  theme(legend.position = "none") +
  scale_colour_manual(breaks = c("historic", "modern"), 
                      values = c("darkgray", "black")) +
  scale_shape_manual(breaks = c("historic", "modern"), 
                     values = c(18, 20)) 

p3 <- meanDensPlot + ULClabel + labs(title = "B              Shady Cove") 
p3

#################################################
# VIOLIN PLOTS
#################################################

dat <- read.csv("./data/bael_histoData.csv", header=TRUE, na.strings="NA")

### Data preparation
### Remove everything unnecessary for histograms of initial data
ini.dat <- droplevels(dat[dat$ini.notes != "angle" & dat$ini.notes != "fuzzy" &
                            dat$ini.notes != "gone" & dat$ini.notes != "nv" &
                            dat$ini.notes != "tentacles", ])

ini.dat <- droplevels(ini.dat[complete.cases(ini.dat$ini.area), ]) # drop NAs

head(ini.dat)

# create new column for Site_era
ini.dat$siteEra <- as.factor(paste(ini.dat$site, ini.dat$time, sep = "_"))
unique(ini.dat$siteEra)
# Reorder levels for plotting
ini.dat$siteEra <- factor(ini.dat$siteEra, rev(c("SC_past", "SC_present", 
                                             "PG_present", "ON_present")))

# violin plot
my_label <- expression(paste("Density (no. ", m^-2, ")"))

ULClabel <- theme(plot.title = element_text(hjust = -0.1, vjust = 1, 
                                            size = rel(1.2)))

sizePlot <- ggplot(ini.dat,  aes(x = siteEra, y = ini.area, fill = time)) +
  ylab(expression(paste("Size (", cm^2, ")"))) + 
  geom_violin(position = position_dodge(1)) + 
  geom_boxplot(width = 0.3, notch = TRUE, color = "black") + 
  scale_x_discrete("", labels = c("SC_past" = "Shady Cove\n1969\nn=164", 
                              "SC_present" = "Shady Cove\n2007\nn=144", 
                              "PG_present" = "Point George\n2007\nn=159",
                              "ON_present" = "O'Neal\n2007\nn=315")) +
  coord_flip() + theme(legend.position = "none") + 
  scale_fill_manual(values = c("darkgray", "white"))
  
p1 <- sizePlot + labs(title = "A") + ULClabel  

##########################################################
# Multi-panel plots
##########################################################

# with subquad densities
pdf("./figs/sizeDensityPlot.pdf", width = 7, height = 3.5)
multiplot(p1, p2, cols = 2)
dev.off()

# with mean quad densities
p3 <- meanDensPlot + ULClabel + 
  labs(title = "B") 

pdf("./figs/sizeMeanDensityPlot.pdf", width = 7, height = 3.5)
multiplot(p1, p3, cols = 2)
dev.off()
