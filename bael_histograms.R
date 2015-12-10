#################################################
# Author: Robin Elahi
# Date: 151209
# Coral histograms for 1969 and 2007
# Figure 1
#################################################

# rm(list=ls(all=TRUE)) 

##### LOAD PACKAGES, DATA #####
library(lme4)
library(ggplot2)
theme_set(theme_classic(base_size = 16))
library(AICcmodavg)
library(dplyr)

dat <- read.csv("./data/bael_histoData.csv", header=TRUE, na.strings="NA")
source("./R/graphicalParams.R")

dat$ini.areaLN <- log(dat$ini.area)
dat$fin.areaLN <- log(dat$fin.area)

##### DATA PREPARATION #####
### Remove everything unnecessary for histograms of initial data
ini.dat <- droplevels(dat[dat$ini.notes != "angle" & dat$ini.notes != "fuzzy" &
                         dat$ini.notes != "gone" & dat$ini.notes != "nv" &
                         dat$ini.notes != "tentacles", ])

ini.dat <- droplevels(ini.dat[complete.cases(ini.dat$ini.area), ]) # drop NAs
write.csv(ini.dat, "./output/ini.dat.csv")
summary(ini.dat)

# Past data
past <- ini.dat[ini.dat$time == "past", ]
summary(past) # 164
(with(past, table(quad, ini.date), ini.area))
(with(past, table(quad, fin.date), fin.area))

quantile(past$ini.area, c(0.95, 0.99))

qplot(quad1, ini.area, data = ini.dat, geom = "boxplot") + coord_flip()

present <- ini.dat[ini.dat$time == "present", ]
(with(present, table(quad, ini.date), ini.area))

ON <- subset(present, site=="ON")
summary(ON) # n = 315
PG <- subset(present, site=="PG")
summary(PG) # n = 159
SC <- subset(present, site=="SC")
summary(SC) # n = 144
quantile(SC$ini.area, c(0.95, 0.99))
scDat <- subset(ini.dat, site == "SC")
dim(scDat)

##### 4-PANEL HISTOGRAM - LINEAR #####
# retrieve proportions, rather than frequencies
past.h <- hist(past$ini.area, breaks = seq(0, 1.8, 0.1))
past.h$density <- past.h$counts/sum(past.h$counts)
plot(past.h, freq=F)

ON.h <- hist(ON$ini.area, breaks = seq(0, 1.8, 0.1))
ON.h$density <- ON.h$counts/sum(ON.h$counts)

PG.h <- hist(PG$ini.area, breaks = seq(0, 1.8, 0.1))
PG.h$density <- PG.h$counts/sum(PG.h$counts)

SC.h <- hist(SC$ini.area, breaks = seq(0, 1.8, 0.1))
SC.h$density <- SC.h$counts/sum(SC.h$counts)

# final plot
pdf("./figs/histo_4panel.pdf", 7, 7)
set_graph_pars(ptype = "panel4")

par(mar=c(1,5,1,1), oma=c(1,1,0,0), cex.axis = 1.2)
plot(past.h, freq=F, main="", xlab="", ylab="", col="darkgray", 
     xlim=c(0,1.7), ylim=c(0, 0.21),  las=1); box()
add_panel_label(ltype = "a")
leg.txtA <- c("Shady Cove", "1969", "n = 164")
legend("topright", leg.txtA, bty="n", adj=c(0, 0))
#abline(v = quantile(past$ini.area, 0.95), lty=3, lwd=2)
abline(v = mean(past$ini.area), lty=2, lwd=2)

plot(SC.h, freq=F, main="", xlab="", ylab="", col="white", 
     xlim=c(0,1.7), ylim=c(0, 0.21),  las=1); box()
add_panel_label(ltype = "b")
leg.txtB <- c("Shady Cove", "2007", "n = 144")
legend("topright", leg.txtB, bty="n", adj=c(0, 0))
#abline(v = quantile(SC$ini.area, 0.95), lty=3, lwd=2)
abline(v = mean(SC$ini.area), lty=2, lwd=2)

plot(ON.h, freq=F, main="", xlab="", ylab="", col="white", 
     xlim=c(0,1.7), ylim=c(0, 0.21),  las=1); box()
add_panel_label(ltype = "c")
leg.txtC <- c("O'Neal", "2007", "n = 315")
legend("topright", leg.txtC, bty="n", adj=c(0, 0))
#abline(v = quantile(ON$ini.area, 0.95), lty=3, lwd=2)
abline(v = mean(ON$ini.area), lty=2, lwd=2)

plot(PG.h, freq=F, main="", xlab="", ylab="", col="white", 
     xlim=c(0,1.7), ylim=c(0, 0.21),  las=1); box()
add_panel_label(ltype = "d")
leg.txtD <- c("Point George", "2007", "n = 159")
legend("topright", leg.txtD, bty="n", adj=c(0, 0))
#abline(v = quantile(PG$ini.area, 0.95), lty=3, lwd=2)
abline(v = mean(PG$ini.area), lty=2, lwd=2)

xlab1 <- expression(paste("Size (", cm^2, ")"))
mtext(xlab1, side = 1, line = 0, outer = TRUE, cex = 1.2)
mtext("Proportion", side = 2, line = -1, outer = TRUE, cex = 1.2)

dev.off()

##### LMER BODY SIZE - FULL DATASET #####
# test whether mean size is differrent among eras
lmerDat <- ini.dat 

Cand.mod <- list()
Cand.mod[[1]] <- lmer(ini.area ~ time +  (1|site) + (1|quad),
							REML = FALSE, data = lmerDat)
Cand.mod[[2]] <- lmer(ini.area ~ 1 + (1|site) + (1|quad),
							REML = FALSE, data = lmerDat)							

#generate AICc table with names
mod_text <- c("Era", "Null model")							
mod.aicctab <- aictab(cand.set= Cand.mod, modnames= mod_text, 
                      sort=TRUE, second.ord=TRUE) # second.ord =TRUE means AICc is used 
print(mod.aicctab, digits=2, LL=TRUE)
#write.csv(mod.aicctab, "sizeFullAIC.csv")

# check normality and homogeneity of variances
bestMod <- lmer(ini.area ~ time + (1|site) +  (1|quad),
							REML = TRUE, data=lmerDat)	

# check normality and homogeneity of variances
par(mfrow = c(1,2))
qqnorm(resid(bestMod)); qqline(resid(bestMod))
plot(resid(bestMod) ~ fitted(bestMod)); abline(h=0)

##### LMER BODY SIZE - TRUNCATED DATASET #####
# test whether mean size is differrent among eras
lmerDatTrunc <- ini.dat[ini.dat$ini.area <= 1, ]

Cand.mod <- list()
Cand.mod[[1]] <- lmer(ini.area ~ time +  (1|site) + (1|quad),
                      REML = FALSE, data = lmerDatTrunc)
Cand.mod[[2]] <- lmer(ini.area ~ 1 + (1|site) + (1|quad),
                      REML = FALSE, data = lmerDatTrunc)							

#generate AICc table with names
mod_text <- c("Era", "Null model")							
mod.aicctab <- aictab(cand.set= Cand.mod, modnames= mod_text, 
                      sort=TRUE, second.ord=TRUE) # second.ord =TRUE means AICc is used 
print(mod.aicctab, digits=2, LL=TRUE)
#write.csv(mod.aicctab, "sizeTruncAIC.csv")

# check normality and homogeneity of variances
bestModTrunc <- lmer(ini.area ~ time + (1|site) +  (1|quad),
                REML = TRUE, data = lmerDatTrunc)	

# check normality and homogeneity of variances
par(mfrow = c(1,2))
qqnorm(resid(bestModTrunc))
qqline(resid(bestModTrunc))
plot(resid(bestModTrunc) ~ fitted(bestModTrunc)); abline(h=0)

##### SUPPLEMENTAL FIGURE - RESIDUALS OF GROWTH MODELS #####
# plot residuals v fitted, and observe v fitted for supplemental

pdf("./figs/histo_4panel_resids.pdf", 7, 7)
set_graph_pars(ptype = "panel4")
plot(resid(bestMod) ~ fitted(bestMod), xlab = "", ylab = "Residuals", col = "black", cex.lab = 1.4, cex.axis = 1.2, las = 1)
abline(h=0, col = "darkgray", lty = 2, lwd = 2)
add_panel_label(ltype = "a")
title(main = "Full dataset", line = 1.5, cex.main = 1.3)

plot(resid(bestModTrunc) ~ fitted(bestModTrunc), xlab = "", ylab = "Residuals", col = "black", cex.lab = 1.4, cex.axis = 1.2, las = 1)
abline(h=0, col = "darkgray", lty = 2, lwd = 2)
add_panel_label(ltype = "b")
title(main = "Truncated dataset", line = 1.5, cex.main = 1.3)

plot(lmerDat$ini.area ~ fitted(bestMod), xlab = "Fitted", ylab = "Observed", col = "black", cex.lab = 1.4, cex.axis = 1.2, las = 1)
abline(0, 1, col = "darkgray", lty = 2, lwd = 2)
add_panel_label(ltype = "c")

plot(lmerDatTrunc$ini.area ~ fitted(bestModTrunc), xlab = "Fitted", ylab = "Observed", col = "black", cex.lab = 1.4, cex.axis = 1.2, las = 1)
abline(0, 1, col = "darkgray", lty = 2, lwd = 2)
add_panel_label(ltype = "d")

dev.off()

##### KOLMOGOROV-SMIRNOV TESTS OF SIZE-FREQUENCY DISTRIBUTIONS #####
ks.test(past$ini.areaLN, SC$ini.areaLN)
# data:  past$ini.area.cm2 and SC$ini.area.cm2
# D = 0.29, p-value = 1.125e-06
# alternative hypothesis: two-sided"

ks.test(past$ini.areaLN, ON$ini.areaLN)
# data:  past$ini.area.cm2 and ON$ini.area.cm2
# D = 0.43, p-value < 2.2e-16
# alternative hypothesis: two-sided

ks.test(past$ini.areaLN, PG$ini.areaLN)
# data:  past$ini.area.cm2 and PG$ini.area.cm2
# D = 0.40, p-value = 3.318e-13
# alternative hypothesis: two-sided

##### ADDITIONAL PLOTS FOR INSPECTION #####
# Plot past histos by quad1

head(arrange(past, desc(ini.area)))
head(past)
head(present)

ggplot(past,  aes(x = ini.area)) +
  geom_density(aes(color = time)) +
  facet_wrap(~ quad1, scales = "fixed", nrow = 4) + 
  xlab("Size (cm2)") + ylab("Probability density") +
  geom_vline(xintercept = 1)

ggplot(past,  aes(x = ini.area)) +
  geom_density(aes(fill = quad1), alpha = 0.5) +
  facet_wrap(~ quad1, scales = "fixed", nrow = 4) + 
  xlab("Size (cm2)") + ylab("Probability density") +
  geom_vline(xintercept = 1)

ggplot(ini.dat,  aes(x = ini.area)) +
  geom_density(aes(fill = time), alpha = 0.5) +
  facet_wrap(~ time + site, scales = "fixed", nrow = 4) + 
  xlab("Size (cm2)") + ylab("Probability density") +
  geom_vline(xintercept = 1)

ggplot(scDat,  aes(x = ini.area)) +
  geom_density(aes(fill = time), alpha = 0.5) +
  # facet_wrap(~ time, scales = "fixed", nrow = 4) + 
  xlab("Size (cm2)") + ylab("Probability density") +
  geom_vline(xintercept = 1)

ggplot(past,  aes(x = ini.area)) +
  geom_histogram(binwidth = 0.1, aes(fill = quad)) +
  facet_wrap(~ quad1, scales = "fixed", nrow = 4) + 
  xlab("Size (cm2)") + ylab("Count") +
  geom_vline(xintercept = 1)

##### DENSITY PLOTS #####
head(scDat)

ggplot(data = scDat, aes(x = ini.area, color = time)) + 
  geom_density()

ggplot(data = scDat, aes(x = time, y = ini.area, color = time)) + 
  geom_violin()


