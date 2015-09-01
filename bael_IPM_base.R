#################################################
# Author: Robin Elahi
# Date: 150825

# Base IPM parameterized by modern data
# Figure 2
#################################################

rm(list=ls(all=TRUE)) # removes all previous material from R's memory
source("./R/baelParamsWA.R")
source("./R/baelParamsCA.R")
source("./R/ipmFunctions.R")
paramsWA
paramsCA

##########################################
##########################################
# MAX AND MIN SIZES
names(hist0710)
range(hist0710$area)

min.size <- .9*min(c(growthDat$size,growthDat$sizeNext), na.rm=T)
max.size <- 1.1*max(c(growthDat$size,growthDat$sizeNext), na.rm=T)

min.size <- .9*min(hist0710$area, na.rm=T)
max.size <- 1.1*max(hist0710$area, na.rm=T)

min.size; max.size

min.size <- 0.02
max.size <- 1.4

binSize <- 0.05
binN <- (max.size - min.size)/binSize
binN

# or, select number of bins
binN <- 50
binSize <- (max.size - min.size)/binN
binSize

##########################################
##########################################
### RUN IPMS
# Basic IPM, original parameters
ipm1 <- bigmatrix(n = binN, params = paramsWA)
res1 <- popF(ipm1, binSize)
res1[1:8]

# This gets the probability of the mean size of the SSD
res1$aboveMean <- with(res1, min(which(cumsum(ssd1) > meanSize)))
res1$belowMean <- with(res1, max(which(cumsum(ssd1) < meanSize)))
res1$meanSizePD <- with(res1, mean(c(ssd2[aboveMean], ssd2[belowMean])))
res1

# Check for eviction
plot(ipm1$meshpts, s.x(ipm1$meshpts, params), xlab = "Size", type = "l", ylab = "Survival probability", lwd = 12) # plot the survival model
points(ipm1$meshpts, apply(ipm1$P, 2, sum), col = "red", lwd = 3, cex = 0.5, pch = 19) # plot the column sums of the survival/growth matrix

# Basic IPM, modified embryo and recruitment parameters
ipm2 <- bigmatrix(n = binN, params = paramsCA)
res2 <- popF(ipm2, binSize)
res2
res2$aboveMean <- with(res2, min(which(cumsum(ssd1) > meanSize)))
res2$belowMean <- with(res2, max(which(cumsum(ssd1) < meanSize)))
res2$meanSizePD <- with(res2, mean(c(ssd2[aboveMean], ssd2[belowMean])))

### Plot for publication

pdf("./figs/ipm_histo_fit.pdf", 7, 3.5)
set_graph_pars(ptype = "panel2")

# Panel A
plot(mx.coral ~ areaR, data = ltDat, xlim = c(0, 1.4), ylim = c(-0.025, 40),
     xlab = xlab2, ylab = "Embryo number", las = 1, type = "n")
abline(mxRegWA, lwd=2, lty=1, col = "darkgray")
abline(mxRegCA, lwd=2, lty=1)
abline(a = 0, b = 0, lty=3, lwd=2, col="darkgray")
points(mx.coral ~ areaR, data = ltDat)
add_panel_label(ltype = "a")
legend("bottomright", "n = 12", cex = 1.1, bty = "n", adj = c(0,-2))

leg.txt <- c("Original", "Modified")
legend("topleft", leg.txt, lwd = 2, bty = "n", 
	col = c("black", "darkgray"), lty = 1, 
	cex = 1, text.col = c("black", "darkgray"))

# Panel B	
hist(hist0710$area, breaks = seq(0.0, 1.4, 0.05), freq = FALSE, 
	xlab = xlab2, ylab="Probability density", col = "gray87", main = "", 
	ylim = c(-0.025, 2), xlim = c(0, 1.4), border = "gray90", las = 1) 
box()
points(ipm1$meshpts , res1$ssd2, type="l", lty=1, lwd = 2, 
		col = "darkgray")
points(ipm2$meshpts , res2$ssd2, type="l", lty=1, lwd = 2, 
		col = "black")

add_panel_label(ltype = "b")

# arrows(res1$max99, 0.6, res1$max99, 0.3, col = "darkgray", 
# 		length = 0.1, lwd = 1.5, angle = 20)
# arrows(res2$max99, 0.6, res2$max99, 0.3, col = "black", 
# 		length = 0.1, lwd = 1.5, angle = 20)		

# Mean size from empirical data
abline(v = mean(hist0710$area), lty=2, lwd=1)

# Mean size from IPM stable size distribution
points(res1$meanSize, -0.05, col = "darkgray", pch = 16, cex = 0.7)
points(res2$meanSize, -0.05, col = "black", pch = 16, cex = 0.7)

# 99% maximum size from IPM stable size distribution
points(res1$max99, -0.05, col = "darkgray", pch = 17, cex = 0.7)
points(res2$max99, -0.05, col = "black", pch = 17, cex = 0.7)

legend("topright", leg.txt, lwd = 2, bty = "n", 
	col = c("black", "darkgray"), lty = 1, 
	cex = 1, text.col = c("black", "darkgray"))

dev.off()
