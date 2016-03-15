#################################################
# Author: Robin Elahi
# Date: 151208
# Base IPM parameterized by modern data
# Figure 2
#################################################

rm(list=ls(all=TRUE)) 

##### LOAD PACKAGES ETC #####
source("./R/bael_params.R")
source("./R/ipmFunctions.R")

# Data from ipmData.csv:
dat <- read.csv("./data/bael_ipmData.csv", header=TRUE, na.strings="NA")
# Select the relevant columns
d <- dat[, which(names(dat) %in% 
                   c("quad", "date", "date.no", "coral.id", "area", 
                     "feret", "code", "sizeOK", "surv", "growth", "recruit"))]

dHIST <- droplevels(d[d$code != "angle" & d$code != "algae" &
                        d$code != "nv" & d$code != "dead", ])
dHIST <- droplevels(dHIST[complete.cases(dHIST$area), ]) # drop NAs

# include only 2007 and 2010: years for IPM
hist07 <- droplevels(dHIST[dHIST$date.no == 39426, ])
hist10 <- droplevels(dHIST[dHIST$date.no == 40515, ])
hist0710 <- rbind(hist07, hist10)
range(hist0710$area)


##### MAX AND MIN SIZES #####
# Will use slightly larger size range for IPM because
# it is otherwise artificially truncated
min.size <- 0.02
max.size <- 1.5

binSize <- 0.02
binN <- (max.size - min.size)/binSize
binN

##### RUN IPMS #####
### Base IPM, original parameters
ipm1 <- bigmatrix(n = binN, params = paramsWA)
res1 <- popF(ipm1, binSize)
res1[1:8]

# Check for eviction
plot(ipm1$meshpts, s.x(ipm1$meshpts, params), xlab = "Size", type = "l", ylab = "Survival probability", lwd = 12) # plot the survival model
points(ipm1$meshpts, apply(ipm1$P, 2, sum), col = "red", lwd = 3, cex = 0.5, pch = 19) # plot the column sums of the survival/growth matrix

### Modified IPM, with modified embryo and recruitment parameters
ipm2 <- bigmatrix(n = binN, params = paramsCA)
res2 <- popF(ipm2, binSize)
res2

##### PLOTTING #####
pdf("./figs/ipm_histo_fit.pdf", 7, 3.5)

set_graph_pars(ptype = "panel2")
xlab2 <- expression(paste("Size (", cm^2, ")"))

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
range(hist0710$area)
hist(hist0710$area, breaks = seq(0.0, 1.5, 0.05), freq = FALSE, 
	xlab = xlab2, ylab="Probability density", col = "gray87", main = "", 
	ylim = c(-0.025, 2), xlim = c(0, 1.4), border = "gray90", las = 1) 
box()
points(ipm1$meshpts , res1$ssd2, type="l", lty=1, lwd = 2, 
		col = "darkgray")
points(ipm2$meshpts , res2$ssd2, type="l", lty=1, lwd = 2, 
		col = "black")

add_panel_label(ltype = "b")

arrows(res1$max99, 0.6, res1$max99, 0.3, col = "darkgray", 
 		length = 0.1, lwd = 1.5, angle = 20)
arrows(res2$max99, 0.6, res2$max99, 0.3, col = "black", 
 		length = 0.1, lwd = 1.5, angle = 20)		

legend("topright", leg.txt, lwd = 2, bty = "n", 
	col = c("black", "darkgray"), lty = 1, 
	cex = 1, text.col = c("black", "darkgray"))

dev.off()

##### INFORMAL TESTS OF IPM SENSITIVITY #####
paramsWA
res2[1:4] # results of the modified IPM using WA parameters
 
paramsTest <- paramsWA; paramsTest$surv.slope = 0.5
paramsTest <- paramsWA; paramsTest$surv.int = 3

paramsTest <- paramsWA; paramsTest$growth.slope = 0.5
paramsTest <- paramsWA; paramsTest$growth.int = 0.1

ipmTest <- bigmatrix(n = binN, params = paramsTest)
resTest <- popF(ipmTest, binSize)
resTest[1:4]

