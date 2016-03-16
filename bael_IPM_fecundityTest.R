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
source("bael_embryos_Ea.R")

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

##### MAX AND MIN SIZES #####
# Will use slightly larger size range for IPM because
# it is otherwise artificially truncated
min.size <- 0.05
max.size <- 1.55

binSize <- 0.05
binN <- (max.size - min.size)/binSize
binN

# breaks = seq(0.02, 1.5, 0.02),

##### IPMs WITH FOR LOOP THROUGH EAS #####
head(loopDat2)

data_in <- loopDat2
N <- dim(data_in)[1]
ipm_out <- result <- vector("list", N) 
res_out <- result <- vector("list", N) 

ipm_out[[1]]


# Initialize parameter dataframe
paramDF <- paramsWA # the cells will be modified in the loop
counter <- 1
EaString <- NULL

for (i in 1:N){
  # Create new paramDF, with modified fecundity parameters
  params.i <- data_in[i, ]
  Ea.i <- as.character(params.i$Ea)
  paramDF$estab.prob.mean <- params.i$estab.prob.mean
  paramDF$estab.prob.sd <- params.i$estab.prob.sd
  paramDF$embryo.int <- params.i$embryo.int
  paramDF$mature.size <- params.i$mature.size
  
  # Now run IPM
  ipm1 <- bigmatrix(n = binN, params = paramDF)
  res1 <- popF(ipm1, binSize)
  
  ipm_out[[i]] <- res1
  res_out[[i]] <- res1
  EaString <- c(EaString, Ea.i)
}

EaString
str(ipm_out[[1]])
str(res_out[[1]])

##### IPMs for WA and CA #####
### Base IPM, original parameters
ipm1 <- bigmatrix(n = binN, params = paramsWA)
res1 <- popF(ipm1, binSize)
res1[1:8]

### Modified IPM, with modified embryo and recruitment parameters
ipm2 <- bigmatrix(n = binN, params = paramsCA)
res2 <- popF(ipm2, binSize)
res2

##### CREATE COLOR RAMP #####
range01 <- function(x)(x - min(x))/diff(range(x))
heatRamp <- function(x) {
  cols <- colorRamp(heat.colors(3)) (range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue = 255))
}
vec_colorsRev <- heatRamp(data_in$Ea)
vec_colors <- rev(vec_colorsRev)

##### PLOTTING #####
# pdf("./figs/ipm_histo_fit.pdf", 7, 3.5)

set_graph_pars(ptype = "panel2")
xlab2 <- expression(paste("Size (", cm^2, ")"))

#### Panel A ####
plot(mx.coral ~ areaR, data = ltDat, xlim = c(0, 1.4), ylim = c(-0.025, 40),
     xlab = xlab2, ylab = "Embryo number", las = 1, type = "n")

add_panel_label(ltype = "a")

for(i in 1:N){
  curve(slopeCA * x + 
          data_in$embryo.int[i], from = 0, to = 1.4, add = TRUE, 
        col = vec_colors[i], lwd = 1)
}

abline(mxRegCA, lwd=2, lty=1)
abline(a = 0, b = 0, lty=3, lwd=2, col="darkgray")
points(mx.coral ~ areaR, data = ltDat)

add_panel_label(ltype = "a")

#### Panel B ####
range(hist0710$area)
hist(hist0710$area, breaks = seq(0.00, 1.55, 0.05), freq = FALSE, 
	xlab = xlab2, ylab="Probability density", col = "gray87", main = "", 
	ylim = c(-0.025, 2.25), xlim = c(0.00, 1.55), border = "gray90", las = 1) 

# # TO MATCH FINER IPM
# hist(hist0710$area, breaks = seq(0.02, 1.5, 0.02), freq = FALSE, 
#      xlab = xlab2, ylab="Probability density", col = "gray87", main = "", 
#      ylim = c(-0.025, 2.5), xlim = c(0, 1.4), border = "gray90", las = 1) 
# 
# box()

for(i in 1:N){
  points(ipm_out[[i]]$meshpts , res_out[[i]]$ssd2, type="l", lty=1, lwd = 1, 
         col = vec_colors[i])
}

points(ipm1$meshpts , res1$ssd2, type="l", lty=1, lwd = 2, 
       col = "darkgray")
points(ipm2$meshpts , res2$ssd2, type="l", lty=1, lwd = 2, 
       col = "black")

add_panel_label(ltype = "b")

# dev.off()

##### CALCULATE LL #####
min.size <- 0.02
max.size <- 1.5

binSize <- 0.02
binN <- (max.size - min.size)/binSize
binN

### Base IPM, original parameters
ipm1 <- bigmatrix(n = binN, params = paramsWA)
res1 <- popF(ipm1, binSize)
res1[1:8]

area.h <- hist(hist0710$area, breaks = seq(0, 1.5, 0.05))
area.h
cnt <- area.h$counts 
dns <- area.h$density
sum(cnt) 
sum(dns); length(dns)
sum(dns * 0.05)

# using density
mod1LL <- - sum(dns * log(res1$ssd2)); mod1LL # LL for empirical data is -1.2


