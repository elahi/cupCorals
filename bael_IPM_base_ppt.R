#################################################
# Author: Robin Elahi
# Date: 150318
# Base IPM parameterized by modern data
# Figure 2
# Plot for ppt
#################################################

rm(list=ls(all=TRUE)) 

##### LOAD PACKAGES ETC #####
source("bael_IPM_fecundityTest.R")

##### PLOT - histo with curve #####

pdf("./figs/ipm_histo_fit_ppt.pdf", 3.5, 3.5)
set_graph_pars(ptype = "panel1")
xlab2 <- expression(paste("Size (", cm^2, ")"))

#### Plot observations with temperature adjusted IPM ####
quantile(hist0710$area, probs = c(0.95, 0.99, 1))
obs99 <- quantile(hist0710$area, probs = 0.99)

hist(hist0710$area, breaks = seq(0.0, 1.5, 0.05), freq = FALSE,
     xlab = xlab2, ylab="Probability density", col = "gray87", main = "",
     ylim = c(-0.025, 2), xlim = c(0, 1.4), border = "gray90", las = 1)

box()

# Optimized Ea
points(ipm_out[[medianEa_position]]$meshpts , 
       res_out[[medianEa_position]]$ssd2, type="l", lty=1, lwd = 2, 
       col = "black")

arrows(res_out[[medianEa_position]]$max99, 0.4, 
       res_out[[medianEa_position]]$max99, 0.13, col = "black", 
       length = 0.1, lwd = 1.5, angle = 20, lty = 1)	

dev.off()

##### PLOT - histo without curve #####
pdf("./figs/ipm_histo_fit_ppt2.pdf", 3.5, 3.5)
set_graph_pars(ptype = "panel1")
xlab2 <- expression(paste("Size (", cm^2, ")"))

#### Plot observations with temperature adjusted IPM ####
quantile(hist0710$area, probs = c(0.95, 0.99, 1))
obs99 <- quantile(hist0710$area, probs = 0.99)

hist(hist0710$area, breaks = seq(0.0, 1.5, 0.05), freq = FALSE,
     xlab = xlab2, ylab="Probability density", col = "gray87", main = "",
     ylim = c(-0.025, 2), xlim = c(0, 1.4), border = "gray90", las = 1)

box()

dev.off()

