#################################################
# Author: Robin Elahi
# Date: 150318
# Base IPM parameterized by modern data
# Figure 2
#################################################

rm(list=ls(all=TRUE)) 

##### LOAD PACKAGES ETC #####
source("bael_IPM_fecundityTest.R")

##### PLOTTING #####
pdf("./figs/ipm_histo_fit.pdf", 7, 3.5)
set_graph_pars(ptype = "panel2")
xlab2 <- expression(paste("Size (", cm^2, ")"))

#### Panel A ####

plot(mx.coral ~ areaR, data = ltDat, xlim = c(0, 1.5), ylim = c(-0.025, 40),
     xlab = xlab2, ylab = "Embryo number", las = 1, type = "n")

add_panel_label(ltype = "a")

abline(mxRegCA, lwd=2, lty=1)
abline(a = 0, b = 0, lty=3, lwd=2, col="darkgray")
points(mx.coral ~ areaR, data = ltDat)

curve(slopeCA * x + data_in$embryo.int[medianEa_position], 
      from = 0, to = 1.5, add = TRUE,
      col = "darkgray", lwd = 2)

add_panel_label(ltype = "a")
legend("bottomright", "n = 12", cex = 1.1, bty = "n", adj = c(0,-2))

leg.txt <- c("Original", "Modified")
legend("topleft", leg.txt, lwd = 2, bty = "n", 
       col = c("black", "darkgray"), lty = 1, 
       cex = 1, text.col = c("black", "darkgray"))

#### Panel B ####
quantile(hist0710$area, probs = c(0.95, 0.99, 1))
obs99 <- quantile(hist0710$area, probs = 0.99)

hist(hist0710$area, breaks = seq(0.0, 1.5, 0.05), freq = FALSE,
     xlab = xlab2, ylab="Probability density", col = "gray87", main = "",
     ylim = c(-0.025, 2), xlim = c(0, 1.4), border = "gray90", las = 1)

box()

# California parameters
points(ipm2$meshpts , res2$ssd2, type="l", lty=1, lwd = 2, 
       col = "black")

# Optimized Ea
points(ipm_out[[medianEa_position]]$meshpts , 
       res_out[[medianEa_position]]$ssd2, type="l", lty=1, lwd = 2, 
       col = "darkgray")

arrows(res2$max99, 0.4, res2$max99, 0.13, col = "black", 
       length = 0.1, lwd = 1.5, angle = 20, lty = 1)

arrows(res_out[[medianEa_position]]$max99, 0.4, 
       res_out[[medianEa_position]]$max99, 0.13, col = "darkgray", 
       length = 0.1, lwd = 1.5, angle = 20, lty = 1)	

add_panel_label(ltype = "b")

legend("topright", leg.txt, lwd = 2, bty = "n", 
       col = c("black", "darkgray"), lty = 1, 
       cex = 1, text.col = c("black", "darkgray"))

dev.off()

##### PLOTTING SIMULATIONS #####

### CREATE COLOR RAMP ###
# range01 <- function(x)(x - min(x))/diff(range(x))
# heatRamp <- function(x) {
#   cols <- colorRamp(heat.colors(3)) (range01(x))
#   apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue = 255))
# }
# vec_colorsRev <- heatRamp(data_in$Ea)
# vec_colors <- rev(vec_colorsRev)

# for(i in 1:N){
#   curve(slopeCA * x + 
#           data_in$embryo.int[i], from = 0, to = 1.4, add = TRUE, 
#         col = vec_colors[i], lwd = 1)
# }


# for(i in 1:N){
#   points(ipm_out[[i]]$meshpts , res_out[[i]]$ssd2, type="l", lty=1, lwd = 1, 
#          col = vec_colors[i])
# }

# # Ea = 0.65 
# points(ipm1$meshpts , res1$ssd2, type="l", lty=1, lwd = 2, 
#        col = "darkgray")




