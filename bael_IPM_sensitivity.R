#################################################
# Author: Robin Elahi
# Date: 151203
# Sensitivity analysis of base IPM
#################################################

rm(list=ls(all=TRUE))

##### LOAD PACKAGES ETC #####
source("./R/bael_params.R") 
params <- paramsWA

source("./R/multiplotF.R") 
source("./R/ipmFunctions.R")

library(grid)

# Upload dataframe of modified vital rates
parDFup <- read.csv("./data/parDFup.csv")
head(parDFup)

##### MODIFYING THE IPM #####
###1.4.4 Make a kernel###
min.size <- 0.02
max.size <- 1.52

n <- 100 # number of cells in the matrix
b <- min.size+c(0:n)*(max.size-min.size)/n # boundary points (edges of cells defining the matrix)
b
y <- 0.5*(b[1:n]+b[2:(n+1)]) # mesh points (i.e., midpoints to be used in numerical integration)
y
h <- y[2]-y[1] # step size (i.e., width of cells)
h

##### SET UP THE FOR LOOP #####
loopDat <- parDFup
# choose specific rows
AllReps <- unique(loopDat$X)
N <- length(AllReps); N

# Set up a matrix for output data 
mat1 <- matrix(nrow = N, ncol = 4)
colnames(mat1) <- c("lambda", "maxSize95", "maxSize99", "meanSize")

# For loop
for(g in 1:N){

  param.g <- data.frame(
  	growth.int = loopDat$growth.int[g], 
  	growth.slope = loopDat$growth.slope[g], 
  	growth.sd = loopDat$growth.sd[g], 
  	embryo.int = loopDat$embryo.int[g], 
  	embryo.slope = loopDat$embryo.slope[g], 
  	mature.size = -loopDat$embryo.int[g]/loopDat$embryo.slope[g],
  	recruit.size.mean = loopDat$recruit.size.mean[g], 
  	recruit.size.sd = loopDat$recruit.size.mean[g], 
  	estab.prob.mean = loopDat$estab.prob.mean[g], 
  	surv.int = loopDat$surv.int[g], 
  	surv.slope = loopDat$surv.slope[g])  	

   	S <- s.x(y, params = param.g) 						# survival 	
  	F <- h*outer(y, y, f.yx, params = param.g)  	# reproduction
  	G <- h*outer(y, y, g.yx, params = param.g) 	# growth
  	P <- G # placeholder; redefine P on the next line
  # fix eviction of offspring
  for(i in 1:(n/2)) {
    G[1,i] <- G[1,i] + 1-sum(G[,i])
    P[,i] <- G[,i]*S[i]
  }
  
  # fix eviction of large adults
  for(i in (n/2+1):n) {
    G[n,i] <- G[n,i] + 1 - sum(G[,i])
    P[,i]<-G[,i]*S[i]
  }
  K <- P + F  			# full matrix
  # calculate lambda
  lam <- Re(eigen(K)$values[1]) # returns dominant eigenvalue, which gives the asymptotic population growth rate
  # calculate max size (99% of stable size distribution)
  w.eigen <- Re(eigen(K)$vectors[,1])
  stable.dist <- w.eigen/sum(w.eigen)
  maxSize95 <- y[min(which(cumsum(stable.dist) > 0.95))] 
  maxSize99 <- y[min(which(cumsum(stable.dist) > 0.99))] 
  meanSize <- sum(y*stable.dist)
# population matrix with continuous variables
  mat1[g,] <- c(lam, maxSize95, maxSize99, meanSize)
}

mat1
mat2 <- cbind(loopDat[1], mat1)
mat2

##### PREPARE OUTPUT FOR PLOTTING #####
mat2$deltaLambda <- mat2$lambda - mat2$lambda[1]
mat2$delta95 <- mat2$maxSize95 - mat2$maxSize95[1]
mat2$delta99 <- mat2$maxSize99 - mat2$maxSize99[1]
mat2$deltaMean <- mat2$meanSize - mat2$meanSize[1]

mat2$perLambda <- with(mat2, 100 * deltaLambda/mat2$lambda[1])
mat2$per95 <- with(mat2, 100 * delta95/mat2$maxSize95[1])
mat2$per99 <- with(mat2, 100 * delta99/mat2$maxSize99[1])
mat2$perMean <- with(mat2, 100 * deltaMean/mat2$meanSize[1])

fullDat <- mat2

# go from wide to long so that i can use facets
# facets - absolute vs percent change; X
names(fullDat)[1] <- "param"

library(reshape2)
datLong <- melt(fullDat, id.vars = c("param"))

datLong$output <- c(rep("IPM", 48), rep("absolute change", 48), 
                    rep("percent change", 48))

### Select the data I want to plot
ggDat <- droplevels(datLong[datLong$output == "percent change" &
				datLong$variable == "perLambda" |
				datLong$variable == "per99" |
				datLong$variable == "perMean", ])

ggDat2 <- droplevels(ggDat[ggDat$param != "Baseline" &
				ggDat$param != "Survival_IS" &
				ggDat$param != "Embryo_IS" &
				ggDat$param != "Growth_IS" , ])

ggDat3 <- ggDat2

# rename levels of variable
ggDat3$variable <- c(rep("Lambda", 8), rep("Maximum size", 8), rep("Mean size", 8))

# Rename data for plotting
plotDat <- ggDat3

### Plotting details

theme_set(theme_classic(base_size = 12))
dodge <- position_dodge(width = 0)
ULClabel <- theme(plot.title = element_text(hjust = -0.15, vjust = 1, size = rel(1.2)))

# Change order of values
paramOrder <- c("Embryo_intercept", "Embryo_slope", "EstabProb", "RecruitSize", 
                "Survival_intercept", "Survival_slope", "Growth_intercept", "Growth_slope")

# Change labels
paramLabels <- c("Fecundity\n(intercept)", "Fecundity\n(slope)", 
                 "Recruitment\nprobability", "Recruit size", "Survival\n(intercept)", 
                 "Survival\n(slope)", "Growth\n(intercept)", "Growth\n(slope)")

unique(plotDat$variable)

# Relevel factors of variable
unique(plotDat$variable)
plotDat$variable <- as.factor(plotDat$variable)
print(levels(plotDat$variable))
plotDat$variable <- with(plotDat, factor(variable, levels(variable)[c(1, 3, 2)]))

##### MULTI-PANEL PLOT #####

ULClabel <- theme(plot.title = element_text(hjust = -0.25, vjust = 1, size = rel(1.2)))

sensitivityPlot <- ggplot(plotDat, aes(param, value)) +
  geom_bar(stat = "identity", color = "black", fill = "darkgray") +
  geom_hline(yintercept = c(0, 10), color = "gray", linetype = "dashed") + 
  # geom_hline(yintercept = 10, lty = 2, color = "darkgray") + 
  theme(legend.justification = c(1, 1), legend.position = c(0.95, 0.95)) +
  scale_fill_discrete(name = "Perturbation") +
  ylab("Percentage change") + xlab("Parameter") +
  scale_x_discrete(limits = rev(paramOrder), labels = rev(paramLabels)) +
  theme(strip.background = element_blank()) +
  # scale_y_continuous(limits = c(-5, 35)) + 
  theme(axis.text.y = element_text(size = 8)) + 
  facet_wrap(~ variable) + coord_flip() + 
  # labs(title = "A") + ULClabel + 
  theme(panel.margin = unit(1.5, "lines"))

sensitivityPlot

ggsave(file = "./figs/ipm_sensitivity.pdf", height = 3.5, width = 7)


