#################################################
# Author: Robin Elahi
# Date: 151203
# Sensitivity analysis of base IPM
# MODFIED
# 160317: Used optimized parameters for fecundity
#################################################

rm(list=ls(all=TRUE))

##### LOAD PACKAGES ETC #####
source("./R/bael_params.R") 

source("./R/multiplotF.R") 
source("./R/ipmFunctions.R")

library(grid)

# Upload dataframe of modified vital rates
params <- paramsOptim

data_out <- params
rownames(data_out) <- "baseline"

paramVec <- names(params)
paramVec[1]
N <- length(paramVec)

counter = 1
i = 1

for(i in 1:N) {
  # Select parameter
  param.i <- paramVec[i]
  # Create new params df that is row-labeled with selected parameter
  paramDF.i <- params 
  rownames(paramDF.i) <- param.i
  # Choose the column with the parameter to be adjusted
  column.i <- which(colnames(paramDF.i) %in% param.i)
  # Increase the column's value by 10 %
  paramDF.i[1, column.i] <- paramDF.i[1, column.i] * 1.1
  # Bind the original dataframe with the new, adjusted dataframe
  data_out <- rbind(data_out, paramDF.i)
}

data_out

##### MODIFYING THE IPM #####
###1.4.4 Make a kernel###
min.size <- 0.02
max.size <- 1.52
binSize <- 0.02

n <- (max.size - min.size)/binSize # number of cells in the matrix
b <- min.size+c(0:n)*(max.size-min.size)/n # boundary points (edges of cells defining the matrix)
b
y <- 0.5*(b[1:n]+b[2:(n+1)]) # mesh points (i.e., midpoints to be used in numerical integration)
y
h <- y[2]-y[1] # step size (i.e., width of cells)
h

##### SET UP THE FOR LOOP #####
loopDat <- data_out
# choose specific rows
AllReps <- unique(rownames(loopDat))
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

mat2 <- cbind(AllReps, data.frame(mat1)) %>% rename(param = AllReps)

##### PREPARE OUTPUT FOR PLOTTING #####
mat2$deltaLambda <- mat2$lambda - mat2$lambda[1]
mat2$delta95 <- mat2$maxSize95 - mat2$maxSize95[1]
mat2$delta99 <- mat2$maxSize99 - mat2$maxSize99[1]
mat2$deltaMean <- mat2$meanSize - mat2$meanSize[1]

mat2$perLambda <- with(mat2, 100 * deltaLambda/mat2$lambda[1])
mat2$per95 <- with(mat2, 100 * delta95/mat2$maxSize95[1])
mat2$per99 <- with(mat2, 100 * delta99/mat2$maxSize99[1])
mat2$perMean <- with(mat2, 100 * deltaMean/mat2$meanSize[1])

fullDat <- mat2 %>%
  filter(param != "recruit.size.sd" & param != "mature.size" &
           param != "growth.sd" & param != "estab.prob.sd" & 
           param != "embryo.sd" & param != "baseline") %>% 
  droplevels()

# Want to plot percentage change, facetted by lambda, mean size, and max size
plotDat <- fullDat %>% select(param, perLambda, per99, perMean) %>% 
  gather(key = variable, value = value, perLambda:perMean) 

# rename levels of variable
length_per_variable <- dim(fullDat)[1]

plotDat$variable <- c(rep("Lambda", length_per_variable), 
                     rep("Maximum size", length_per_variable), 
                     rep("Mean size", length_per_variable))

### Plotting details

theme_set(theme_classic(base_size = 12))
dodge <- position_dodge(width = 0)
ULClabel <- theme(plot.title = element_text(hjust = -0.15, vjust = 1, size = rel(1.2)))

unique(plotDat$param)

# Change order of values
paramOrder <- c("Embryo_intercept", "Embryo_slope", "EstabProb", "RecruitSize", 
                "Survival_intercept", "Survival_slope", "Growth_intercept", "Growth_slope")

paramOrder <- c("embryo.int", "embryo.slope", "estab.prob.mean", 
                "recruit.size.mean", "surv.int", "surv.slope", 
                "growth.int", "growth.slope")

# Change labels
paramLabels <- c("Fecundity\n(intercept)", "Fecundity\n(slope)", 
                 "Recruitment\nprobability", "Recruit size", 
                 "Survival\n(intercept)", "Survival\n(slope)", 
                 "Growth\n(intercept)", "Growth\n(slope)")

unique(plotDat$variable)

# Relevel factors of variable
unique(plotDat$variable)
plotDat$variable <- as.factor(plotDat$variable)
print(levels(plotDat$variable))
plotDat$variable <- with(plotDat, factor(variable, levels(variable)[c(1, 3, 2)]))

##### MULTI-PANEL PLOT #####

ULClabel <- theme(plot.title = element_text(hjust = -0.25, vjust = 1, size = rel(1.2)))

ggplot(plotDat, aes(param, value)) +
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

ggsave(file = "./figs/ipm_sensitivity.pdf", height = 3.5, width = 7)


