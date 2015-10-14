#################################################
# Author: Robin Elahi
# Date: 151013

# Sensitivity analysis of base IPM
# Figure 3
#################################################

rm(list=ls(all=TRUE))
source("./R/baelParamsWA.R") 
source("./R/multiplotF.R") 

params

################################################
#### MODIFYING THE IPM
#### CHANGE EACH PARAMETER BY 10%
#### CALCULATE NEW MAX 95, 99, MEAN, LAMBDA
################################################
# params
# # what will I actually change?
# parDF <- params
# parDF
# 
# # remove mature size
# parDF$mature.size <- NULL
# parDF
# write.csv(parDF, "parDF.csv")

# upload relevant scenarios
parDFup <- read.csv("./data/parDFup.csv")

###1.4.4 Make a kernel###
min.size <- 0.02
max.size <- 1.5

binSize <- 0.02
n <- (max.size - min.size)/binSize # number of cells in the matrix

b <- min.size+c(0:n)*(max.size-min.size)/n # boundary points (edges of cells defining the matrix)
b
y <- 0.5*(b[1:n]+b[2:(n+1)]) # mesh points (i.e., midpoints to be used in numerical integration)
y
h <- y[2]-y[1] # step size (i.e., width of cells)
h

### Now calculate lambda and max size of the stable size distribution after running the ipm
# as a function of different temps
head(parDFup)
# for loop
loopDat <- parDFup

dim(loopDat)

# choose specific rows
AllReps <- unique(loopDat$X)
AllReps
N <- length(AllReps); N

######################################
######################################
### CHANGE GROWTH INT AND SLOPE
######################################
######################################

# matrix
mat1 <- matrix(nrow = N, ncol = 4)
colnames(mat1) <- c("lambda", "maxSize95", "maxSize99", "meanSize")
mat1
head(loopDat)

# Loop over each row in loopDat; get lambda and maxSize
for(g in 1:N){
  row.g <- AllReps[g]
  
  param.g <- data.frame(
  	growth.int = loopDat$growth.int[row.g], 
  	growth.slope = loopDat$growth.slope[row.g], 
  	growth.sd = loopDat$growth.sd[row.g], 
  	embryo.int = loopDat$embryo.int[row.g], 
  	embryo.slope = loopDat$embryo.slope[row.g], 
  	mature.size = -loopDat$embryo.int[row.g]/loopDat$embryo.slope[row.g],
  	recruit.size.mean = loopDat$recruit.size.mean[row.g], 
  	recruit.size.sd = loopDat$recruit.size.mean[row.g], 
  	estab.prob.mean = loopDat$estab.prob.mean[row.g], 
  	surv.int = loopDat$surv.int[row.g], 
  	surv.slope = loopDat$surv.slope[row.g])  	

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

mat2$deltaLambda <- mat2$lambda - mat2$lambda[1]
mat2$delta95 <- mat2$maxSize95 - mat2$maxSize95[1]
mat2$delta99 <- mat2$maxSize99 - mat2$maxSize99[1]
mat2$deltaMean <- mat2$meanSize - mat2$meanSize[1]

mat2$perLambda <- with(mat2, 100 * deltaLambda/mat2$lambda[1])
mat2$per95 <- with(mat2, 100 * delta95/mat2$maxSize95[1])
mat2$per99 <- with(mat2, 100 * delta99/mat2$maxSize99[1])
mat2$perMean <- with(mat2, 100 * deltaMean/mat2$meanSize[1])

mat2

fullDat <- mat2
names(fullDat)
# go from wide to long so that i can use facets
# facets - absolute vs percent change; X
names(fullDat)[1] <- "param"

library(reshape2)
datLong <- melt(fullDat, id.vars = c("param"))

head(datLong)
dim(datLong)

datLong$output <- c(rep("IPM", 48), rep("absolute change", 48), rep("percent change", 48))

datLong
summary(datLong)
unique(datLong$variable)
unique(datLong$param)

#write.csv(datLong, "sensitivityDatLong_150423.csv")

####
# read in sensitivity data

#####
ggDat <- droplevels(datLong[datLong$output == "percent change" &
				datLong$variable == "perLambda" |
				datLong$variable == "per99" |
				datLong$variable == "perMean", ])
dim(ggDat)
ggDat2 <- droplevels(ggDat[ggDat$param != "Baseline" &
				ggDat$param != "Survival_IS" &
				ggDat$param != "Embryo_IS" &
				ggDat$param != "Growth_IS" , ])
dim(ggDat2)

unique(ggDat2$param)
unique(ggDat2$variable)
unique(ggDat2$output)

ggDat3 <- ggDat2
ggDat3
dim(ggDat3)
# rename levels of variable
ggDat3$variable <- c(rep("Lambda", 8), rep("Maximum size", 8), rep("Mean size", 8))
ggDat3


####
plotDat <- ggDat3
plotDat

dodge <- position_dodge(width = 0)

theme_set(theme_classic(base_size = 12))
plotDat$param

# want to change order of values
paramOrder <- c("Embryo_intercept", "Embryo_slope", "EstabProb", "RecruitSize", 
                "Survival_intercept", "Survival_slope", "Growth_intercept", "Growth_slope")

paramLabels <- c("Fecundity\n(intercept)", "Fecundity\n(slope)", 
                 "Recruitment\nprobability", "Recruit size", "Survival\n(intercept)", 
                 "Survival\n(slope)", "Growth\n(intercept)", "Growth\n(slope)")

ULClabel <- theme(plot.title = element_text(hjust = -0.15, vjust = 1, size = rel(1.2)))

##########################################################
# Single plot, maximum size
##########################################################
head(plotDat)
plotMaxSize <- ggplot(plotDat[plotDat$variable == "Maximum size", ], 
	aes(param, value)) +
	geom_bar(stat = "identity", color = "black", fill = "darkgray") +
	geom_hline(yintercept = 0, color = "gray", linetype = "dashed") + 
	theme(legend.justification = c(1, 1), legend.position = c(0.95, 0.95)) +
	scale_fill_discrete(name = "Perturbation") +
	ylab("Percentage change") + xlab("Parameter") +
	scale_x_discrete(limits = rev(paramOrder), labels = rev(paramLabels)) +
	theme(strip.background = element_blank()) +
  scale_y_continuous(limits = c(-3, 33)) + 
  theme(axis.text.y = element_text(size = 8))
	
plotMaxSize + coord_flip()
ggsave(file = "./figs/ipm_sensitivity.pdf", height = 3.5, width = 3.5)
