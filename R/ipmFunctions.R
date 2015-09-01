# RE 150521
#IPM functions, after Merow, Easterling

##########################################
##########################################
###1.4.3 Define vital rate functions to describe life history
# 1. survival probability function
s.x <-function(x, params) {
  u = exp(params$surv.int + params$surv.slope*x)
  return(u/(1+u))
}

# 2. growth function
g.yx <- function(xp, x, params) {
  dnorm(xp, mean = params$growth.int + params$growth.slope*x,
        sd  = params$growth.sd)
}

# 3. reproduction function
f.yx <- function(xp, x, params) {
  embryos = ifelse(x < params$mature.size, 0, 
        x*params$embryo.slope + params$embryo.int)
  embryos * params$estab.prob.mean *
  dnorm(xp, mean = params$recruit.size.mean, 
          sd = params$recruit.size.sd)
}

##########################################
##########################################
### Function to create the 'bigmatrix' of size n x n, with paramsX
bigmatrix <- function(n, params = paramsX) {
b <- min.size+c(0:n)*(max.size-min.size)/n # boundary points (edges of cells defining the matrix)
b
y <- 0.5*(b[1:n]+b[2:(n+1)]) # mesh points (i.e., midpoints to be used in numerical integration)
y
h <- y[2]-y[1] # step size (i.e., width of cells)
h

# now make the kernel 
G <- h*outer(y, y, g.yx, params = params)  # growth matrix
dim(G)
S <- s.x(y, params = params)              	# survival
P <- G				# placeholder; redefine P on the next line
# fix eviction of offspring
for(i in 1:(n/2)) {
  G[1,i] <- G[1,i] + 1-sum(G[,i])
  P[,i] <- G[,i]*S[i]
}

# fix eviction of large adults and small recruits
for(i in (n/2+1):n) {
  G[n,i] <- G[n,i] + 1 - sum(G[,i])
  P[,i] <- G[,i]*S[i]
}

F <- h*outer(y, y, f.yx, params = params) 	# reproduction matrix

K <- P + F				# full matrix

return(list(matrix = K, meshpts = y, P = P)); 

}

##########################################
##########################################
# FUNCTION TO GET RELEVANT POPULATION PARAMETERS
popF <- function(bigmatrix, binSize) {
	lam <- Re(eigen(bigmatrix$matrix)$values[1]) # dominant eigenvalue		
	w.eigen <- Re(eigen(bigmatrix$matrix)$vectors[,1]) # right eigenvector
	ssdSUM <- w.eigen/sum(w.eigen) # right eigen normalized to sum
	ssdNORM <- ssdSUM * 1/binSize
	meanSize <- sum(ssdSUM * bigmatrix$meshpts)
	maxSize95 <- bigmatrix$meshpts[min(which(cumsum(ssdSUM) > 0.95))] 	
	maxSize99 <- bigmatrix$meshpts[min(which(cumsum(ssdSUM) > 0.99))] 	
	v.eigen <- Re(eigen(t(bigmatrix$matrix))$vectors[,1]) # left eigenvector
	repro.val <- v.eigen/v.eigen[1] # normalized eigenvector
	y <- bigmatrix$meshpts
	h <- y[2] - y[1] # step size (width of cells)
	v.dot.w <- sum(ssdSUM*repro.val)*h
	sens <- outer(repro.val, ssdSUM)/v.dot.w
	elas <- matrix(as.vector(sens)*as.vector(bigmatrix$matrix)/lam, 
		nrow = length(bigmatrix$meshpts))
	#
	return(list(lambda = lam, meanSize = meanSize, 
	            max95 = maxSize95, max99 = maxSize99,
		ssd1 = ssdSUM, ssd2 = ssdNORM, meshpts = y, 
		reproVal = repro.val, 
		sensitivity = sens, elasticity = elas))
}
