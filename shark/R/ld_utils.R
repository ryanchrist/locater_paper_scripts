
r2 <- function(X, idx = 1L){
  # X should be a matrix of haplotypes with L (number of variants) rows and N (number of haplotypes) columns 
  # idx should be an index for the row / variant in X against which we will calculate r^2 for all variants
  
  if(!is.matrix(X) || nrow(X) < 2){
    stop("X must be a matrix with at least 2 rows")
  }
  
  orig.idx <- idx
  idx <- as.integer(idx)
  
  if(orig.idx != idx || length(idx) != 1 || idx < 1 || idx > nrow(X)){
    stop("idx must be an integer index corresponding to a row of X")
  }
  
  daf <- c(rowMeans(X))
  
  ((c(X %*% c(X[idx,]))/ncol(X) - daf  * daf[idx])^2) / (daf * (1-daf) * daf[idx] * (1-daf[idx]))
}

D.prime <- function(X, idx = 1L){
  # X should be a matrix of haplotypes with L (number of variants) rows and N (number of haplotypes) columns 
  # idx should be an index for the row / variant in X against which we will calculate D' for all variants
  
  if(!is.matrix(X) || nrow(X) < 2){
    stop("X must be a matrix with at least 2 rows")
  }
  
  orig.idx <- idx
  idx <- as.integer(idx)
  
  if(orig.idx != idx || length(idx) != 1 || idx < 1 || idx > nrow(X)){
    stop("idx must be an integer index corresponding to a row of X")
  }
  

  
  EX <- c(rowMeans(X))
  Ey <- EX[idx]
  EX.Ey <- EX * Ey
  
  D <- c(X %*% c(X[idx,]))/ncol(X) - EX.Ey
  
  denom <- ifelse(D >= 0, pmin(EX, Ey), pmax(0,Ey + EX - 1)) - EX.Ey
  
  D / denom
}



Jaccard <- function(X, idx = 1L){
  # X should be a matrix of haplotypes with L (number of variants) rows and N (number of haplotypes) columns 
  # idx should be an index for the row / variant in X against which we will calculate the Jaccard index for all variants
  
  if(!is.matrix(X) || nrow(X) < 2){
    stop("X must be a matrix with at least 2 rows")
  }
  
  orig.idx <- idx
  idx <- as.integer(idx)
  
  if(orig.idx != idx || length(idx) != 1 || idx < 1 || idx > nrow(X)){
    stop("idx must be an integer index corresponding to a row of X")
  }
  
  y <- c(X[idx,])
  
  X.and.y <- c(X %*% y)
  X.union.y <- c(rowSums(X)) + sum(y) - X.and.y
  
  X.and.y / X.union.y
}

Overlap <- function(X, idx = 1L){
  # X should be a matrix of haplotypes with L (number of variants) rows and N (number of haplotypes) columns 
  # idx should be an index for the row / variant in X against which we will calculate the nestedness index for all variants
  
  if(!is.matrix(X) || nrow(X) < 2){
    stop("X must be a matrix with at least 2 rows")
  }
  
  orig.idx <- idx
  idx <- as.integer(idx)
  
  if(orig.idx != idx || length(idx) != 1 || idx < 1 || idx > nrow(X)){
    stop("idx must be an integer index corresponding to a row of X")
  }
  
  dac <- c(rowSums(X))
  
  c(X %*% c(X[idx,])) / ifelse(dac >= dac[idx],dac[idx],dac)
}
# 
# require(kalis)
# X <- SmallHaps
# lead.var.idx <- 250
# daf <- c(rowMeans(X))
# layout(matrix(1:4,2))
# 
# yy <- r2(QueryCache(),lead.var.idx)
# plot(yy,ylim = c(0,1),bty="n",las=1,col="darkgrey", main = "LD with red target variant: r^2",ylab="r^2",xlab = "pos")
# points(lead.var.idx,yy[lead.var.idx],col="red",pch=19)
# 
# yy <- D.prime(QueryCache(),lead.var.idx)
# plot(yy,ylim = c(0,1),bty="n",las=1,col="darkgrey", main = "LD with red target variant: D'",ylab="D'",xlab = "pos")
# points(lead.var.idx,yy[lead.var.idx],col="red",pch=19)
# 
# yy <- Jaccard(QueryCache(),lead.var.idx)
# plot(yy,ylim = c(0,1),bty="n",las=1,col="darkgrey", main = "LD with red target variant: Jaccard",ylab="Jaccard index",xlab = "pos")
# points(lead.var.idx,yy[lead.var.idx],col="red",pch=19)
# 
# yy <- Overlap(QueryCache(),lead.var.idx)
# plot(yy,ylim = c(0,1),bty="n",las=1,col="darkgrey", main = "LD with red target variant: Overlap",ylab="Overlap coefficient",xlab = "pos")
# points(lead.var.idx,yy[lead.var.idx],col="red",pch=19)
# 
# 
# X <- matrix(0,nrow=8298,ncol=2)
# X[1:1151,] <- 1
# X[1152:(1151+4794),1] <- 1
# X[(1151+4794+1):(1151+4794+70),2] <- 1
# 
# X <- t(X)
# 
# r2(X,1)
# D.prime(X,1)
# Jaccard(X,1)
# Overlap(X,1)

