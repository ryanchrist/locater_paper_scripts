
TestLocus.Univariate <- function(mat,y,Ey,SDy,k=c(20,200),neg.log10.cutoff=5, eval.for.popstruc = FALSE, lower.tail = TRUE, return.bounds = TRUE){ #, project.struc = NULL){

  # if(is.null(project.struc)){
  mat <- sweep(mat*SDy,2,SDy,"*")
  y.residuals <- (y - Ey)/SDy
  obs <- sum(y.residuals*(mat%*%y.residuals))
  sumM2 <- sum(mat^2)
  trM <- sum(diag(mat))

  #
  # }else{
  #   xx <- y - Ey - c(project.struc$Q%*%(project.struc$Z%*%(y-Ey)))
  #   obs <- sum(xx*(mat%*%xx))
  # }
  #

  CalcBounds.Univariate(args = list(mat),
                        obs = obs,
                        trace1 = trM,
                        upper_bound_on_trace2 = sumM2)
}



CalcBounds.Univariate <- function(f = function(k,args){eigs_sym(args[[1]],k)},
                                  args,
                                  obs,
                                  trace1,
                                  upper_bound_on_trace2,
                                  k=c(20,200),
                                  neg.log10.cutoff=5,
                                  eval.for.popstruc = FALSE,
                                  lower.tail = TRUE){

  res <- c(1,NA,NA,NA) # start as potentially interesting

  if(eval.for.popstruc){ k <- max(k)} # eval.for.popstruc means just go directly to evaluating the largest k

  for(j in 1:length(k)){

    e <- f(k[j], args)

    R.sum.etasq <- upper_bound_on_trace2 - sum(e$values^2)
    E.R <-  trace1 - sum(e$values)

    gauss.tcdf <- QFGauss(e$values,sigma = sqrt(2*R.sum.etasq))

    gauss.approxfullcdf <- function(x, density = FALSE, lower.tail = TRUE, log.p = FALSE) gauss.tcdf(x-E.R, density = density, lower.tail = lower.tail, log.p = log.p)
    attr(gauss.approxfullcdf,"mu") <- attr(gauss.tcdf,"mu") + E.R
    attr(gauss.approxfullcdf,"Q.sd") <- attr(gauss.tcdf,"Q.sd")

    tcdf <- QFGauss(e$values)

    if(lower.tail){
      res[2] <-  -gauss.approxfullcdf(obs,log.p=T)/log(10)
      res[3:4] <- -log10(QFGaussBounds(tcdf,"identity",min(abs(e$values)),E.R, R.sum.etasq)(obs)[1,1:2])
      if(res[3] < neg.log10.cutoff){
        res[1] <- 0;
        break
      }

    }else{
      res[2] <-  -gauss.approxfullcdf(obs,lower.tail = F, log.p=T)/log(10)
      res[3:4] <- -log10(QFGaussBounds(tcdf,"identity",min(abs(e$values)),E.R, R.sum.etasq)(obs)[1,3:4])

      if(res[4] < neg.log10.cutoff){
        res[1] <- 0;
        break
      }

    }
  }
  # first entry: indicator of still potentially interesting
  # second entry: -log approx p-value
  # third entry: lower bound on p-value
  # fourth entry: upper bound on p-value
  res
}


Fast.One.Eig.Screen <- function(first.eig = NULL,
                                obs,
                                trace1,
                                upper_bound_on_trace2){

  R.sum.etasq <- upper_bound_on_trace2 - first.eig$values^2
  E.R <-  trace1 - first.eig$values

  tcdf <- suppressWarnings(QFGauss(first.eig$values))

  # return lower bound on p-value in the lower tail
  c(QFGaussBounds(tcdf,"identity",abs(first.eig$values),E.R, R.sum.etasq)(obs)[1,1])
}

project.matrix <- function(M,Q,sigma = NULL){
  X <- M%*%Q
  temp <- tcrossprod(X,Q)
  M <- M - temp - t(temp) + tcrossprod( Q%*%(crossprod(Q,X)), Q)
  rm(X,temp); gc()
  if(is.null(sigma)){
    return(M)
  }else{
    return(sweep(M*sigma,2,sigma,"*"))
  }
}

# Below is general version for univariate and multivariate testing
CalcBounds <- function(f = function(k,args){eigs_sym(args[[1]],k)},
                       args,
                       obs,
                       traces,
                       k=c(20,200),
                       neg.log10.cutoff=5,
                       tau = 1, # may be a vector as long as obs
                       delta = 0, # may be a vector as long as obs
                       only.point.est = FALSE, # if TRUE, then the first element of k is used for calculated trunc part
                       eval.for.popstruc = FALSE,
                       unfinished = NULL,
                       lower.tail = TRUE){

  res <- matrix(NA_real_,nrow=5,ncol=length(obs))
  res[1,] <- 1

  one.inflation.setting <- length(unique(tau)) == 1 & length(unique(delta)) == 1

  if(only.point.est){
    e <- f(k[1], args)
    res[2,] <- k[1]

    if(one.inflation.setting){  # calculate function only once
      calc.func <- CalcBounds2(traces, e, tau[1], delta[1], only.point.est)
    }

    for(p in 1:length(obs)){
      if(!one.inflation.setting){  # calculate function separately for each phenotype
        calc.func <- CalcBounds2(traces, e, tau[p], delta[p], only.point.est)
      }
      res[5,p] <- calc.func(obs[p],lower.tail)
    }

    return(res)
  }


  if(eval.for.popstruc){ k <- max(k)} # eval.for.popstruc means just go directly to evaluating the largest k

  if(is.null(unfinished)){unfinished <- rep(TRUE,length(obs))}
  j <- 0

  while(any(unfinished)){

    j <- j+1 # advance to next k
    e <- f(k[j], args)

    if(one.inflation.setting){ # calculate function only once
      calc.func <- CalcBounds2(traces, e, tau[1], delta[1])
    }

    for(p in which(unfinished)){
      res[2,p] <- k[j]

      if(!one.inflation.setting){ # calculate function separately for each phenotype
        calc.func <- CalcBounds2(traces, e, tau[p], delta[p])
      }

      res[3:5,p] <- calc.func(obs[p],lower.tail)

      if( res[ if(lower.tail){3}else{4}, p ] < neg.log10.cutoff ){
        res[1,p] <- 0
        unfinished[p] <- FALSE
      }
    }

    if(j==length(k)){unfinished <- FALSE} # we've run out of k to evaluate for now
  }

  # first row: still interesting indicator
  # second row: final k used to calculate this particular bound and point estimate
  # third row: -log10 lower bound
  # fourth row: -log10 upper bound
  # fifth row: -log10 pvalue point estimate

  res
}


CalcBounds2 <- function(traces, e = NULL, tau = 1, delta = 0, only.point.est = FALSE){

  # this returns a function that is NOT vectorized (it can only take observed value at a time)
  # the returned function returns a vector of three numbers:
  # first entry: -log10 lower bound
  # second entry: -log10 upper bound (NULL if e is NULL)
  # third entry: -log10 pvalue point estimate (NULL if e is NULL)

  if(tau < 0 | delta < 0 | length(tau) != 1 | length(delta) != 1){
    stop("tau and delta must be greater than or equal to 0 and have length 1.")
  }

  delta2 <- delta^2


  if(is.null(e)){ # Use fully Gaussian approximation
    mu <- tau * (1 + delta2) * traces$trace
    sigma <- tau * sqrt((2 + 4*delta2) * traces$hsnorm2)

    return(function(obs, lower.tail = FALSE){
      if(length(obs) != 1){stop("while this function could be easily parallelized,
                              to keep the code simple, for now it only takes one obs at a time.")}
      -pnorm(obs, mean = mu, sd = sigma,
             lower.tail = lower.tail, log.p = TRUE)/log(10)
    })
  }

  if(length(e$values)==length(traces$diag) |
     (traces$hsnorm2 - sum(e$values^2)) <= 0 ){
    # bounds are not needed because either we have all of the eigenvalues or the eigenvalues that
    # are in the remainder are non-zero due to numerical imprecision.

    return(function(obs, lower.tail = FALSE){
      gauss.tcdf <- QFGauss(tau * e$values)
      if(length(obs) != 1){stop("while this function could be easily parallelized,
                              to keep the code simple, for now it only takes one obs at a time.")}
      rep(-gauss.tcdf(obs, lower.tail = lower.tail, log.p = TRUE)/log(10),3)
    })
  }


  # Calc Required Traces
  R.max.abs.eta <- tau * min(abs(e$values))
  R.sum.eta <-  tau * (traces$trace - sum(e$values))
  R.sum.etasq <- tau^2 * (traces$hsnorm2 - sum(e$values^2))
  if(R.sum.etasq < 0){stop("Squared HS norm of matrix is less than sum of eigenvalues squared -- check that the provided traces are from the same matrix as the eigenvalues")}
  R.sum.eta.deltasq <-  delta2 * R.sum.eta
  R.sum.etasq.deltasq <- delta2 * R.sum.etasq


  # Point Estimate Function
  gauss.tcdf <- suppressWarnings(QFGauss(tau * e$values, sigma = sqrt((2 + 4*delta2)*R.sum.etasq)))
  E.R <- (1 + delta2) * R.sum.eta
  gauss.approxfullcdf <- function(x, density = FALSE, lower.tail = TRUE, log.p = FALSE) gauss.tcdf(x-E.R, density = density, lower.tail = lower.tail, log.p = log.p)
  attr(gauss.approxfullcdf,"mu") <- attr(gauss.tcdf,"mu") + E.R
  attr(gauss.approxfullcdf,"Q.sd") <- attr(gauss.tcdf,"Q.sd")

  if(only.point.est){

    function(obs, lower.tail = FALSE){
      if(length(obs) != 1){stop("while this function could be easily parallelized,
                              to keep the code simple, for now it only takes one obs at a time.")}
      -gauss.approxfullcdf(obs, lower.tail = lower.tail, log.p=T)/log(10)
    }

  } else {

    # Bound Function
    tcdf <- suppressWarnings(QFGauss(e$values))
    bound.func <- QFGaussBounds(tcdf,"identity", R.max.abs.eta, R.sum.eta, R.sum.etasq, R.sum.eta.deltasq, R.sum.etasq.deltasq)

    function(obs, lower.tail = FALSE){
      if(length(obs) != 1){stop("while this function could be easily parallelized,
                              to keep the code simple, for now it only takes one obs at a time.")}
      -c(log(bound.func(obs)[1,1:2 + if(lower.tail){0}else{2}]),gauss.approxfullcdf(obs,lower.tail = lower.tail,log.p=T))/log(10)
    }

  }
}




# Takes in matrix x
calc_obs <- function(M,x, fwd, bck,pars, Q, tQ, from_recipient = 1, nthreads = 1) {
  x <- x - Q %*% (tQ %*% x)
  res <- kalis:::DedipPar(M, x, fwd, bck, pars, from_recipient = from_recipient, nthreads = nthreads)
  list("min_min"  = sum(x*res$min_min),
       "min_mean" = sum(x*res$min_mean),
       "min_max"  = sum(x*res$min_max),
       "min2nd"   = sum(x*res$min2nd))
}


calc_traces <- function(M, Q, sym.M = FALSE,
                        from_recipient = 1, nthreads = 1){
  if(sym.M){M <- 0.5*(M + t(M))}
  J <- crossprod(Q, M)
  tX <- t((Q %*% (J%*%Q)) - (M %*% Q))

  kalis:::CalcTraces(M,tX,t(Q),J,from_recipient,nthreads)
}



QFBinary <- function(Ey,M){
  # quick implementation of diagonal correction for binary random variables

  # Ey is the expected value of the vector y \in {0,1}
  # M a symmetric N x N  matrix

  M.diag <- diag(M)

  u <- M.diag*(2*Ey-1)/sqrt(4*Ey*(1-Ey))

  diag(M) <- 0
  e <- eigen(M,symmetric = T)

  epsilon <- c(crossprod(e$vectors,u))

  raw.cdf <- QForm::QFGauss(e$values,epsilon/e$values)

  C <- sum( (epsilon^2)/e$values - M.diag )

  function(x, density = FALSE, lower.tail = TRUE, log.p = FALSE){
    raw.cdf(x + C ,density,lower.tail,log.p)}
}

CalcNNWeights <- function(n, # number of haplotypes in data
                          maf.cutoff, # desired maf cut-off for weighting alleles
                          type = "coal",
                          x = NULL, # allele count (x and y encode AFS)
                          y = NULL # number of times that allele count is present
){

  max.allele.count <- ceiling(n * maf.cutoff) # largest allele count we're looking to weight

  if(type == "inverse"){
    res <- rev(1/1:(max.allele.count-1))
  } else if(type == "unif"){
    res <- cumsum(rev(1/1:(max.allele.count-1)))
  } else if(type == "skat"){
    res <- cumsum(rev(dbeta(2:max.allele.count/n,1,25)/(1:(max.allele.count-1))))
  } else if(type == "coal"){
    res <- cumsum(rev(1/(2:max.allele.count * 1:(max.allele.count-1))))
  } else if(type == "afs"){
    if(is.null(x) | is.null(y)){ stop("both x and y must be supplied if it is desired to calc weights using allele counts.")}
    # line below ensures this will extrapolate out to doubletons and max.allele.count even if none in data
    f.hat <- approx(max.allele.count - rev(x),cummax(rev(y)),0:c(max.allele.count-2),
                    method="constant",
                    rule=2L)
    f.hat$x <- max.allele.count - f.hat$x
    res <- cumsum(f.hat$y/(f.hat$x-1))
  } else {
    stop("type must be coal, unif, skat, or afs")
  }

  # res[i] is the weight to give to NN size i.  Weights normalized to 1e3
  rev(res)*(1e3/sum(res))
}


CalcNNWeights2 <- function(x, afs = NULL, mu = 1e-8, n = NULL, m = 1, weight.by.eff.num.vars = TRUE, eta = rep(1,length(x))){

  # this code assumes that the distance from a recipient to itself is encoded in x as an NA and that that is the only NA!

  if(sum(is.na(x))!=1){stop("x should have one NA entry that corresponds to the recipient's distance to itself.")}

  N = length(x)

  if(is.null(afs)){
    pi <- c(0, 1/2:(N-1), 0) # standard coalescent
    pi <- pi / sum(pi)
  } else {
    if(length(afs)!=N){stop("afs must have the same length as N, the number of haplotypes.  afs[N] should be zero.")}
    afs[1] <- 0
    afs[N] <- 0
    pi <- afs / sum(afs)
  }

  r <- data.table::frank(x, na.last = FALSE, ties.method = "random")
  d <- x[order(r)]

  d[1] <- d[2]
  Z <- d[N] - d[1]
  theta <- c(diff(d),0)
  theta <- theta / sum(theta)

  # if(is.null(n)){ n <- -Z/log(mu) }
  #
  # phi <- theta^n * pi^m
  #theta <- theta
  #theta <- theta / sum(theta)

  phi <- n * theta + m * pi

  #phi <- phi / sum(phi)

  w <- c(0,rev(cumsum(rev(eta[-1] * phi[-1] / 1:(N-1)))))

  if(weight.by.eff.num.vars){ w <- w * (n+m) }

  w[r]
}

#x <- c(4,2,8,NA,3,1,4,6)
#CalcNNWeights2(x)
phi2 <- function(x,z){

  el1 <- exp(z[1] + z[2] * (x>0) + z[3] * x )
  1-1/(1+el1)
}
phi3 <- function(x,z){
  el1 <- exp(z[1] + z[2] * (x>0) + z[3] * x)
  el2 <- exp(z[4] + z[5] * (x>0) + z[6] * x)
  1-1/( (1+el1) * (1+el2) )
}

phi2.logdac <- function(x,logdac,z){
  el1 <- exp(z[1] + z[4] * logdac + (z[2] + z[5] * logdac) * (x>0) + (z[3] + z[6] * logdac) * x)
  1-1/(1+el1)
}

phi3.logdac <- function(x,logdac,z){
  el1 <- exp(z[1] + z[4] * logdac + (z[2] + z[5] * logdac) * (x>0) + (z[3] + z[6] * logdac) * x)
  el2 <- exp(z[7] + z[10] * logdac + (z[8] + z[11] * logdac) * (x>0) + (z[9] + z[12] * logdac) * x)
  1-1/( (1+el1) * (1+el2) )
}

CalcNNWeights3 <- function(x = NULL, # a column of distances
                           r = NULL, # length(r) = length(x) = N, rank
                           dd = NULL, # difference of ordered distances. length(dd) = N-1
                           daf.range = c(0,1),
                           include.theta = TRUE,
                           theta.cutoff = 0,
                           dd.theta.ctl = FALSE,
                           logistic.regression.phi = FALSE,
                           phi.pars = NULL,
                           da.cdf = NULL, # da.cdf[j] = proportion of derived alleles present in j copies.  Thus, da.cdf is N-1 long.
                           alpha = 0,
                           beta = 0,
                           eta = NULL, # eta overrides maf.cutoff, should be N-1 long, with
                           # eta[j] = the probability that a variant present in j copies is causal given it is present
                           precalc.logdac = if(!is.null(x)){log(1:(length(x)-1))}else{log(1:(length(r)-1))},
                           zero.dd.zero.phi = TRUE
){
  # this code assumes that the distance from a recipient to itself is encoded in x as an NA and that that is the only NA!
  if(!(is.null(r) | is.null(dd))){
    N <- length(r)
  } else if(!is.null(x)){
    if(sum(is.na(x))!=1){stop("x should have one NA entry that corresponds to the recipient's distance to itself.")}
    N <- length(x)

    r <- data.table::frank(x, na.last = FALSE, ties.method = "random")
    dd <- diff(x[order(r)]) # x is reordered and dd was calculated. length(dd) = N-1
    dd[1] <- 0 # for now all singletons are excluded and this gets rid of the leading NA.
  } else {
    stop("Either x needs to be specified OR {r, dd} need to be specified.")
  }
  if(include.theta){
    if (dd.theta.ctl){
      phi <- dd
    } else{ if (logistic.regression.phi){
      if(phi.pars$model.family == "binomial"){
        if(phi.pars$formula.type == "logdac"){
          phi <- phi2.logdac(dd,precalc.logdac,phi.pars$pars)
        } else {
          phi <- phi2(dd,phi.pars$pars)
        }
      } else {
        if(phi.pars$formula.type == "logdac"){
          phi <- phi3.logdac(dd,precalc.logdac,phi.pars$pars)
        } else {
          phi <- phi3(dd,phi.pars$pars)
        }
      }
    } else {
      phi <- as.integer(dd > theta.cutoff)
      }
    }

    if(zero.dd.zero.phi){ phi[dd==0] <- 0 }
    # if(mac.cutoff<=(N-1)){
    #   phi[mac.cutoff:(N-1)] <- 0 # don't count any mutations as potentially causal that we assume cannot be b/c they're above the theshold
    # }
  } else {

    phi <- rep(0,N-1)
  }


  if(!is.null(da.cdf) & (alpha!=0 | beta!=0)){

    if(length(da.cdf) != N - 1){stop("Length of da.cdf must be length(x)-1.")}
    # Instead of renormalizing da.cdf like below, build cut.off into eta and thereby give less weight to
    # recipients who are less likely to carry a rare variant.
    # da.cdf <- diff(da.cdf[1:ceiling(N * maf.cutoff)])
    # da.cdf <- da.cdf / sum(da.cdf)
    # da.cdf <- c(0,cumsum(da.cdf),rep(1,N -ceiling(N * maf.cutoff)))
    # indexing for s is +1 def of s in the tex: s[0] in latex = s[1] in R. s has length N.
    # s <- rep(0,N)
    # s[N] <- N-1
    jump.points <- which(dd > 0) # dd should always >= 0 because it is calculated on sorted matrix
    if(length(jump.points)){
    pi <- rep(0,N-1)
    # Add on a jump at the end of the end of the CDF to assign weight to potential variants above the oldest jump
    if(jump.points[length(jump.points)]!=(N-1)){jump.points <- c(jump.points,N-1)}
    pi[jump.points]<- c(da.cdf[jump.points[1]],diff(da.cdf[jump.points]))
    } else {
      pi <- rep(1/(N-1),N-1)
    }

    # if(mac.cutoff<=(N-1)){
    #   pi[mac.cutoff:(N-1)] <- 0 # don't account for any hidden mutations above the causal variant threshold
    # }
    # diff(da.cdf[c(1,jump.points+1)])
    # jump.points <- jump.points[jump.points < ceiling(N * maf.cutoff)]
    # s[1 + jump.points] <- jump.points
    # s <- cummax(s)
    # if(ceiling(N * maf.cutoff)<=(N-1)){s[N] <- s[N-1]}

    phi <- phi + (alpha + beta * sum(phi)) * pi #(da.cdf[s[-1]+1]-da.cdf[s[-N]+1])

  }


  if(!is.null(eta)){ # eta overrides maf.cutoff
    phi <- eta * phi

  } else {

    if(!is.null(daf.range)){

      # we KEEP weights on variants that are present in a number of copies C such that
      # maf.range[1] * N <= C <= maf.range[2] * N

      lower.mac.cutoff <- ceiling(daf.range[1]*N) - 1
      if(lower.mac.cutoff >= 1){
        phi[1:lower.mac.cutoff] <- 0
      }

      upper.mac.cutoff <- floor(daf.range[2]*N) + 1
      if(upper.mac.cutoff<=(N-1)){
        phi[upper.mac.cutoff:(N-1)] <- 0
      }

    }
  }

  # be default, all neighborhood sides are evenly weighted.
  w <- c(0,rev(cumsum(rev(phi[-1] / 1:(N-2)))),0) # any variants where we are dividing by N-1 would be monomorphic so it always gets 0 weight

  w[r]
}
#

#x <- sort(c(1,10,4.2,5.8,5.8,4.5,7,3.9,3.9,9,4,NA,4,2),na.last = FALSE)
# plot(x,CalcNNWeights3(x,da.cdf = seq(0,1,len=length(x)-1),include.theta = TRUE,dd.theta.ctl = TRUE)) #,da.cdf = seq(0,1,len=length(x)))#,alpha=1,beta=1)
# plot(x,CalcNNWeights3(x,da.cdf = seq(0,1,len=length(x)-1),include.theta = TRUE,dd.theta.ctl = TRUE,daf.range =c(0.2,0.6))) #,da.cdf = seq(0,1,len=length(x)))#,alpha=1,beta=1)
# plot(x,CalcNNWeights3(x,da.cdf = seq(0,1,len=length(x)-1),include.theta = TRUE,theta.cutoff=0.32)) #,da.cdf = seq(0,1,len=length(x)))#,alpha=1,beta=1)
# plot(x,CalcNNWeights3(x,da.cdf = seq(0,1,len=length(x)-1),include.theta = TRUE,theta.cutoff=0.32,daf.range =c(0.2,0.6))) #,da.cdf = seq(0,1,len=length(x)))#,alpha=1,beta=1)



HapDAFS <- function(pop.labels = NULL){

  # by default, everyone belongs to their own population
  # we can linearly interpolate the density between non-zero points to get a somewhat smooth CDF.


  if(is.null(N()) | is.null(L())){stop("Unable to calc DAFS: no haplotypes found loaded in cache.")}

  # Calc derived allele counts (in blocks of loci to reduce striding and preserve memory)
  loci.blocks <- seq(1,L(),by=1e3)
  loci.blocks <- cbind(loci.blocks,c(loci.blocks[-1]-1,L()))
  dac <- rep(0L,L()) # derived allele count
  for(i in 1:nrow(loci.blocks)){ dac[loci.blocks[i,1]:loci.blocks[i,2]] <- c(rowSums(QueryCache(loci.idx = loci.blocks[i,1]:loci.blocks[i,2])))}


  if(is.null(pop.labels)){

    calc_dafs <- function(x){
      y <- tabulate(x,nbins=N())
      y / sum(y)
    }

    hap.blocks <- seq(1,N(),by=1e3)
    hap.blocks <- cbind(hap.blocks,c(hap.blocks[-1]-1,N()))
    dafs <- matrix(0,N(),N())
    for(i in 1:nrow(hap.blocks)){
      dafs[,hap.blocks[i,1]:hap.blocks[i,2]] <- apply(X = QueryCache(hap.idx = hap.blocks[i,1]:hap.blocks[i,2]) * dac,
                                                      MARGIN = 2, FUN = calc_dafs)
    }

  } else {
    pops <- unique(pop.labels)
    dafs <- matrix(0,N(),length(pops))
    colnames(dafs) <- pops
    for(i in 1:length(pops)){
      slab <- QueryCache(hap.idx = which(pop.labels == pops[i]))
      num.derived.carried <- colSums(slab)

      order.dac <- order(dac)
      dac.ordered <- dac[order.dac]
      weight.per.variant <- c(rowMeans(sweep(slab,2,num.derived.carried,"/")))[order.dac]
      weight.per.dac <- tapply(weight.per.variant,dac.ordered,sum)
      temp.dafs.col <- approx(as.numeric(names(weight.per.dac)),cumsum(as.numeric(weight.per.dac)), xout=1:N(),rule=1:2)$y
      temp.dafs.col[is.na(temp.dafs.col)] <- 0
      dafs[,i] <- c(temp.dafs.col[1],diff(temp.dafs.col))
    }

  }

  dafs
}


plot.daf.matrix <- function(mat,pop.labels = rep(1,ncol(mat))){
  unique.pop.labels <- unique(pop.labels)
  temp.colors <- RColorBrewer::brewer.pal(8,"Dark2")[1:length(unique.pop.labels)]

  selected <- c()
  for(i in unique.pop.labels){
    candidates <- which(pop.labels==i)
    if(length(candidates)>1){
      selected <- c(selected,sample(candidates,size = min(length(candidates),100)))
    } else {
      selected <- c(selected,candidates)
    }
  }

  plot(cumsum(mat[,selected[1]]),type="l",col=temp.colors[which(unique.pop.labels ==pop.labels[selected[1]])],
       main="Hap Specific Derived Allele Count Distributions",
       las=1,bty="n",ylab="ECDF",xlab="Derived Allele Count",ylim=c(0,1))
  for(i in selected[-1]){lines(cumsum(mat[,i]),col= temp.colors[which(unique.pop.labels ==pop.labels[i])])}
  legend("topleft",legend = unique.pop.labels, fill=temp.colors,border=temp.colors, bty="n")
  # xx <- c(0,1/(2:ncol(mat)))
  # lines(cumsum(xx) /sum(xx),lwd=1.8,lty=2)
  #xx <- c(0,rep(1,ncol(mat)-1))
  #lines(cumsum(xx) /sum(xx),lwd=1.8,lty=2)
  segments(0,0,N(),1,lty=2,lwd=1.8)
}


find_sprigs_old <- function(x, inclusive = FALSE){
  # x here is a list that's N long st x[[i]] gives the indices of the (tied) nearest neighbors of i
  roster <- rep(NA_integer_,length(x))

  if(inclusive){
    last.label <- 0L

    for(i in 1:length(x)){
      temp <- roster[x[[i]]]

      # combine the NAs
      temp2 <- na.omit(unique(temp))

      if(length(temp2)==1){
        label <- temp2
      } else {
        last.label <- last.label + 1L
        label <- last.label
      }
      roster[c(i,x[[i]][is.na(temp)])] <- label

      if(length(temp2)>1){
        roster[roster %in% temp2] <- label
      }
    }

  } else {

    label <- 0L
    # add self to own neighborhood
    x <- mapply(c,x,1:length(x))

    done <- rep(FALSE,length(x))
    to.prune <- rep(NA_integer_,length(x))

    for(i in 1:length(x)){
      if(!done[i]){

        # pulling out clique of the graph that are fully connected bi-directionally:
        # if i is in a clique, rather trivially, this will return the full clique
        # Note, we require i %in% proposed.set to prevent called cliques from being broken up later in the for loop:
        # if i is not in a clique, then it's still possible for a partial clique to be returned that doesn't include i if i projects onto
        # a superset or subset of a clique.  If i superseeds a clique member and projects
        # onto a subset, this clique subset will be overwritten later by the larger clique.  However, it would still be possible for a i that comes
        # after all of the clade members in our for loop to break up the clique by projecting onto a subset of them.
        # Enforcing i %in% proposed.set avoids that possibility.

        # we also require that length(proposed.set) > 1 so that we don't end up with solo cliques being called that are just i by itself.

        proposed.set <- Reduce(intersect,x[x[[i]]])
        if(length(proposed.set) > 1 && i %in% proposed.set){
          label <- label + 1L
          roster[proposed.set] <- label
          done[proposed.set] <- TRUE
          to.prune[proposed.set[sapply(x[proposed.set], setequal, y = proposed.set)]] <- length(proposed.set)
        }
      }
    }

    # individuals that are not part of a fully connected clique are left with NA_integer_ on the roster
  }

  # Size frequency spectrum: table(table(roster))

  attr(roster,"to.prune") <- to.prune
  roster
}


find_sprigs <- function(x){
  # x here is a list that's N long st x[[i]] gives the indices of the (tied) nearest neighbors of i
  roster <- rep(NA_integer_,length(x))

  label <- 0L
  # add self to own neighborhood
  x <- mapply(c,x,1:length(x))

  done <- rep(FALSE,length(x))
  to.prune <- rep(NA_integer_,length(x))

  # the randomness in indices here is not really essential but
  for(i in sample.int(length(x))){
    if(!done[i]){

      # pulling out cliques in the graph that are fully connected bi-directionally:
      # if i is in a clique, rather trivially, this will return the full clique
      # Note, we require i %in% proposed.set to prevent called cliques from being broken up later in the for loop:
      # if i is not in a clique, then it's still possible for a partial clique to be returned that doesn't include i if i projects onto
      # a superset or subset of a clique.  If i supercedes a clique member and projects
      # onto a subset, this clique subset will be overwritten later by the larger clique.  However, it would still be possible for a i that comes
      # after all of the clade members in our for loop to break up the clique by projecting onto a subset of them.
      # Enforcing i %in% proposed.set avoids that possibility.

      # we also require that length(proposed.set) > 1 so that we don't end up with solo cliques being called that are just i by itself.

      #missing_sprig_6 <- c(6103,1804, 6015, 4726, 4752, 807,3118,3991,6466,6068,  10,1250, 3669, 3658, 1997, 1399, 1116, 3738, 5015)
      proposed.set <- Reduce(intersect,x[x[[i]]])
      # in case the neighborhood of i overshoots into previously established cliques, it has a BIG effect in real data
      proposed.set <- proposed.set[!done[proposed.set]]

      if(length(proposed.set) > 1 && i %in% proposed.set){

        label <- label + 1L
        # this repeated intersection step has truly a small effect but
        # helps guard us against the case where i might have erroneously added some candidates in its neighborhood that
        # do not belong in the clade and do not include some clade members.  This steps helps us recover those clade members
        proposed.set <- Reduce(intersect,x[proposed.set])
        proposed.set <- proposed.set[!done[proposed.set]]

        # if(!all(is.na(roster[missing_sprig_6])) && !all(roster[missing_sprig_6]==6L)){
        #   print(i)
        #   print(label)
        #   print(roster[missing_sprig_6])
        #   browser()
        # }
        roster[proposed.set] <- label
        done[proposed.set] <- TRUE
      }
    }

    # individuals that are not part of a fully connected clique are left with NA_integer_ on the roster
  }

  # Size frequency spectrum: table(table(roster))

  attr(roster,"n.sprigs") <- label
  roster
}
#
# start.old <- proc.time()
# old.sprigs <- find_sprigs(nM)
# proc.time() - start.old
# #
#
# start.new1 <- proc.time()
# new.sprigs.1 <- find_sprigs_test_ver(nM)
# proc.time() - start.new1
# #
#
# start.new2 <- proc.time()
# new.sprigs.2 <- find_sprigs_test_ver(nM,TRUE)
# proc.time() - start.new2
# #
#
# old.sprigs[is.na(old.sprigs)] <- 0
# new.sprigs.1[is.na(new.sprigs.1)] <- 0
# new.sprigs.2[is.na(new.sprigs.2)] <- 0
#
# told <- table(old.sprigs) # missing clade 6
# t1 <- table(new.sprigs.1) # has 19 members in clade 6
# t2 <- table(new.sprigs.2) # has 29 members in clade 6
#
#
# plot(old.sprigs,new.sprigs.1) # new sprigs 1 does not allow reassignment of folks who are already in a sprig
# abline(0,1)
# told
# t1
# t2
#
# plot(new.sprigs.1,new.sprigs.2)
# abline(0,1)
# plot(t1,t2)
# abline(0,1)

# sprigs1 <- find_sprigs(nM)
# sprigs2 <- find_sprigs(nM)
# sprigs3 <- find_sprigs(nM)
# table(sapply(tapply(sprigs,sprigs2,unique),length))
# table(sapply(tapply(sprigs2,sprigs3,unique),length))

#Testing find_sprigs
# find_sprigs(list(
#   1:5,
#   3:7,
#   1:10,
#   1:10,
#   1:10,
#   5:11
# ))

find_irreducible_components <- function(m, keep.isolated.nodes = TRUE){
  # m is the adjacency matrix of the graph of 1s and 0s or TRUE and FALSE.

  cliques <- list()
  done <- rep(FALSE,nrow(m))

  for(i in 1:nrow(m)){
    if(!done[i]){

      proposed.set <- boundary.set <- i

      repeat{

        # update boundary set
        if(length(boundary.set)==1){
          boundary.set <- setdiff(which(m[boundary.set,]>0),proposed.set)
        } else {
          boundary.set <- setdiff(which(colSums(m[boundary.set,])>0),proposed.set)
        }

        # update proposed set
        if(length(boundary.set)){
          proposed.set <- c(boundary.set,proposed.set)
        } else {
          break
        }

      }

      done[proposed.set] <- TRUE

      if(!keep.isolated.nodes){ if(length(proposed.set)==1){next} }

      cliques[[length(cliques)+1]] <- proposed.set
    }
  }
  cliques
}


#testin find_irreducible_components
# demo.m <- matrix(0,10,10)
# demo.m[sample(1:10,4),sample(1:10,4)] <- 1
# demo.m <- demo.m + t(demo.m)
# demo.m[demo.m > 1] <-1
# demo.m[1:3,1:3] <- 1
# demo.m[9:10,9:10] <- 1
# find_irreducible_components(demo.m)

haps2genotypes <- function(haps, # p x N matrix of 1s and 0s
                           ploidy = 1L,
                           inheritance.model = "additive"){
  ploidy.in <- ploidy
  ploidy <- as.integer(ploidy)
  if(ploidy!=ploidy.in || ploidy <= 0){stop("ploidy must be a positive integer.")}

  N <- ncol(haps)

  if(ploidy > 1){
    genotypes <- haps[,seq(1,N,by = ploidy)]
    for(j in 2:ploidy){
      if(inheritance.model == "additive"){
        genotypes <- genotypes + haps[,seq(j,N,by = ploidy)]
      } else if(inheritance.model == "dominant"){
        genotypes <- genotypes | haps[,seq(j,N,by = ploidy)]
      } else {
        stop("inheritance.model mis-specified: for now, must be additive or dominant")
      }
    }
    storage.mode(genotypes) <- "integer"
  } else {
    genotypes <- haps
  }

  if(! "matrix" %in% class(genotypes)){
    genotypes <- matrix(genotypes,ncol=N/ploidy)
  }

  genotypes # a p x n matrix
}


new_renyi <- function(x){
  signs <- sign(x)
  x <- x^2
  x.rank <- rank(x, ties.method = "random")

  x <- qexp(pchisq(x[order(x.rank)], df = 1, log.p = TRUE), log.p = TRUE)
  x <- c(x[1],diff(x)) * seq.int(length(x),1)
  r <- list("n" = length(x), "signs" = signs, "ranks" = x.rank, "exps" = x)
  class(r) <- c("renyi","list")
  return(r)
}

print.renyi <- function(x, digits = getOption("digits")){
  if(is.null(digits)){
    digits <- getOption("digits")
  }
  cat("", "Renyi Representation\n", "  signs :", sep = "")
  utils::str(x$signs)
  cat("", "  ranks :", sep = "")
  utils::str(x$ranks, digits.d = digits)
  cat("", "  i.i.d. exponential r.v.s :", sep = "")
  utils::str(x$exps, digits.d = digits)
  invisible(x)
}



# If exclude.idx was indexing the Renyi instead of the target Gaussian vector
# for excluding some to.exclude subset:
# use exclude.idx = r$ranks[to.exclude] (it must be the indices of the exponential r.v.s to exclude)
# for excluding the top k
# use exclude.idx = seq.int(length(r$exps) - k + 1,length(r$exps))



inverse.renyi <- function(x, exclude.idx = NULL){
  # exclude.idx are the indices in the target Gaussian vector that we want to remove
  # that means that for excluding the top k:
  # use exclude.idx = which(x$ranks >= x$n - k + 1)
  # In order to exclude some set based on indices corresponding to the iid exponentials use:
  # exclude.idx = which(x$ranks %in% exps.to.exclude)

  x <- validate_renyi(x)

  if(!is.null(exclude.idx) && length(exclude.idx)){
    x$exps <- x$exps[-x$ranks[exclude.idx]]
    x$ranks[exclude.idx] <- NA
    x$ranks <- rank(x$ranks, na.last = "keep")
  }

  x$signs * sqrt(qchisq(pexp(cumsum(x$exps / seq.int(length(x$exps),1)), log.p = TRUE), df = 1, log.p = TRUE)[x$ranks])
}

validate_renyi <- function(x){
  # check validity of input renyi representation
  if(length(x$signs) != length(x$exps) | length(x$ranks) != length(x$exps)){
    stop("Length of all 3 vectors in a Renyi representation X must be equal.",call. = FALSE)}
  x
}

renyi <- function(x){

  if(!is.numeric(x)){stop("x must be numeric")}
  if(any(!is.finite(x))){stop("all elements of x must be finite with no NAs")}

  r <- new_renyi(x)
  validate_renyi(r)
}


renyi.outlier.test <- function(x, k = 2^(0:6), log.p = FALSE){
  pgamma(cumsum(rev(x$exps)[1:max(k)])[k],shape = k,lower.tail = FALSE, log.p = log.p)}
#
# # testing
# z0 <- rnorm(1e4)
# r <- renyi(z0)
#
# z1 <- inverse.renyi(r)
# all.equal(z0,z1)
#
# z2 <- inverse.renyi(r, exclude.idx = sample.int(r$n,size = 9e3))
# plot(z0,z2); abline(0,1)
#
# y <- rnorm(1e3)
# r <- renyi(y)
# k <- 2^(0:6)
# signif.trace <- -renyi.outlier.test(r,k,log.p = TRUE)/log(10)
# k.max <- k[which.max(signif.trace)]
# updated.y <- inverse.renyi(r, exclude.idx = which(r$ranks >= r$n - k.max + 1))
# plot(y,updated.y); abline(0,1)


rank2gauss <- function(x){
  # x is a vector of ranks with no ties
  # See ?quantile Type 9 to help justify below transformation
  qnorm((seq_len(length(x)) - 3/8)/(length(x) + 1/4))[x]
}

new_clade_assignments <- function(clades,sprigs.1,sprigs.2,inactive.sprigs, overwrite.single.copy.carriers = TRUE){

  inactive.sprigs.1.ind <- sprigs.1 %in% inactive.sprigs
  inactive.sprigs.2.ind <- sprigs.2 %in% inactive.sprigs

  if(overwrite.single.copy.carriers){
    res <- list("assign.to.sprigs.1" = which(inactive.sprigs.2.ind & !inactive.sprigs.1.ind),
                "assign.to.sprigs.2" = which(inactive.sprigs.1.ind & !inactive.sprigs.2.ind))
  } else {
    res <- list("assign.to.sprigs.1" = which(inactive.sprigs.2.ind & !inactive.sprigs.1.ind & !is.na(sprigs.1)),
                "assign.to.sprigs.2" = which(inactive.sprigs.1.ind & !inactive.sprigs.2.ind & !is.na(sprigs.2)))

    dual.btw.inactive.idx <- which(inactive.sprigs.1.ind & inactive.sprigs.2.ind)

    if(length(dual.btw.inactive.idx)){
      assign.vec <- sample(c(FALSE,TRUE),length(dual.btw.inactive.idx),replace = TRUE)

      temp.assign.to.sprigs.1 <- c()
      temp.assign.to.sprigs.2 <- c()

      for(i in 1:length(dual.btw.inactive.idx)){
        if(assign.vec[i]){
          temp.assign.to.sprigs.1 <- c(temp.assign.to.sprigs.1,dual.btw.inactive.idx[i])
        } else {
          temp.assign.to.sprigs.2 <- c(temp.assign.to.sprigs.2,dual.btw.inactive.idx[i])
        }
      }
      res$assign.to.sprigs.1 <- c(res$assign.to.sprigs.1,temp.assign.to.sprigs.1)
      res$assign.to.sprigs.2 <- c(res$assign.to.sprigs.2,temp.assign.to.sprigs.2)
    }
  }
  res
}

solve_active_nodes <- function(x, weights){

  if(nrow(x) >= 16){
    # as a fail safe for now in case the graph component is large, we just set all nodes to active
    # and then inactivate them one at a time from highest degree to lowest degree nodes until we find a solution
    config <- rep(1L,nrow(x))
    priority.to.remove <- nrow(x) - rank(colSums(x),ties.method = "random") + 1

    for(i in 1:nrow(x)){
      config[which(priority.to.remove == i)] <- 0L
      if(sum(config * (x%*%config)) == 0){break}
    }
    return(as.logical(config))
  }


  x <- as.matrix(x)
  storage.mode(x) <- "integer"
  configs <- as.matrix(expand.grid(replicate(n = nrow(x),c(0L,1L),simplify = FALSE)))

  valid.config.ind <- rowSums(configs * (configs%*%x)) == 0

  weights <- (1:nrow(x))^4

  config.scores <- configs[valid.config.ind,] %*% weights

  as.logical(configs[which(valid.config.ind)[which.max(config.scores)],])
}

calc_sprig_phenotypes <- function(y.resids,sprigs,n.unique.sprigs,ploidy){

  if(all(is.na(sprigs))){stop("all provided sprigs are NA")}

  if(ploidy == 1){

    sprigs.to.remove.sizes <- tabulate(sprigs, nbins = n.unique.sprigs)
    sprig.coefficient <- 1/sqrt(sprigs.to.remove.sizes)


    return(list("skip.renyi" = FALSE,
                "y.sprigs" = c(tapply(y.resids,
                                      sprigs,
                                      function(x){sum(x)/sqrt(length(x))})),
                "renyi.sprigs" = sprigs,
                "non.renyi.sprigs" = integer(),
                "sprigs.to.remove.sizes" = sprigs.to.remove.sizes,
                "sprig.coefficient" = sprig.coefficient))

  }

  if(ploidy!=2){
    stop("currently only ploidy = 1 or 2 supported")
  }

  ###########################################
  ###########################################
  # Make Unambiguous Clade Assignments
  ###########################################
  ###########################################


  sprigs.1 <- sprigs[seq.int(1,length(sprigs),2)]
  sprigs.2 <- sprigs[seq.int(2,length(sprigs),2)]

  # build clade assignments based on unambiguous sprig assignment:
  # either both alleles are assigned to the same sprig or only one allele is assigned to a sprig
  renyi.sprigs <- pmax(sprigs.1,sprigs.2,na.rm = T)
  dual.samples.idx <- which(sprigs.1 != sprigs.2)
  renyi.sprigs[dual.samples.idx] <- NA_integer_

  # the resulting renyi.sprigs here may be missing some nodes/sprigs, meaning that
  # there are no samples unambiguously assigned to those nodes/sprigs (these are nodes in our initial network with weight 0)
  # instead of trying to keep those nodes in the Renyi test by randomly assigning some of their edges to them, at least for now
  # we just remove those nodes.
  # If an edge in our graph connects a zero weight node to a positively weighted node, we assign that edge (that dual sample)
  # to the positively weighted node. In the rare case that an edge in our graph connects two zero weight nodes, we simply
  # say that edge doesn't belong to a clade for Renyi testing and test those nodes/clades/sprigs in the quadratic form
  # This procedure ensures that the resulting graph will only consist of positively weighted nodes that are still connected via
  # edges to some other nodes



  # these are the sprigs/nodes without anyone unambiguously assigned to them which will need to be tested in the quadratic form
  zero.weighted.sprigs <- which(tabulate(renyi.sprigs,n.unique.sprigs) == 0)

  if(length(zero.weighted.sprigs)==n.unique.sprigs){
    warning("All sprigs had zero weight")
    return(list("skip.renyi" = TRUE))
  }



  temp.assign.list <- new_clade_assignments(renyi.sprigs,sprigs.1,sprigs.2,inactive.sprigs = zero.weighted.sprigs)

  renyi.sprigs[temp.assign.list$assign.to.sprigs.1] <- sprigs.1[temp.assign.list$assign.to.sprigs.1]
  renyi.sprigs[temp.assign.list$assign.to.sprigs.2] <- sprigs.2[temp.assign.list$assign.to.sprigs.2]

  # Remove samples from dual.samples.idx that have now been unambiguously assigned
  dual.samples.idx <- setdiff(dual.samples.idx,union(temp.assign.list$assign.to.sprigs.1,
                                                     temp.assign.list$assign.to.sprigs.2))
  #length(dual.samples.idx)
  #dual.samples.idx[(sprigs.1[dual.samples.idx] %in% zero.weighted.sprigs & sprigs.2[dual.samples.idx] %in% zero.weighted.sprigs)]

  # Remove samples from dual.samples.idx that were edges between two zero.weighted.sprigs
  dual.samples.idx <- dual.samples.idx[!(sprigs.1[dual.samples.idx] %in% zero.weighted.sprigs & sprigs.2[dual.samples.idx] %in% zero.weighted.sprigs)]
  #length(dual.samples.idx)

  # dual.samples.idx now indicates samples that form edges between two positively weighted nodes and need to be assigned/resolved

  # In summary:
  # Q: what do you do when a sprig is in zero.weighted.sprigs (ie: all samples under a sprig are ambiguously assigned) ?
  # A: we assign NA to that sprig's putative.sprig.y, immediately remove it from Renyi testing, and test it only in the quadratic form
  # Q: what do you do when a doubleton has one homozygous sample under it?
  # A: we just keep it and test it in the Renyi


  # this should always be true: renyi.sprigs at this stage includes all sprigs that will considered as candidates from Renyi screening
  #all.equal(sort(union(na.omit(unique(renyi.sprigs)),zero.weighted.sprigs)),1:200)



  ###########################################
  ###########################################
  # Sort Dual Samples Into Clades
  ###########################################
  ###########################################

  # generate graph
  sprig.dag <- as(Matrix::sparseMatrix(
    sprigs.1[dual.samples.idx],
    sprigs.2[dual.samples.idx],dims = rep(n.unique.sprigs,2),symmetric = T),"ngCMatrix")


  # calculate node weights
  sample.weights <- rep(1L,length(renyi.sprigs))
  sample.weights[which(sprigs.1 == sprigs.2)] <- 2L

  temp.sprig.sd <- c(tapply(sample.weights^2, renyi.sprigs, function(x){sqrt(sum(x))}))
  putative.sprig.sd <- rep(NA_real_,n.unique.sprigs)
  putative.sprig.sd[as.integer(names(temp.sprig.sd))] <- temp.sprig.sd

  temp.sprig.y <- c(tapply(y.resids * sample.weights, renyi.sprigs, sum))
  putative.sprig.y <- rep(NA_real_,n.unique.sprigs)
  putative.sprig.y[as.integer(names(temp.sprig.y))] <- temp.sprig.y
  putative.sprig.y <- putative.sprig.y / putative.sprig.sd


  if(!setequal(zero.weighted.sprigs, which(is.na(putative.sprig.y)))){
    stop("error: zero weighted sprigs did not line up with those missing putative sprig phenotypes")
  }

  message(paste("# of putative sprigs was ",length(putative.sprig.y)))
  message(paste("# of zero weighted sprigs: ",length(zero.weighted.sprigs)))


  w <- rep(NA_real_,n.unique.sprigs)
  if(length(zero.weighted.sprigs)){
    r <- renyi(putative.sprig.y[-zero.weighted.sprigs])
    message(paste("# of ranks in initial Renyi transform was ",length(r$ranks)))
    w[-zero.weighted.sprigs] <- rank2gauss(rank(r$signs * r$ranks))^4
  } else {
    r <- renyi(putative.sprig.y)
    message(paste("# of ranks in initial Renyi transform was ",length(r$ranks)))
    w<- rank2gauss(rank(r$signs * r$ranks))^4
  }

  old.y.sprigs <- putative.sprig.y
  old.y.sprigs[zero.weighted.sprigs] <- NA_real_
  #plot(putative.sprig.y,w)


  components <- find_irreducible_components(sprig.dag,keep.isolated.nodes = FALSE)
  inactive.sprigs <- active.sprigs <- replicate(length(components),integer(),simplify=FALSE)

  # Loop through components, solving for optimal assignments, making those assignments and recording nodes/sprigs left for the quadratic form to test
  for(i in 1:length(components)){
    if(length(components[[i]])==1){next}
    optimal.config <- solve_active_nodes(sprig.dag[components[[i]],components[[i]]], w[components[[i]]])
    active.sprigs[[i]] <- components[[i]][optimal.config]
    inactive.sprigs[[i]] <- components[[i]][!optimal.config]
  }

  active.sprigs <- unlist(active.sprigs)
  inactive.sprigs <- unlist(inactive.sprigs)
  non.renyi.sprigs <- c(inactive.sprigs,zero.weighted.sprigs)



  temp.assign.list <- new_clade_assignments(renyi.sprigs, sprigs.1, sprigs.2, inactive.sprigs = inactive.sprigs, overwrite.single.copy.carriers = FALSE)
  renyi.sprigs[temp.assign.list$assign.to.sprigs.1] <- sprigs.1[temp.assign.list$assign.to.sprigs.1]
  renyi.sprigs[temp.assign.list$assign.to.sprigs.2] <- sprigs.2[temp.assign.list$assign.to.sprigs.2]

  dual.samples.idx <- setdiff(dual.samples.idx,union(temp.assign.list$assign.to.sprigs.1,
                                                     temp.assign.list$assign.to.sprigs.2))

  dual.samples.idx <- dual.samples.idx[!(sprigs.1[dual.samples.idx] %in% inactive.sprigs &
                                           sprigs.2[dual.samples.idx] %in% inactive.sprigs)]

  if(length(dual.samples.idx)){stop("Some dual samples that should be assigned/resolved remain unassigned/unresolved.")}



  temp.assign.list <- new_clade_assignments(renyi.sprigs, sprigs.1, sprigs.2, inactive.sprigs = inactive.sprigs, overwrite.single.copy.carriers = TRUE)
  new.assign.renyi.sprigs <- rep(NA_integer_,length(renyi.sprigs))
  new.assign.renyi.sprigs[temp.assign.list$assign.to.sprigs.1] <- sprigs.1[temp.assign.list$assign.to.sprigs.1]
  new.assign.renyi.sprigs[temp.assign.list$assign.to.sprigs.2] <- sprigs.2[temp.assign.list$assign.to.sprigs.2]



  ###########################################
  ###########################################
  # Update sprig-level phenotypes
  # (removing sprigs that are now inactive
  # and adding newly assigned clade mamebers)
  ###########################################
  ###########################################


  # caution: we performed the Renyi on putative.sprig.y[-zero.weighted.sprigs] but all definitions of
  # inactive.sprigs were calculated in the original reference frame (same sprig labels as sprigs.1 etc.)
  # so we need to convert the indexing before inverting the Renyi



  updated.putative.sprig.y <- rep(NA_real_,n.unique.sprigs)
  temp.exclude.ind <- rep(FALSE,n.unique.sprigs)
  temp.exclude.ind[inactive.sprigs] <- TRUE

  if(length(zero.weighted.sprigs)){
    temp.exclude.ind <- temp.exclude.ind[-zero.weighted.sprigs]
    updated.putative.sprig.y[-zero.weighted.sprigs] <- inverse.renyi(r, exclude.idx = which(temp.exclude.ind))
  } else {
    updated.putative.sprig.y <- inverse.renyi(r, exclude.idx = which(temp.exclude.ind))
  }


  inactive.putative.sprig.y <- rep(NA_real_,n.unique.sprigs)
  temp.exclude.ind <- rep(FALSE,n.unique.sprigs)
  temp.exclude.ind[active.sprigs] <- TRUE

  if(length(zero.weighted.sprigs)){
    temp.exclude.ind <- temp.exclude.ind[-zero.weighted.sprigs]
    inactive.putative.sprig.y[-zero.weighted.sprigs] <- inverse.renyi(r, exclude.idx = which(temp.exclude.ind))
  } else {
    inactive.putative.sprig.y <- inverse.renyi(r, exclude.idx = which(temp.exclude.ind))
  }

  updated.putative.sprig.y <- pmax(updated.putative.sprig.y,inactive.putative.sprig.y,na.rm = TRUE)



  if(any(is.na(updated.putative.sprig.y[-c(zero.weighted.sprigs)]))){
    stop("There is at least one NA in updated.putative.sprig.y that does not correspond to a sprig listed in zero.weighted.sprigs!")
  }

  # show that the Renyi procedure made on subtle changes to putative.sprig.y
  # plot(putative.sprig.y,updated.putative.sprig.y); abline(0,1)

  new.putative.sprig.sd <- rep(NA_real_,n.unique.sprigs)
  new.temp.sprig.sd <- c(tapply(sample.weights^2, new.assign.renyi.sprigs, function(x){sqrt(sum(x))}))
  new.putative.sprig.sd[as.integer(names(new.temp.sprig.sd))] <- new.temp.sprig.sd

  new.putative.sprig.y <- rep(NA_real_,n.unique.sprigs)
  new.temp.sprig.y <- c(tapply(y.resids * sample.weights, new.assign.renyi.sprigs, sum))
  new.putative.sprig.y[as.integer(names(new.temp.sprig.y))] <- new.temp.sprig.y
  new.putative.sprig.y <- new.putative.sprig.y / new.putative.sprig.sd

  if(!setequal(which(!is.na(new.putative.sprig.y)),na.omit(unique(new.assign.renyi.sprigs))) ){
    stop("The sprigs for which newly-assigned sprig level phenotypes were calculated do not match the sprigs present in new.assign.renyi.sprigs.")
  }


  new.putative.sprig.sd[is.na(new.putative.sprig.sd)] <- 0
  new.putative.sprig.y[is.na(new.putative.sprig.y)] <- 0


  y.sprigs <- (putative.sprig.sd * updated.putative.sprig.y + new.putative.sprig.sd * new.putative.sprig.y) /
    sqrt(putative.sprig.sd^2 + new.putative.sprig.sd^2)

  if(!setequal(which(is.na(y.sprigs)), c(zero.weighted.sprigs))){
    stop("There is at least one NA in sprig.y that does not correspond to a sprig listed in zero.weighted.sprigs!")
  }



  sprigs.to.remove.sizes <- tabulate(sprigs, nbins = n.unique.sprigs)
  sprigs.to.remove.sizes[non.renyi.sprigs] <- NA_integer_


  sprig.sd <- rep(NA_real_,n.unique.sprigs)
  temp.sprig.sd <- c(tapply(sample.weights^2, renyi.sprigs, function(x){sqrt(sum(x))}))
  sprig.sd[as.integer(names(temp.sprig.sd))] <- temp.sprig.sd

  sprig.weight <- rep(NA_real_,n.unique.sprigs)
  temp.sprig.weight <- c(tapply(sample.weights, renyi.sprigs, sum))
  sprig.weight[as.integer(names(temp.sprig.weight))] <- temp.sprig.weight

  sprig.coefficient <- sprig.sd / sprig.weight



  return(list("skip.renyi" = FALSE,
              "y.sprigs" = y.sprigs,
              "renyi.sprigs" = renyi.sprigs,
              "old.y.sprigs" = old.y.sprigs,
              "non.renyi.sprigs" = non.renyi.sprigs,
              "sprigs.to.remove.sizes" = sprigs.to.remove.sizes,
              "sprig.coefficient" = sprig.coefficient))

  # demo:
  #plot(putative.sprig.y,y.sprigs); abline(0,1)
  #qqnorm(y.sprigs);abline(0,1)

}

make_factor <- function(x,select){
  outlier.clades <- gl(1,length(x),labels=NA_character_)
  outlier.clades.set <- which(x %in% select)
  sprigs.outlier.clades.set <- x[outlier.clades.set]
  levels(outlier.clades) <- c(NA_character_,unique(sprigs.outlier.clades.set))
  outlier.clades[outlier.clades.set] <- sprigs.outlier.clades.set
  outlier.clades
}



rotate <- function(x, w = rep(1,length(x))){
  k <- length(x)
  if(k!=length(w)){stop("length of x and w must be equal.")}
  w_k <- w[k]
  C <- sqrt(sum(w^2))
  beta <- C / abs(w_k)
  alpha <- (beta * (beta+1) * (w_k / w)^2)[-k]
  # the last element of alpha, alpha_k, here is not necessary & is dropped but kept for making
  # the expression below simple
  y_1 <- sum(x * w)
  c(y_1, y_1 - x[-k] * alpha * w[-k] + x[length(x)] * beta * w_k) / c(C,abs(w[-k]) * alpha)
}


inverse.rotate <- function(x, w = rep(1,length(x))){
  k <- length(x)
  if(k!=length(w)){stop("length of x and w must be equal.")}
  w_k <- w[k]
  C <- sqrt(sum(w^2))
  beta <- C / abs(w_k)
  alpha <- (beta * (beta+1) * (w_k / w)^2)[-k]
  # the last element of alpha, alpha_k, here is not necessary & is dropped but kept for making
  # the expression below simple

  x <- x/(c(C,abs(w[-k])*alpha))
  x.sum <- sum(x)
  x.sum * w + c(-x[-1] * w[-k] * alpha, w_k * beta * (x.sum - x[1]))
}
#
# x <- 1:10
# w <- rep(1,10)
# z <- rotate(x,w)
# inverse.rotate(z,w)
#
# z2 <- z
# z2[1] <- z2[1] + 1
# inverse.rotate(z2,w)
#
#
# shapiro.test(rotate(rnorm(1e3)))
#
#
# w <- rnorm(1e4)
# x <- rnorm(1e4)
# plot(x,inverse.rotate(rotate(x,w),w))
# abline(0,1)
#
# x <- rnorm(1e6)
# w <- rnorm(1e6)
#
# start <- proc.time()
# x <- inverse.rotate(rotate(x,w),w)
# proc.time() - start
# #

low_mem_STAAR <- function (genotype, obj_nullmodel, annotation_phred = NULL, rare_maf_cutoff = 0.01,
                           rv_num_cutoff = 2,
                           w_1_beta_shape_params = c(1,25),
                           w_2_beta_shape_params = c(1,1))
{
  if (!inherits(genotype, "matrix") && !inherits(genotype,
                                                 "Matrix")) {
    stop("genotype is not a matrix!")
  }
  if (dim(genotype)[2] == 1) {
    stop(paste0("Number of rare variant in the set is less than 2!"))
  }
  annotation_phred <- as.data.frame(annotation_phred)
  if (dim(annotation_phred)[1] != 0 & dim(genotype)[2] != dim(annotation_phred)[1]) {
    stop(paste0("Dimensions don't match for genotype and annotation!"))
  }

  ######################################
  ######################################
  # Begin RC edits (edit block 1 of 2)
  # here we make an additional assumption which STAAR does not make:
  # we assume that all variants have already been imputed so all entries
  # of genotype are in {0,1,2}

  csG <- Matrix::colSums(genotype)
  AF <- csG / (2*nrow(genotype))

  # recode genotypes so that they all encode the minor allele
  if(any(AF > 0.5)){
    genotype[, AF > 0.5, drop=FALSE] <- 2 - genotype[, AF > 0.5, drop=FALSE]
  }

  MAF <- pmin(AF,1-AF)
  RV_label <- as.vector((MAF < rare_maf_cutoff) & (MAF > 0))

  if(!all(RV_label)){
    Geno_rare <- genotype[, RV_label]
  } else {
    Geno_rare <- genotype
  }
  # only other edit besides this code chunk is that
  # we added the STAAR::: prefix to all internal
  # STAAR package function calls below
  # End RC edits (edit block 1 of 2)
  ######################################
  ######################################
  rm(genotype)
  gc()

  annotation_phred <- annotation_phred[RV_label, , drop = FALSE]
  if (sum(RV_label) >= rv_num_cutoff) {
    G <- as(Geno_rare, "dgCMatrix")
    MAF <- MAF[RV_label]
    rm(Geno_rare)
    gc()
    annotation_rank <- 1 - 10^(-annotation_phred/10)

    ######################################
    ######################################
    # Begin RC edits (edit block 2 of 2)
    # Here we allow for manual overriding of beta shape parameters
    #' *WARNING: with this override, the output is still labeled 1,25 and 1,1*
    #' *so w_1_beta_shape_params results will be called 1,25 and*
    #' *w_2_beta_shape_params results will be called 1,1*
    w_1 <- dbeta(MAF, w_1_beta_shape_params[1], w_1_beta_shape_params[2]) # 1,25
    w_2 <- dbeta(MAF, w_2_beta_shape_params[1], w_2_beta_shape_params[2]) # 1,1
    # End RC edits (edit block 2 of 2)
    ######################################
    ######################################

    if (dim(annotation_phred)[2] == 0) {
      w_B <- w_S <- as.matrix(cbind(w_1, w_2))
      w_A <- as.matrix(cbind(w_1^2/dbeta(MAF, 0.5, 0.5)^2,
                             w_2^2/dbeta(MAF, 0.5, 0.5)^2))
    }
    else {
      w_B_1 <- annotation_rank * w_1
      w_B_1 <- cbind(w_1, w_B_1)
      w_B_2 <- annotation_rank * w_2
      w_B_2 <- cbind(w_2, w_B_2)
      w_B <- cbind(w_B_1, w_B_2)
      w_B <- as.matrix(w_B)
      w_S_1 <- sqrt(annotation_rank) * w_1
      w_S_1 <- cbind(w_1, w_S_1)
      w_S_2 <- sqrt(annotation_rank) * w_2
      w_S_2 <- cbind(w_2, w_S_2)
      w_S <- cbind(w_S_1, w_S_2)
      w_S <- as.matrix(w_S)
      w_A_1 <- annotation_rank * w_1^2/dbeta(MAF, 0.5,
                                             0.5)^2
      w_A_1 <- cbind(w_1^2/dbeta(MAF, 0.5, 0.5)^2, w_A_1)
      w_A_2 <- annotation_rank * w_2^2/dbeta(MAF, 0.5,
                                             0.5)^2
      w_A_2 <- cbind(w_2^2/dbeta(MAF, 0.5, 0.5)^2, w_A_2)
      w_A <- cbind(w_A_1, w_A_2)
      w_A <- as.matrix(w_A)
    }
    if (obj_nullmodel$relatedness) {
      if (!obj_nullmodel$sparse_kins) {
        P <- obj_nullmodel$P
        residuals.phenotype <- obj_nullmodel$scaled.residuals
        pvalues <- STAAR:::STAAR_O_SMMAT(G, P, residuals.phenotype,
                                         weights_B = w_B, weights_S = w_S, weights_A = w_A,
                                         mac = as.integer(round(MAF * 2 * dim(G)[1])))
      }
      else {
        Sigma_i <- obj_nullmodel$Sigma_i
        Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
        cov <- obj_nullmodel$cov
        residuals.phenotype <- obj_nullmodel$scaled.residuals
        pvalues <- STAAR:::STAAR_O_SMMAT_sparse(G, Sigma_i, Sigma_iX,
                                                cov, residuals.phenotype, weights_B = w_B,
                                                weights_S = w_S, weights_A = w_A, mac = as.integer(round(MAF *
                                                                                                           2 * dim(G)[1])))
      }
    }
    else {
      X <- model.matrix(obj_nullmodel)
      working <- obj_nullmodel$weights
      sigma <- sqrt(summary(obj_nullmodel)$dispersion)
      if (obj_nullmodel$family[1] == "binomial") {
        fam <- 1
      }
      else if (obj_nullmodel$family[1] == "gaussian") {
        fam <- 0
      }
      residuals.phenotype <- obj_nullmodel$y - obj_nullmodel$fitted.values
      pvalues <- STAAR:::STAAR_O(G, X, working, sigma, fam, residuals.phenotype,
                                 weights_B = w_B, weights_S = w_S, weights_A = w_A,
                                 mac = as.integer(round(MAF * 2 * dim(G)[1])))
    }
    num_variant <- sum(RV_label)
    cMAC <- sum(G)
    num_annotation <- dim(annotation_phred)[2] + 1
    results_STAAR_O <- STAAR:::CCT(pvalues)
    results_ACAT_O <- STAAR:::CCT(pvalues[c(1, num_annotation + 1,
                                            2 * num_annotation + 1, 3 * num_annotation + 1, 4 *
                                              num_annotation + 1, 5 * num_annotation + 1)])
    pvalues_STAAR_S_1_25 <- STAAR:::CCT(pvalues[1:num_annotation])
    pvalues_STAAR_S_1_1 <- STAAR:::CCT(pvalues[(num_annotation +
                                                  1):(2 * num_annotation)])
    pvalues_STAAR_B_1_25 <- STAAR:::CCT(pvalues[(2 * num_annotation +
                                                   1):(3 * num_annotation)])
    pvalues_STAAR_B_1_1 <- STAAR:::CCT(pvalues[(3 * num_annotation +
                                                  1):(4 * num_annotation)])
    pvalues_STAAR_A_1_25 <- STAAR:::CCT(pvalues[(4 * num_annotation +
                                                   1):(5 * num_annotation)])
    pvalues_STAAR_A_1_1 <- STAAR:::CCT(pvalues[(5 * num_annotation +
                                                  1):(6 * num_annotation)])
    results_STAAR_S_1_25 <- c(pvalues[1:num_annotation],
                              pvalues_STAAR_S_1_25)
    results_STAAR_S_1_25 <- data.frame(t(results_STAAR_S_1_25))
    results_STAAR_S_1_1 <- c(pvalues[(num_annotation + 1):(2 *
                                                             num_annotation)], pvalues_STAAR_S_1_1)
    results_STAAR_S_1_1 <- data.frame(t(results_STAAR_S_1_1))
    results_STAAR_B_1_25 <- c(pvalues[(2 * num_annotation +
                                         1):(3 * num_annotation)], pvalues_STAAR_B_1_25)
    results_STAAR_B_1_25 <- data.frame(t(results_STAAR_B_1_25))
    results_STAAR_B_1_1 <- c(pvalues[(3 * num_annotation +
                                        1):(4 * num_annotation)], pvalues_STAAR_B_1_1)
    results_STAAR_B_1_1 <- data.frame(t(results_STAAR_B_1_1))
    results_STAAR_A_1_25 <- c(pvalues[(4 * num_annotation +
                                         1):(5 * num_annotation)], pvalues_STAAR_A_1_25)
    results_STAAR_A_1_25 <- data.frame(t(results_STAAR_A_1_25))
    results_STAAR_A_1_1 <- c(pvalues[(5 * num_annotation +
                                        1):(6 * num_annotation)], pvalues_STAAR_A_1_1)
    results_STAAR_A_1_1 <- data.frame(t(results_STAAR_A_1_1))
    if (dim(annotation_phred)[2] == 0) {
      colnames(results_STAAR_S_1_25) <- c("SKAT(1,25)",
                                          "STAAR-S(1,25)")
      colnames(results_STAAR_S_1_1) <- c("SKAT(1,1)", "STAAR-S(1,1)")
      colnames(results_STAAR_B_1_25) <- c("Burden(1,25)",
                                          "STAAR-B(1,25)")
      colnames(results_STAAR_B_1_1) <- c("Burden(1,1)",
                                         "STAAR-B(1,1)")
      colnames(results_STAAR_A_1_25) <- c("ACAT-V(1,25)",
                                          "STAAR-A(1,25)")
      colnames(results_STAAR_A_1_1) <- c("ACAT-V(1,1)",
                                         "STAAR-A(1,1)")
    }
    else {
      colnames(results_STAAR_S_1_25) <- c("SKAT(1,25)",
                                          paste0("SKAT(1,25)-", colnames(annotation_phred)),
                                          "STAAR-S(1,25)")
      colnames(results_STAAR_S_1_1) <- c("SKAT(1,1)", paste0("SKAT(1,1)-",
                                                             colnames(annotation_phred)), "STAAR-S(1,1)")
      colnames(results_STAAR_B_1_25) <- c("Burden(1,25)",
                                          paste0("Burden(1,25)-", colnames(annotation_phred)),
                                          "STAAR-B(1,25)")
      colnames(results_STAAR_B_1_1) <- c("Burden(1,1)",
                                         paste0("Burden(1,1)-", colnames(annotation_phred)),
                                         "STAAR-B(1,1)")
      colnames(results_STAAR_A_1_25) <- c("ACAT-V(1,25)",
                                          paste0("ACAT-V(1,25)-", colnames(annotation_phred)),
                                          "STAAR-A(1,25)")
      colnames(results_STAAR_A_1_1) <- c("ACAT-V(1,1)",
                                         paste0("ACAT-V(1,1)-", colnames(annotation_phred)),
                                         "STAAR-A(1,1)")
    }
    return(list(num_variant = num_variant, cMAC = cMAC, RV_label = RV_label,
                results_STAAR_O = results_STAAR_O, results_ACAT_O = results_ACAT_O,
                results_STAAR_S_1_25 = results_STAAR_S_1_25, results_STAAR_S_1_1 = results_STAAR_S_1_1,
                results_STAAR_B_1_25 = results_STAAR_B_1_25, results_STAAR_B_1_1 = results_STAAR_B_1_1,
                results_STAAR_A_1_25 = results_STAAR_A_1_25, results_STAAR_A_1_1 = results_STAAR_A_1_1))
  }
  else {
    stop(paste0("Number of rare variant in the set is less than ",
                rv_num_cutoff, "!"))
  }
}






