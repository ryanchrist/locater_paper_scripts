
# Returns a function that for a given observation returns a vector of three numbers
# first entry: -log10 lower bound 
# second entry: -log10 upper bound
# third entry: -log10 pvalue point estimate





CalcBoundsAttemptHigherOrderApprox <- function(traces, e = NULL, delta2 = 0, higher.order = FALSE, num.bg.covs = NULL){
  
  # Returns a function that provides point estimate approximation in the case where we divide YPMPY by YPIPY to get 
  # the scale of a quantitative trait
  
  # traces must be the traces of PMP even when using higher order approximation
  
  if(higher.order){
    if(is.null(num.bg.covs){stop("num.bg.covs required if higher.order = TRUE")})
    mat.trace <- 0
    mat.hsnorm2 <- traces$hsnorm2 - ((traces$trace)^2)/(length(traces$diag) - num.bg.covs) 
  } else {
    mat.trace <- traces$trace
    mat.hsnorm2 <- traces$hsnorm^2
  }
  
  
  if(is.null(e)){ # Use fully Gaussian approximation

    mu <- traces$trace  
    sigma <- sqrt((traces$hsnorm2 - if(higher.order){(traces$trace^2) / (length(traces$diag) - num.bg.covs) }else{0} ) * (2+4*delta2) / (1+delta2)^2)
    
    cdf.approx <- function(x, lower.tail = TRUE, log.p = FALSE){
      pnorm(x, mean = mu, sd = sigma,
            lower.tail = lower.tail, log.p = log.p)}
    attr(cdf.approx,"mu") <- mu
    attr(cdf.approx,"Q.sd") <- sigma
    
  } else {
    
    if(!higher.order){
      
    E.R <- (1+delta2) * (traces$trace - sum(e$values))

    temp.cdf <- QFGauss(f.eta = e$values,
                        delta = rep(sqrt(delta2),length(e$values)),
                        sigma = sqrt( (2+4*delta2) * (traces$hsnorm2 - sum(e$values^2)) ))
    cdf.approx <- function(x, lower.tail = TRUE, log.p = FALSE){ temp.cdf(x * (1+delta2) - E.R, lower.tail = lower.tail, log.p = log.p) }
    attr(cdf.approx,"mu") <- traces$trace  
    attr(cdf.approx,"Q.sd") <- sqrt((traces$hsnorm2) * (2+4*delta2) / (1+delta2)^2)
    
    } else {
      
      E.R <- (1+delta2) * (0 - sum(e$values))
      
      temp.cdf <- QFGauss(f.eta = e$values,
                          delta = rep(sqrt(delta2),length(e$values)),
                          sigma = sqrt( (2+4*delta2) * (traces$hsnorm2 - (traces$trace^2) / (length(traces$diag) - num.bg.covs) - sum(e$values^2)) ))
      cdf.approx <- function(x, lower.tail = TRUE, log.p = FALSE){ temp.cdf( (x - traces$trace) * (1+delta2) - E.R, lower.tail = lower.tail, log.p = log.p) }
      attr(cdf.approx,"mu") <- traces$trace
      attr(cdf.approx,"Q.sd") <- sqrt((traces$hsnorm2 - (traces$trace^2) / (length(traces$diag) - num.bg.covs) ) * (2+4*delta2) / (1+delta2)^2)
      
      
    }
  }

}
  
  


# Calculate standard CDF approximations

EvalCDFApprox <- function(obs,
                          samps, # vector of true samples sorted from lowest to highest
                          traces, # traces as returned by calc traces
                          ekurt, # a scalar or vector of excess kurtosis (for the phenotype residuals)
                          e = NULL, # an eigendecomposition returned by eg: eigen or eigs 
                          higher.order = FALSE # if TRUE use higher order approximation to quotient of quadratic forms
){
  
  if(is.null(e)){ # Use fully Gaussian approximation
    per.true.var.in.trunc.part <- 0 
    
    cdf.approx <- function(x, lower.tail = TRUE, log.p = FALSE){
      pnorm(x, mean = traces$trace, sd = sqrt(2*traces$hsnorm2),
            lower.tail = lower.tail, log.p = log.p)}
    attr(cdf.approx,"mu") <- traces$trace
    attr(cdf.approx,"Q.sd") <- sqrt(2*traces$hsnorm2)
    
  } else {
    
    true.var <- sum(ekurt * traces$diag^2) + 2*traces$hsnorm2
    per.true.var.in.trunc.part <- sum(e$values^2) / true.var
    
    E.R <- traces$trace - sum(e$values)
    
    if(include.ekurt){ 
      diag.part.of.var <- sum(ekurt * (traces$diag - (e$vectors^2) %*% e$values)^2 ) 
    } else {
      diag.part.of.var <- 0
    }
    
    temp.cdf <- QFGauss(f.eta = e$values,sigma = sqrt(diag.part.of.var + 2 * (traces$hsnorm2 - sum(e$values^2)) ))
    cdf.approx <- function(x, lower.tail = TRUE, log.p = FALSE){ temp.cdf(x - E.R, lower.tail = lower.tail, log.p = log.p) }
    attr(cdf.approx,"mu") <- attr(temp.cdf,"mu") + E.R
    attr(cdf.approx,"Q.sd") <- attr(temp.cdf,"Q.sd")
  }
  
  
  test <- goftest::ad.test(samps, null = cdf.approx, estimated=FALSE)
  An <- as.numeric(test$statistic) # as.numeric to remove name
  AD.neglog10pval <- -log10(test$p.value)
  
  desired.quantiles <- qbeta(ppoints(50),0.4,0.4)
  yy <- ceiling(length(samps)*desired.quantiles) # desired ranks
  xx <- samps[yy]
  approx.probs <- cdf.approx(xx)
  
  log.lower.tail <- pbinom(yy,length(samps),approx.probs,TRUE,TRUE)
  
  log.upper.tail <- dbinom(yy,length(samps),approx.probs,TRUE)
  log.upper.tail <- log1p(exp(pbinom(yy,length(samps),approx.probs,FALSE,log.p = TRUE) 
                              - log.upper.tail)) + log.upper.tail 
  
  nl.pval <- -pmin(log.lower.tail,log.upper.tail)/log(10)
  
  max.departure.neglog10pval <- max(nl.pval)
  
  which.max.nl.pval <- which(nl.pval==max.departure.neglog10pval)-25
  which.max.nl.pval <- which.min(abs(which.max.nl.pval))
  
  quantile.of.max.departure <- desired.quantiles[which.max.nl.pval]
  
  ecdf.neglog10pval <- -log10(mean(samps<=obs))
  approx.neglog10pval <- -cdf.approx(obs,log.p=TRUE)/log(10)
  
  list(
    "per.true.var.in.trunc.part" = per.true.var.in.trunc.part,
    "An" = An, 
    "AD.neglog10pval" = AD.neglog10pval,
    "quantile.of.max.departure" = quantile.of.max.departure,
    "max.departure.neglog10pval" = max.departure.neglog10pval,
    "obs" = obs,
    "ecdf.neglog10pval" = ecdf.neglog10pval,
    "approx.neglog10pval" = approx.neglog10pval
  )
}




EvalCDFBinary <- function(obs,
                          samps, # vector of true samples
                          traces, # traces calculated on hollow matrix
                          p, # probability of equalling 1
                          hollow.e = NULL # an eigendecomposition of hollow matrix returned by eg: eigen or eigs 
){
  
  mu <- 2 * p - 1
  sigma <- 2 * sqrt(p*(1-p))
  ekurt <- 1 / (p * (1-p)) - 6
  
  true.var <- sum(ekurt * traces$diag^2) + 2*traces$hsnorm2
  
  if(is.null(hollow.e)){ # Use fully Gaussian approximation
    per.true.var.in.trunc.part <- 0 
    
    cdf.approx <- function(x, lower.tail = TRUE, log.p = FALSE){
      pnorm(x, mean = traces$trace, sd = sqrt(true.var),
            lower.tail = lower.tail, log.p = log.p)}
    attr(cdf.approx,"mu") <- traces$trace
    attr(cdf.approx,"Q.sd") <- sqrt(true.var)
    
  } else {
    
    epsilon <- c(crossprod(hollow.e$vectors,traces$diag*mu/sigma))
    
    trunc.var <- 2*sum(hollow.e$values^2) + 4*sum(epsilon^2)
    
    per.true.var.in.trunc.part <- trunc.var / true.var
    
    raw.cdf <- QFGauss(f.eta = hollow.e$values,
                       delta = epsilon/hollow.e$values,
                       sigma = sqrt(true.var - trunc.var))
    
    mu0 <- traces$trace - sum(hollow.e$values + (epsilon^2)/hollow.e$values)
    
    cdf.approx <- function(x, density = FALSE, lower.tail = TRUE, log.p = FALSE){
      raw.cdf(x - mu0,density,lower.tail,log.p)
    }
    
  }
  
  test <- goftest::ad.test(samps, null = cdf.approx, estimated=FALSE)
  An <- as.numeric(test$statistic) # as.numeric to remove name
  AD.neglog10pval <- -log10(test$p.value)
  
  
  desired.quantiles <- qbeta(ppoints(50),0.4,0.4)
  yy <- ceiling(length(samps)*desired.quantiles) # desired ranks
  xx <- samps[yy]
  approx.probs <- cdf.approx(xx)
  
  log.lower.tail <- pbinom(yy,length(samps),approx.probs,TRUE,TRUE)
  
  log.upper.tail <- dbinom(yy,length(samps),approx.probs,TRUE)
  log.upper.tail <- log1p(exp(pbinom(yy,length(samps),approx.probs,FALSE,log.p = TRUE) 
                              - log.upper.tail)) + log.upper.tail 
  
  nl.pval <- -pmin(log.lower.tail,log.upper.tail)/log(10)
  
  max.departure.neglog10pval <- max(nl.pval)
  
  which.max.nl.pval <- which(nl.pval==max.departure.neglog10pval)-25
  which.max.nl.pval <- which.min(abs(which.max.nl.pval))
  
  quantile.of.max.departure <- desired.quantiles[which.max.nl.pval]
  
  ecdf.neglog10pval <- -log10(mean(samps<=obs))
  approx.neglog10pval <- -cdf.approx(obs,log.p=TRUE)/log(10)
  
  list(
    "per.true.var.in.trunc.part" = per.true.var.in.trunc.part,
    "An" = An, 
    "AD.neglog10pval" = AD.neglog10pval,
    "quantile.of.max.departure" = quantile.of.max.departure,
    "max.departure.neglog10pval" = max.departure.neglog10pval,
    "obs" = obs,
    "ecdf.neglog10pval" = ecdf.neglog10pval,
    "approx.neglog10pval" = approx.neglog10pval
  )
}




LinearModelResidBootQForm<- function(n.empirical.samps,# number of samples requested
                                     y,
                                     y.hat,
                                     M,
                                     Q,
                                     Z){
  
  # Draw samples from quadratic form true distribution in batches of 1000 and sort the resulting samples
  # from lowest to highest 
  
  epsilon <- y - y.hat
  
  n <- nrow(M)
  
  p <- ncol(Q) # number of parameters in the null model
  
  n.sets <- n.empirical.samps%/%1e3
  if(n.sets > 0){
    samps <- as.list(1:n.sets)
    for(ii in 1:n.sets){
      samp.resids <- y.hat + matrix(sample(epsilon, n*1e3,replace = TRUE),n)
      samp.resids <- samp.resids - Q %*% (Z %*% samp.resids)
      samps[[ii]] <- (n-p) * c(colSums(samp.resids * (M %*% samp.resids))) / c(colSums(samp.resids^2))
    }
  }
  
  # Add on any remainder of samples here.
  n.remainder <- n.empirical.samps%%1e3
  
  if(n.remainder==0){
    return(sort(c(simplify2array(samps))))
  } else {
    samp.resids <- y.hat + matrix(sample(epsilon, n*n.remainder,replace = TRUE),n)
    samp.resids <- samp.resids - Q %*% (Z %*% samp.resids)
    temp <- (n-p) * c(colSums(samp.resids * (M %*% samp.resids))) / c(colSums(samp.resids^2))
    return(sort(c(c(simplify2array(samps)),temp)))
  }
}




HaploidSKAT <- function(y.resids,Q,G,maf,SDy.for.logistic = NULL){
  
  if(!is.null(SDy.for.logistic)){
    G <- G * SDy.for.logistic
  }
  
  # G is a matrix of 0s and 1s that is N x p for p variants
  
  G.tilde <- (G - Q %*%crossprod(Q,G) ) * dbeta(maf,1,25)
  
  eta <- (svd(G.tilde,0,0)$d)^2
  
  obs <- sum(c(crossprod(G.tilde,y.resids))^2)
  
  -QForm::QFGauss(eta)(obs,lower.tail = FALSE,log.p=T)/log(10)
}






SampleCaseControlSKATCalcTestStat <- function(y,G,A){
  
  # Calculate maf for observed G (do they calculate maf based on cases and controls in their method or just based on controls?)
  daf <- colMeans(G)/2   # calculate this before generating haplotype level matrix
  maf <- pmin(daf,1-daf)
  
  # QR decompose A to get Q
  # a <- qr(A)
  # Q <- qr.Q(a) # Q is n x p
  # rm(a)
  # 
  
  m1 <- glm(y ~ 0 + A, family = "binomial")
  y.resids <- y - m1$fitted.values
  
  G.tilde <- sweep(G,2,dbeta(maf,1,25),"*")
  
  sum(c(crossprod(G.tilde,y.resids))^2)
}




SampleCaseControlSKATQForm <- function(x,
                                       y_pop,
                                       G_pop,
                                       A_pop,
                                       n.cases,
                                       n.controls){
  
  # Sample cases and controls from y_pop 
  cases.idx <- which(y_pop == 1)
  controls.idx <- which(y_pop == 0)
  
  # enable super sampling
  subset.case <- sample(cases.idx, size = n.cases, replace = T)
  subset.control <- sample(controls.idx, size = n.controls, replace = T)
  subset <- c(subset.case,subset.control)
  
  # Subset y_pop, G_pop and A_pop accordingly
  y <- y_pop[subset]
  G  <- G_pop[subset,]
  A <- A_pop[subset,]
  
  
  SampleCaseControlSKATCalcTestStat(y,G,A)
}




# 
# 
# obs <- SampleCaseControlSKATCalcTestStat(y,G,A)
# 
# require(parallel)
# nthreads <- 4
# samps <- as.list(1:n.samps)
# samps <- unlist(mclapply(samps, SampleCaseControlSKATQForm,
#          y_pop = y_pop, G_pop = G_pop, A_pop = A_pop, n.cases = n.cases, n.controls = n.controls,
#          mc.cores = nthreads))
# 
# p.value <- mean(samps > obs)
# 




#test <- SampleBinaryQForm(5e3+500,rbeta(1e3,1,1),matrix(rchisq(1e6,1),1e3),rep(1,1e3),rep(1,1e3))
#test <- BootstrapQForm(5e3+500,rnorm(1e3),matrix(rchisq(1e6,1),1e3),rep(1,1e3),rep(1,1e3))
