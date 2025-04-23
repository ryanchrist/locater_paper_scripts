splot <- function(x, y, z,z.name){
  # set axis label
  axx <- list(title = "neglog10Ne")
  axy <- list(title = "neglog10mu")
  axz <- list(title = z.name)
  
  A <- akima::interp(x= x,y = y,z = z,
                     xo = seq(min(x), max(x), length = 100),
                     yo = seq(min(y), max(y), length = 100),
                     linear = TRUE, extrap = TRUE)
  A$z <- t(A$z)
  names(A) <- c("int","slope","neglog10pval")
  
  p <-  plotly::plot_ly() %>%
    plotly::add_surface(x= A$int , y= A$slope, z= A$neglog10pval) %>%
    plotly::add_markers(x= x , y= y, z= z) %>%
    plotly::layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz, # set axis lable
                                # set the aspectmode to fix the ratio of 3 axes 
                                aspectmode='cube'))
  
  p
}



# project.matrix <- function(M,Q,Z,sigma = NULL){
#   X <- M%*%Q
#   temp <- X%*%Z
#   M <- M - temp - t(temp) + (Q%*%(Z%*%X))%*%Z
#   rm(X,temp); gc()
#   if(is.null(sigma)){
#     return(M)
#   }else{
#     return(sweep(M*sigma,2,sigma,"*"))
#   }
# }

MakeTuneMatMul <- function(M,Q,Z,sigma = NULL){
  if(is.null(sigma)){
    out <- function(x, args = NULL) {
      x <- x - Q %*% (Z %*% x)
      x <- M%*%x
      x - Q %*% (Z %*% x)
    }
  }else{
    out <- function(x, args = NULL) {
      x <- sigma * x
      x <- x - Q %*% (Z %*% x)
      x <- M%*%x
      c(x - Q %*% (Z %*% x)) * sigma
    }
  }
}



bckfuturefunc <- function(mat, y, Ey, current.target, k, pre.std = TRUE){

  allele.f <- mean(y)
  sigma <- sqrt(Ey*(1-Ey))
  y.residuals <- (y - Ey)/sigma


  if(pre.std){
    diag(mat) <- NA
    mat <- scale(mat)
  }else{
    mat <- sweep(mat,MARGIN = 2,STATS=log(colSums(exp(-mat))),FUN="+")
  }

  diag(mat) <- 0
  mat <- (mat + t(mat))/2


  # Initialize output vector
  out <- rep(0,9)


  out[1] <- TruncApproxSignificance(sweep(x=mat*sigma,MARGIN = 2,STATS = sigma,FUN = "*"),y.residuals, k)

  c(current.target, allele.f, out) # returns a vector with 11 elements
}
# Make the downstream changes




CalcApproxPvalue <- function(evalues,
                             obs,
                             trace1,
                             upper_bound_on_trace2,
                             lower.tail = TRUE){

  if(any(!is.numeric(evalues)) || any(!is.finite(evalues))){return(NA_real_)}

  R.sum.etasq <- upper_bound_on_trace2 - sum(evalues^2)
  E.R <-  trace1 - sum(evalues)

  gauss.tcdf <- QFGauss(evalues,sigma = sqrt(2*R.sum.etasq))

  gauss.approxfullcdf <- function(x, density = FALSE, lower.tail = TRUE, log.p = FALSE) gauss.tcdf(x-E.R, density = density, lower.tail = lower.tail, log.p = log.p)
  attr(gauss.approxfullcdf,"mu") <- attr(gauss.tcdf,"mu") + E.R
  attr(gauss.approxfullcdf,"Q.sd") <- attr(gauss.tcdf,"Q.sd")

  -gauss.approxfullcdf(obs,lower.tail = lower.tail, log.p=T)/log(10)
}


VariantResolution <- function(y, fwd, bck, pars, standardize = FALSE, symmetrize = TRUE, k = 100,
                              M = NULL, A = NULL, robust = FALSE, from_recipient = 1, nthreads = 1){

  # let y.tilde = (y - mu )/ sqrt(mu * (1-mu)), the standardized variant at a core locus l.
  # let fwd$l = l - 1 and  bck$l = l + 1.
  # For its test statistic, this function uses y.tilde' P D M D P y where
  # M is our core distance / probability matrix coming from kalis
  # D = Diag( sqrt(mu * (1-mu) ) )
  # P = I - QQ' where Q is from the QR decomposition of D %*% A
  # A is the background PCs / covariates


  if(is.null(A)){
    A <- matrix(1,length(y),1)
    scale <- FALSE
  } else {
    scale <- TRUE
  }

  if(is.null(M)){M <- matrix(0,length(y),length(y))}

  m1 <- glm(y ~ 0 + A, family="binomial")
  Ey <- m1$fitted.values
  SDy <- sqrt(m1$fitted.values * (1-m1$fitted.values))
  y.resids <- (y - Ey) / SDy
  if(scale){ A <- A * SDy }

  # QR decomposition
  a <- qr(A); Q <- qr.Q(a); rm(a)
  tQ <- t(Q)


  z <- c(y.resids - Q %*% (tQ %*% y.resids))
  if(scale){ z <- z * SDy}

  rho.list <- kalis:::input_checks_for_probs_and_dist_mat(fwd,bck,beta.theta.opts =
                                                            list("pars" = pars, "bias" = 0.5))

  if(robust){
    M <- sweep((fwd$alpha),2,fwd$alpha.f,"/")*sweep((bck$beta),2,bck$beta.g*(length(bck$beta.g)-1),"/")
    M[,which(fwd$alpha.f==0 | bck$beta.g==0)] <- 0
    M <- -log(M)
    M[!is.finite(M)] <- 744.4400719213812180897
    M <- 0.5 * (M + t(M))
    diag(M) <- 0
    obs <- sum(z* (M%*%z))

  } else {
    obs <- sum(z * kalis:::MatAndMulBtwVar(M,fwd,bck,z,standardize, FALSE, FALSE, rho.list$rho.fwd, rho.list$rho.bck, from_recipient,nthreads))

    if(symmetrize){
      # explicitly symmetrize M to get the HS norm exactly, otherwise, the following just yields an upper bound
      M <- 0.5 * (M + t(M))
    }

  }

  if(scale){
    M <- M * SDy
    M <- sweep(M,2,SDy,"*")
  }

  if(k == nrow(M)){

    evals <- eigen(M,TRUE,TRUE)$values

    return(-QForm::QFGauss(evals)(obs,lower.tail = TRUE,log.p=TRUE)/log(10))

  }else{

  # calculate traces
  J <- tQ %*% M
  tX <- t((Q %*% (J%*%Q)) - (M %*% Q))
  traces <- kalis:::CalcTraces(M,tX,tQ,tQ,J,from_recipient,nthreads)

  evals <- RSpectra::eigs_sym(
    function(x, args = list("M" = M, "Q" = Q, "SDy" = SDy, "from_recipient" = from_recipient, "nthreads" = nthreads)) {
      x <- x - args$Q %*% (crossprod(args$Q, x))
      x <- kalis:::MatOnlyMul(args$M, x, args$from_recipient, args$nthreads)
      x - args$Q %*% (crossprod(args$Q, x))
    },
    k = k,
    n = length(y),
    args = list("M" = M, "Q" = Q, "SDy" = SDy,
                                          "from_recipient" = from_recipient, "nthreads" = nthreads),
                              opts = list("retvec" = FALSE,   "ncv" = min(length(y), max( 4*((2*k+1)%/%4+1), 20))))$values


  return(CalcApproxPvalue(evals,
                   obs,
                   traces$trace,
                   traces$hsnorm2))
  }
}



# deviance: -2*loglik
multinom.df <- function(theta,y,X,lambda,alpha){
  
  end.of.first.layer <- length(theta)/2
  subset1 <- 2:end.of.first.layer
  subset2 <- (end.of.first.layer+2):length(theta)
  
  l1 <- theta[1] + c(X %*% theta[subset1])
  l2 <- theta[end.of.first.layer+1] + c(X %*% theta[subset2])
  el1 <- exp(l1)
  el2 <- exp(l2)
  
  -2 * sum( y[,1] * l1 + y[,2] * l2 - (1L-y[,1]) * log1p(el2) - log1p(el1) ) +
    lambda * (alpha*sum(abs(theta[c(subset1,subset2)])) + (1-alpha)*sum(theta[c(subset1,subset2)]^2))
}

# gradient of deviance function
multinom.sf <- function(theta,y,X,lambda,alpha){
  
  end.of.first.layer <- length(theta)/2
  subset1 <- 2:end.of.first.layer
  subset2 <- (end.of.first.layer+2):length(theta)
  
  l1 <- theta[1] + c(X %*% theta[subset1])
  l2 <- theta[end.of.first.layer+1] + c(X %*% theta[subset2])
  el1 <- exp(l1)
  el2 <- exp(l2)
  
  r1 <- y[,1] - el1/(1+el1)
  r2 <- y[,2] - (1-y[,1]) * el2/(1+el2)
  
  theta.sign <- rep(0L,length(theta))
  theta.sign[theta>0] <- 1L
  theta.sign[theta<0] <- -1L
  
  -2 * c(sum(r1),
         c(colSums(X*r1)),
         sum(r2),
         c(colSums(X*r2))) + 
    lambda * (alpha * c(0,theta.sign[subset1],0,theta.sign[subset2]) + (1-alpha) * 2 * c(0,theta[subset1],0,theta[subset2]))
}

multinom.percent.deviance.explained <- function(theta,y,X){
  
  end.of.first.layer <- length(theta)/2
  subset1 <- 2:end.of.first.layer
  subset2 <- (end.of.first.layer+2):length(theta)
  
  l1 <- theta[1] + c(X %*% theta[subset1])
  l2 <- theta[end.of.first.layer+1] + c(X %*% theta[subset2])
  el1 <- exp(l1)
  el2 <- exp(l2)
  
  prop <- c(colMeans(y))
  l1_null <- log(prop[1]) - log1p(-prop[1])
  l2_null <- log(prop[2]) - log1p(-prop[1]-prop[2]) # prop 1 included here because prop 2 is defined as prob of being hidden conditional on not observed.
  el1_null <- exp(l1_null)
  el2_null <- exp(l2_null)
  
  
  sum( y[,1] * l1 + y[,2] * l2 - (1L-y[,1]) * log1p(el2) - log1p(el1) ) /
    sum( y[,1] * l1_null + y[,2] * l2_null - (1L-y[,1]) * log1p(el2_null) - log1p(el1_null) )
}


fit.phi.multinom <- function(y,X){
  # X here is the design matrix of covariates WITHOUT an intercept that will be used in both layers of the multinomial
  if(any(rowSums(y) > 1) | all(rowSums(y)>0)){warning("y is not valid")}
  
  opt.res <- optim(par = runif(2*(ncol(X)+1)),fn = multinom.df,gr = multinom.sf,
                   y = y, X = X, lambda = 1e-8, alpha = 0.1,
                   method = "BFGS", control=list("reltol" = 1e-12))
  
  # check convergence
  if(opt.res$convergence){warning("BFGS had trouble converging!")}
  
  list(
    "percent.dev.explained" = multinom.percent.deviance.explained(opt.res$par,y,X),# % deviance explained (objective for tuning!)
    "pars" = unname(opt.res$par / rep(c(1,attr(X,"scaled:scale")),2)) # Calculate optimal parameters on the original scale of x
  )
}


fit_logistic_model_internal <- function(formula, x) {
  
  pops <- unique(x$pop)
  
  # Warning: for now formula must be of the form y ~ dist.jump
  classifier.list <- as.list(1:length(pops))
  names(classifier.list) <- pops
  
  coeffs <- matrix(0,length(pops),2)
  rownames(coeffs) <- pops
  
  for(i in 1:length(pops)){
    
    data.subset <- filter(x, pop == pops[i])
    
    null.cutoff <- quantile(data.subset$dist.jump[data.subset$type=="null"],0.8)
    data.subset <- filter(data.subset, (type != "null" | dist.jump <= null.cutoff))
    
    obs.cutoff <- quantile(data.subset$dist.jump[data.subset$type=="obs"],0.2)
    data.subset <- filter(data.subset, (type != "obs" | dist.jump >= obs.cutoff))
    
    hid.cutoff <- quantile(data.subset$dist.jump[data.subset$type=="hid"],0.2)
    data.subset <- filter(data.subset, (type != "hid" | dist.jump >= hid.cutoff))
    
    m1 <- glm(formula, data = data.subset)
    if(!m1$converged){warning(paste("glm failed to converge when fitting binomial logistic classifier for population",pops[i]))}
    coeffs[i,] <- m1$coefficients
  }

  res <- function(dist.jump, pop){ 
    logistic(coeffs[which(pops == pop),1] + coeffs[which(pops == pop),2] * dist.jump)
    }
  
  attr(res,"coeffs") <- coeffs
  res
}

fit_multinomial_model_internal <- function(formula, x) {
  # Warning: for now formula must be of the form y ~ dist.jump
  
  pops <- unique(x$pop)
  
  predict.multinomial <- function(dist.jump,theta){
    el1 <- exp(theta[1] + theta[2] * dist.jump)
    el2 <- exp(theta[3] + theta[4] * dist.jump)
    1-1/( (1+el1) * (1+el2) )
  }
  
  classifier.list <- as.list(1:length(pops))
  names(classifier.list) <- pops
  
  coeffs <- matrix(0,length(pops),4)
  rownames(coeffs) <- pops
  
  for(i in 1:length(pops)){
    data.subset <- filter(x, pop == pops[i])
    
    null.cutoff <- quantile(data.subset$dist.jump[data.subset$type=="null"],0.8)
    data.subset <- filter(data.subset, (type != "null" | dist.jump <= null.cutoff))
    
    obs.cutoff <- quantile(data.subset$dist.jump[data.subset$type=="obs"],0.2)
    data.subset <- filter(data.subset, (type != "obs" | dist.jump >= obs.cutoff))
    
    hid.cutoff <- quantile(data.subset$dist.jump[data.subset$type=="hid"],0.2)
    data.subset <- filter(data.subset, (type != "hid" | dist.jump >= hid.cutoff))
                          
    X <- scale(model.matrix(formula, data = data.subset)[,-1], center = FALSE) # here we use [,-1] to remove the intercept
    m1 <- fit.phi.multinom(cbind(data.subset$obs.var.ind,data.subset$hid.var.ind),X,lambda = 0)
    if(!m1$converged){warning(paste("optim failed to converge when fitting multinomial classifier for population",pops[i]))}
    coeffs[i,] <- m1$pars
  }
  
  res <- function(dist.jump, pop){ 
    predict.multinomial(dist.jump,coeffs[which(pops == pop),])
  }
  
  attr(res,"coeffs") <- coeffs
  res
}


make_fit_classifiers_internal <- function(x, model.name){

  if (model.name == "linear") {
    
    res <- function(dist.jump, pop){ dist.jump }
    
  } else if (model.name == "heavyside"){
    
    res <- function(dist.jump, pop){as.integer(dist.jump>0)}
    
  } else if (model.name == "obs.vs.null") {
    
    res <- fit_logistic_model_internal(obs.var.ind ~ dist.jump, filter(x, hid.var.ind == 0))
    
  } else if (model.name == "hid.vs.null") {
    
    res <- fit_logistic_model_internal(hid.var.ind ~ dist.jump, filter(x, obs.var.ind == 0))
    
  } else if (model.name == "obs.or.hid.vs.null.binomial") {
    
    res <- fit_logistic_model_internal(I(obs.var.ind | hid.var.ind) ~ dist.jump, x)
    
  } else if (model.name == "obs.or.hid.vs.null.multinomial") {
    
    res <- fit_multinomial_model_internal(obs.var.ind ~ dist.jump, x)
    
  } else {
    stop("model.name provided is not a recognized option")
  }
  
  res
}


make_fit_classifiers <- function(x, model.names = c("linear")){
  # x is a data.table with the samples/variables needed to classify variant clades 
  
  res <- as.list(1:length(model.names))
  names(res) <- model.names
  
  for(i in 1:length(model.names)){ res[[i]] <- make_fit_classifiers_internal(x, model.names[[i]]) }
  
  res
}



