


f <- function(x){rbinom(1,1,min(x,1))} # must take arguments x >= 0
n <- 1e3
#target.p <- c(0.05,0.25,0.9)
target.p <- c(0.1,0.75,0.95)


logistic.df <- function(theta,y,x,lambda,alpha){
  
  l1 <- -theta[1] + x * theta[2]
  el1 <- exp(l1)
  
  -2 * sum( y * l1 - log1p(el1) ) +
    lambda * (alpha*abs(theta[2]) + (1-alpha)*theta[2]^2)
}

# gradient of deviance function
logistic.sf <- function(theta,y,x,lambda,alpha){
  
  l1 <- -theta[1] + x * theta[2]
  el1 <- exp(l1)
  
  r1 <- y - el1/(1+el1)
  
  -2 * c(-sum(r1),
         sum(x*r1)) + 
    lambda * (alpha * c(0,sign(theta[2])) + (1-alpha) * 2 * c(0,theta[2]))
}




fit.logistic <- function(X){
  opt.res <- optim(par = c(2*log(1/(target.p[1]) - 1),1e-2),fn = logistic.df,gr = logistic.sf,
                   y = X[,1], x = X[,2], lambda = 0.01, alpha = 0.5,
                   method = "L-BFGS-B",lower=c(log(1/(target.p[1]*.98) - 1),0)) 
  # these lower bounds ensure that we'll never propose a negative effect size
  # check convergence
  if(opt.res$convergence){warning("L-BFGS-B had trouble converging!")}
  unname(opt.res$par)
}

# let's just run with this but instead of trying to fit all three points, we just focus on fitting 2 points at a time.
# using the previous points to give us a head start.



# Initialize starting data
#############################

# start x0 at 0.5 and doubling till we find the first 1 x0
# then run that doubling again this time starting at x1 = 2*x0 until we find first 1.
n.init <- 2
y <- x <- matrix(0,n+n.init,3)

y0 <- y1 <- 0
x0 <- 0.25
while(!y0){
  x0 <- 2 * x0
  y0 <- f(x0^2)
}
x1 <- 2 * x0
while(!y1){
  x1 <- 2 * x1
  y1 <- f(x1^2)
}

# now start with that data of three 0s at 0 and three 1s and begin fitting a logistic model.
y[1,] <- c(0,0,1)
x[1,] <- c(0,0,x0)
y[2,] <- c(0,1,1)
x[2,] <- c(0,x1,x1)


# Iterate
########################

for(i in 1:n){
  temp.pars <- fit.logistic(na.omit(cbind(c(y[1:(i+n.init-1),]),c(x[1:(i+n.init-1),]))))
  
  #temp.pars <- fit.logistic(cbind(c(y[1:(i+n.init-1),]),c(x[1:(i+n.init-1),])))
  
  #poly.f <- function(x){temp.pars[1] + temp.pars[2] * x + temp.pars[3] * 0.5 * (3*x^2-1) + temp.pars[4] * 0.5 * (5*x^3 - 3*x) + log(1/target.p - 1)}
  #x[i+n.init] <- uniroot(poly.f,interval = c(0,1),extendInt = "upX")$root
  
  x[i+n.init,] <- ( temp.pars[1] - log(1/target.p - 1) )/temp.pars[2]
  
  
  
  # y[i+n.init,] <- c(f(x[i+n.init,1]),
  #                   f(x[i+n.init,2]),
  #                   if(i <= 25){f(x[i+n.init,3])}else{NA})
  #                   #if(i <= n.burnin){sample.f(x[i+n.init,3])}else{NA_integer_})

  y[i+n.init,] <- c(if(i <= 20
                       ){f(x[i+n.init,1])}else{NA},
                    f(x[i+n.init,2]),
                    f(x[i+n.init,3]))
  
  
  print(i)
}

plot(x[,1],type="l",ylim=c(0,1))
lines(x[,2],col="red")
lines(x[,3],col="blue")
# abline(h=0.05,lty=2)
# abline(h=0.25,lty=2)
# abline(h=0.9,lty=2)
abline(h=0.1,lty=2)
abline(h=0.75,lty=2)
abline(h=0.95,lty=2)

# model y ~ sqrt(effect.size)


