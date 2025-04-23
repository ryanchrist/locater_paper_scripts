

SMT <- function(region.index, # a scalar index
                region.mat, # a matrix of non-overlapping segment cores
                y, # list (typically raw) phenotypes
                use.logistic, # a vector the same length as y indicating that a variable is binary
                A, # matrix of background covariates including an intercept
                null.deviance # deviance of the null model, note that this is simply the raw, unscaled RSS for linear models
){

  n <- nrow(A)

  target.loci <- region.mat[region.index,1]:region.mat[region.index,2]

  chisq_pval <- function(x,deg_freedom){-pchisq(x,df=deg_freedom,lower.tail=F,log.p=T)/log(10)}

  test_stat <- matrix(NA_real_,nrow=c(region.mat[region.index,2]-region.mat[region.index,1]+1),ncol=4*length(y) + 1)

  for(i in 1:length(target.loci)){

    t <- target.loci[i]
    phased_variant <- c(QueryCache(loci.idx = t))
    geno <- phased_variant[seq(1,2*n-1,by = 2)] + phased_variant[seq(2,2*n,by = 2)]

    test_stat[i,1] <- t
    #SNP Methods
    for(j in 1:length(y)){

      yy <- y[[j]]

      if(use.logistic[j]){
        #Dominant
        if(any(geno==0)){
          test_stat[i,(2+4*(j-1))] <- chisq_pval(null.deviance[j]-glm(yy ~ 0 + A + I(geno>0),family="binomial")$deviance,1)
        }
        # Recessive
        if(any(geno==2)){
          test_stat[i,(3+4*(j-1))] <- chisq_pval(null.deviance[j]-glm(yy ~ 0 + A + I(geno==2),family="binomial")$deviance,1)
        }
        # Additive
        if(!all(geno==1)){
          test_stat[i,(4+4*(j-1))] <- chisq_pval(null.deviance[j]-glm(yy ~ 0 + A + geno,family="binomial")$deviance,1)
        }
        # General
        if(all(0:2 %in% geno)){
          gen_null <- glm(yy ~ 0 + A + I(geno==1) + I(geno==2),family="binomial")
          test_stat[i,(5+4*(j-1))] <- chisq_pval(null.deviance[j]-gen_null$deviance,2)
        }

      } else {

        # In general we have an F test where the statistic F = (RSS_null / RSS_alt - 1) * (n - p_alt) / (p_alt - p_null)
        # which is ~ F(df1 = p_alt - p_null , df2 = n - p_alt)

        p_alt <- ncol(A) + 1

        #Dominant
        if(any(geno==0)){
          test_stat[i,(2+4*(j-1))] <- -pf((null.deviance[j] / sum( lm(yy ~ 0 + A + I(geno>0) )$residuals^2 ) - 1)*(n-p_alt),
                                          1,n-p_alt,lower.tail= FALSE, log.p = TRUE)/log(10)
        }
        # Recessive
        if(any(geno==2)){
          test_stat[i,(3+4*(j-1))] <- -pf((null.deviance[j] / sum( lm(yy ~ 0 + A + I(geno==2))$residuals^2 ) - 1)*(n-p_alt),
                                          1,n-p_alt,lower.tail= FALSE, log.p = TRUE)/log(10)
        }
        # Additive
        if(!all(geno==1)){
          test_stat[i,(4+4*(j-1))] <- -pf((null.deviance[j] / sum( lm(yy ~ 0 + A + geno)$residuals^2 ) - 1)*(n-p_alt),
                                          1,n-p_alt,lower.tail= FALSE, log.p = TRUE)/log(10)
        }
        # General
        if(all(0:2 %in% geno)){
          test_stat[i,(5+4*(j-1))] <- -pf((null.deviance[j] / sum( lm(yy ~ 0 + A + I(geno==1) + I(geno==2))$residuals^2 ) - 1)*(n-p_alt-1)/2,
                                          2,n-p_alt-1,lower.tail= FALSE, log.p = TRUE)/log(10)
        }
      }
    }
  }
  return(test_stat)
}




HaploidSMT <- function(region.index, # a scalar index
                region.mat, # a matrix of non-overlapping segment cores
                y, # list (typically raw) phenotypes
                use.logistic, # a vector the same length as y indicating that a variable is binary
                A, # matrix of background covariates including an intercept
                null.deviance # deviance of the null model, note that this is simply the raw, unscaled RSS for linear models
){

  n <- nrow(A)

  target.loci <- region.mat[region.index,1]:region.mat[region.index,2]

  chisq_pval <- function(x,deg_freedom){-pchisq(x,df=deg_freedom,lower.tail=F,log.p=T)/log(10)}

  test_stat <- matrix(NA_real_,nrow=c(region.mat[region.index,2]-region.mat[region.index,1]+1),ncol=length(y) + 1)

  for(i in 1:length(target.loci)){

    t <- target.loci[i]
    phased_variant <- c(QueryCache(loci.idx = t))
    geno <- phased_variant

    test_stat[i,1] <- t
    #SNP Methods
    for(j in 1:length(y)){

      yy <- y[[j]]

      if(use.logistic[j]){
          test_stat[i,j+1] <- chisq_pval(null.deviance[j]-glm(yy ~ 0 + A + geno,family="binomial")$deviance,1)
      } else {

        # In general we have an F test where the statistic F = (RSS_null / RSS_alt - 1) * (n - p_alt) / (p_alt - p_null)
        # which is ~ F(df1 = p_alt - p_null , df2 = n - p_alt)

        p_alt <- ncol(A) + 1

        test_stat[i,j+1] <- -pf((null.deviance[j] / sum( lm(yy ~ 0 + A + geno)$residuals^2 ) - 1)*(n-p_alt),
                                          1,n-p_alt,lower.tail= FALSE, log.p = TRUE)/log(10)
      }
    }
  }
  return(test_stat)
}



SingleHaploidSMT <- function(variant.index, # a row out of haps to test
                             haps, # a L x N matrix of haplotypes {0,1}
                             y, # list (typically raw) phenotypes
                             use.logistic, # a vector the same length as y indicating that a variable is binary
                             A, # matrix of background covariates including an intercept
                             null.deviance # vector of deviances of the null model (one entry per phenotype), note that this is simply the raw, unscaled RSS for linear models
){
  
  geno <- c(haps[variant.index,])
  
  n <- nrow(A)
  
  chisq_pval <- function(x,deg_freedom){-pchisq(x,df=deg_freedom,lower.tail=F,log.p=T)/log(10)}
  
  pvals <- rep(NA_real_,length(y))
  
  #SNP Methods
  for(j in 1:length(y)){
    
    yy <- y[[j]]
    
    if(use.logistic[j]){
      
      pvals[j] <- chisq_pval(null.deviance[j]-glm(yy ~ 0 + A + geno,family="binomial")$deviance,1)
      
    } else {
      
      # In general we have an F test where the statistic F = (RSS_null / RSS_alt - 1) * (n - p_alt) / (p_alt - p_null)
      # which is ~ F(df1 = p_alt - p_null , df2 = n - p_alt)
      p_alt <- ncol(A) + 1
      pvals[j] <- -pf((null.deviance[j] / sum( lm(yy ~ 0 + A + geno)$residuals^2 ) - 1)*(n-p_alt),
                      1,n-p_alt,lower.tail= FALSE, log.p = TRUE)/log(10)
      
    }
  }
  return(pvals)
}


