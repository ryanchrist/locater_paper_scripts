# NULL SIMS

# Load Simulation Packages and Parameters

# PLATFORM-SPECIFIC DIRECTORIES TO DECLARE:
run_group <- "locater_sims"
run_name <- "locater12"
#shark.dir <- "/Users/ryan/Dropbox/shark/" # FIXME
shark.dir <- "/rchrist/shark/"
user <- "rchrist"

#res.storage.path <- paste0("/Users/ryan/Documents/lab/hall_group/temp/res/",run_group,"/",run_name,"/") # FIXME
res.storage.path <- paste0("/",user,"/data/simulation/",run_group,"/",run_name,"/")
if(!dir.exists(res.storage.path)){dir.create(res.storage.path,recursive = TRUE)}

#sim.output.dir <- paste0("/Users/ryan/Documents/lab/hall_group/temp/sims/",run_group,"/",run_name,"/") # FIXME
sim.output.dir <- paste0("/",user,"/data/simulation/temp/",run_group,"/",run_name,"/")
if(!dir.exists(sim.output.dir)){dir.create(sim.output.dir,recursive = TRUE)}

# Assuming that nthreads is the first commandArg and job.index is the second as usual
nthreads <- as.integer(commandArgs(TRUE)[1])
job.index <- commandArgs(TRUE)[2]
#job.index <- 76L # FIXME
#nthreads <- 1L # FIXME
.libPaths(c("/usr/local/lib/R/site-library",.libPaths()))

############################################################
# BELOW HERE SHOULD NOT NEED TO BE EDITED FOR REPLICATION
############################################################

# set seed based on job_index
set.seed(job.index)

#### Load Libraries ####
library(rdistill)
library(kalis)
library(data.table)
library(locater)
library(Matrix)

###  Load Shark Functions  ###
source(paste0(shark.dir,"simulation/shark_simulation_header.R"))


# declare local utility functions
neg_log10_simple_anova <- function(rss1,p1,rss2,p2,n){
  # p1 is number of predictors in model 1 yielding residual sum of squares 1 (rss1)
  # p2 is number of predictors in model 2 yielding residual sum of squares 2 (rss2)
  if(n <= p1 | n <= p2){stop("n = length(y) must be greater than p1 and p2 to run F test")}
  if(p1 >= p2){stop("p2 must be greater than p1 to run F test")}
  -pf( ((rss1 - rss2)/(p2-p1)) / (rss2/(n-p2)), p2-p1, n-p2, lower.tail = FALSE,log.p = TRUE)/log(10)
}


### Declare Run Parameters ###

n_reps <- 1000L # number of phenotype reps

#N.haps <- rep(3e2,3)# FIXME
N.haps <- rep(2e4,3)
n <- sum(N.haps)/2
chr_length <- 1e6L

neglog10Ne <- 4
neglog10mu <- 16

use.forking <- FALSE

# set simulation parameters governing X
x_key <- CJ("n" = n,
            "p_active" = 0L, # c(2,4,8,16)
            "gene_size" = 1e4L) # c(5e3,2e4,5e4)
x_key$id <- seq.int(nrow(x_key))

# set observed genotyping parameters
obs_key <- data.table::CJ("prop_hidden" = c(0)) #c(0,0.5,1)
obs_key$id <- seq.int(nrow(obs_key))

# set simulation parameters governing y
y_key <- data.table::CJ("noise" = c("std"), # c("proj","std")
                        "target_signal" = 0,
                        "rep" = 1:n_reps)
y_key$id <- seq.int(nrow(y_key))


# initialize storage for results
res_list <- as.list(1:(nrow(x_key)*nrow(obs_key)))



# Simulate Genomic Region & background covariates

region <- sim_genomic_region(job.index,
                             shark.dir,
                             champs.path = "champs",
                             sim.output.dir,
                             N = N.haps,
                             chr.length = chr_length)


pos <- region$pos

p_chr <- length(pos)

A <- cbind(1,as(t(fac2sparse(rep(1:3,times=N.haps/2))),"denseMatrix"),rbinom(n,1,0.5),rnorm(n))
rsA <- rowSums(A)
qr.A <- qr(A)


# LOOP OVER X_KEY (ARCHITECTURES)
#####################################
for(i in 1:nrow(x_key)){

  #i <- 1 # FIXME
  # SIMULATE PHENOTYPES
  #############################################

  y <- rsA + matrix(rnorm(n*n_reps),n,n_reps)

  # set.seed(43)
  # res <- rdistill::rdistill_pivot_par(y = y, x = X, Q = Q, max_num_causal = 16)
  # set.seed(41)
  # res2 <- rdistill::rdistill_pivot_par(y = y, x = X, Q = Q, max_num_causal = 16)
  #
  # plot(-log10(res$p_value),-log10(res2$p_value))
  # abline(0,1)


  # LOOP OVER OBS_KEY (WHETHER OR NOT CAUSAL VARIANTS OBSERVED)
  ###################################################################

  for(j in 1:nrow(obs_key)){
    # Reload cache depending on whether causal variants should be hidden before passing to iterate over target loci

    #j <- 1 # FIXME

    CacheHaplotypes(region$haps,
                    transpose = TRUE,
                    hdf5.pkg = "rhdf5",
                    warn.singletons = TRUE)

    # Run Tests
    ############

    # Identify Target Loci
    # smt.res <- locater::TestCachedMarkers(y, A = A)
    #
    # # Take at most 100 target loci / phenotype
    # smt_thresh <- pmax(4,unname(apply(smt.res,2,quantile,probs=1-50/nrow(smt.res),type=1)))
    #
    # # for smt.res, rows are loci, cols are phenotypes
    # target.loci.ind <- rep(FALSE,nrow(smt.res))
    # target.mat <- matrix(FALSE,nrow(smt.res),ncol(smt.res))
    # for(ii in 1:ncol(smt.res)){
    #   target.mat[,ii] <- smt.res[,ii] > smt_thresh[ii]
    #   target.loci.ind <- target.loci.ind | target.mat[,ii]
    # }

    target.loci <- floor(L()/2)  # which(target.loci.ind)


    # Run LOCATER over target loci
    #################################

    print(paste("Number of Target Loci:",length(target.loci)))

    pars <- Parameters(CalcRho(diff(region$map),s = 10^-neglog10Ne),
                       mu = 10^-neglog10mu,
                       use.speidel = TRUE)

    #target.loci <- target.loci[100:105] # FIXME

    #set.seed(100)
    start1 <- proc.time()[3]
    tryCatch({
      res <- TestLoci(y, pars, target.loci = target.loci,
                      A = A,
                      test.opts = data.frame(
                        "smt.noise" = c("raw"),#,TRUE,FALSE,FALSE), # Clade-free testing options (eg: SMT, might be more complex)
                        "thresh" = 0.2,
                        "old.sprigs" = c(FALSE), # Clade calling options
                        "max1var" = c(TRUE),#,FALSE,TRUE,FALSE),
                        "max.k" = 512,
                        "sw.thresh" = 0,
                        "eig.thresh" = 0,
                        "calc.obs.T" = TRUE),
                      verbose = TRUE,
                      num.ckpts = 0L,
                      ckpt.first.locus = FALSE,
                      use.forking = use.forking,
                      nthreads = nthreads)


      # Annotate res
      res$p_chr <- p_chr
      res$x_id = x_key$id[i]
      res$obs_id <- obs_key$id[j]
      res$y_id = rep(y_key$id,length(target.loci))

      res_list[[j + (i-1)*nrow(obs_key)]] <- res

      # Select set of target loci for follow up
      # new.target.loci <- as.list(1:max(approx.res$test.config))
      # for(i in 1:length(new.target.loci)){
      #   new.target.loci[[i]] <- as.integer(
      #     approx.res[test.config==i,locus.idx[which.thresh.middle(tot,screening.thresh)],phenotype]$V1)}
      # new.target.loci <- sort(unique(unlist(new.target.loci)))
      #
      # Run TestLoci with sw.approx = FALSE to revisit and exactly evaluate candidate loci
      # exact.res <- TestLoci(y, pars, target.loci = new.target.loci,
      #                       A = A,
      #                       sw.approx = FALSE,
      #                       test.opts = test.opts,
      #                       verbose = TRUE,
      #                       num.ckpts = num.ckpts,
      #                       ckpt.first.locus = FALSE,
      #                       use.forking = use.forking,
      #                       nthreads = nthreads)
      #
    }, error = function(e){print(e)})

    print(paste("TestLoci call",j + (i-1)*nrow(obs_key),"out of",nrow(obs_key)*nrow(x_key),
                "covering",length(target.loci),"target loci done in",proc.time()[3] - start1,"seconds."))


    #res_orig <- res  # TestLoci returns same answer if set.seed is set identical ahead of it.
    # plot(res_orig$rd,res$rd); abline(0,1)
    #
    # plot(res_orig$num.sprigs,res$num.sprigs)
    #
    #plot(res_orig$rd,pmax(res_orig$rd,res$rd)-log10(2)); abline(0,1)
    #
    # plot(ecdf(pmax(res_orig$rd,res$rd)-log10(2) - res_orig$rd))
    # abline(v=0)
    #
    # plot(ecdf(pmax(res_orig$tot,res$tot)-log10(2) - res_orig$tot))
    # abline(v=0)
    #
    # plot(ecdf(abs(res_orig$rd-res$rd)))
    #plot(res$tot,res_orig$tot-res$tot,xlim=c(0,8))

    #all.equal(res_orig,res)



    # Transform target.mat into a new  are phenotypes # build this in as new column in res

  }
}


# Store output
saveRDS(rbindlist(res_list, use.names = TRUE), paste0(res.storage.path,"res_",job.index,".rds"))
details <- lapply(res_list,function(x){attr(x,"details")})
saveRDS(details, paste0(res.storage.path,"details_",job.index,".rds"))
print("ALL SIMS COMPLETED.")





