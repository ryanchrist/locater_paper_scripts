# locater14

# This is an extension based on the pilot run of locater13
# that runs sims for three allele classes: any, 2, 150-749

# Load Simulation Packages and Parameters

# PLATFORM-SPECIFIC DIRECTORIES TO DECLARE:
run_group <- "locater_sims"
run_name <- "locater14"
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

select_causal_vars <- function(G, p){
  # G is a sparse genotype matrix
  # p is an active number of causal variants desired (a positive integer)
  # function returns columns of G that have been selected as causal variants,
  # breaking ties among repeated causal variants at random.

  # place duplicated variants into the same cluster
  R <- crossprod(G)
  R <- R == diag(R)
  R <- R & t(R)

  cluster_list <- list()
  idx_to_visit <- seq_len(nrow(R))
  ii <- 0L
  while(length(idx_to_visit)){
    ii <- ii + 1L
    cluster_list[[ii]] <- which(R[,idx_to_visit[1],drop=FALSE])
    idx_to_visit <- setdiff(idx_to_visit,cluster_list[[ii]])
  }

  if(length(cluster_list)<p){return(NULL)}

  # sample causal cluster (of the same variants according to the number of identical
  # variants in each cluster)
  cluster_list <- cluster_list[sample.int(length(cluster_list),
                                          size = p,
                                          replace = FALSE,
                                          prob = sapply(cluster_list,length))]

  # now that we've selected our unique clusters, select a representative variant
  # from each to be causal
  sapply(cluster_list,function(x){x[sample.int(length(x),size=1)]})
}

# G <- Matrix::sparseMatrix(5,1,dims=c(1000,1))
# select_causal_vars(G,1)
# select_causal_vars(G,2)
#
# G <- Matrix::sparseMatrix(c(5,5,2,5,6,7,10,2,5,6),c(1,2,3,3,3,4,5,6,6,6),dims=c(1000,6),x = c(2,2,1,2,1))
# select_causal_vars(G,1)
# select_causal_vars(G,2)


### Declare Run Parameters ###

n_reps <- 10L # number of phenotype reps

#N.haps <- rep(3e2,3)# FIXME
N.haps <- rep(2e4,3)
n <- sum(N.haps)/2
chr_length <- 1e6L

neglog10Ne <- 4
neglog10mu <- 16

use.forking <- FALSE

# set simulation parameters governing X

x_key <- CJ("n" = n,
            "p_active" = c(3,9,15),
            "gene_size" = c(1e4), # c(5e3,2e4,5e4)
            "ac_range"=c("any","2","150-749") # "2-5","6-29","30-149"
)

x_key$id <- seq.int(nrow(x_key))

x_key <- x_key[(as.integer(job.index)-1L)%%9+1L,]

# set observed genotyping parameters
obs_key <- data.table::CJ("prop_hidden" = c(0)) # c(0,1)  c(0,0.5,1)
obs_key$id <- seq.int(nrow(obs_key))

# set simulation parameters governing y
y_key <- data.table::CJ("noise" = c("proj"), # c("proj","std")
                        "target_signal" = seq(0,150,by=3),
                        "rep" = 1:n_reps)
y_key$id <- seq.int(nrow(y_key))


# initialize storage for results
res_list <- as.list(1:(nrow(x_key)*nrow(obs_key)))





# LOOP OVER X_KEY (ARCHITECTURES)
#####################################
for(i in 1:nrow(x_key)){

  #i <- 1 # FIXME


  # SIMULATE GENOMIC REGION
  #############################################

  p_active <- x_key$p_active[i]

  repeat({

    # Simulate Genomic Region
    region <- sim_genomic_region(job.index,
                                 shark.dir,
                                 champs.path = "champs",
                                 sim.output.dir,
                                 N = N.haps,
                                 chr.length = chr_length)

    pos <- region$pos
    p_chr <- length(pos)


    # define and load functional window
    functional_window <- floor((chr_length + x_key$gene_size[i] * c(-1,1))/2)
    loci.idx <- which(pos >= functional_window[1] & pos <= functional_window[2])

    p_gene <- length(loci.idx) # number of variants in functional window

    X <- ReadHaplotypes(region$haps,
                        loci.idx = loci.idx,
                        transpose = TRUE)$haps
    X <- t(as(X[,seq(1,2*n,by=2)] + X[,seq(2,2*n,by=2)],"sparseMatrix"))

    # select causal variants

    # enforce optional allele count restriction
    if(x_key$"ac_range"[i]!="any"){
      ac_range <- as.integer(strsplit(x_key$"ac_range"[i],"-")[[1]])

      ac <- colSums(X)

      if(length(ac_range)==1){
        ac_index <- which(ac==ac_range)
      } else {
        ac_index <- which(ac >= ac_range[1] & ac <= ac_range[2])
      }
      # move on to next architecture if too few variants in ac range in gene in window
      if(length(ac_index) < p_active){next}

      loci.idx <- loci.idx[ac_index]
      X <- X[,ac_index,drop=FALSE]
    }

    rep_var <- select_causal_vars(X, p_active)

    if(!is.null(rep_var)){
      # indicating that there were enough unique variants to get
      # p_active variants, keep the chromosome simulation we have and move on
      break}
  })


  active_predictors_idx <- loci.idx[rep_var]
  X <- as.matrix(X[,rep_var,drop=FALSE])



  # SIMULATE PHENOTYPE
  ######################################

  A <- cbind(1,as(t(fac2sparse(rep(1:3,times=N.haps/2))),"denseMatrix"),rbinom(n,1,0.5),rnorm(n))
  rsA <- rowSums(A)
  qr.A <- qr(A)

  X_active_norm <- qr.resid(qr.A, X)
  X_active_norm <- scale(X_active_norm, center = FALSE, scale = sqrt(colSums(X_active_norm^2)))

  qr.X_active_norm <- qr(X_active_norm)
  v <- as.vector(forwardsolve(qr.R(qr.X_active_norm),matrix(1,p_active,1),transpose = T))
  v <- v/sqrt(sum(v^2))

  effect_vec <- qr.Q(qr.X_active_norm) %*% v

  z <- matrix(rnorm(n*n_reps),n,n_reps)
  z_resid <- qr.resid(qr.X_active_norm,z)

  y <- matrix(0,nrow=n,ncol=nrow(y_key))
  for(ii in 1:nrow(y_key)){
    y[,ii] <- rsA +
      sqrt(qchisq(-log(10)*y_key$target_signal[ii],df = p_active,lower.tail = FALSE,log.p = TRUE)) * effect_vec +
      if(y_key$noise[ii]=="proj"){z_resid[,y_key$rep[ii]]}else{z[,y_key$rep[ii]]}
  }
  rm(z); rm(z_resid)


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


    if(obs_key$prop_hidden[j]){
      locus_idx_subset <- seq_len(length(pos))[
        -active_predictors_idx[sample.int(p_active,
                                          size = ceiling(p_active * obs_key$prop_hidden[j]))]]
      temp_pos <- pos[locus_idx_subset]
    }else{
      locus_idx_subset <- NULL
      temp_pos <- pos
    }

    CacheHaplotypes(region$haps,
                    transpose = TRUE,
                    hdf5.pkg = "rhdf5",
                    warn.singletons = TRUE,
                    loci.idx = locus_idx_subset)

    # Run Tests
    ############

    # Identify Target Loci
    smt.res <- locater::TestCachedMarkers(y, A = A)

    # Take at most 100 target loci / phenotype
    smt_thresh <- pmax(4,unname(apply(smt.res,2,quantile,probs=1-50/nrow(smt.res),type=1)))

    # for smt.res, rows are loci, cols are phenotypes
    target.loci.ind <- rep(FALSE,nrow(smt.res))
    target.mat <- matrix(FALSE,nrow(smt.res),ncol(smt.res))
    for(ii in 1:ncol(smt.res)){
      target.mat[,ii] <- smt.res[,ii] > smt_thresh[ii]
      target.loci.ind <- target.loci.ind | target.mat[,ii]
    }

    target.loci <- which(target.loci.ind)

    # Run Methods @ Target Loci
    # Could consider dropping windows down around these target loci
    # and making that a comparison
    # for(t in target_loci){
    #   # Cauchy
    #   # gdistill
    #   # STAAR
    # }

    # RUN ORACLE
    oracle_signal <- neg_log10_simple_anova(colSums(qr.resid(qr.A,y)^2),
                                            ncol(A),
                                            colSums(qr.resid(qr(cbind(A,X)),y)^2),
                                            ncol(A)+p_active,
                                            n)

    # Run LOCATER over target loci
    #################################

    print(paste("Number of Target Loci:",length(target.loci)))

    pars <- Parameters(
      CalcRho(cM = diff(if(is.null(locus_idx_subset)){
        region$map
      }else{
        region$map[locus_idx_subset]
      }), s = 10^-neglog10Ne),
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
                        #"max.k" = 128, # FIXME
                        "sw.thresh" = 0,
                        "eig.thresh" = 0,
                        "calc.obs.T" = FALSE),
                      verbose = TRUE,
                      num.ckpts = 3L,
                      ckpt.first.locus = FALSE,
                      use.forking = use.forking,
                      nthreads = nthreads)

      res$targeted <- c(t(target.mat[target.loci,]))

      # Annotate res
      res$p_chr <- p_chr
      res$p_gene <- p_gene
      res$x_id = x_key$id[i]
      res$obs_id <- obs_key$id[j]
      res$y_id = rep(y_key$id,length(target.loci))
      res$smt_thresh <- rep(smt_thresh,length(target.loci))
      res$oracle <- rep(oracle_signal,length(target.loci))

      res$pos <- temp_pos[res$locus.idx]

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
saveRDS(rbindlist(res_list, use.names = TRUE,idcol="details_idx"), paste0(res.storage.path,"res_",job.index,".rds"))
details <- lapply(res_list,function(x){attr(x,"details")})

saveRDS(details, paste0(res.storage.path,"details_",job.index,".rds"))
print("ALL SIMS COMPLETED.")





