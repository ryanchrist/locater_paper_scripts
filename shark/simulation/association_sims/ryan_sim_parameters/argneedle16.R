# argneedle14

# This is an exact copy of locater14 that saves the simulated ARG and then runs ARG-Needle on it

# Load Simulation Packages and Parameters

# PLATFORM-SPECIFIC DIRECTORIES TO DECLARE:
run_group <- "locater_sims"
run_name <- "argneedle16"
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
            "p_active" = c(9),
            "gene_size" = c(1e5), # c(5e3,2e4,5e4)
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
                                 chr.length = chr_length,
                                 save_tree = TRUE)

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

  # Below here we start making changes for running ARG-NEEDLE

  # Make Q available for passing residuals of y to ARG-NEEDLE
  Q <- qr.Q(qr.A)

  for(j in 1:nrow(obs_key)){
    # Reload cache depending on whether causal variants should be hidden before passing to iterate over target loci

    start1 <- proc.time()[3]

    # Write ARG Needle Sample file

    sample_path <- paste0(sim.output.dir,"job_",job.index,".sample")

    fwrite(x = data.table("ID_1" = 0:n,
                          "ID_2" = 0:n,
                          "missing" = integer(n+1)),
           file = sample_path,
           sep = " ")

    pheno_path <- paste0(sim.output.dir,"job_",job.index,".pheno")

    out_path <- paste0(sim.output.dir,"results_",job.index)

    res <- data.table("chisq" = rep_len(NA_real_,nrow(y_key)))
    res$p_value <- NA_real_
    res$avg_abs_dist <- NA_real_
    res$abs_avg_dist <- NA_real_
    res$avg_MAC <- NA_integer_
    res$num_tied_edges <- NA_integer_

    res3 <- res2 <- res

    res$sampling_rate <- 1e-3
    res2$sampling_rate <- 1e-5
    res3$sampling_rate <- 1e-7



    for(ii in 1:nrow(y_key)){

      # Write phenotypes to text file formatted as expected by arg_association
      fwrite(x = data.table("ID_1" = seq_len(n),
                            "ID_2" = seq_len(n),
                            "pheno" = y[,ii] - Q %*% crossprod(Q,y[,ii])),
             file = pheno_path,
             sep = "\t",
             col.names = FALSE)


      # Declare base arguments
      base.args <- c(paste("--arg_path",region$ts),
                     paste("--arg_sample_path",sample_path),
                     paste("--arg_id",paste0("chr1.chunk",job.index)),
                     paste("--residualised_pheno_path",pheno_path),
                     paste("--out_path",out_path),
                     paste("--random_seed",job.index),
                     "--min_mac 2",
                     "--sampling_rate 1e-3")

      tryCatch({

        # run arg-needle
        system2("python3", c(paste0(shark.dir,"arg-needle/arg_needle_association_without_logging.py"),
                             base.args))
        #system2("ls",c("-lthra","/rchrist/data/simulation/temp/locater_sims/argneedle14/"))

        # collect and store results
        system2("gunzip", c("-f",paste0(out_path,".tab.gz")))

        an_res <- fread(paste0(out_path,".tab"))
        an_res[,idx:=.I]

        temp_max <- an_res[,max(CHISQ)]

        max_carriers <- an_res[CHISQ == temp_max,idx]
        res$num_tied_edges[ii] <- length(max_carriers)
        res$chisq[ii] <- an_res[max_carriers,CHISQ[1]]
        res$p_value[ii] <- an_res[max_carriers,P[1]]
        res$avg_MAC[ii] <- an_res[max_carriers,mean(MAC)]

        # These two ways of averaging will agree in cases when START_BP and END_BP
        # are both on the same side of 5e5.  However, they will disagree when that's not true.
        res$abs_avg_dist[ii] <- mean(an_res[max_carriers,abs(5e5 - ((START_BP+1L) + (END_BP+1L))/2)])

        temp_mean_abs <- double(length(max_carriers))
        for(jj in 1:length(max_carriers)){
          temp_mean_abs[jj] <- an_res[max_carriers[jj],mean(abs(5e5 - (START_BP+1L):(END_BP+1L)))]
        }
        res$avg_abs_dist[ii] <- mean(temp_mean_abs)
        res$y_id <- y_key$id
        res$x_id <- x_key$id[i]
        res$obs_id <- obs_key$id[j]

        # Redo for different cutoff!

        temp_max <- an_res[MU_STAR <= 1e-5, max(CHISQ)]

        max_carriers <- an_res[MU_STAR <= 1e-5 & CHISQ == temp_max,idx]
        res2$num_tied_edges[ii] <- length(max_carriers)
        res2$chisq[ii] <- an_res[max_carriers,CHISQ[1]]
        res2$p_value[ii] <- an_res[max_carriers,P[1]]
        res2$avg_MAC[ii] <- an_res[max_carriers,mean(MAC)]

        # These two ways of averaging will agree in cases when START_BP and END_BP
        # are both on the same side of 5e5.  However, they will disagree when that's not true.
        res2$abs_avg_dist[ii] <- mean(an_res[max_carriers,abs(5e5 - ((START_BP+1L) + (END_BP+1L))/2)])

        temp_mean_abs <- double(length(max_carriers))
        for(jj in 1:length(max_carriers)){
          temp_mean_abs[jj] <- an_res[max_carriers[jj],mean(abs(5e5 - (START_BP+1L):(END_BP+1L)))]
        }
        res2$avg_abs_dist[ii] <- mean(temp_mean_abs)
        res2$y_id <- y_key$id
        res2$x_id <- x_key$id[i]
        res2$obs_id <- obs_key$id[j]


        # Redo for different cutoff!

        temp_max <- an_res[MU_STAR <= 1e-7, max(CHISQ)]

        max_carriers <- an_res[MU_STAR <= 1e-7 & CHISQ == temp_max,idx]
        res3$num_tied_edges[ii] <- length(max_carriers)
        res3$chisq[ii] <- an_res[max_carriers,CHISQ[1]]
        res3$p_value[ii] <- an_res[max_carriers,P[1]]
        res3$avg_MAC[ii] <- an_res[max_carriers,mean(MAC)]

        # These two ways of averaging will agree in cases when START_BP and END_BP
        # are both on the same side of 5e5.  However, they will disagree when that's not true.
        res3$abs_avg_dist[ii] <- mean(an_res[max_carriers,abs(5e5 - ((START_BP+1L) + (END_BP+1L))/2)])

        temp_mean_abs <- double(length(max_carriers))
        for(jj in 1:length(max_carriers)){
          temp_mean_abs[jj] <- an_res[max_carriers[jj],mean(abs(5e5 - (START_BP+1L):(END_BP+1L)))]
        }
        res3$avg_abs_dist[ii] <- mean(temp_mean_abs)
        res3$y_id <- y_key$id
        res3$x_id <- x_key$id[i]
        res3$obs_id <- obs_key$id[j]


        res_list[[j + (i-1)*nrow(obs_key)]] <- rbindlist(list(res,res2,res3))

      }, error = function(e){print(e)})

    }

    print(paste("ARG NEEDLE call",j + (i-1)*nrow(obs_key),"out of",nrow(obs_key)*nrow(x_key),
                "done in",proc.time()[3] - start1,"seconds."))

  }
}


# Store output
saveRDS(rbindlist(res_list, use.names = TRUE,idcol=TRUE), paste0(res.storage.path,"res_",job.index,".rds"))
print("ALL SIMS COMPLETED.")





