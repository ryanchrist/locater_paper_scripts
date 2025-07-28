# Load data.table and set the number of threads it can use
require(data.table)
setDTthreads(4L)

# Declare Paths
# (these paths will need to be modified for reproducibility)
###############################################################

# path to locater14 results folder (containing .rds files containing simulation output based on the locater12.R script)
results_path <- "~/Dropbox/Hall_group/locater_sims/locater14/"


# path to where a csv.gz with all results collected should be stored
collected_results_path <- "~/Dropbox/Hall_group/locater_sims/collected_results/"



#  BELOW HERE SHOULD NOT NEED TO BE MODIFIED FOR REPLICABILITY
##################################################################

n.samps <- 900
samp.ind <- rep(FALSE,n.samps)
res.list <- as.list(1:n.samps)
for(i in 1:n.samps){
  target.file <- paste0(results_path,"res_",i,".rds")

  if(!file.exists(target.file)){next}
  tryCatch({
    res.list[[i]] <- readRDS(target.file)
    samp.ind[i] <- TRUE
  },error=function(e){print(e)})

  # Remove redundant/un-needed variables to help conserve memory
  res.list[[i]][,c("obs.qform.T",
                   "precise",
                   "sw.thresh",
                   "eig.thresh",
                   "smt.noise",
                   "thresh",
                   "max1var",
                   "details_idx",
                   "old.sprigs",
                   "max.k",
                   "calc.obs.T",
                   "exit.status",
                   "k",
                   "test.config") := NULL]


  print(i)
}

res <- rbindlist(res.list[samp.ind],idcol = TRUE)
rm(res.list);gc()

data.table::fwrite(res,paste0(collected_results_path,"locater14.csv.gz"))
