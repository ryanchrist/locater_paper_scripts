# Load data.table and set the number of threads it can use
require(data.table)
setDTthreads(4L)

# Declare Paths
# (these paths will need to be modified for reproducibility)
###############################################################

# path to locater12 results folder (containing .rds files containing simulation output based on the locater12.R script)
results_path <- "~/Dropbox/Hall_group/locater_sims/locater12/"


# path to where a csv.gz with all results collected should be stored
collected_results_path <- "~/Dropbox/Hall_group/locater_sims/collected_results/"



#  BELOW HERE SHOULD NOT NEED TO BE MODIFIED FOR REPLICABILITY
##################################################################

n.samps <- 1000
samp.ind <- rep(FALSE,n.samps)
res.list <- details.list<- as.list(1:n.samps)

for(i in 1:n.samps){
  target.file <- paste0(results_path,"res_",i,".rds")
  target.file2 <- paste0(results_path,"details_",i,".rds")

  if(!file.exists(target.file)){next}
  tryCatch({
    res.list[[i]] <- readRDS(target.file)
    details.list[[i]] <- readRDS(target.file2)
    samp.ind[i] <- TRUE
  },error=function(e){print(e)})
}

res.list <- res.list[samp.ind]
res <- rbindlist(res.list,idcol = TRUE)
rm(res.list);gc()

details.list <- details.list[samp.ind]
rm(details.list);gc() # remove these unless needed for further diagnostics

# Were all replicates retreived from directory?
if(!all(samp.ind)){warning("Not all replicates were found")}

data.table::fwrite(res,paste0(collected_results_path,"locater12.csv.gz"))

# Check upon reading that storage classes are aligned
# res2 <- data.table::fread(paste0(collected_results_path,"locater12.csv.gz"))
# cbind(sapply(res,class),sapply(res2,class))
