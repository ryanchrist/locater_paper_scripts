# locater15 log timing
#
# Script for reproducing benchmark timing results reported in the
# main text.
#
###################################################################

# Load Packages
require(data.table)
setDTthreads(4L)
require(stringr)

# path to where a directory where the csv.gz files with our simulation results are stored
collected_results_path <- "~/Dropbox/Hall_group/locater_sims/collected_results/"


##################################################################
##################################################################
##################################################################
#  BELOW HERE SHOULD NOT NEED TO BE MODIFIED FOR REPLICABILITY
##################################################################
##################################################################
##################################################################


d <- fread(paste0(collected_results_path,"locater15_log_statistics.txt"),header = FALSE)

d <- as.data.table(transpose(str_split(d[[2]],pattern = " in ")))

names(d) <- c("num_targets","seconds")

d$num_targets <- as.integer(transpose(str_split(transpose(str_split(d[[1]],pattern = " covering "))[[2]],pattern=" target "))[[1]])
d$seconds <- as.numeric(transpose(str_split(d[[2]],pattern = " seconds."))[[1]])

# For 30k samples it took LOCATER built on top of kalis
paste(round(sum(d$seconds)/sum(d$num_targets)/60,2),"minutes per target locus.")
#' *6.42 minutes / locus*



d <- fread(paste0(collected_results_path,"locater15_locater_only_log_statistics.txt"),header = FALSE)

s <- as.integer(transpose(str_split(transpose(str_split(d[[2]]," took "))[[2]]," seconds."))[[1]])

plot(ecdf(s))
mean(s)/60
#'*4.00 minutes*

sd(s)/60
#'*0.99 minute*

# but really what we need to benchmark is the SMT routine
# the SD routine (given X preclaculated so testing only)
# the QForm routine (given Omega precalculated so testing only)

d <- fread(paste0(collected_results_path,"locater15_testsprigs_only_log_statistics.txt"),header = FALSE)

s <- as.integer(transpose(str_split(transpose(str_split(d[[2]]," took "))[[2]]," seconds."))[[1]])

plot(ecdf(s))
mean(s)
#'*19.14 seconds*
sd(s)
#'*2.88 seconds*


d <- fread(paste0(collected_results_path,"locater15_testclademat_only_log_statistics.txt"),header = FALSE)

s <- as.integer(transpose(str_split(transpose(str_split(d[[2]]," took "))[[2]]," seconds."))[[1]])

plot(ecdf(s))
mean(s)/60
#'*3.07 minutes*
sd(s)/60
#'*0.98 minutes*


