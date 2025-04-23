require(locater)
require(data.table)
require(dplyr)
require(stringr)
shark.dir <- "~/Dropbox/shark/"
#plot.dir <- "~/Dropbox/Hall_group/locater_sims/plots/locater14/"
plot.dir <- "~/Desktop/"
if(!dir.exists(plot.dir)){dir.create(plot.dir,recursive = TRUE)}

# Declare Plotting Resources
###############################
source(paste0(shark.dir,"simulation/association_sims/ryan_viz/locater_paper_fig_utils.R"))

##################################################################
##################################################################
##################################################################
#  BELOW HERE SHOULD NOT NEED TO BE MODIFIED FOR REPLICABILITY
##################################################################
##################################################################
##################################################################


n.samps <- 900
samp.ind <- rep(FALSE,n.samps)
res.list <- as.list(1:n.samps)
for(i in 1:n.samps){
  target.file <- paste0("~/Dropbox/Hall_group/locater_sims/locater17/res_",i,".rds")

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

# check we have all of our simulations
all(samp.ind[(1:900-1L)%%9+1L <= 3]) # if FALSE we have missing jobs!
which((1:900-1L)%%9+1L <= 3)[!samp.ind[(1:900-1L)%%9+1L <= 3]]


# LOAD ARG-NEEDLE SIMS
n.samps <- 900
samp.ind <- rep(FALSE,n.samps)
res.list <- as.list(1:n.samps)
for(i in 1:n.samps){
  target.file <- paste0("~/Dropbox/Hall_group/locater_sims/argneedle16/res_",i,".rds")

  if(!file.exists(target.file)){next}
  tryCatch({
    res.list[[i]] <- readRDS(target.file)
    samp.ind[i] <- TRUE
  },error=function(e){print(e)})

  res.list[[i]][,.id:=NULL]

  print(i)
}

res_an <- rbindlist(res.list[samp.ind],idcol = TRUE)
rm(res.list);gc()



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

x_key$x_id <- seq.int(nrow(x_key))

# set observed genotyping parameters
obs_key <- data.table::CJ("prop_hidden" = c(1)) # c(0,1)  c(0,0.5,1)
obs_key$obs_id <- seq.int(nrow(obs_key))

# set simulation parameters governing y
y_key <- data.table::CJ("noise" = c("proj"), # c("proj","std")
                        "target_signal" = seq(0,150,by=3),
                        "rep" = 1:n_reps)
y_key$y_id <- seq.int(nrow(y_key))


res <- merge(res,x_key,by="x_id")
res <- merge(res,obs_key,by="obs_id")
res <- merge(res,y_key,by="y_id")

res_an <- merge(res_an,x_key,by="x_id")
res_an <- merge(res_an,obs_key,by="obs_id")
res_an <- merge(res_an,y_key,by="y_id")

#recode
res[,ac_range:= case_match(ac_range,
                           "any" ~ "all variants",
                           "150-749" ~ "rare variants only",
                           "2" ~ "doubletons only")]

res_an[,ac_range:= case_match(ac_range,
                              "any" ~ "all variants",
                              "150-749" ~ "rare variants only",
                              "2" ~ "doubletons only")]


x_key[,ac_range:= case_match(ac_range,
                             "any" ~ "all variants",
                             "150-749" ~ "rare variants only",
                             "2" ~ "doubletons only")]

res100 <- res

# Incorporate locater14 results:
# same sims with 9 causal vars but
# 10kb rather than 100kb causal window
############################################
# res10 <- readRDS("~/Dropbox/Hall_group/locater_sims/locater14_9_active_vars_res.rds")
# res <- rbind(res10,res100)
# rm(res10);rm(res100)

cutoffs <- readRDS("~/Dropbox/Hall_group/locater_sims/locater14_derived_genome_wide_cutoffs.rds")



# Confirm we have all 1000 simulations / point for power curves:
##################################################################
pa <- 9
j <- 3
  res[p_active==pa & ac_range == ac_names[j], max(smt), list(.id,rep,target_signal)][
    ,.N,target_signal]

  res[p_active==pa & ac_range == ac_names[j] & targeted == TRUE, max(tot), list(.id,rep,target_signal)][
    ,.N,target_signal] # it's expected for LOCATER to have less than 1000 obs here for smaller effect sizes because those are sims where no
  # variants met the threshold to be targeted, so if we get up to 1000 observed jobs for the larger effect sizes,
  # we know we have 1000 jobs/observations going towards estimating the smaller

  res_an[p_active==pa & ac_range == ac_names[j] & sampling_rate==1e-3, p_value, list(.id,rep,target_signal)][
    ,.N,target_signal]
  res_an[p_active==pa & ac_range == ac_names[j] & sampling_rate==1e-5, p_value, list(.id,rep,target_signal)][
    ,.N,target_signal]
  res_an[p_active==pa & ac_range == ac_names[j] & sampling_rate==1e-7, p_value, list(.id,rep,target_signal)][
    ,.N,target_signal]




# Calculate Power Curves
######################################################

pa <- 9
cairo_pdf(paste0(plot.dir,"/power_9_active_vars_hidden_vars_100kb_causal_window_v2.pdf"),width = 10,height=8)

par(oma = c(1, 1, 1, 1), mar = c(2.5, 2.5, 0.5, 0.5))

plot(0,0,type="n",ylim=c(0,1),xlim=c(0,150),las=1,bty="n",
     xlab="",ylab="",
     main="", xaxt="n",yaxt="n")
axis(2,seq(0,1,by=0.2),pos=0,las=1,cex.axis=1.4)
axis(1,seq(0,150,by=30),pos=0,las=1,cex.axis=1.4)
mtext("Signal Strength",side = 1,line = 2,cex=1.4)
mtext("Power",side = 2,line = 2,cex=1.4)

for(j in 1:length(ac_names)){
  # lines(res[x_id==i, oracle[1], list(.id,rep,target_signal)][#.id!=38
  #   ,mean(V1 > cutoff),target_signal],col=colpal[i],lwd=3,lty=1)
  lines(res[p_active==pa & ac_range == ac_names[j], max(smt), list(.id,rep,target_signal)][#.id!=38
    ,mean(V1 > cutoffs$smt),target_signal],col=colpal[1],lwd=2,lty=j)
  lines(res[p_active==pa & ac_range == ac_names[j] & targeted == TRUE, max(tot), list(.id,rep,target_signal)][#.id!=38
    ,mean(V1 > cutoffs$locater),target_signal],col=colpal[2],lwd=2,lty=j)

  lines(res_an[p_active==pa & ac_range == ac_names[j] & sampling_rate==1e-3, mean(p_value < 10^(-cutoffs$an3)), target_signal],
        col=colpal[3],lwd=2,lty=j)

  lines(res_an[p_active==pa & ac_range == ac_names[j] & sampling_rate==1e-5, mean(p_value < 10^(-cutoffs$an5)), target_signal],
        col=colpal[4],lwd=2,lty=j)

  lines(res_an[p_active==pa & ac_range == ac_names[j] & sampling_rate==1e-7, mean(p_value < 10^(-cutoffs$an7)), target_signal],
        col=colpal[5],lwd=2,lty=j)

}
legend(x = 125,y = 0.75,legend = c("SMT","LOCATER","AN3","AN5","AN7"),fill=colpal[1:5],
       border = NA,bty="n")
legend(x = 122,y = 0.55,y.intersp = 2,legend = c("all variants",expression(atop("DAC [150,750)","variants only")),"doubletons only"), lty=c(1,2,3),lwd=2, col=colpal[9], border = NA,bty="n")
dev.off()




# Calculate Power Curves with STAAR, LOCATER, and SMT
#########################################################
# cairo_pdf(paste0(plot.dir,"/power_9_active_vars_observed_vars_10kb_vs_100kb_causal_window_with_acat_oracle.pdf"),width = 10,height=8)
#
# par(oma = c(1, 1, 1, 1), mar = c(2.5, 2.5, 0.5, 0.5))
#
# plot(0,0,type="n",ylim=c(0,1),xlim=c(0,150),las=1,bty="n",
#      xlab="",ylab="",
#      main="", xaxt="n",yaxt="n")
# axis(2,seq(0,1,by=0.2),pos=0,las=1,cex.axis=1.4)
# axis(1,seq(0,150,by=30),pos=0,las=1,cex.axis=1.4)
# mtext("Signal Strength",side = 1,line = 2,cex=1.4)
# mtext("Power",side = 2,line = 2,cex=1.4)
#
# for(j in 1:length(ac_names)){
#   lines(res[p_active==9 & ac_range == ac_names[j] & targeted == TRUE & gene_size==1e4, max(tot), list(.id,rep,target_signal)][#.id!=38
#     ,mean(V1 > cutoff),target_signal],col=colpal[j],lwd=3,lty=1)
#   lines(res[p_active==9 & ac_range == ac_names[j] & targeted == TRUE & gene_size==1e5, max(tot), list(.id,rep,target_signal)][#.id!=38
#     ,mean(V1 > cutoff),target_signal],col=colpal[j],lwd=3,lty=3)
#
#   lines(res_staar[p_active==9 & ac_range == ac_names[j], acat_with_cutoff, list(.id,rep,target_signal)][
#     ,mean(acat_with_cutoff > staar_discovery_cutoff),target_signal],col=colpal[j],lwd=3,lty=5)
#
#   lines(res_staar[p_active==9 & ac_range == ac_names[j], acat_no_cutoff, list(.id,rep,target_signal)][
#     ,mean(acat_no_cutoff > staar_discovery_cutoff),target_signal],col=colpal[j],lwd=3,lty=6)
#
# }
# legend(x = 100,y = 0.6,legend = c("LOCATER 10kb Window","LOCATER 100kb Window","OA - rare vars - 100kb",
#                                   "OA - all vars - 100kb"),lty=c(1,3,5,6),lwd=3,border = NA,bty="n")
# legend(x = 100,y = 0.4,legend = ac_names,fill = colpal[1:3],border = NA,bty="n")

#dev.off()



