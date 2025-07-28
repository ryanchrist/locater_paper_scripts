# locater16
#
# Script for reproducing power, localization, and source plots
# based on simulations where the causal variants are OBSERVED
# and all within a 100kb causal region
#
# Should be run AFTER running locater14_viz_revision1.R
#
###################################################################

# Load Packages
require(data.table)
setDTthreads(4L)
require(dplyr)
require(stringr)

# Declare path to shark directory (top of the locater_paper_scripts repository)
shark.dir <- "~/Dropbox/shark/"

# Declare path for putting plots
plot.dir <- "~/Desktop/"
if(!dir.exists(plot.dir)){dir.create(plot.dir,recursive = TRUE)}

# path to where a directory where the csv.gz files with our simulation results are stored
collected_results_path <- "~/Dropbox/Hall_group/locater_sims/collected_results/"

# path_to_temp_storage . Must match definition of path_to_temp_storage
# used in locater14_viz_revision_1.R since that is where locater14_viz_revision_1.R
# will have stored the significance cutoffs to be used
path_to_temp_storage <- "~/Dropbox/Hall_group/locater_sims/"

# Import cutoffs
cutoffs <- readRDS(paste0(path_to_temp_storage,"locater14_derived_genome_wide_cutoffs.rds"))




##################################################################
##################################################################
##################################################################
#  BELOW HERE SHOULD NOT NEED TO BE MODIFIED FOR REPLICABILITY
##################################################################
##################################################################
##################################################################


# Load Plotting Resources
###############################
source(paste0(shark.dir,"simulation/association_sims/ryan_viz/locater_paper_fig_utils.R"))


# LOAD DATA
##################################################################

# LOAD LOCATER SIMS
res <- fread(paste0(collected_results_path,"locater16.csv.gz"))

# LOAD ARG-NEEDLE SIMS
res_an <- fread(paste0(collected_results_path,"argneedle16.csv.gz"))



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
obs_key <- data.table::CJ("prop_hidden" = c(0)) # c(0,1)  c(0,0.5,1)
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



# Calculate Power Curves
######################################################

pa <- 9
cairo_pdf(paste0(plot.dir,"/power_9_active_vars_observed_vars_100kb_causal_window_v2.pdf"),width = 10,height=8)

par(oma = c(1, 1, 1, 1), mar = c(2.5, 2.5, 0.5, 0.5))

plot(0,0,type="n",ylim=c(0,1),xlim=c(0,150),las=1,bty="n",
     xlab="",ylab="",
     main="", xaxt="n",yaxt="n")
axis(2,seq(0,1,by=0.2),pos=0,las=1,cex.axis=1.4)
axis(1,seq(0,150,by=50),pos=0,las=1,cex.axis=1.4)
axis(1,seq(0,150,by=10),pos=0,las=1,labels = FALSE)
mtext("Total Association Signal Strength",side = 1,line = 2,cex=1.4)
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



