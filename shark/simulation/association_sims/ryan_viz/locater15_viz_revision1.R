# locater15
#
# Script for reproducing power, localization, and source plots
# based on simulations where the causal variants are HIDDEN
# and all within a 10kb causal region
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
res <- fread(paste0(collected_results_path,"locater15.csv.gz"))

# LOAD ARG-NEEDLE SIMS
# here we load argneedle14 because argneedle yields the same testing
# results whether or not the causal variants were hidden in the
# original haplotype dataset since we just
# give argneedle the true underlying ARG.
# argneelde just places new mutations on that ARG according to
# a Poisson process.
res_an <- fread(paste0(collected_results_path,"argneedle14.csv.gz"))





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

# res_staar <- merge(res_staar,x_key,by="x_id")
# res_staar <- merge(res_staar,obs_key,by="obs_id")
# res_staar <- merge(res_staar,y_key,by="y_id")

res_an <- merge(res_an,x_key,by="x_id")
res_an <- merge(res_an,obs_key,by="obs_id")
res_an <- merge(res_an,y_key,by="y_id")


#recode
res[,ac_range:= case_match(ac_range,
                           "any" ~ "all variants",
                           "150-749" ~ "rare variants only",
                           "2" ~ "doubletons only")]

# res_staar[,ac_range:= case_match(ac_range,
#                                  "any" ~ "all variants",
#                                  "150-749" ~ "rare variants only",
#                                  "2" ~ "doubletons only")]


res_an[,ac_range:= case_match(ac_range,
                              "any" ~ "all variants",
                              "150-749" ~ "rare variants only",
                              "2" ~ "doubletons only")]


x_key[,ac_range:= case_match(ac_range,
                             "any" ~ "all variants",
                             "150-749" ~ "rare variants only",
                             "2" ~ "doubletons only")]




# Calculate Power Curves with SMT, LOCATER and ARG-NEEDLE
###########################################################

for(pa in unique(x_key$p_active)){

  cairo_pdf(paste0(plot.dir,"/power_",pa,"_active_vars_hidden_vars_v2.pdf"),width = 10,height=8)

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
}




# CALC 80% POWER SIGNAL CUTOFFS
#################################

x_key$smt_80p_signal <- x_key$locater_80p_signal <-
  x_key$an3_80p_signal <- x_key$an5_80p_signal <- x_key$an7_80p_signal <-
  x_key$acat_all <- x_key$acat_rare <- 0


x_key$smt_localization <- x_key$locater_localization <-
  x_key$an3_localization <- x_key$an5_localization <- x_key$an7_localization <- 0


for(pa in unique(x_key$p_active)){
  for(j in 1:length(ac_names)){
    x_key[p_active==pa & ac_range == ac_names[j], "smt_80p_signal"] <- res[
      p_active==pa & ac_range == ac_names[j], max(smt), list(.id,rep,target_signal)
    ][
      ,mean(V1 > cutoffs$smt),target_signal
    ][
      ,approx(V1,target_signal,xout = 0.8,method = "linear", ties = "ordered")$y
    ]

    x_key[p_active==pa & ac_range == ac_names[j], "locater_80p_signal"] <- res[
      p_active==pa & ac_range == ac_names[j] & targeted == TRUE, max(tot), list(.id,rep,target_signal)
    ][
      ,mean(V1 > cutoffs$locater),target_signal
    ][
      ,approx(V1,target_signal,xout = 0.8,method = "linear", ties = "ordered")$y
    ]

    x_key[p_active==pa & ac_range == ac_names[j], "an3_80p_signal"] <- res_an[
      p_active==pa & ac_range == ac_names[j] & sampling_rate==1e-3,
      mean(p_value < 10^(-cutoffs$an3)),
      target_signal
    ][,approx(V1,target_signal,xout = 0.8,method = "linear", ties = "ordered")$y]

    x_key[p_active==pa & ac_range == ac_names[j], "an5_80p_signal"] <- res_an[
      p_active==pa & ac_range == ac_names[j] & sampling_rate==1e-5,
      mean(p_value < 10^(-cutoffs$an5)),
      target_signal
    ][,approx(V1,target_signal,xout = 0.8,method = "linear", ties = "ordered")$y]

    x_key[p_active==pa & ac_range == ac_names[j], "an7_80p_signal"] <- res_an[
      p_active==pa & ac_range == ac_names[j] & sampling_rate==1e-7,
      mean(p_value < 10^(-cutoffs$an7)),
      target_signal
    ][,approx(V1,target_signal,xout = 0.8,method = "linear", ties = "ordered")$y]


    x_key[p_active==pa & ac_range == ac_names[j], "smt_localization"] <-
      res[p_active==pa & ac_range == ac_names[j] & smt >= cutoff,
          2*mean(abs(5e5-pos[smt==max(smt)]))/1e3,
          list(.id,rep,target_signal)
      ][,quantile(V1,0.8),target_signal
      ][target_signal >= x_key[p_active==pa & ac_range == ac_names[j],smt_80p_signal],]$V1[1]

    x_key[p_active==pa & ac_range == ac_names[j], "locater_localization"] <-
      res[p_active==pa & ac_range == ac_names[j] & targeted == TRUE & tot >= cutoff,
          2*mean(abs(5e5-pos[tot==max(tot)]))/1e3,
          list(.id,rep,target_signal)
      ][,quantile(V1,0.8),target_signal
      ][target_signal >= x_key[p_active==pa & ac_range == ac_names[j],locater_80p_signal],]$V1[1]

    x_key[p_active==pa & ac_range == ac_names[j], "an3_localization"] <-
      res_an[p_active==pa & ac_range == ac_names[j] & sampling_rate==1e-3,
             2*avg_abs_dist/1e3,
             list(.id,rep,target_signal)
      ][,quantile(V1,0.8,na.rm=TRUE),target_signal
      ][target_signal >= x_key[p_active==pa & ac_range == ac_names[j],an3_80p_signal],]$V1[1]

    x_key[p_active==pa & ac_range == ac_names[j], "an5_localization"] <-
      res_an[p_active==pa & ac_range == ac_names[j] & sampling_rate==1e-5,
             2*avg_abs_dist/1e3,
             list(.id,rep,target_signal)
      ][,quantile(V1,0.8,na.rm=TRUE),target_signal
      ][target_signal >= x_key[p_active==pa & ac_range == ac_names[j],an5_80p_signal],]$V1[1]

    x_key[p_active==pa & ac_range == ac_names[j], "an7_localization"] <-
      res_an[p_active==pa & ac_range == ac_names[j] & sampling_rate==1e-7,
             2*avg_abs_dist/1e3,
             list(.id,rep,target_signal)
      ][,quantile(V1,0.8,na.rm=TRUE),target_signal
      ][target_signal >= x_key[p_active==pa & ac_range == ac_names[j],an7_80p_signal],]$V1[1]

  }
}


# MAKE POWER DOTPLOT
############################
custom_dotplot <- function(x, x_max, col = "#000000"){

  n <- nrow(x)
  p <- ncol(x)

  # color blind palette as recommended: https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible

  plot_bg <- paste0(col,c("00","00","AA","AA","AA"))

  op <- par(oma = c(0,1,0,0),
            mar = c(4.1, 6.1, 3.1, 1.1)
            #mar = c(4.1, 5.1, 2.1, 6.1)
  )
  on.exit(par(op))
  #"#000000" "#D55E00" "#56B4E9" "#F0E442" "#009E73" "#CC79A7" "#E69F00"

  plot(0, type="n",
       xlim=c(0,ceiling(max(x_max))),ylim=c(0,n),
       bty="n",ylab="",xlab="",xaxt="n",yaxt="n"
  )
  segments(rep(0,nrow(x)),1:nrow(x),rep(x_max,nrow(x)),1:nrow(x),
           col = gray.colors(3,0.7,0.9),lwd=2,lty=1)
  for(i in 1:nrow(x)){
    points(as.numeric(x[i,c("smt_80p_signal",
                            "locater_80p_signal",
                            "an3_80p_signal",
                            "an5_80p_signal",
                            "an7_80p_signal"
    )]),rep(n-i+1, 5),
    col=col, bg=plot_bg, cex = c(2.25,2.25,2,2,2),pch=c(2,6,3,4,1),
    lwd=c(3,3,3,3,3), xpd=TRUE)}
  # segments(x$ORACLE,(nrow(x):1)-0.1,x$ORACLE,(nrow(x):1)+0.1,
  #          col = "#000000",lwd=4)
  axis(1,at=pretty(seq(0,x_max,len=5)),pos=0,lwd=2)
  axis(1,at=seq(0,150,by=10),pos=0,lwd=1,labels = FALSE)

  y_labels <- paste0(x$ac_range,c("         ",
                                  "         ",
                                  "       ",
                                  "   ",
                                  "   ",
                                  " ",
                                  "   ",
                                  "   ",
                                  " "),
                     x$p_active)

  for(i in 1:nrow(x)){
    mtext(y_labels[i],2,at=n-i+1,las=1)
    #text(x = -25, y=n-i+1, labels = y_labels[i],las=1)
  }
  mtext("causal variant",2,at=n+1.6,las=1)
  mtext("type        #",2,at=n+1,las=1)
  segments(-100,n+1/2,-4,n+1/2,xpd=T)


  legend(0,-1.2,
         legend = c("SMT","LOCATER","AN3","AN5","AN7"),
         pt.cex = c(2.25,2.25,2,2,2),pch=c(2,6,3,4,1),
         col=col,
         pt.bg=plot_bg,pt.lwd=c(3,3,3,3,3),
         horiz = T,xpd=T,bty = "n",
         text.width = c(20,35,20,20,20))
}



x <- copy(x_key)
setkey(x,ac_range,p_active)
x[,ac_range:=case_match(ac_range,"all variants" ~ "any",
                        "doubletons only" ~ "doubletons",
                        "rare variants only" ~ "DAC [150,750)")]

cairo_pdf(paste0(plot.dir,"/power_summary_dotplot_hidden_vars_v2.pdf"),width = 7,height=5)
custom_dotplot(x,
               x_max = 150,
               col = colpal[1:5])
dev.off()




# SHOW SOURCE OF SIGNAL
############################

for(pa in unique(x_key$p_active)){

  cairo_pdf(paste0(plot.dir,"/source_",pa,"_active_vars_hidden_vars_v2.pdf"),width = 10,height=8)

  par(oma = c(1, 1, 1, 1), mar = c(2.5, 2.5, 0.5, 0.5))

  plot(0,0,type="n",ylim=c(0,1),xlim=c(0,150),las=1,bty="n",
       xlab="",ylab="",
       main="", xaxt="n",yaxt="n")
  axis(2,seq(0,1,by=0.2),pos=0,las=1,cex.axis=1.4)
  axis(1,seq(0,150,by=50),pos=0,las=1,cex.axis=1.4)
  axis(1,seq(0,150,by=10),pos=0,las=1,labels = FALSE)
  mtext("Total Association Signal Strength",side = 1,line = 2,cex=1.4)
  mtext("Proportion of LOCATER Signal",side = 2,line = 2,cex=1.4)

  for(j in 1:length(ac_names)){
    # lines(res[x_id==i, oracle[1], list(.id,rep,target_signal)][#.id!=38
    #   ,mean(V1 > cutoff),target_signal],col=colpal[i],lwd=3,lty=1)

    lines(res[p_active==pa & ac_range == ac_names[j] & targeted == TRUE, (smt/(smt+rd+qform))[which.max(tot)], list(.id,rep,target_signal)][,mean(V1),target_signal],
          col=colpal[1],lwd=3,lty=j)
    lines(res[p_active==pa & ac_range == ac_names[j] & targeted == TRUE, (rd/(smt+rd+qform))[which.max(tot)], list(.id,rep,target_signal)][,mean(V1),target_signal],
          col=colpal[2],lwd=3,lty=j)
    lines(res[p_active==pa & ac_range == ac_names[j] & targeted == TRUE, (qform/(smt+rd+qform))[which.max(tot)], list(.id,rep,target_signal)][,mean(V1),target_signal],
          col=colpal[3],lwd=3,lty=j)
  }
  legend(x = 75,y = 1,legend = expression(italic(SMT),italic(SD),italic(Q)),fill = colpal[1:3],border = NA,bty="n")
  legend(x = 100,y = 1,legend = c("all variants","DAC [150,750) variants only","doubletons only"),lty=1:3,col = colpal[9], lwd=3,border = NA,bty="n")
  dev.off()
}





# SHOW LOCALIZATION OF SIGNAL
########################################

# given the signal is strong enough that we have 80% power,
# how many kb large would a CI need to be to cover the middle of the gene 50% of the time?
# note, this is equivalent to 2 x the median distance to the center of the gene in kb.
# we can generalize this and ask how wide a 80% CI needs to be.
# Only look at points where the power is at least 50% so that these confidence intervals may be estimated with some precision.

for(pa in unique(x_key$p_active)){

  cairo_pdf(paste0(plot.dir,"/localization_",pa,"_active_vars_hidden_vars_v2.pdf"),width = 10,height=8)

  par(oma = c(1, 1, 1, 1), mar = c(2.5, 2.5, 0.5, 0.5))

  plot(0,0,type="n",ylim=c(0,1e3),xlim=c(0,150),las=1,bty="n",
       xlab="",ylab="",
       main="", xaxt="n",yaxt="n")
  axis(2,seq(0,1000,by=200),pos=0,las=1,cex.axis=1.4)
  axis(1,seq(0,150,by=50),pos=0,las=1,cex.axis=1.4)
  axis(1,seq(0,150,by=10),pos=0,las=1,labels = FALSE)
  mtext("Total Association Signal Strength",side = 1,line = 2,cex=1.4)
  mtext("80% CI Width (kb)",side = 2,line = 2,cex=1.4)


  for(j in 1:length(ac_names)){
    # lines(res[x_id==i, oracle[1], list(.id,rep,target_signal)][#.id!=38
    #   ,mean(V1 > cutoff),target_signal],col=colpal[i],lwd=3,lty=1)

    lines(res[p_active==pa & ac_range == ac_names[j] & smt >= cutoffs$smt,
              #2*min(abs(5e5-pos[smt==max(smt)]))/1e3,
              2*mean(abs(5e5-pos[smt==max(smt)]))/1e3,
              list(.id,rep,target_signal)][,quantile(V1,0.8),target_signal]
          [target_signal >= x_key[p_active==pa & ac_range == ac_names[j],smt_80p_signal],],
          col=colpal[1],lwd=2,lty=j)

    lines(res[p_active==pa & ac_range == ac_names[j] & targeted == TRUE & tot >= cutoffs$locater,
              #2*min(abs(5e5-pos[tot==max(tot)]))/1e3,
              2*mean(abs(5e5-pos[tot==max(tot)]))/1e3,
              list(.id,rep,target_signal)][,quantile(V1,0.8),target_signal]
          [target_signal >= x_key[p_active==pa & ac_range == ac_names[j],locater_80p_signal],],
          col=colpal[2],lwd=2,lty=j)


    lines(res_an[p_active==pa & ac_range == ac_names[j] & sampling_rate==1e-3 & p_value < 10^(-cutoffs$an3),
                 2*avg_abs_dist/1e3,
                 list(.id,rep,target_signal)]
          [,quantile(V1,0.8,na.rm=TRUE),target_signal]
          [target_signal >= x_key[p_active==pa & ac_range == ac_names[j],an3_80p_signal],],
          col=colpal[3],lwd=2,lty=j)

    lines(res_an[p_active==pa & ac_range == ac_names[j] & sampling_rate==1e-5 & p_value < 10^(-cutoffs$an5),
                 2*avg_abs_dist/1e3,
                 list(.id,rep,target_signal)]
          [,quantile(V1,0.8,na.rm=TRUE),target_signal]
          [target_signal >= x_key[p_active==pa & ac_range == ac_names[j],an5_80p_signal],],
          col=colpal[4],lwd=2,lty=j)

    lines(res_an[p_active==pa & ac_range == ac_names[j] & sampling_rate==1e-7 & p_value < 10^(-cutoffs$an7),
                 2*avg_abs_dist/1e3,
                 list(.id,rep,target_signal)]
          [,quantile(V1,0.8,na.rm=TRUE),target_signal]
          [target_signal >= x_key[p_active==pa & ac_range == ac_names[j],an5_80p_signal],],
          col=colpal[5],lwd=2,lty=j)

  }

  legend(x = 70,y = 1e3,legend = c("SMT","LOCATER","AN3","AN5","AN7"),fill=colpal[1:5],
         border = NA,bty="n")
  legend(x = 105,y = 1.05e3,y.intersp = 1.2,legend = c("all variants",expression(atop("DAC [150,750)","variants only")),"doubletons only"),
         lty=c(1,2,3),lwd=2, col=colpal[9], border = NA,bty="n",xpd=TRUE)
  dev.off()
}




