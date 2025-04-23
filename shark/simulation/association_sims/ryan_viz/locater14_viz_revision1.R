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



# LOAD DATA
##################################################################


# LOAD LOCATER SIMS
n.samps <- 900
samp.ind <- rep(FALSE,n.samps)
res.list <- as.list(1:n.samps)
for(i in 1:n.samps){
  target.file <- paste0("~/Dropbox/Hall_group/locater_sims/locater14/res_",i,".rds")

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


# LOAD STAAR SIMS
n.samps <- 900
samp.ind <- rep(FALSE,n.samps)
res.list <- as.list(1:n.samps)
for(i in 1:n.samps){
  target.file <- paste0("~/Dropbox/Hall_group/locater_sims/staar14/res_",i,".rds")

  if(!file.exists(target.file)){next}
  tryCatch({
    res.list[[i]] <- readRDS(target.file)
    samp.ind[i] <- TRUE
  },error=function(e){print(e)})


  print(i)
}

res_staar <- rbindlist(res.list[samp.ind],idcol = TRUE)
rm(res.list);gc()
# quality control checks
any(!is.finite(res_staar$ACAT_O_no_cutoff))



# LOAD ARG-NEEDLE SIMS
n.samps <- 900
samp.ind <- rep(FALSE,n.samps)
res.list <- as.list(1:n.samps)
for(i in 1:n.samps){
  target.file <- paste0("~/Dropbox/Hall_group/locater_sims/argneedle14/res_",i,".rds")

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
            "p_active" = c(3,9,15),
            "gene_size" = c(1e4), # c(5e3,2e4,5e4)
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

res_staar <- merge(res_staar,x_key,by="x_id")
res_staar <- merge(res_staar,obs_key,by="obs_id")
res_staar <- merge(res_staar,y_key,by="y_id")

res_an <- merge(res_an,x_key,by="x_id")
res_an <- merge(res_an,obs_key,by="obs_id")
res_an <- merge(res_an,y_key,by="y_id")

#recode
res[,ac_range:= case_match(ac_range,
                           "any" ~ "all variants",
                           "150-749" ~ "rare variants only",
                           "2" ~ "doubletons only")]

res_staar[,ac_range:= case_match(ac_range,
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


# Quality Check Results
###############################

# check STAAR sim null distribution

qqexp <- function(x,...){
  plot(qexp(ppoints(length(x)),rate=log(10)),sort(x),
       xlab="Expected Quantiles",
       ylab="Observed Quantiles",
       bty="n",las=1,col="navy",...)
  abline(0,1,col="darkorange",lwd=1.5)
}

qqexp(-log10(runif(1e4)))
qqexp(res_staar[target_signal==0,ACAT_O_no_cutoff])


#
# # Look back to understand why we could be getting a max LOCATER signal of 6.7 over a chromosome
# # if we indeed have target signal of 150
# res[target_signal==150 & p_active==9 & ac_range == "all", max(smt), list(.id,rep)][,summary(V1)]
# temp_res <- res[target_signal==150 & p_active==9 & ac_range == "all" & targeted == TRUE, max(tot), list(.id,rep)]
# head(temp_res[order(temp_res$V1),],20)
# # looks like we need to investigate the details of simulation 725
# td <- readRDS("~/Dropbox/Hall_group/locater_sims/locater14/details_725.rds")
#
# res[target_signal==150 & p_active==9 & ac_range == "all" & .id==725, max(smt), list(rep)]
#
# res[target_signal==150 & p_active==9 & ac_range == "all" & .id==725 & rep==1 & targeted==TRUE,]
#
#
#
# # see sims that are yielding 0 for the oracle
# res[target_signal==150,oracle[1],list(.id,rep)][V1==0,]
# # no discernable pattern:
# unique(res[target_signal==150,oracle[1],list(.id,rep)][V1==0,]$.id)%%9


# Calculate cutoffs
###########################

# corresponds to testing entire genome with 10kb sliding windows with 5kb overlap
staar_discovery_cutoff <- -log10(0.05 / (3e9/5e3)) #2.3*10^-7 # 7 # 0.05 / (3e9/5e3) = 8.33e-8 < 1e-7

res_an_cutoffs <- res_an[target_signal==0,-log10(quantile(p_value,0.05)/3e3),sampling_rate]

# Number of simulations where LOCATER did not test any variants because none of the
# variants achieved the 10^-4 threshold
sum(!is.finite(res[target_signal==0,suppressWarnings(max(tot[targeted == TRUE])),list(.id,rep,x_id)]$V1))

cutoffs <- list(
  "smt" = unname(res[target_signal==0,max(smt),list(.id,rep,x_id)][,quantile(V1,0.95)+log10(3e3)]),
  # Out of our 9000 independent null simulations (target_signal == 0), close to half of them had no LOCATER signal tested (all variants have targeted==FALSE)
  # because none of them met the SMT threshold of 10^-4. Since these corresponding to genomic regions that were not tested by LOCATER, we can just assign
  # a max locater signal to each of these simulations of 0.  This slightly reduces the genome-wide LOCATER cutoff (from about 8.6 to about 8.4).
  # Since max(tot[targeted == TRUE]) returns -Inf and a warning when tot[targeted == TRUE] has length 0, we supress warnings and remove infinite results in the line below
  "locater" = unname(res[target_signal==0,suppressWarnings(max(tot[targeted == TRUE])),list(.id,rep,x_id)][,quantile(ifelse(is.finite(V1),V1,0),0.95)+log10(3e3)]),
  "an3" = res_an_cutoffs[sampling_rate==1e-3,V1],
  "an5" = res_an_cutoffs[sampling_rate==1e-5,V1],
  "an7" = res_an_cutoffs[sampling_rate==1e-7,V1]
)

cutoffs$acat_all <- res_staar[target_signal==0,quantile(ACAT_O_no_cutoff,0.95)+log10(3e9/5e3)]
cutoffs$acat_rare <- res_staar[target_signal==0,quantile(ACAT_O_with_cutoff,0.95)+log10(3e9/5e3)]

xtable::xtable(as.data.frame(cutoffs))


# Calculate Power Curves with SMT, LOCATER and ARG-NEEDLE
###########################################################

for(pa in unique(x_key$p_active)){

  cairo_pdf(paste0(plot.dir,"/power_",pa,"_active_vars_observed_vars_v2.pdf"),width = 10,height=8)

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
}



# Calculate Power Curves with STAAR, LOCATER, and SMT
#########################################################
#
# for(pa in unique(x_key$p_active)){
#
#   cairo_pdf(paste0(plot.dir,"/power_",pa,"_active_vars_observed_vars_with_acat_v2.pdf"),width = 10,height=8)
#
#   par(oma = c(1, 1, 1, 1), mar = c(2.5, 2.5, 0.5, 0.5))
#
#   plot(0,0,type="n",ylim=c(0,1),xlim=c(0,150),las=1,bty="n",
#        xlab="",ylab="",
#        main="", xaxt="n",yaxt="n")
#   axis(2,seq(0,1,by=0.2),pos=0,las=1,cex.axis=1.4)
#   axis(1,seq(0,150,by=30),pos=0,las=1,cex.axis=1.4)
#   mtext("Signal Strength",side = 1,line = 2,cex=1.4)
#   mtext("Power",side = 2,line = 2,cex=1.4)
#
#   for(j in 1:length(ac_names)){
#     # lines(res[x_id==i, oracle[1], list(.id,rep,target_signal)][#.id!=38
#     #   ,mean(V1 > cutoff),target_signal],col=colpal[i],lwd=3,lty=1)
#     lines(res[p_active==pa & ac_range == ac_names[j], max(smt), list(.id,rep,target_signal)][#.id!=38
#       ,mean(V1 > cutoffs$smt),target_signal],col=colpal[1],lwd=2,lty=j)
#     lines(res[p_active==pa & ac_range == ac_names[j] & targeted == TRUE, max(tot), list(.id,rep,target_signal)][#.id!=38
#       ,mean(V1 > cutoffs$locater),target_signal],col=colpal[2],lwd=2,lty=j)
#
#
#     lines(res_staar[p_active==pa & ac_range == ac_names[j], ACAT_O_with_cutoff, list(.id,rep,target_signal)][
#       ,mean(ACAT_O_with_cutoff > cutoffs$acat_rare),target_signal],
#       col=colpal[6],lwd=2,lty=j)
#
#     lines(res_staar[p_active==pa & ac_range == ac_names[j], ACAT_O_no_cutoff, list(.id,rep,target_signal)][
#       ,mean(ACAT_O_no_cutoff > cutoffs$acat_all),target_signal],
#       col=colpal[8],lwd=2,lty=j)
#
#   }
#   legend(x = 125,y = 0.75,legend = c("SMT","LOCATER","ACAT-O (rare vars)","ACAT-O (all vars)"),fill=colpal[c(1,2,6,8)],
#          border = NA,bty="n")
#   legend(x = 122,y = 0.55,y.intersp = 2,legend = c("all variants",expression(atop("intermediate","variants only")),"doubletons only"), lty=c(1,2,3),lwd=2, col=colpal[9], border = NA,bty="n")
#   dev.off()
# }


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


    x_key[p_active==pa & ac_range == ac_names[j], "acat_rare_80p_signal"] <- res_staar[
      p_active==pa & ac_range == ac_names[j], ACAT_O_with_cutoff, list(.id,rep,target_signal)
    ][
      ,mean(ACAT_O_with_cutoff > cutoffs$acat_rare),target_signal
    ][
      ,approx(V1,target_signal,xout = 0.8,method = "linear", ties = "ordered")$y
    ]

    x_key[p_active==pa & ac_range == ac_names[j], "acat_all_80p_signal"] <- res_staar[
      p_active==pa & ac_range == ac_names[j], ACAT_O_no_cutoff, list(.id,rep,target_signal)
    ][
      ,mean(ACAT_O_no_cutoff > cutoffs$acat_all),target_signal
    ][
      ,approx(V1,target_signal,xout = 0.8,method = "linear", ties = "ordered")$y
    ]

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

# x_key[,locater_gain:=pmax(0,smt_80p_signal-locater_80p_signal)]
# x_key[,smt_gain:=pmax(0,locater_80p_signal-smt_80p_signal)]
#
# cairo_pdf(paste0(plot.dir,"/power_summary_barplot_observed_vars.pdf"),width = 8,height=8)
# make_locater_summary_boxplot(
#   x_key[c(3,6,9,2,5,8,1,4,7),list(locater_80p_signal,locater_gain,smt_gain)],
#   "Signal Strength Required to Achieve 80% Power")
# dev.off()




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

cairo_pdf(paste0(plot.dir,"/power_summary_dotplot_observed_vars_v2.pdf"),width = 7,height=5)
custom_dotplot(x,
               x_max = 150,
               col = colpal[1:5])
dev.off()



# MAKE POWER DOTPLOT COMPARING TO STAAR (JUST IN THE ALL VARIANT CASE)
#######################################################################

custom_dotplot2 <- function(x, x_max, col = "#000000"){

  # x_max = 120
  # col = colpal[1:4]

  n <- 3
  p <- ncol(x)

  # color blind palette as recommended: https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible

  plot_bg <- paste0(col,c("00","00","FF","00"))

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
  segments(rep(0,3),1:3,rep(x_max,3),1:3,
           col = gray.colors(3,0.7,0.9),lwd=2,lty=1)

  for(i in 1:3){
    points(as.numeric(x[i,c("smt_80p_signal","locater_80p_signal",
                            "acat_rare_80p_signal",
                            "acat_all_80p_signal")]),rep(n-i+1, 4),
           col=col, bg=plot_bg, cex = c(2.25,2.25,2,2),pch=c(2,6,23,5),
           lwd=c(3,3,3,3), xpd=TRUE)}
  # segments(x$ORACLE,(nrow(x):1)-0.1,x$ORACLE,(nrow(x):1)+0.1,
  #          col = "#000000",lwd=4)
  axis(1,at=pretty(seq(0,x_max,len=5)),pos=0,lwd=2)

  y_labels <- paste0(x$ac_range[1:3],c("         ",
                                  "         ",
                                  "       "),
                     x$p_active[1:3])

  for(i in 1:3){
    mtext(y_labels[i],2,at=n-i+1,las=1)
    #text(x = -25, y=n-i+1, labels = y_labels[i],las=1)
  }
  mtext("causal variant",2,at=n+1.6,las=1)
  mtext("type        #",2,at=n+1,las=1)
  segments(-100,n+1/2,-4,n+1/2,xpd=T)


  legend(-1,-1.2,
         legend = c("SMT","LOCATER","ACAT-O (rare)","ACAT-O (all)"),
         col=col,
         pt.bg=plot_bg,pt.cex = c(2.25,2.25,2,2),pch=c(2,6,23,5),
         pt.lwd=c(3,3,3,3),
         horiz = T,xpd=T,bty = "n",
         text.width = c(13,24,30,28))
}



x <- copy(x_key)
setkey(x,ac_range,p_active)
x[,ac_range:=case_match(ac_range,"all variants" ~ "any",
                        "doubletons only" ~ "doubletons",
                        "rare variants only" ~ "DAC [150,750)")]

cairo_pdf(paste0(plot.dir,"/power_summary_dotplot_observed_vars_with_acat_v2.pdf"),width = 7,height=5/2)
custom_dotplot2(x,
               x_max = 120,
               col = colpal[c(1,2,6,8)])
dev.off()

#
# custom_dotplot3 <- function(x, x_max, col = "#000000"){
#
#   n <- nrow(x)
#   p <- ncol(x)
#
#   # color blind palette as recommended: https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
#
#   plot_bg <- rev(paste0(col,c("00","AA")))
#
#   op <- par(oma = c(0,1,0,0),
#             mar = c(4.1, 6.1, 3.1, 1.1)
#             #mar = c(4.1, 5.1, 2.1, 6.1)
#   )
#   on.exit(par(op))
#   #"#000000" "#D55E00" "#56B4E9" "#F0E442" "#009E73" "#CC79A7" "#E69F00"
#
#   plot(0, type="n",
#        xlim=c(0,ceiling(max(x_max))),ylim=c(0,n),
#        bty="n",ylab="",xlab="",xaxt="n",yaxt="n"
#   )
#   segments(rep(0,nrow(x)),1:nrow(x),rep(x_max,nrow(x)),1:nrow(x),
#            col = gray.colors(3,0.7,0.9),lwd=2,lty=1)
#   for(i in 1:nrow(x)){
#     points(as.numeric(x[i,c("smt_80p_signal","locater_80p_signal",
#                             "A_half_half_no_cutoff_80p_signal",
#                             "S_half_half_no_cutoff_80p_signal")]),rep(n-i+1, 4),
#            col=col, bg=plot_bg, cex = 2,pch=23,
#            lwd=c(3,0,0,0),xpd=T)}
#   # segments(x$ORACLE,(nrow(x):1)-0.1,x$ORACLE,(nrow(x):1)+0.1,
#   #          col = "#000000",lwd=4)
#   axis(1,at=pretty(seq(0,x_max,len=5)),pos=0,lwd=2)
#
#   y_labels <- paste0(x$ac_range,c("         ",
#                                   "         ",
#                                   "       ",
#                                   "   ",
#                                   "   ",
#                                   " ",
#                                   "   ",
#                                   "   ",
#                                   " "),
#                      x$p_active)
#
#   for(i in 1:nrow(x)){
#     mtext(y_labels[i],2,at=n-i+1,las=1)
#     #text(x = -25, y=n-i+1, labels = y_labels[i],las=1)
#   }
#   mtext("causal variant",2,at=n+1.6,las=1)
#   mtext("type        #",2,at=n+1,las=1)
#   segments(-100,n+1/2,-4,n+1/2,xpd=T)
#
#
#   legend(35,-1.2,
#          legend = c("SMT","LOCATER","ACAT-V","SKAT"),
#          pch=23,
#          col=col,
#          pt.bg=plot_bg,pt.cex=2,pt.lwd=c(3,0),
#          horiz = T,xpd=T,bty = "n")
# }
#
#
# x <- copy(x_key)
# setkey(x,ac_range,p_active)
# x[,ac_range:=case_match(ac_range,"all variants" ~ "any",
#                         "doubletons only" ~ "doubletons",
#                         "rare variants only" ~ "intermediate")]
#
# dotchart(t(as.matrix(x[,list(S_half_half_no_cutoff_80p_signal,A_half_half_no_cutoff_80p_signal,smt_80p_signal,locater_80p_signal)])))
#
# cairo_pdf(paste0(plot.dir,"/power_summary_dotplot_observed_vars_with_skat_v2.pdf"),width = 7,height=5/2)
# custom_dotplot3(x,
#                 x_max = 120,
#                 col = colpal[1:4])
# dev.off()
#





# SHOW SOURCE OF SIGNAL
############################

for(pa in unique(x_key$p_active)){

  cairo_pdf(paste0(plot.dir,"/source_",pa,"_active_vars_observed_vars_v2.pdf"),width = 10,height=8)

  par(oma = c(1, 1, 1, 1), mar = c(2.5, 2.5, 0.5, 0.5))

  plot(0,0,type="n",ylim=c(0,1),xlim=c(0,150),las=1,bty="n",
       xlab="",ylab="",
       main="", xaxt="n",yaxt="n")
  axis(2,seq(0,1,by=0.2),pos=0,las=1,cex.axis=1.4)
  axis(1,seq(0,150,by=30),pos=0,las=1,cex.axis=1.4)
  mtext("Signal Strength",side = 1,line = 2,cex=1.4)
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

  cairo_pdf(paste0(plot.dir,"/localization_",pa,"_active_vars_observed_vars_v2.pdf"),width = 10,height=8)

  par(oma = c(1, 1, 1, 1), mar = c(2.5, 2.5, 0.5, 0.5))

  plot(0,0,type="n",ylim=c(0,1e3),xlim=c(0,150),las=1,bty="n",
       xlab="",ylab="",
       main="", xaxt="n",yaxt="n")
  axis(2,seq(0,1000,by=200),pos=0,las=1,cex.axis=1.4)
  axis(1,seq(0,150,by=30),pos=0,las=1,cex.axis=1.4)
  mtext("Signal Strength",side = 1,line = 2,cex=1.4)
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


# x_key[,min_localization:=pmin(locater_localization,smt_localization)]
# x_key[,locater_localization_gain:=pmax(0,locater_localization-smt_localization)]
# x_key[,smt_localization_gain:=pmax(0,smt_localization-locater_localization)]
#
#
# cairo_pdf(paste0(plot.dir,"/localization_summary_barplot_observed_vars.pdf"),width = 8,height=8)
# make_locater_summary_boxplot(
#   x_key[c(3,6,9,2,5,8,1,4,7),list(min_localization,locater_localization_gain,smt_localization_gain)],
#   "80% CI Width (kb)")
# dev.off()

#####
# Save results for processing alongside locater16 results in locater16_viz.R
################################################################################

saveRDS(res[p_active==9,],"~/Dropbox/Hall_group/locater_sims/locater14_9_active_vars_res.rds")
saveRDS(cutoffs,"~/Dropbox/Hall_group/locater_sims/locater14_derived_genome_wide_cutoffs.rds")





