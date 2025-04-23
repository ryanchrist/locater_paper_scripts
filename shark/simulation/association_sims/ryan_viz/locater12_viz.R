require(locater)
require(data.table)
require(dplyr)
require(stringr)
shark.dir <- "~/Dropbox/shark/"
plot.dir <- "~/Dropbox/Hall_group/locater_sims/plots/locater12/"
if(!dir.exists(plot.dir)){dir.create(plot.dir,recursive = TRUE)}

#  BELOW HERE SHOULD NOT NEED TO BE MODIFIED FOR REPLICABILITY
##################################################################

n.samps <- 1000
samp.ind <- rep(FALSE,n.samps)
res.list <- details.list<- as.list(1:n.samps)

for(i in 1:n.samps){
  target.file <- paste0("~/Dropbox/Hall_group/locater_sims/locater12/res_",i,".rds")
  target.file2 <- paste0("~/Dropbox/Hall_group/locater_sims/locater12/details_",i,".rds")

  if(!file.exists(target.file)){next}
  tryCatch({
    res.list[[i]] <- readRDS(target.file)
    details.list[[i]] <- readRDS(target.file2)
    samp.ind[i] <- TRUE
  },error=function(e){print(e)})
}

res.list <- res.list[samp.ind]
res <- rbindlist(res.list,idcol = TRUE)
details.list <- details.list[samp.ind]
rm(res.list);gc()

# Were all replicates retreived from directory?
all(samp.ind)

# Declare helper functions
check_tail <- function(x,cutoff = 0.001,p = NULL){

  sx <- sum(x<cutoff)
  binom_signal <-max(
    -pbinom(sx,size = if(is.null(p)){length(x)}else{p},prob = cutoff,lower.tail = TRUE, log.p = TRUE)/log(10),
    -pbinom(sx,size = if(is.null(p)){length(x)}else{p},prob = cutoff,lower.tail = FALSE,log.p = TRUE)/log(10)
  ) - log(2)/log(10)

  z <- qnorm(x[x<cutoff]/cutoff)

  print(round(c(
    binom_signal,
    -log10(shapiro.test(z)$p.value),
    -log10(t.test(z)$p.value)
  ),digits=4))
  qqnorm(z)
  qqline(z)
}

qq <- function(x,cutoff = NULL,...){
  xx <- qexp(ppoints(length(x)),log(10))
  yy <- sort(x)
  if(!is.null(cutoff)){
    q <- floor(length(x)*cutoff)
    xx <- tail(xx,q)
    yy <- tail(yy,q)
  }
  plot(xx,yy,bty="n",las=1,
       xlab="Expected Quantiles",ylab="Observed Quantiles",...);abline(0,1)}


# LOOK AT ALL SIMS TOGETHER
cutoff <- 0.001
check_tail(10^(-res$smt),cutoff = cutoff)
check_tail(10^(-res$rd),cutoff = cutoff)

check_tail(10^(-res$qform),cutoff = cutoff,p = length(res$qform))
check_tail(10^(-res$tot),cutoff = cutoff)

plot_idx <- c(seq(1,1e6-1e4,by=1e2),989902:1e6)

cairo_pdf(paste0(plot.dir,"/locater_smt_subtest_1million_null_sims_qqplot.pdf"),width = 10,height=8)

par(oma = c(1, 1, 1, 1), mar = c(2.5, 2.5, 0.5, 0.5))

plot(0,0,type="n",ylim=c(0,7),xlim=c(0,7),las=1,bty="n",
     xlab="",ylab="",
     main="", xaxt="n",yaxt="n")
axis(2,seq(0,7,by=1),pos=0,las=1,cex.axis=1.4)
axis(1,seq(0,7,by=1),pos=0,las=1,cex.axis=1.4)
mtext("Expected Quantiles",side = 1,line = 2,cex=1.4)
mtext("Observed Quantiles",side = 2,line = 2,cex=1.4)

points(qexp(ppoints(1e6)[plot_idx],rate = log(10)),sort(res$smt)[plot_idx])
abline(0,1,col="red")

dev.off()



cairo_pdf(paste0(plot.dir,"/locater_sd_subtest_1million_null_sims_qqplot.pdf"),width = 10,height=8)

par(oma = c(1, 1, 1, 1), mar = c(2.5, 2.5, 0.5, 0.5))

plot(0,0,type="n",ylim=c(0,7),xlim=c(0,7),las=1,bty="n",
     xlab="",ylab="",
     main="", xaxt="n",yaxt="n")
axis(2,seq(0,7,by=1),pos=0,las=1,cex.axis=1.4)
axis(1,seq(0,7,by=1),pos=0,las=1,cex.axis=1.4)
mtext("Expected Quantiles",side = 1,line = 2,cex=1.4)
mtext("Observed Quantiles",side = 2,line = 2,cex=1.4)

points(qexp(ppoints(1e6)[plot_idx],rate = log(10)),sort(res$rd)[plot_idx])
abline(0,1,col="red")

dev.off()



cairo_pdf(paste0(plot.dir,"/locater_qform_subtest_1million_null_sims_qqplot.pdf"),width = 10,height=8)

par(oma = c(1, 1, 1, 1), mar = c(2.5, 2.5, 0.5, 0.5))

plot(0,0,type="n",ylim=c(0,7),xlim=c(0,7),las=1,bty="n",
     xlab="",ylab="",
     main="", xaxt="n",yaxt="n")
axis(2,seq(0,7,by=1),pos=0,las=1,cex.axis=1.4)
axis(1,seq(0,7,by=1),pos=0,las=1,cex.axis=1.4)
mtext("Expected Quantiles",side = 1,line = 2,cex=1.4)
mtext("Observed Quantiles",side = 2,line = 2,cex=1.4)

points(qexp(ppoints(1e6)[plot_idx],rate = log(10)),sort(res$qform)[plot_idx])
abline(0,1,col="red")

dev.off()


cairo_pdf(paste0(plot.dir,"/locater_combined_subtest_1million_null_sims_qqplot.pdf"),width = 10,height=8)

par(oma = c(1, 1, 1, 1), mar = c(2.5, 2.5, 0.5, 0.5))

plot(0,0,type="n",ylim=c(0,7),xlim=c(0,7),las=1,bty="n",
     xlab="",ylab="",
     main="", xaxt="n",yaxt="n")
axis(2,seq(0,7,by=1),pos=0,las=1,cex.axis=1.4)
axis(1,seq(0,7,by=1),pos=0,las=1,cex.axis=1.4)
mtext("Expected Quantiles",side = 1,line = 2,cex=1.4)
mtext("Observed Quantiles",side = 2,line = 2,cex=1.4)

points(qexp(ppoints(1e6)[plot_idx],rate = log(10)),sort(res$tot)[plot_idx])
abline(0,1,col="red")

dev.off()


plot(qexp(ppoints(1000),rate = log(10)),quantile(res$tot,probs = ppoints(1000),na.rm=T),
     main="LOCATER (combined)",
     bty="n",las=1, xlab="Expected Quantiles",ylab="Observed Quantiles"); abline(0,1)



cutoff <- 0.01
qq(res$smt,cutoff = cutoff, main="SMT",ylim=c(0,6),xlim=c(0,6))
qq(res$rd,cutoff = cutoff,main = "Stable Distillation",ylim=c(0,6),xlim=c(0,6))
qq(res$qform,cutoff = cutoff, main = "QForm",ylim=c(0,6),xlim=c(0,6))
qq(res$tot,cutoff=cutoff, main = "LOCATER (Combined)",ylim=c(0,6),xlim=c(0,6))
