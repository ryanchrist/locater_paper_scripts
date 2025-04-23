
# Declare Plotting Resources
###############################
# first 8 colors taken from color blind palette proposed by
# Wong, B. Points of view: Color blindness. Nat Methods 8, 441 (2011). https://doi.org/10.1038/nmeth.1618
# with gray "#AAAAAA as our added 9th color

# original:
#colpal <- c("#000000", "#D55E00","#56B4E9","#F0E442","#009E73","#CC79A7","#E69F00","#0072B2","#AAAAAA")

#swapped:
colpal <- c("#000000", "#D55E00","#56B4E9","#CC79A7","#009E73","#F0E442","#E69F00","#0072B2","#AAAAAA")


# "#FFE442" "#F0E442"
# "#FFFF00" "#F0E442"
# "#B59410"

# mean(res$p_chr) =  14848
# Bonferroni Threshold for this experiment:
#0.05/(3300*14848)
# so here we go with 10^-9 cutoff (somewhat typical for modern large sequencing studies)
cutoff <- 8.5

#ac_names <- c("all","common","doubleton")
ac_names <- c("all variants",
              "rare variants only", # "150-749"
              "doubletons only")



make_locater_summary_boxplot <- function(x,ylab=ylab){

  left_scale <- pretty(c(0,rowSums(x)))
  scale_max <- max(left_scale)

  barplot(t(x * 120/scale_max),
          col=c("#009E73","#F0E442","#0072B2"),
          width=0.8,
          space = c(0.25),
          las=1,border = NA,
          ylab=ylab,cex.lab=1.2,yaxt="n",ylim=c(0,120))

  axis(2, at = left_scale*120/scale_max, labels = left_scale,las=1)

  text(x = seq(0.8/2+0.8*(0.25),by=1,length.out=9),
       y = -4,
       labels = rep(c(3,9,15),3),
       xpd=TRUE)
  text(-0.5,-4,labels = "# causal:",xpd=TRUE)

  segments(x0 = 0.2, y0 = -10, x1 = 3, y1 = -10,xpd=TRUE)
  segments(x0 = 0.2, y0 = -9, x1 = 0.2, y1 = -11,xpd=TRUE)
  segments(x0 = 3, y0 = -9, x1 = 3, y1 = -11,xpd=TRUE)

  segments(x0 = 3.2, y0 = -10, x1 = 6, y1 = -10,xpd=TRUE)
  segments(x0 = 3.2, y0 = -9, x1 = 3.2, y1 = -11,xpd=TRUE)
  segments(x0 = 6, y0 = -9, x1 = 6, y1 = -11,xpd=TRUE)

  segments(x0 = 6.2, y0 = -10, x1 = 9, y1 = -10,xpd=TRUE)
  segments(x0 = 6.2, y0 = -9, x1 = 6.2, y1 = -11,xpd=TRUE)
  segments(x0 = 9, y0 = -9, x1 = 9, y1 = -11,xpd=TRUE)

  text(x = c(1.6,4.6,7.6),
       y = -12.5,
       labels = c("all", "doubletons only", "rare variants only"), #"DAF in [0.0025,0.0125) only"),
       xpd=TRUE)
  text(-0.5,-9.5,labels = " causal",xpd=TRUE)
  text(-0.5,-12.5,labels = "variants:",xpd=TRUE)

  legend("topright",legend=c("SMT","LOCATER","OVERLAP"),fill=c("#F0E442","#0072B2","#009E73"),bty="n",
         border=c("#F0E442","#0072B2","#009E73"))
}
