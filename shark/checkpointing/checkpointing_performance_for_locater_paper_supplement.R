# Demonstration of Optimal Checkpointing Efficiency
# For Supplement of Clade Distillation paper
# partly dervied from checkpointing_trial1.R

# load libraries
require(RColorBrewer)
require(kalis) # must be kalis 2.0

# Declare output plot directory
plot_dir <- "~/Desktop/"



#  BELOW HERE SHOULD NOT NEED TO BE MODIFIED FOR REPLICABILITY
##################################################################

# The following call calculates the cost table of using optimal
# checkpointing for 1 million equally spaced variants
# using up to 10 checkpoints
# depending on your system, this may take about half an hour

max.num.checkpoints <- 10
cost.table <- kalis:::calc_tables(1:1e6,max.num.checkpoints)$cost


# Make checkpointing cost figure for supplement
cairo_pdf(paste0(plot_dir,"uniform_ckpt_cost_as_function_of_K.pdf"))
par(mar=c(5,5,2,1))
cols <- brewer.pal(3,"Set1")
l <- c(1e6,1e5,1e4)
#l <- c(1e4)
for(i in 1:length(l)){
  if(i == 1){
    plot(0:max.num.checkpoints,log10(cost.table[l[i],])-log10(l[i]),col = cols[i], type="l",ylab = "log10 Distance - log10 L", xlab = "Number of Checkpoints Used",
         bty="n",ylim=c(0,6),xlim = c(0,10),las=1)
  }else{
    lines(0:max.num.checkpoints,log10(cost.table[l[i],])-log10(l[i]),col = cols[i])
  }
  legend("topright",legend = c("1,000,000","100,000","10,000"),fill=cols)
  abline(h = 1,lty=2,col="lightgray")
  abline(h = 2,lty=2,col="lightgray")
}
dev.off()


# For validation, the much slower pure R implementation
# of the kalis:::calc_tables function above is given
# commented out below
#####################################

#
# start <- proc.time()
# max.n <- 1e6
# max.num.checkpoints <- 10
#
# s <- 1:max.n
#
# # the first row corresponds to solving a 0 locus problem
# cost.table <- matrix(0L,nrow=max.n + 1,ncol= max.num.checkpoints + 1)
# index.table <- matrix(0L,nrow=max.n,ncol= max.num.checkpoints)
#
# #rev.cost.table <- cost.table
# cost.table[,1] <- c(0,s*(s+1)/2)
# #rev.cost.table[,1] <- rev(cost.table[,1])
#
#
# for(k in 1:max.num.checkpoints){
#   for(n in 1:max.n){
#     # now solving a n long problem with k checkpoints
#     v <- cost.table[1:n,k + 1] + s[1:n] + cost.table[n:1,k]
#     x <- which.min(v)
#     index.table[n, k] <- x
#     cost.table[n + 1,k + 1] <- v[x]
#   }
#   print(paste(k,"done at",c(proc.time() - start)[3]/3600,"hours from start."))
# }
#
# cost.table <- cost.table[-1,]
# finish <- proc.time() - start


