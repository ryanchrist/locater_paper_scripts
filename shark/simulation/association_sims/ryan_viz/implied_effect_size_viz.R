# Displaying Effect Sizes Underlying Power Simulations
# This script produces two plots and a table for the supplement
# of our LOCATER manuscript providing estimates of the effect sizes
# used in our power simulations.

plot.dir <- "~/Desktop/"

f <- function(x,df = 3){
  sqrt(qchisq(-log(10)*x,df = df,lower.tail = FALSE,log.p = TRUE))/sqrt(df)
}

# verified as the inverse!
-pchisq(q = 9*(f(1:50,9)^2),df = 9,lower.tail = FALSE,log.p = TRUE)/log(10)


cairo_pdf(paste0(plot.dir,"rescaled_effect_size.pdf"),width = 7,height = 7)
xx <- seq(0,150,len=1e3)
par(oma=double(4),mar=c(3,3,0.5,0.5))
plot(xx,f(xx),bty="n",las=1,type="l",xlab="",
     ylab="",lty=3,xaxt="n",yaxt="n")
lines(xx,f(xx,9),bty="n",las=1,lty=2)
lines(xx,f(xx,15),bty="n",las=1,lty=1)
legend(5,15,legend=c(3,9,15),lty=c(3,2,1),
       bty="n",title = "# Causal Variants")
axis(1,at=pretty(xx),pos=0,lwd=2)
axis(1,at=seq(0,150,by=10),pos=0,lwd=1,labels = FALSE)
axis(2,at=seq(0,15,by=5),pos=0,lwd=2,las=1)
mtext(expression(paste("Standardized Effect Size ",symbol(beta),"*(s,a)")),2,1.5,cex=1.2)
mtext("Total Association Signal Strength (s)",1,1.5,cex=1.2)
dev.off()

g <- function(x,N = 2e4){
  p <- x/N
  1/sqrt(N*p*(1-p)*(1-1/N-p/(N-1)))
}


cairo_pdf(paste0(plot.dir,"inverse_length.pdf"),width = 7,height = 7)
xxx <- 1:200
par(oma=double(4),mar=c(3,3,0.5,0.5))
plot(xxx,g(xxx),bty="n",las=1,xlab="",
     ylab="",lty=3,xaxt="n",yaxt="n",
     ylim=c(0,1),xlim=c(0,200),col="darkblue")
points(c(2,20,200),g(c(2,20,200)),
       col="darkorange",pch=20,cex=1.5)

legend(170,.98,legend=c("Yes","No"),pch=c(20,1),col=c("darkorange","darkblue"),
       pt.cex = c(1.5,1),bty="n",title = "In Table")

axis(1,at=pretty(xxx),pos=0,lwd=2)
axis(1,at=seq(0,200,by=10),pos=0,lwd=1,labels = FALSE)
axis(2,at=seq(0,1,by=0.2),pos=0,lwd=2,las=1)
mtext("Inverse Length",2,1.5,cex=1.2)
mtext("Minor Allele Count (m)",1,1.5,cex=1.2)
dev.off()


# generate table

xx <- seq(10,130,by=10)
m <- cbind(f(xx,df=3),f(xx,df=9),f(xx,df=15))
m <- cbind(m*g(2),m*g(20),m*g(200))
rownames(m) <- as.character(xx)
colnames(m) <- c(t(outer(c(2,20,200),c(3,9,15),FUN = paste,sep="_")))
xtable::xtable(t(m))


# Check claim about projection
n <- 3000
A <- cbind(1,
           c(double(n/3),rep_len(1,n/3),double(n/3)),
           c(double(n/3),double(n/3),rep_len(1,n/3)))
Q <- qr.Q(qr(A))

x <- c(rep_len(1,324),double(n-324))

xp <- x - Q%*%crossprod(Q,x)
plot(xp)













