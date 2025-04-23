PlotManhattan <- function(x,y,...){
  
  # x is a vector of position / x-axis coordinates
  # y is a vector or a list, each element a vector the same length as x
  # if y is a named list, it's names are used for a legend
  # additional arguments can be used to add an x-axis label, a title, etc.
  
  if(!is.list(y)){ y <- list(y)}
  if(!all(sapply(y,is.numeric))){
    stop("y must be a  numeric vector of a list of numeric vectors")
  }
  
  if(!all(sapply(y,length)==length(x))){
    stop("length of x does not match length of entries in y")
  }
  
  
  if(length(y) > length(colors())){stop(paste("y cannot be a list with more than",length(colors()),"entries : there are not enough colors() in base R."))}
  
  custom.palette <- c("black","blue2","darkorange2","chartreuse4","brown2","darkgoldenrod2", "darkmagenta","coral4","cornsilk4","aquamarine3")
  
  if(length(y) <= length(custom.palette)){
    current.palette <- custom.palette[1:length(y)]
  }else{
    current.palette <- c(custom.palette,sample(colors()[!(colors() %in% custom.palette)],size = length(y) - length(custom.palette),replace=FALSE))
  }
  
  
  y.lim = c(0,max(unlist(y),na.rm = T))
  
  yy <- pretty(seq(y.lim[1],y.lim[2],len=5))
  xx <- pretty(x)
  plot(x,y[[1]],ylab=expression(paste(-log[10],"  p-value")), col = current.palette[1], bty="n",las=1,yaxt="n",xaxt="n",xlim=range(xx),ylim=c(min(yy)-0.1,max(yy)),...)
  axis(2,at=yy, las=1, pos = (min(xx) - 0.01*diff(xx)[1]) )
  axis(1,at=xx, las=1, pos = (min(yy) - 0.01*diff(yy)[1]) )
  
  if(length(y) > 1){
    for(i in 2:length(y)){
      points(x,y[[i]], col = current.palette[i],...)
    }
  }
  
  if(any(is.null(names(y)))){names(y) <- paste("Series",1:length(y))}
  legend("topright",legend=names(y),fill = current.palette)
}



InflationPlot <- function(x,title = "Inflation Plot", zoom = 2,...){

  # x should be -log10 p-values
  # title should be the overall title
  # zoom is a factor giving how much we should zoom in on the bottom left corner in the second plot (adjust as necessary)

  # the blue line is the x = y line
  # the organge line is OLS run on the 0.25 to 0.5 quantiles of the observed values,
  # yielding an estimate of the inflation present

  # Note that inflation observed here may be true to proximal contamination from nearby true signals
  # Also, if we are plotting the distribution of approximate p-values (point estimates),
  # then there could be a bias introduced via that approximation.

  InflationPlot.internal <- function(x,y,zz){
    plot(x,y,ylab="",
         xlab = "",
         bty="n",las=1,pch="*",yaxt="n",xaxt="n",ylim=range(zz),
         xlim = range(zz))
    axis(2,at=zz, las=1, pos = 0 )
    axis(1,at=zz, las=1, pos = 0 )
    mtext("Observed Quantiles",2,2.5,cex = 1.25)
    mtext("Expected Quantiles",1,2,cex=1.25)
    segments(0,0,max(zz),max(zz),lty=2,col="blue")
  }

  old.par <- par(no.readonly = T)

  par(mfrow=c(1,2),oma = c(0, 0, 3, 0),mar = c(4.1,4.1,0.1,1.1))

  y <- quantile(x,probs = ppoints(1e6))
  if(max(y) < max(x)){
    extreme.obs <- x[x>max(y)]
    extreme.probs <- ecdf(x)(extreme.obs)
    extreme.probs[extreme.probs==1] <- 1-0.5/length(x)
    y <- c(y,extreme.obs)
    x <- qexp(c(ppoints(1e6),extreme.probs),rate=log(10))
  }else{
    x <- qexp(ppoints(1e6),rate=log(10))
  }
  zz <- pretty(seq(0,max(x,y),len=10))
  
  y2 <- y[floor(length(y)/4):floor(length(y)/2)]
  x2 <- x[floor(length(y)/4):floor(length(y)/2)]
  
  m1 <- lm(y2 ~ x2)
  
  InflationPlot.internal(x,y,zz)
  segments(0,predict(m1,data.frame(x2=0)),max(zz),predict(m1,data.frame(x2=max(zz))),lty=2,col="darkorange")
  
  y <- y[1:floor(length(y)/zoom)]
  x <- x[1:floor(length(x)/zoom)]
  zz <- pretty(seq(0,max(x,y),len=10))
  
  InflationPlot.internal(x,y,zz)
  segments(0,predict(m1,data.frame(x2=0)),max(zz),predict(m1,data.frame(x2=max(zz))),lty=2,col="darkorange")
  
  mtext(title, outer = TRUE, cex = 1.5,line=1.5)
  mtext(paste0("(intercept: ",signif(m1$coefficients[1],4),"   slope: ",signif(m1$coefficients[2],4),")"),
        outer = TRUE,line = 0.25 )
  par(old.par)
  
  ans <- m1$coefficients
  names(ans) <- c("intercept","slope")
  return(ans)
}






