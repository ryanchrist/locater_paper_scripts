rwcombn <- function(weights,k){
  # sample a random combination of size k from the integers 1:length(weights)
  # where the combination is drawn with probability proportional
  # to the product of the weights of its components

  # this function uses a recursion to avoid the high memory footprint required
  # by the straightforward approach using utils::combn to generate and store
  # all possible k-combinations

  # weights must be a vector of non-negative weights
  # k must be a positive integer <=length(weights)

  # input checks
  #####################
  if( !is.vector(weights) || any(!is.numeric(weights)) || any(weights < 0)){stop("weights must be a numeric vector with non-negative entries.")}

  if( !is.vector(k) || (length(k) > 1) || !is.numeric(k) ){stop("k must be a positive integer (with length 1).")}
  k.new <- as.integer(k)
  if(k.new != k || k.new <= 0L ){stop("k must be a positive integer.")}

  N <- as.integer(length(weights))


  # define required recursion
  #####################
  g <- function(i,K){

    if(K == 0){return(1)}
    if(K > 0 & i > N){return(0)}

    idx.list <- as.list(i:N)
    mask.length <- N-i+1L
    mask <- rep(0,mask.length)

    for(j in i:N){
      ii <- j-i+1L
      temp <- g(j+1L,K-1L)
      mask[ii] <- temp
      if(temp>0){ idx.list[[ii]] <- c(j,attr(temp,"idx")) }
    }

    if(all(mask==0)){return(0)}

    prob <- mask * weights[i:N]
    temp <- sum(prob)
    attr(temp,"idx") <- idx.list[[sample.int(mask.length,1L,prob = prob)]]
    temp
  }

  # here is a simplified version of the above recursion function
  # if we were just aiming to calculate the normalizing constant and
  # not worry about sampling a k-combination.
  # f <- function(i,K){
  #   if(K == 0){return(1)}
  #   if(K > 0 & i > N){return(0)}
  #   tot <- 0
  #   for(j in i:N){
  #     tot <- tot + w[j]*f(j+1,K-1)
  #   }
  #   tot
  # }

  # run recursion
  #####################
  ans <- g(1L,k.new)

  list("idx" = attr(ans,"idx"),
       "sample.weight" = prod(weights[attr(ans,"idx")]),
       "norm.const" = as.numeric(ans))
}

# Example:
# w <- runif(52)
# k <- 4
# choose(length(w),k) # number of combinations we're about to iterate over
# system.time({res <- rwcombn(w,k)})
# res$idx # combination sampled
# signif(res$sample.weight) # weight of the sampled combination
# signif(res$norm.const) # normalizing constant of weights (summed over all possible combinations)
# signif(res$sample.weight/res$norm.const) # probability of selecting the choosen combination

#champs.path = "/Library/Frameworks/Python.framework/Versions/3.8/bin/champs"
#systems2(champs.path)



########
# Estimate a realistic maf density for array data from UK biobank segment
########

# HERE IS COMMENTED CODE USED TO GENERATE THE FUNCTION WITH THE INBUILT LOOK UP TABLE BELOW
# array.AC.file="~/Dropbox/Hall_group/kalis_power_simulations/data/ukb_snp_chr22_v2.acount"
# array.AC <- read.csv(array.AC.file,sep = "\t")
# # 976754 haps in UK biobank from our analysis and the UKBiobank paper
# num.haps.ukbb <- 976754
# array.AC <- array.AC[array.AC$OBS_CT > 97e4,] # remove variants with lot of missingness
# array.AC$AAF <- array.AC$ALT_CTS / num.haps.ukbb
# maf <- pmin(array.AC$AAF,1-array.AC$AAF)
#
# maf <- maf[maf>0]
# lower.maf.array <- maf[maf <= 0.01]
# upper.maf.array <- maf[maf > 0.01]
# lower.den <- density(lower.maf.array,bw="SJ", n = 128,from = -0.001,to=0.011)
# upper.den <- density(upper.maf.array,bw="SJ", n = 128,from = 0, to = 0.52)
# mixture.proportion <- length(lower.maf.array)/(length(lower.maf.array) + length(upper.maf.array))
#
# paste(lower.den$y,collapse = ",")
# paste(upper.den$y,collapse = ",")
# print(mixture.proportion,digits=22)


calc_array_maf_density <- function(){
  lower.maf.density <- approxfun(
    x = seq(-0.001,0.011,len=128),
    y = c(2.00564974131786e-09,9.38648433234171e-08,3.21009213593481e-06,7.75726117295077e-05,0.00119236693532329,0.0134172015154073,0.103468436288655,0.544648859021534,2.10902643152082,6.09596523242826,14.5621053561795,32.7634086159006,70.5070296176484,136.438532316045,225.004907928636,312.532059288732,376.343279379251,410.361802817355,425.512019744952,432.399719898765,432.755277904341,423.36805235063,403.661695965921,376.459929236942,343.25460848061,305.164857052373,267.286046767356,235.277705736229,211.461902601383,195.420654994318,185.696043974716,179.22422202715,171.552729304097,161.037209269545,150.587549161449,143.777957526034,140.162759279178,135.815902434597,128.089878291347,118.498988703778,110.873935317945,107.032974672339,104.957589884928,100.656962317706,92.8280528118216,84.2500136407172,78.402273515012,75.0685469425902,71.3111004589491,65.6048016005027,59.837291424508,56.8399396345087,56.5711179291444,55.353003423713,50.0732610780989,42.3460118954953,36.7687259215623,36.2917806381774,39.017154776558,40.7702167257885,40.2135932819505,40.0527967470966,41.7246736774501,42.639023283613,40.3062587208412,36.2893028247743,34.0653461715554,34.2031206872749,34.3253941586221,33.11717210902,32.8461796525425,36.0503609126386,41.1111227617081,43.2332004463749,39.8677388014063,33.5589575510033,29.0402391596416,28.73914869411,31.5170137494894,34.5453873577481,35.7780300834306,36.0217363710426,38.0159557190393,42.7337774417493,47.3841524895419,48.2054073216291,43.9309677701916,37.2643894029782,32.1797775764097,31.1337664199359,33.2039732044951,35.6960977068686,36.4512549182912,35.0866641858246,32.7951368880135,32.0152854330412,34.0225039200591,37.0420757052993,38.060472764207,36.7151277741502,34.7841545623895,33.2842573331131,32.3375795043983,32.569517248856,33.7939001577476,34.5551034095876,34.4138794008715,34.3314506274407,34.7864212779999,35.4421916103514,35.7267132668698,35.1396081838707,33.7416283261848,32.6816539033618,32.8699552195441,32.5693427913415,28.1023407156908,18.954491179485,9.38911480920395,3.31224155915431,0.811359570643546,0.142798468862798,0.0169877339788183,0.0013862226947281,8.34399951473368e-05,3.19878896691355e-06,8.74287017100397e-08,1.75909381467735e-09),
    method = "linear",rule = 2)

  upper.maf.density <- approxfun(
    x = seq(0,0.52,len=128),
    y = c(0.0106237791843794,0.328576467312649,2.81262991704202,7.93823134521883,11.5971527346151,12.8364602456208,12.0393369031773,10.4225045674638,8.76365486837091,7.40042222782672,6.37894724918335,5.66091620523063,5.1736033343666,4.78510329847798,4.45987653737903,3.96135490859821,3.69543388519855,3.42970637395796,3.22429461837295,3.05200442801456,2.74140342024274,2.68596604655819,2.6553035686577,2.60461153461981,2.60677829418859,2.45502728348925,2.2987565179282,1.98421253835978,1.6998796855768,1.55840635526189,1.56979793795303,1.57012014536264,1.51731663765083,1.47607925537791,1.52016887158938,1.60033590826038,1.43945674998076,1.38115876540773,1.49434854935489,1.62304388576153,1.79904896894004,1.71290623207706,1.40723810351915,1.3337515594508,1.36926257104884,1.35816698738644,1.40355317708767,1.36068338186844,1.11132874468276,1.08331021240552,1.21689062515942,1.35628738810133,1.37169023102031,1.23114541109678,1.01778704146528,1.01194654106321,1.10598747208803,1.21384363451636,1.27061700118153,1.25526621560854,1.29565216957121,1.20864362271513,1.05532339177471,1.0352232288658,1.16477632172595,1.18946649283725,1.26769155688158,1.32929516145679,1.24619769503278,1.2248481837852,1.25299906541151,1.09165668416429,0.984071437224724,0.988348512863761,0.979470592890605,1.01277828631465,0.990713595987226,1.018209512138,0.986117195865378,0.919728606140354,0.921687554695844,0.893939392369599,0.923819729812912,0.919853682870172,0.97808058942348,1.05190215929353,0.942785995048442,0.860206335978903,0.894728800467327,0.899764893368233,0.97681261029377,1.03271606417706,0.96610940257794,0.936093913956895,0.958285799406,0.954954654647557,0.895678313092632,0.947320732448659,0.996278604261303,0.872959016889425,0.882901560644202,0.992719262370136,1.01168209060488,0.968526762925748,0.956274986133653,1.02672093735073,0.983529074708035,0.740487458878662,0.672672899368676,0.834304575898841,1.01338052754489,1.02640951858102,1.0177172756171,1.02685340973358,0.964346991725961,0.818242770295206,0.781946430910657,0.873476680050878,0.876833382040123,0.933443558533416,0.982092245085582,0.88227498320666,0.532729333199611,0.133668633897931,0.00987422223127303,0.000172382202698358,6.17137478687726e-07,6.49385541038937e-10),
    method = "linear",rule = 2)

  mixture.proportion <- 0.15782417161463321853

  function(x){
    y <- lower.maf.density(x) * mixture.proportion + upper.maf.density(x) * (1-mixture.proportion)
    y[ x< 0 | x > 0.5] <- 0
    y
  }
}

# array_maf_density <- calc_array_maf_density()
# xxx <- seq(0,0.5,len=1e5)
# plot(xxx,array_maf_density(xxx),type="l",ylab="Density",las=1,bty="n",xlab = "maf",main="Estimated UK Biobank Chr22 Array Data MAF Density")
# sum(array_maf_density(xxx))*mean(diff(xxx))




make_reject_sampler <- function(maf){

  # Get estimate MAF density: f.seq
  maf <- maf[maf>0]
  lower.maf.array <- maf[maf <= 0.01]
  upper.maf.array <- maf[maf > 0.01]
  lower.den <- density(lower.maf.array,bw="nrd", n = 128,from = -0.001,to=0.011)
  upper.den <- density(upper.maf.array,bw="nrd", n = 128,from = 0, to = 0.52)
  mixture.proportion <- length(lower.maf.array)/(length(lower.maf.array) + length(upper.maf.array))

  lower.maf.density <- approxfun(lower.den, method="linear", rule=2)
  upper.maf.density <- approxfun(upper.den, method="linear", rule=2)

  f.seq <- function(x){
    y <- lower.maf.density(x) * mixture.proportion + upper.maf.density(x) * (1-mixture.proportion)
    y[ x< 0 | x > 0.5] <- 0
    y
  }

  f.array <- calc_array_maf_density()


  # Calculate Dominance Constant for Rejection Sampler
  #########################################################
  # Need to check for dominance at all knots used for KDE inside the domain [0,0.5] and at the end points of the domain 0 and 0.5.  With linear interpolation,
  # dominance at all of these points implies dominance everywhere over the domain.
  xxx <- c(0,0.5,seq(-0.001,0.011,len=128),seq(0,0.52,len=128))
  xxx <- xxx[xxx >= 0 & xxx <= 0.5]
  xxx <- sort(unique(xxx))

  f <- function(x) min(f.seq(xxx)*x-f.array(xxx))
  B <- uniroot(f,c(0,100),extendInt = "upX")$root
  # correct for any numerical tolerance issues by making sure we indeed have dominance
  while(any(B * f.seq(xxx) < f.array(xxx))){ B <- B * 1.001}

  # return function that, for a set of samples x, returns TRUE for accepted samples and FALSE for rejected samples
  function(x){ runif(length(x)) < (f.array(x) / (B*f.seq(x))) }
}



sim_genomic_region <- function(job.index,
                               shark.dir,
                               champs.path,
                               sim.output.dir,
                               N,
                               chr.length = 1e6L,
                               causal.window.size = 1e4L,
                               num.causal.vars = 0L,
                               count.range.causal.vars = c(1,sum(N)),
                               fixed.map = NULL,
                               remove_singletons = TRUE,
                               diagnostics = FALSE,
                               max.iter = Inf,
                               save_tree = FALSE
){

  # parse options
  ########################

  if(num.causal.vars >= 1){
    if(as.integer(num.causal.vars) != num.causal.vars){stop("num.causal.vars must be an integer if >= 1")}
    num.causal.vars <- as.integer(num.causal.vars)
  } else if(num.causal.vars > 0){
    num.causal.vars <- as.double(num.causal.vars) # takes care of case when
  } else if(num.causal.vars == 0){
    num.causal.vars <- 0L
  } else {
    stop("num.causal.vars must be >= 0.")
  }
  # From here, if num.causal.vars is an integer, it will be taken as a count,
  # if it is a double, it will be taken as a lower bound on the proportion of causal variants in the gene.


  haps.path <- paste0(sim.output.dir,"msprime_haps","_",job.index,".h5")
  pos.path <- paste0(sim.output.dir,"positions","_",job.index,".txt")
  map.path <- paste0(sim.output.dir,"map","_",job.index,".txt")
  carriers.path <- paste0(sim.output.dir,"carriers_class_0_run_",job.index,".txt")

  causal.window <- (chr.length + c(-causal.window.size,causal.window.size))/2

  # Simulate Haplotypes
  counter <- 0L
  ready <- FALSE
  while(!ready){

    if(counter == max.iter){stop(paste("Failed to simulate a haplotype dataset that met the specified clade constraint in first",max.iter,"iterations."))}
    counter <- counter + 1L

    if(file.exists(haps.path)){file.remove(haps.path)}
    if(file.exists(pos.path)){file.remove(pos.path)}
    if(file.exists(map.path)){file.remove(map.path)}
    if(file.exists(carriers.path)){file.remove(carriers.path)}


    base.args <- c(
      paste("-N",N[1],N[2],N[3]),
      paste("-o",sim.output.dir),
      paste("-i",job.index))

    if (is.null(fixed.map)){
      print("Using random map")
      base.args <- c(base.args,
                     paste("-L",as.integer(chr.length)),
                     "-M any")
    } else {
      print(paste("Using fixed map from",fixed.map))
      base.args <- c(base.args,
                     paste("-M",fixed.map))
    }

    if(remove_singletons){
      base.args <- c(base.args,"-remove singletons")
    } else {
      base.args <- c(base.args,"-remove none")
    }

    if(save_tree){
      base.args <- c(base.args,"-save_tree")
      tree.path <- paste0(sim.output.dir,"ts","_",job.index,".tsz")
    } else {
      tree.path <- NA_character_
    }



    if(is.integer(num.causal.vars) && num.causal.vars > 0){

      # if num.causal.vars is integer, for efficiency try to enforce the clade constraint inside python
      system2(champs.path,c(base.args,
                            paste("-c",num.causal.vars,"obs","any",
                                  count.range.causal.vars[1]/sum(N),
                                  (count.range.causal.vars[2]+0.5)/sum(N)),
                            paste("-w",as.integer(causal.window.size))))
    } else {
      system2(champs.path,base.args)
    }

    if(diagnostics){print(paste(counter,"calls to champs have been made."))}

    pos <- c(read.table(pos.path)[,1])
    map <- c(read.table(map.path)[,2])

    causal.window.idx <- which(pos >= causal.window[1] & pos <= causal.window[2])   # variants in gene / causal region

    if(length(causal.window.idx) < num.causal.vars){ next }

    # Test if we meet clade constraint
    ######################################
    genotypes <- t(ReadHaplotypes(haps.path,
                                  loci.idx = causal.window.idx,
                                  transpose=T)$haps)
    counts <- c(colSums(genotypes))
    candidate_causal_vars <- which(count.range.causal.vars[1] <= counts & counts <= count.range.causal.vars[2])

    # remove perfectly linked causal vars as candidates
    candidate_causal_vars <- candidate_causal_vars[!duplicated(t(genotypes[,candidate_causal_vars]))]

    if(length(candidate_causal_vars) < num.causal.vars){ next }

    ready <- TRUE
  }

  if(diagnostics){print(paste(length(candidate_causal_vars),"unique vars meet frequency constraint."))}


  #subsample candidate_causal_vars and genotypes to match count or proportion given by num.causal.vars
  select_causal_vars <- sample(candidate_causal_vars,
                                  size = if(is.integer(num.causal.vars)){
                                    num.causal.vars
                                  }else{
                                    ceiling(num.causal.vars * length(candidate_causal_vars))
                                  },
                                  replace = FALSE)

  genotypes <- as.matrix(genotypes[,select_causal_vars])
  causal.var.idx <- causal.window.idx[select_causal_vars]

  if(!file.exists(haps.path)){stop("haplotype simulation via champs may have failed b/c haplotype file haps.path does not exist")}

  list(
    "sim.output.dir" = sim.output.dir,
    "haps" = haps.path,
    "pos" = pos,
    "map" = map,
    "chr.length" = chr.length,
    "chr.start.pos" = 1L,
    "N.haps" = sum(N),
    "pop.labels" = as.factor(c(rep("YRI",N[1]),rep("CEU",N[2]),rep("CHB",N[3]))),
    # pop.labels must be declared for each haplotype, not each individual
    "causal.vars" = genotypes, # this is used only in oracles
    "causal.window" = causal.window,
    "causal.window.idx" = causal.window.idx,
    "count.range.causal.vars" = count.range.causal.vars,
    "causal.var.idx" = causal.var.idx,
    "candidate.causal.var.idx" = candidate_causal_vars + causal.window.idx[1] - 1L,
    "ts" = tree.path
  )
}



calc_var_exclude_idx <- function(x,
                                 exclude.vars = "none",
                                 max.iter = Inf
){

  if(exclude.vars == "none"){
    return(integer())
  } else if(exclude.vars == "causal"){
    return(x$causal.var.idx)
  } else if(exclude.vars!="array"){
    stop("exclude.vars must be none, causal, or array when calculating variant filter")
  }

  mu.nbinom = 2e-5 * (log(x$N.haps-1) - digamma(1)) * x$chr.length # expected number of mutations in chr.length
  # From Waterson's Estimator: log(n-1) - digamma(1) ~ the n-1 harmonic sum since digamma(1) = Euler-Mascheroni constant
  # We took the variant density of the uk biobank array data
  # and used Waterson's Estimator to calculate some mutation parameter theta which can be thought of the mutation rate of array data
  # since randomly filtering a Poisson process is still a Poisson process.  This yielded an estimate of theta that was 2e-5 / base pair.
  # We noticed that the resulting distribution was still massively overdispersed
  # (the number of mutations should be the sum of geometric rvs with different rates so that the
  # variance should equal the mean but the variance was much higher when taking 1MB subsets of the data).
  # So we settled on using a negative binomial distribution for modelling the number of mutations observed in array data.
  # Setting the neg.binom parameter size=10 allowed us to recover the 1st and 3rd quartiles observed in a histogram of the number of
  # variants in many 1MB subsets of the UK biobank data.
  # When testing out these parameters for smaller sample sizes, the negative binomial yielded reasonable estimates for the number of
  # observed variants.

  range_of_acceptable_num_vars <- qnbinom(c(0.25,0.75),
                                          mu = mu.nbinom,
                                          size = 10)
                                          # 2e-5 and size = 10 here both derived from uk biobank chip data

  maf <- rowMeans(ReadHaplotypes(x$haps,
                                 transpose = TRUE)$haps)
  maf <- pmin(maf, 1-maf)

  keep <- make_reject_sampler(maf)

  # perform rejection sampling
  counter <- 0L
  ready <- FALSE
  while(!ready){
    if(counter == max.iter){warning(paste("Rejection Sampler rejected too many variants even after",max.iter,"attempts at sub-sampling: check MAF density of input dataset.")); return(NULL)}
    counter <- counter + 1L

    vars.to.keep <- keep(maf)
    if(sum(vars.to.keep) > range_of_acceptable_num_vars[2]){ ready <- TRUE}
  }


  # further subsample variants to get the variant density to match that of array data
  # (based on 1st and 3rd quartiles of many 1MB segments extracted from array data)
  num.vars.to.keep <- sample(range_of_acceptable_num_vars[1]:range_of_acceptable_num_vars[2],
                             size = 1,
                             prob = dnbinom(range_of_acceptable_num_vars[1]:range_of_acceptable_num_vars[2],
                                            mu = mu.nbinom,
                                            size = 10))

  vars.to.keep[sample(which(vars.to.keep), size =  sum(vars.to.keep) - num.vars.to.keep, replace = FALSE)] <- FALSE

  which(!vars.to.keep)
}




sim_phenotypes <- function(x, # a simulated genomic region with causal variant annotations
                           type = "quantitative", # for now, quantitative or binary
                           ploidy = 1L,
                           inheritance = "additive",
                           f = function(maf){1}, # a function of maf that returns the effect size to be attributed to each causal variant
                           neg.effect.rate = 0.5,
                           null.sim.second.mode = NULL,
                           A = NULL,
                           skip.null.model.fitting = FALSE,
                           option.list = list("disease.prevalence" = 0.01) # for now, only relevant for binary traits
){

  n <- x$N.haps / ploidy

  if(is.null(A)){
    if (length(levels(x$pop.labels))<2){
      A <- cbind(matrix(rep(1,n),ncol=1),
                 rbinom(n,1,0.5),
                 rnorm(n))
    } else {
      A <- cbind(model.matrix(1:n ~ x$pop.labels[seq(1,x$N.haps,by=ploidy)]),
                 rbinom(n,1,0.5),
                 rnorm(n))
    }
  }

  daf <- colMeans(x$causal.vars)
  #maf <- pmin(daf,1-daf)

  sign.vec <- rep(-1,floor(neg.effect.rate*ncol(x$causal.vars)))
  sign.vec <- c(sign.vec,rep(1,ncol(x$causal.vars)-length(sign.vec)))
  sign.vec <- sample(sign.vec,length(sign.vec),FALSE)

  beta <- sign.vec * f(daf)

  if(ploidy > 1){
    genotypes <- x$causal.vars[seq(1,x$N.haps,by = ploidy),]
    for(j in 2:ploidy){
      if(inheritance == "additive"){
        genotypes <- genotypes + x$causal.vars[seq(j,x$N.haps,by = ploidy),]
      } else if(inheritance == "dominant"){
        genotypes <- genotypes | x$causal.vars[seq(j,x$N.haps,by = ploidy),]
      } else {
        stop("inheritance mis-specified: for now, must be additive or dominant")
      }
    }
    storage.mode(genotypes) <- "integer"
  } else {
    genotypes <- x$causal.vars
  }

  genotypes <- as.matrix(genotypes)

  if(type == "quantitative"){

    genetic.effects <- c(genotypes %*% beta)

    if(is.null(null.sim.second.mode)){
      y <- c(rowSums(A)) + genetic.effects + rnorm(n)
    } else {
      omega <- genetic.effects / sqrt(sum(genetic.effects^2))
      zz <- rnorm(n)
      fz <- sum(omega * zz)
      xx <- rnorm(1) + if(runif(1)>0.5){0}else{null.sim.second.mode}

      log.importance.weight <- log(2)-log1p(exp(dnorm(xx,mean = null.sim.second.mode,log = T)-dnorm(xx,log = T)))
      y <- c(rowSums(A)) + zz + omega * (xx - fz)
    }

    if(!skip.null.model.fitting){
      # Fit null model
      m1 <- lm(y ~ 0 + A)
      Ey <- m1$fitted.values
      SDy <- sqrt(sum(m1$residuals^2)/(n - ncol(A)))
    }

  } else if(type == "pop.binary"){

    initial.logit.Ey <- c(genotypes %*% beta) + c(rowSums(A)) # population membership + random effect + G * beta
    case.num <- option.list$disease.prevalence * n  # how many cases are we expecting

    f <- function(a){ ee <- exp(-(initial.logit.Ey+a)); return(sum(1/(1+ee)) - case.num)}

    alpha_0 <- uniroot(f,lower=-100,upper=100,extendInt = "yes")$root  # find the root of the function, find alpla_0
    # alpha_0 is how much we need to add to the intercept of 1 to get the desired exected disease prevalence.
    logit.Ey <- c(genotypes %*% beta) + c(rowSums(A)) + alpha_0

    Ey <- 1/(1+exp(-logit.Ey))

    y <- as.numeric(runif(n) < Ey)

    if(!skip.null.model.fitting){
      m1 <- glm(y ~ 0 + A, family="binomial")
      Ey <- m1$fitted.values
      SDy <- sqrt(m1$fitted.values * (1-m1$fitted.values))
    }

  } else if(type == "case.control.binary"){
    stop("case.control.binary phenotype simulation not yet implemented.")
  } else {
    stop("Invalid phenotype simulation type requested of sim_phenotypes.")
  }

  list("y" = y,
       "Ey" = if(skip.null.model.fitting){NA}else{Ey},
       "SDy" = if(skip.null.model.fitting){NA}else{SDy},
       "A" = A,
       "ploidy" = ploidy,
       "n" = n,
       "inheritance" = inheritance,
       "log.importance.weight" = if(is.null(null.sim.second.mode)){NULL}else{log.importance.weight})
}


sim_par_phenotypes <- function(x, a, n.effect.sizes, A, genotypes, counts, neg.effect.rate = 0, ploidy=2L){ #, use.forking, nthreads){

  # returns a vector indicating which architectures (in the dataframe a) genomic region x is valid for

  included <- rep(FALSE,n.effect.sizes*nrow(a))

  y <- replicate(n.effect.sizes*nrow(a),NULL,simplify = FALSE)

  #if(!use.forking){

  for(i in 1:nrow(a)){

    # Test if sampled genomic region meets constraints of architecture i
    ###########################################################################

    if(length(x$causal.window.idx) < a$num.causal.vars[i]){ next }

    candidate_causal_vars <- which(a$count.range.min[i] <= counts & counts <= a$count.range.max[i])

    # remove perfectly linked causal vars as candidates
    candidate_causal_vars <- candidate_causal_vars[!duplicated(t(genotypes[,candidate_causal_vars]))]

    if(length(candidate_causal_vars) < a$num.causal.vars[i]){ next }

    # By this point, we know sampled genomic region meets constraints of architecture i
    ######################################################################################

    # simulate phenotype vectors under architecture i for seq and array data
    ######################################################################################

    select_causal_vars <- sample(candidate_causal_vars,
                                 size = a$num.causal.vars[i],
                                 replace = FALSE)

    x$causal.vars <- as.matrix(genotypes[,select_causal_vars])
    x$causal.var.idx <- x$causal.window.idx[select_causal_vars]

    # simulate phenotype for y and y.array

    for(j in 1:n.effect.sizes){

      temp.effect.size.scale <- getElement(a,paste0("effect.size.scale.",j))[i]

      if(!is.na(temp.effect.size.scale)){

        included[(i-1)*n.effect.sizes + j] <- TRUE

        y[[(i-1)*n.effect.sizes + j]] <- sim_phenotypes(x,
                                           type = "quantitative", # for now, quantitative or binary
                                           ploidy = ploidy,
                                           inheritance = a$inheritance[i],
                                           f = switch(a$func.shape[i],
                                                      "1" = function(x){temp.effect.size.scale},
                                                      "-log(x)" = function(x){-log(x) * temp.effect.size.scale},
                                                      "1/x" = function(x){temp.effect.size.scale/x},
                                                      "1/x^2" = function(x){temp.effect.size.scale/x^2},
                                                      "1/sqrt(x)" = function(x){temp.effect.size.scale/sqrt(x)}), # a function of maf that returns the effect size to be attributed to each causal variant
                                           neg.effect.rate = neg.effect.rate,
                                           null.sim.second.mode = NULL,
                                           A = A,
                                           skip.null.model.fitting = TRUE,
                                           option.list = list("disease.prevalence" = 0.01))$y
      }

    }

  }

  return(list("included" = included,
              "y" = do.call(cbind,y)))
}





make_tuning_pseudo_phenotype <- function(x,ploidy){

  n <- x$N.haps / ploidy
  if (ploidy==1){
    y <- rep(0L,x$N.haps)
    y[rowSums(x$causal.vars) > 0] <- 1L
  } else {
    if (ploidy==2){
      y <- sample(c(0,1), replace=TRUE, size=n)
      print("Ploidy=2 tuning mode, testing only.")
    } else {
      stop("only ploidy==1 or 2 supported")
    }
  }
   # this construction of y allows us to do testing either at a core locus (to be removed from the data during testing) or
  # to actually remove a whole set of variants and do testing somewhere in the middle

  # only intercept and population membership indicator

  if (length(levels(x$pop.labels))<2){
    A <- matrix(rep(1,n),ncol=1)
  } else {
    A <- model.matrix(1:n ~ x$pop.labels)
  }
  if (!is.null(x$bg.vecs)){
    A <- cbind(A, x$bg.vecs)
  }

  m1 <- glm(y ~ 0 + A, family="binomial")

  Ey <- m1$fitted.values
  SDy <- sqrt(m1$fitted.values * (1-m1$fitted.values))


  list("y" = y,
       "Ey" = Ey,
       "SDy" = SDy,
       "A" = A,
       "ploidy" = ploidy,
       "inheritance" = "additive")

}



sim_background_cov <- function(x, ploidy = 2L){
  n <- x$N.haps / ploidy

  if (length(levels(x$pop.labels))<2){
    cbind(matrix(rep(1,n),ncol=1),
          rbinom(n,1,0.5),
          rnorm(n))
  } else {
    cbind(model.matrix(1:n ~ x$pop.labels[seq(1,x$N.haps,by=ploidy)]),
          rbinom(n,1,0.5),
          rnorm(n))
  }
}


