
create_genomic_region_real_data <- function(#job.index,
                               #shark.dir,
                               #champs.path,
                               sim.output.dir,
                               #N,
                                file.name,
                               # chr.length = 1e6L,
                               chr.start.pos=1L,
                               # causal.window.size = 1e4L,
                               # num.causal.vars = 0L,
                               # 
                               # fixed.map = NULL,
                               
                               diagnostics = FALSE
){
  
  # parse options
  ########################
  force(sim.output.dir)
  force(file.name)
  force(chr.start.pos)
  force(diagnostics)
  #require(rhdf5)

  haps.path <- paste0(sim.output.dir,file.name)
  pos <- c(h5read(haps.path,name="POS")) 
  map <- c(h5read(haps.path,name="map")) 
  chr.length <- max(pos)-min(pos)+1
  one_var <- t(ReadHaplotypes(haps.path,
                              loci.idx = 1,
                              transpose=T)$haps)
  N.haps <- length(one_var)
  pop.labels <- rep("FIN",N.haps)
  causal.var.idx <- NA
  
  list(
    "sim.output.dir" = sim.output.dir,
    "haps" = haps.path,
    "pos" = pos,
    "map" = map,
    "chr.length" = chr.length,
    "chr.start.pos" = min(pos),
    "N.haps" = N.haps,
    "pop.labels" = pop.labels, 
    # pop.labels must be declared for each haplotype, not each individual
   # "causal.vars" = genotypes, # this is used only in oracles
   # "causal.window" = causal.window,
   # "causal.window.idx" = causal.window.idx,
  #  "count.range.causal.vars" = count.range.causal.vars,
    "causal.var.idx" = causal.var.idx#,
  #  "candidate.causal.var.idx" = candidate_causal_vars + causal.window.idx[1] - 1L
  )
}






######################################
create_genomic_region_WashU_CCDG <- function(
  genome.input.dir,
  file.prefix,
  propagation.window.size = NULL,
  diagnostics = FALSE){
  # parse options
  ########################
  force(genome.input.dir)
  force(file.prefix)
  #force(chr.start.pos)
  force(propagation.window.size)
  force(diagnostics)
  #require(rhdf5)
  
  haps.path <- paste0(genome.input.dir,file.prefix,".phased_multi_merged_unique_ID_ancestral_split_sample_AC_var.hap.gz")
  #pos <- c(h5read(haps.path,name="POS")) 
  #map <- c(h5read(haps.path,name="map")) 
  
  map.dir <- "/xinxin/data/WashU_CCDG/recomb_map/"
  map.file <- read.table(paste0(map.dir,file.prefix,"_AF_weighted_cM.txt"),header=TRUE)
  map <- map.file$cM
  
  pos <- map.file$POS
  
  pass <- map.file$REF != "N" & nchar(map.file$REF) + nchar(map.file$ALT) <= 2
  
  chr.length <- max(pos)-min(pos)+1
  # one_var <- t(ReadHaplotypes(haps.path,
  #                             loci.idx = 1,
  #                             transpose=T)$haps)
  
  legend.file <- paste0(genome.input.dir,file.prefix,".phased_multi_merged_unique_ID_ancestral_split_sample_AC_var.samples")
  legend.table <- read.table(legend.file,header=TRUE)
  
  N.haps <- nrow(legend.table) * 2
  
  pop.labels <- rep("NA",N.haps)
  causal.var.idx <- NA
  
  if (diagnostics & !is.null(propagation.window.size)){
    propagation.window  <- (chr.length + c(-propagation.window.size,propagation.window.size))/2 + min(pos)
    propagation.window.idx <- which(pos >= propagation.window[1] & pos <= propagation.window[2] & pass) 
    print("Diagnostics mode on. not propagating all pass vars")
  } else {
    propagation.window.idx <- pass
  }
   
  
  list(
    "genome.input.dir" = genome.input.dir,
    "haps" = haps.path,
    "pos" = pos,
    "map" = map,
    "chr.length" = chr.length,
    "chr.start.pos" = min(pos),
    "N.haps" = N.haps,
    "pop.labels" = pop.labels, 
    "pass" = pass,
    "propagation.window.idx" = propagation.window.idx,
    # pop.labels must be declared for each haplotype, not each individual
    # "causal.vars" = genotypes, # this is used only in oracles
    # "causal.window" = causal.window,
    # "causal.window.idx" = causal.window.idx,
    #  "count.range.causal.vars" = count.range.causal.vars,
    "causal.var.idx" = causal.var.idx#,
    #  "candidate.causal.var.idx" = candidate_causal_vars + causal.window.idx[1] - 1L
  )
  
}




