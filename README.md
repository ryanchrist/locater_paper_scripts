# Introduction
This repository provides instructions, scripts, and source code needed to replicate the simulations and figures found in 

```Christ R. R., Wang X., Aslett L. J. M., Steinsaltz D., Hall, I. M. (2025). "Clade Distillation for Genome-wide Association Studies." Genetics, to appear.```

An early pre-print of this paper is available here (https://doi.org/10.1101/2024.09.30.615852):

```Christ R. R., Wang X., Aslett L. J. M., Steinsaltz D., Hall, I. M. (2025). "Clade Distillation for Genome-wide Association Studies." bioRxiv. doi:10.1101/2024.09.30.615852.```

This paper introduces a the LOCATER methodology for testing testing tree structures (represented via a set of discrete clades and/or a relatedness matrix) for association with traits of interest.  Installation instructions, vignettes, documentation, and source code for the latest version of our R package implementing this methodology, `locater`, can be found here: https://github.com/ryanchrist/locater.

Guidance for replicating our simulations and figures is divided into three sections below. The first is for investigators aiming to run the version of `locater` used in our paper (rather than the latest version of `locater`) on their own data. The second provides the code and simulation data tables needed for recreate the figures presented in the `locater` paper. For investigators aiming to reproduce our simulation data tables from scratch, the third section provides instructions, code, and a link to the docker image we built for running the simulations on our cluster. A Zenodo record containing a snapshot of this repo at its first release (v1.0.0) is available at [doi:10.5281/zenodo.16543318](www.doi.org/10.5281/zenodo.16543318).

# Running locater 1.0.0
For investigators aiming to run their own experiments using the precise version of `locater` used in our paper, the `locater_v1` directory provides a Dockerfile for building a minimalist version-specific Docker image. The image imports a rocker R 4.4 base image. For all R package dependencies that are stored on CRAN or Bioconductor, the appropriate time-stamped versions are imported. For the package dependencies we maintain and `locater`, versioned tar.gz files are stored in the `locater_v1` folder and installed by the Dockerfile. The version of `locater` provided as a tar.gz corresponds to release 1.0.0 (https://github.com/ryanchrist/locater/releases/tag/v1.0.0). The version of the ancestry inference engine `kalis` provided as a tar.gz corresponds to release 2.0.2 (https://github.com/louisaslett/kalis/releases/tag/v2.0.2). The resulting Docker image takes up only 1.18 Gb. Investigators are encouraged to build onto this Docker image to facilitate their own experiments.

# Recreating Figures from Simulation Data Tables
Data tables storing our simulation results (as csv.gz files) are available on Zenodo at [doi:10.5281/zenodo.16423205](www.doi.org/10.5281/zenodo.16423205). Each csv.gz corresponds to a set of experiments as delineated below. Each csv.gz can be directly loaded by the R script with the corresponding name in the `shark/simulation/association_sims/ryan_viz/` directory to produce the corresponding paper figures. Each R script requires that the paths for inputs and outputs be declared at the top of the script, but otherwise no modification of the scripts should be necessary to produce the figures from the csv.gz files. For instance, the path to the csv.gz files downloaded from Zenodo should be provided as the `collected_results_path` at the top of each script. The `locater14_viz_revision1.R` script should be run before running all other visualization scripts (since it calculates the genome wide discovery thresholds to be used), otherwise, the scripts can be run in any order.

- `locater14_viz_revision1.R` produces the main text dotplot figures and supplemental plots based on power simulations where all causal variants were observed and fell within a 10kb window
- `locater15_viz_revision1.R` produces the main text dotplot figure and supplemental plots based on power simulations where all causal variants were hidden and fell within a 10kb window
- `locater16_viz_revision1.R` produces supplemental plot based on power simulations where all causal variants were observed and fell within a 100kb window
- `locater17_viz_revision1.R` produces supplemental plot based on power simulations where all causal variants were hidden and fell within a 100kb window
- `locater12_viz.R` produces the QQ-plots in the supplement resulting from our null simulations
- `implied_effect_size_viz.R` produces the figures and table in the supplement reporting the effect sizes implicit in our simulations
- `locater15_log_timing_viz.R` recalculates the timing benchmark results reported in the main text from the four corresponding log files


# Rerunning Simulations

All scripts required to run the simulations presented in our pre-print https://www.biorxiv.org/content/10.1101/2024.09.30.615852 are available within the `shark` sub-directory. All simulations were run under the `mini-shark` docker image (see below). All simulations were run using the IBM LSF job management system on an academic university cluster using a job submission template found in `shark/simulation/association_sims/ryan_run_locater_sim.sh`. That script can be used to call to any of the seven simulation R scripts in `shark/simulation/ryan_sims_parameters` described below.

- `locater12.R` was used to run our null simulations
- `locater14.R` was used to run our power simulations assuming all causal variants were observed and fell within a 10kb window
- `locater15.R` was used to run our power simulations assuming all causal variants were hidden and fell within a 10kb window
- `argneedle14.R` was used to re-run the power simulations assuming all causal variants fell within a 10kb window (observed or hidden) using `ARG-NEEDLE` rather than `locater`
- `staar14.R` was used to re-run the power simulations assuming all causal variants fell within a 10kb window were observed using `STAAR` rather than locater
- `locater16.R` was used to run our power simulations assuming all causal variants were observed and fell within a 100kb window
- `locater17.R` was used to run our power simulations assuming all causal variants were hidden and fell within a 100kb window
- `argneedle16.R` was used to re-run the power simulations assuming all causal variants fell within a 100kb window (observed or hidden) using ARG-NEEDLE rather than locater

The results of these scripts were then written to intermediary files (`.rds` files). The corresponding scripts we used for collecting these `.rds` files into the csv.gz data tables described in the previous section can be found in `shark/simulation/association_sims/ryan_collect_sims/`. Please see the *Recreating Figures from Simulation Data Tables* section above to produce our paper figures from these csv.gz files.

As described in our paper, we used [msprime](https://tskit.dev/msprime/docs/stable/installation.html) to run our haplotype simulations. Our public docker image containing an installation of `msprime` and other python dependencies required to directly run our simulation scripts from `R` is available on DockerHub: https://hub.docker.com/repositories/rchrist7 . The version of the image used to run our simulations is pinned under tag `v1.0.0`. On a system with Docker installed, R can be launched in a interactive session under this image by running

```
docker run -it rchrist7/mini-shark:v1.0.0 /bin/bash R
```
