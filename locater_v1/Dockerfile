FROM rocker/r-ver:4.4

RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
    git \
    zlib1g-dev

ENV WORKING_DIR=/tmp
RUN cd $WORKING_DIR

# Set up R for package versioning 
#########################################
RUN R -e "install.packages('remotes')"
RUN R -e "install.packages('groundhog')"
RUN R -e "groundhog::set.groundhog.folder('/root/R_groundhog/')" 
ENV CRAN_DATE=2025-04-09
# Bioconductor version 3.20 corresponds to the CRAN_DATE of 2025-04-09
ENV BIOC_VER=3.20

# Install kalis
#########################################

# install kalis rhdf5 dependency from Bioconductor 
# Since rhdf5 doesn't have any dependencies outside base R, 
# remotes won't pull in any other dependencies here that 
# require version control to be consistent with CRAN_DATE
RUN R -e "remotes::install_bioc('$BIOC_VER/rhdf5')"

# install kalis CRAN dependencies 
RUN R -e "groundhog::groundhog.library('dplyr','$CRAN_DATE')"
RUN R -e "groundhog::groundhog.library('glue','$CRAN_DATE')"
RUN R -e "groundhog::groundhog.library('rlang','$CRAN_DATE')"
RUN R -e "groundhog::groundhog.library('digest','$CRAN_DATE')"
RUN R -e "groundhog::groundhog.library('checkmate','$CRAN_DATE')"
RUN R -e "groundhog::groundhog.library('RColorBrewer','$CRAN_DATE')"
RUN R -e "groundhog::groundhog.library('lobstr','$CRAN_DATE')"
RUN R -e "groundhog::groundhog.library('prettyunits','$CRAN_DATE')"

# set kalis installation flags
RUN mkdir ~/.R && \
    echo 'CFLAGS= -march=native -mtune=native -O3 -Wall -pipe -pedantic -pthread' > ~/.R/Makevars

# install kalis from source
COPY kalis.tar.gz $WORKING_DIR
RUN R -e "setwd('$WORKING_DIR'); install.packages('kalis.tar.gz', repos = NULL, type = 'source', build_vignettes = TRUE, verbose = TRUE); stopifnot(require(kalis))"


# Install locater
#########################################

# Install locater dependencies 
RUN R -e "groundhog::groundhog.library('RcppRoll','$CRAN_DATE')"
RUN R -e "groundhog::groundhog.library('RSpectra','$CRAN_DATE')"
RUN R -e "groundhog::groundhog.library('data.table','$CRAN_DATE')"
COPY QForm.tar.gz $WORKING_DIR
RUN R -e "setwd('$WORKING_DIR'); install.packages('QForm.tar.gz', repos = NULL, type = 'source', build_vignettes = TRUE, verbose = TRUE); stopifnot(require(QForm))"
COPY renyi.tar.gz $WORKING_DIR
RUN R -e "setwd('$WORKING_DIR'); install.packages('renyi.tar.gz', repos = NULL, type = 'source', build_vignettes = TRUE, verbose = TRUE); stopifnot(require(renyi))"

# Install locater
COPY locater.tar.gz $WORKING_DIR
RUN R -e "setwd('$WORKING_DIR'); install.packages('locater.tar.gz', repos = NULL, type = 'source', build_vignettes = TRUE, verbose = TRUE); stopifnot(require(locater))"



ENTRYPOINT \bin\bash