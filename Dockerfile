FROM openanalytics/r-base
#FROM wg99526/micloudtest

#MAINTAINER Tobias Verbeke "tobias.verbeke@openanalytics.eu"
MAINTAINER Won Gu "rndnjs526@gmail.com"

# system libraries of general use
RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libssl1.1 \
    libxml2-dev \
    build-essential \
    r-base-dev \
    pkg-config \
    cmake \
    && rm -rf /var/lib/apt/lists/*

# system library dependency for the euler app
RUN apt-get update && apt-get install -y \
    libmpfr-dev \
    && rm -rf /var/lib/apt/lists/*

# basic shiny functionality
RUN R -e "install.packages(c('shiny', 'rmarkdown'), repos='https://cloud.r-project.org/')"

# install dependencies of the MiCloud app
RUN R -e "install.packages(c('seqinr', 'shinydashboard', 'proxy', 'dashboardthemes', 'tidyverse', 'plotly', 'shinyWidgets', 'shinyjs', 'googleVis', 'xtable'), repos='https://cloud.r-project.org/')"
RUN R -e "install.packages(c('DT', 'htmltools', 'phangorn', 'bios2mds', 'zip', 'zCompositions', 'dplyr', 'forestplot', 'quantreg', 'fossil', 'picante' ), repos='https://cloud.r-project.org/')"
RUN R -e "install.packages(c('entropart', 'lme4', 'lmerTest', 'dirmult', 'robustbase', 'robCompositions', 'BiasedUrn'), repos='https://cloud.r-project.org/')"
RUN R -e "install.packages(c('CompQuadForm', 'GUniFrac', 'ecodist', 'MiRKAT', 'gridExtra', 'ggplot2', 'patchwork', 'ggthemes', 'erer', 'DiagrammeR', 'stringr'), repos='https://cloud.r-project.org/')"
RUN R -e "install.packages(c('devtools', 'betareg', 'remotes'), repos='https://cloud.r-project.org/')"

RUN R -e "remotes::install_github('joey711/phyloseq')"
RUN R -e "remotes::install_github('joey711/biomformat')"

RUN R -e "remotes::install_github('hk1785/GLMM-MiRKAT')"
RUN R -e "remotes::install_github('nyiuab/NBZIMM')"
RUN R -e "install.packages(c('gridGraphics', 'compositions'), repos = 'https://cloud.r-project.org/')"
RUN R -e "install.packages(c('rgl', 'vegan3d', 'pca3d', 'jpeg', 'splitTools', 'survival', 'survminer', 'coin'), repos = 'https://cloud.r-project.org/')"
RUN R -e "install.packages(c('randomForestSRC', 'kableExtra', 'caret', 'randomForest', 'glmnet'), repos = 'https://cloud.r-project.org/')"

# copy the app to the image
RUN mkdir /root/app
COPY app /root/app
COPY Rprofile.site /usr/lib/R/etc/

COPY app/Data/MiCloud_S_Ex.Rdata /root/app

COPY app/MiDataProc.Data.Upload.R /root/app
COPY app/MiDataProc.Alpha.Cross.Sectional.R /root/app
COPY app/MiDataProc.Beta.Cross.Sectional.R /root/app

COPY app/MiDataProc.Surv.Model1.R /root/app
COPY app/MiDataProc.Surv.Model3.Alpha.R /root/app
COPY app/MiDataProc.Surv.Model3.Beta.R /root/app
COPY app/MiDataProc.Surv.Model3.Taxa.R /root/app
COPY app/MiDataProc.Surv.Model4.R /root/app
COPY app/www/MiSurv_workflow.png /root/app

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/root/app')"]
