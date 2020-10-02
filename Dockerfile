FROM bioconductor/bioconductor_docker:RELEASE_3_11
LABEL authors="Kaur Alasoo" \
      description="Docker image containing all requirements for the eQTL Catalogue colocalisation workflow"

RUN R -e "BiocManager::install(c('dplyr','tidyr','assertthat','devtools','GenomicRanges','readr','coloc','optparse','Rsamtools','readr','stringr'))"
RUN R -e "devtools::install_github('mrcieu/gwasvcf')"

