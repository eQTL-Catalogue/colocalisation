FROM bioconductor/bioconductor_docker:RELEASE_3_11
LABEL authors="Kaur Alasoo" \
      description="Docker image containing all requirements for the eQTL Catalogue colocalisation workflow"

RUN R -e "BiocManager::install(c('Biostrings','remotes', 'dplyr', 'optparse', 'readr', 'tidyr','assertthat', 'seqminer', 'qvalue', 'coloc','GenomicRanges','SummarizedExperiment','VariantAnnotation'))"
RUN R -e "remotes::install_github('mrcieu/gwasvcf')"

