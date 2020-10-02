FROM nfcore/base
LABEL authors="Kaur Alasoo" \
      description="Docker image containing all requirements for the eQTL Catalogue colocalisation workflow"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/colocalisation-1.0dev/bin:$PATH
