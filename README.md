# eQTL-Catalogue/colocalisation
**This pipeline runs colocalisation analysis using summary statistics from the eQTL Catalogue (or any summary statistics in the same format) and GWAS summary statistics in VCF format downloaded from the [MRC OpenGWAS database](https://gwas.mrcieu.ac.uk/datasets/).**


### To run in University of Tartu HPC

1. Clone the repo
```
git clone https://github.com/eQTL-Catalogue/colocalisation.git
```

2. Go to [nextflow.config](https://github.com/eQTL-Catalogue/colocalisation/blob/master/nextflow.config) and set the parameters as you need
   
  - _gwas_ss_tsv = "${baseDir}/testdata_coloc/gwas_sumstats_all.tsv"_<br /> // path to the GWASs summary stats files. [see example here](https://github.com/eQTL-Catalogue/colocalisation/blob/master/testdata/gwas_sumstats_test.tsv)<br />  
  - _qtl_ss_tsv = "${baseDir}/testdata_coloc/eqtl_sumstats_tx.tsv"_<br />  // path to the QTL summary statistics [see example here](https://github.com/eQTL-Catalogue/colocalisation/blob/master/testdata/eqtl_sumstats_permuted.tsv)<br />  
  - _gwas_lift_chain = "/gpfs/hpc/projects/eQTLCatalogue/GRCh37_to_GRCh38/GRCh37_to_GRCh38.chain.gz"_<br />  // This is a chain file to lift up the version of GWAS variants from GRCh37 to GRCh38. Don't change it if you don't know what you are doing.<br />
  - _hg38_ref_genome = "/gpfs/hpc/projects/genomic_references/annotations/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"_<br />  // This is a reference genome. Don't change it if you don't know what you are doing.
  - _outdir = './results_coloc_tx'_<br />  // Output directory of the pipeline. The results will be here
  - _use_permutation = false_<br />  // A flag to inform pipeline if you are using lead_var_pairs (credible sets) to do colocalisation in specific pairs of molecular_trait_id and variant_id. If **true** the 4th column of **qtl_ss_tsv** should be permutation run result file of the qtl_subset. If **false** the same column should be the lead_var_pairs file [see example](https://github.com/eQTL-Catalogue/colocalisation/blob/master/testdata_coloc/lead_pairs_rnaseq_2.tsv)
  - _cis_window = 200000_<br />  // cis window where you wanna perform colocalisation. defaults is +-200000 basepairs
  - _n_batches = 10_<br />  // Number of batches needed to split the lead_var_pairs file in processing. I.E. [] has 42445 pairs in it. so if **n_batches**=10 in each job pipeline will process 4245 pairs.

3. start a screen session
```
screen -S coloc_nf
```
4. ssh to stage1
```
ssh stage1
```
5. load the needed modules
```
module load java-1.8.0_40
module load singularity/3.5.3
module load nextflow
```
6. Change directory to where you cloned the repo
```
cd colocalisation/
```

7. run the pipeline
```
nextflow run main.nf -profile tartu_hpc
```

### Credits
[eQTL-Catalogue/colocalisation](https://github.com/eQTL-Catalogue/colocalisation) was originally written by [Nurlan Kerimov](https://github.com/kerimoff) under supervision of [Kaur Alasoo](https://github.com/kauralasoo)
