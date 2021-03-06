# Convert GWAS summary statistics to VCF format with [gwas2vcf](https://github.com/MRCIEU/gwas2vcf) tool

### Step 1: Clone the gwas2vcf repository

```bash
git clone git@github.com:MRCIEU/gwas2vcf.git
cd gwas2vcf
```

### Step 2: Build the singularity container from Quay.io:

```bash
module load singularity/3.5.3
singularity build gwas2vcf.img docker://quay.io/eqtlcatalogue/gwas2vcf:v20.11.1
```

### Step 3: Run gwas2vcf using the sinularity container

```bash
singularity exec -B /gpfs/:/gpfs/ gwas2vcf.img python main.py\
    --out ../LDLC/test.vcf.gz\
    --data ../continuous-LDLC-both_sexes-medadj_irnt.tsv\
    --ref human_g1k_v37.fasta\
    --id LDLC\
    --cohort_controls 457526\
    --json params.json
```

Example params.json file for the Pan-UKBB data can be found [here](params.json). 
