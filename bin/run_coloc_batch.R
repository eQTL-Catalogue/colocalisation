message(" ## Loading libraries: dplyr, readr, coloc, GenomicRanges, Rsamtools, optparse")
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("coloc"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("Rsamtools"))
suppressPackageStartupMessages(library("optparse"))

#Parse command-line options
option_list <- list(
  make_option(c("--gwas_sumstats"), type="character", default=NULL,
              help="GWAS summary statistics in VCF format (version GRCh38)", metavar = "type"),
  make_option(c("--gwas_id"), type="character", default="gwas_id",
              help="GWAS identifier", metavar = "type"),
  make_option(c("--qtl_sumstats"), type="character", default=NULL,
              help="eQTL Catalogue summary statistics path ", metavar = "type"),
  make_option(c("--qtl_subset"), type="character", default="qtl_subset",
              help="QTL subset identifier", metavar = "type"),
  make_option(c("--lead_pairs"), type="character", default=NULL,
              help="TSV file containing gene_variant pairs, chromosome and position of variant.", metavar = "type"),
  make_option(c("--window_coloc"), type="integer", default=200000,
              help="Colocalisation window from each side of the variant position. [default \"%default\"]", metavar = "type"),
  make_option(c("--chunk"), type="character", default="1 1", 
              help="Perform analysis in chunks. Eg value 5 10 would indicate that phenotypes are split into 10 chunks and the 5th one of those will be processed. [default \"%default\"]", metavar = "type"),
  make_option(c("--output_prefix"), type="character", default=NULL,
              help="Prefix of the output files.", metavar = "type"),
  make_option(c("--outdir"), type="character", default="./coloc_results/",
              help="Path to the output directory. [default \"%default\"]", metavar = "type"))

message(" ## Parsing options")
opt <- optparse::parse_args(OptionParser(option_list=option_list))

#Deubgging
if(FALSE){
  opt = list(
    gwas_sumstats = "ebi-a-GCST004599.GRCh38.sorted.vcf.gz",
    gwas_id = "gwas_id",
    qtl_sumstats = "Alasoo_2018.macrophage_naive_ge.nominal.sorted.tsv.gz",
    qtl_subset ="qtl_subset",
    lead_pairs = "lead_pairs_rnaseq.tsv",
    window_coloc = 200000,
    chunk = "1 100",
    output_prefix = "ebi-a-GCST004599_Alasoo_2018.macrophage_naive_ge_1_100.tsv",
    outdir = "./coloc_results/")
}

gwas_sumstats = opt$gwas_sumstats
gwas_id = opt$gwas_id
qtl_sumstats = opt$qtl_sumstats
qtl_subset = opt$qtl_subset
lead_pairs = opt$lead_pairs
coloc_window = opt$window_coloc
chunk = opt$chunk
output_prefix = opt$output_prefix
outdir = opt$outdir

message("__________ OPTIONS __________")
message("gwas_sumstats: ", gwas_sumstats)
message("gwas_id: ", gwas_id)
message("qtl_sumstats: ", qtl_sumstats)
message("qtl_subset: ", qtl_subset)
message("lead_pairs_file: ", lead_pairs)
message("window_coloc: ", coloc_window)
message("chunk: ", chunk)
message("output_prefix: ", output_prefix)
message("output_file_path: ", outdir)

# Define a small helper function to quickly read regions from tabix-indexed summary statistics files into R.
#' A general function to quickly import tabix indexed tab-separated files into data_frame
#'
#' @param tabix_file Path to tabix-indexed text file
#' @param param An instance of GRanges, RangedData, or RangesList
#' provide the sequence names and regions to be parsed. Passed onto Rsamtools::scanTabix()
#' @param ... Additional parameters to be passed on to readr::read_delim()
#'
#' @return List of data_frames, one for each entry in the param GRanges object.
#' @export
scanTabixDataFrame <- function(tabix_file, param, ...){
  tabix_list = Rsamtools::scanTabix(tabix_file, param = param)
  df_list = lapply(tabix_list, function(x,...){
    if(length(x) > 0){
      if(length(x) == 1){
        #Hack to make sure that it also works for data frames with only one row
        #Adds an empty row and then removes it
        result = paste(paste(x, collapse = "\n"),"\n",sep = "")
        result = readr::read_delim(result, delim = "\t", ...)[1,]
      }else{
        result = paste(x, collapse = "\n")
        result = readr::read_delim(result, delim = "\t", ...)
      }
    } else{
      #Return NULL if the nothing is returned from tabix file
      result = NULL
    }
    return(result)
  }, ...)
  return(df_list)
}

# In eQTL Catalogue, **variants with multiple rsids are split over multiple rows** in the summary statistics files. 
# Thus, we first want to retain only one unique record per variant. 
# To simplify colocalisation analysis, we also want to exclude multi-allelic variants. 
# The following function imports summary statistics from a tabix-index TSV file and performs necessary filtering.
import_eQTLCatalogue <- function(ftp_path, region, selected_molecular_trait_id, column_names, verbose = TRUE){
  if(verbose){
    print(ftp_path)
  }
  
  #Fetch summary statistics with Rsamtools
  summary_stats = scanTabixDataFrame(ftp_path, region, col_names = column_names)[[1]] 
  if (is.null(summary_stats)) {
    return(NULL)
  }
  summary_stats = summary_stats %>% 
    dplyr::filter(molecular_trait_id == selected_molecular_trait_id)
  
  #Remove rsid duplicates and multi-allelic variant
  summary_stats = dplyr::select(summary_stats, -rsid) %>% 
    dplyr::distinct() %>% #rsid duplicates
    dplyr::mutate(id = paste(chromosome, position, sep = ":")) %>% 
    dplyr::group_by(id) %>% 
    dplyr::mutate(row_count = n()) %>% dplyr::ungroup() %>% 
    dplyr::filter(row_count == 1) #Multialllics
  
  return(summary_stats)
}

run_coloc <- function(eqtl_sumstats, gwas_sumstats){
  eQTL_dataset = list(varbeta = eqtl_sumstats$se^2, 
                      N = (eqtl_sumstats$an)[1]/2, # Samples size is allele number (AN) dvided by 2
                      MAF = eqtl_sumstats$maf, 
                      type = "quant", 
                      beta = eqtl_sumstats$beta,
                      snp = eqtl_sumstats$id)
  gwas_dataset = list(beta = gwas_sumstats$ES,
                      varbeta = gwas_sumstats$SE^2, 
                      type = "quant", 
                      snp = gwas_sumstats$id,
                      MAF = gwas_sumstats$MAF, 
                      N = gwas_sumstats$SS)
  coloc_res = coloc::coloc.abf(dataset1 = eQTL_dataset, dataset2 = gwas_dataset,p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  res_formatted = dplyr::as_tibble(t(as.data.frame(coloc_res$summary)))
  return(res_formatted)
}

splitIntoBatches <- function(n, batch_size){
  n_batches = ceiling(n/batch_size)
  batch_ids = rep(seq(1:n_batches), each = batch_size)[1:n]
  return(batch_ids)
}

splitIntoChunks <- function(chunk_number, n_chunks, n_total){
  chunk_size = floor(n_total/(n_chunks))
  batches = splitIntoBatches(n_total,chunk_size)
  batches[batches > n_chunks] = n_chunks
  selected_batch = batches == chunk_number
  return(selected_batch)
}

col_names_eqtl <- c("molecular_trait_id","chromosome","position","ref","alt","variant","ma_samples","ac","an","maf","pvalue","beta","se","molecular_trait_object_id","gene_id","median_tpm","r2","type","rsid")

# Define a function which performs coloc between variant and phenotype in given window around the variant position
#' @param pair One row of dataframe with [molecular_trait_id,variant,chromosome,position] columns
#' @param gwas_ss GWAS vcf file. Should have a tabix index file in the same path
#' @param eqtl_ss eQTL summary statistics in archived tsv format. Should have a tabix index file in the same path
#' @param coloc_window Window size around given variant to perform colocalisation
#'
#' @return one row of tibble which will be rbinded after apply function finishes the process.
#' @export
coloc_in_region <- function(pair, gwas_ss, eqtl_ss, coloc_window, col_names_eqtl, gwas_id, qtl_subset){
  print(pair)
  var_id = as.character(pair["variant"]) %>% trimws()
  var_pos = as.numeric(pair["position"])
  var_chrom = as.character(pair["chromosome"]) %>% trimws()
  molecular_trait_id = as.character(pair["molecular_trait_id"])
  region_granges = GenomicRanges::GRanges(
    seqnames = var_chrom, 
    ranges = IRanges::IRanges(start = max(0, var_pos - coloc_window), end = var_pos + coloc_window), 
    strand = "*")
  
  eqtl_ss_df <- import_eQTLCatalogue(eqtl_ss, region_granges, selected_molecular_trait_id = molecular_trait_id, col_names_eqtl)
  if (is.null(eqtl_ss_df) || nrow(eqtl_ss_df) == 0) {
    return(data.frame())
  }
  
  gwas_stats = gwasvcf::query_gwas(gwas_ss, chrompos = paste0(var_chrom, ":", max(0, var_pos - coloc_window), "-", var_pos + coloc_window))
  gwas_stats = gwasvcf::vcf_to_granges(gwas_stats) %>% 
    keepSeqlevels(var_chrom)
  
  gwas_stats <- gwas_stats %>%
    dplyr::as_tibble() %>%
    dplyr::transmute(chromosome = seqnames, position = start, AF, ES, SE, LP, SS) %>%
    dplyr::mutate(id = paste(chromosome, position, sep = ":")) %>%
    dplyr::mutate(MAF = pmin(AF, 1-AF)) %>% #Calculate MAF
    dplyr::group_by(id) %>% #Keep bi-alleilic variants
    dplyr::mutate(row_count = n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(row_count == 1) 
  
  results = run_coloc(eqtl_ss_df, gwas_stats)
  results = results %>% mutate(variant = var_id, 
                     molecular_trait_id = molecular_trait_id, 
                     chromosome = var_chrom,
                     position=var_pos,
                     gwas_id = gwas_id,
                     qtl_subset = qtl_subset)

  return(results)
}

# Read molecular_trait_id variant pairs
pheno_var_df <- read_tsv(lead_pairs)

#Split phenotype list into chunks
chunk_vector = strsplit(chunk, split = " ") %>% unlist() %>% as.numeric()
chunk_id = chunk_vector[1]
n_chunks = chunk_vector[2]
selected_chunk = splitIntoChunks(chunk_id, n_chunks, nrow(pheno_var_df))
selected_pairs = pheno_var_df[selected_chunk,]

coloc_results <- do.call("rbind", apply(selected_pairs, 1, coloc_in_region, 
                                        gwas_ss = gwas_sumstats, 
                                        eqtl_ss = qtl_sumstats, 
                                        coloc_window = coloc_window, 
                                        col_names_eqtl=col_names_eqtl,
                                        gwas_id = gwas_id,
                                        qtl_subset = qtl_subset))

# If there are colocalisation results then write it into the file
if(!is.na(coloc_results) && nrow(coloc_results) > 0){
  if (!dir.exists(outdir)) dir.create(outdir)
  if (!is.null(output_prefix)) {
    file_name = file.path(outdir, output_prefix)
  } else {
    file_name = file.path(outdir, paste(gwas_id, qtl_subset, chunk_id, n_chunks, sep = "_")) %>% paste0(".tsv")
  }
  message(" ## write colocalisation results to ", file_name )
  coloc_results = coloc_results %>% select(c("gwas_id", "qtl_subset", "variant", "molecular_trait_id", "chromosome", "position"), everything())
  utils::write.table(coloc_results, file = file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
} else {
  message("Coloc results are empty or null!")
}

























