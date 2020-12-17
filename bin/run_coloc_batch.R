message(" ## Loading libraries: dplyr, readr, coloc, GenomicRanges, Rsamtools, optparse, gwasvcf")
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("coloc"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("seqminer"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("gwasvcf"))

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
    gwas_sumstats = "ebi-a-GCST004133-af.GRCh38.sorted.vcf.gz",
    gwas_id = "ebi-a-GCST004133",
    qtl_sumstats = "Alasoo_2018.macrophage_naive_txrev.nominal.sorted.tsv.gz",
    qtl_subset ="Alasoo_2018.macrophage_naive_txrev",
    lead_pairs = "Alasoo_2018.macrophage_naive_txrev.leadpairs.tsv",
    window_coloc = 200000,
    chunk = "4 10",
    output_prefix = "ebi-a-GCST002318_Alasoo_2018.macrophage_naive_tx_4_10.tsv",
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

# In eQTL Catalogue, **variants with multiple rsids are split over multiple rows** in the summary statistics files. 
# Thus, we first want to retain only one unique record per variant. 
# To simplify colocalisation analysis, we also want to exclude multi-allelic variants. 
# The following function imports summary statistics from a tabix-index TSV file and performs necessary filtering.
import_eQTLCatalogue <- function(ftp_path, region, selected_molecular_trait_id, column_names, verbose = TRUE){
  if(verbose){
    print(ftp_path)
  }
  
  #Fetch summary statistics with seqminer
  fetch_table = seqminer::tabix.read.table(tabixFile = ftp_path, tabixRange = region, stringsAsFactors = FALSE) %>%
    dplyr::as_tibble()
  colnames(fetch_table) = column_names
  
  if (is.null(fetch_table)) {
    return(NULL)
  }
  summary_stats = fetch_table %>% 
    dplyr::filter(molecular_trait_id == selected_molecular_trait_id)
  
  #Remove rsid duplicates and multi-allelic variant
  summary_stats = dplyr::select(summary_stats, -rsid) %>% 
    dplyr::distinct() %>% #rsid duplicates
    dplyr::mutate(id = paste(chromosome, position, sep = ":")) %>% 
    dplyr::group_by(id) %>% 
    dplyr::mutate(row_count = n()) %>% dplyr::ungroup() %>% 
    dplyr::filter(row_count == 1) %>% #Multialllics
    dplyr::filter(!is.nan(se)) %>% # remove variants with unknown SE
    dplyr::filter(!is.na(se)) %>% 
    dplyr::select(-row_count) # remove added row_count columns

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
                      MAF = gwas_sumstats$maf, 
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
  
  #Construct character string for the region
  region_string = paste0(var_chrom, ":", as.integer(max(0, var_pos - coloc_window)), "-", as.integer(var_pos + coloc_window))
  
  eqtl_ss_df <- import_eQTLCatalogue(ftp_path = eqtl_ss, 
                                     region = region_string, 
                                     selected_molecular_trait_id = molecular_trait_id, 
                                     column_names = col_names_eqtl)
  if (is.null(eqtl_ss_df) || nrow(eqtl_ss_df) == 0) {
    return(data.frame())
  }
  eqtl_maf_df = dplyr::select(eqtl_ss_df, id, maf)
  
  gwas_stats = gwasvcf::query_gwas(gwas_ss, chrompos = region_string)
  # return empty dataframe if there is no any variants in the given region of GWAS SS
  if(nrow(gwas_stats) == 0){
    return(data.frame())
  }
  gwas_stats = gwasvcf::vcf_to_granges(gwas_stats) %>% 
    keepSeqlevels(var_chrom)
  
  gwas_stats <- gwas_stats %>%
    dplyr::as_tibble() %>%
    dplyr::transmute(chromosome = seqnames, position = start, ES, SE, LP, SS) %>%
    dplyr::mutate(id = paste(chromosome, position, sep = ":")) %>%
    dplyr::group_by(id) %>% #Keep bi-alleilic variants
    dplyr::mutate(row_count = n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(row_count == 1) %>%
    dplyr::left_join(eqtl_maf_df, by = "id") %>%
    dplyr::filter(!is.na(maf))

  # return an empty dataframe if GWAS and QTL SummStats does not share any variant.
  if (sum(gwas_stats$id %in% eqtl_ss_df$id) == 0) {
    return(data.frame())
  }
  
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

#Read column names from the eQTL summary file
qtl_header = read.table(qtl_sumstats, nrows = 1, header = T, sep = "\t")
col_names_eqtl <- colnames(qtl_header)

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
                                        col_names_eqtl = col_names_eqtl,
                                        gwas_id = gwas_id,
                                        qtl_subset = qtl_subset))

if (!dir.exists(outdir)) dir.create(outdir)

if (!is.null(output_prefix)) {
  file_name = file.path(outdir, output_prefix)
} else {
  file_name = file.path(outdir, paste(gwas_id, qtl_subset, chunk_id, n_chunks, sep = "_")) %>% paste0(".tsv")
}

# If there are colocalisation results then write it into the file, if not create an empty file
if(!is.na(coloc_results) && nrow(coloc_results) > 0){
  message(" ## write colocalisation results to ", file_name )
  coloc_results = coloc_results %>% select(c("gwas_id", "qtl_subset", "variant", "molecular_trait_id", "chromosome", "position"), everything())
  utils::write.table(coloc_results, file = file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
} else {
  file.create(file_name)
  message("Coloc results are empty or null!")
}
