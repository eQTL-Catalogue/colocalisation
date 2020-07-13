#!/usr/bin/env nextflow

Channel.fromPath(params.gwas_ss_tsv)
    .ifEmpty { error "Cannot find gwas_vcf: ${params.gwas_vcf}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ row.gwas_id, file(row.gwas_vcf)]}
    .set { gwas_vcf_ch }

Channel.fromPath(params.gwas_lift_chain)
    .ifEmpty { error "Cannot find gwas_vcf: ${params.gwas_lift_chain}" }
    .set { gwas_lift_chain_ch }

Channel.fromPath(params.hg38_ref_genome)
    .ifEmpty { error "Cannot find gwas_vcf: ${params.hg38_ref_genome}" }
    .set { hg38_ref_genome_ch }

Channel.fromPath(params.gene_variant_list)
    .ifEmpty { error "Cannot find gene_variant_list: ${params.gene_variant_list}" }
    .set { gene_variant_list_ch }

Channel.fromPath(params.qtl_ss_tsv)
    .ifEmpty { error "Cannot find gene_variant_list: ${params.eqtl_summ_stats_path}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ row.qtl_subset, file(row.qtl_ss), file(row.qtl_ss_index)]}
    .set { eqtl_summ_stats_ch }

process lift_to_GRCh38{
    tag "${gwas_id}"
    publishDir "${params.outdir}/GRCh38_conv/", mode: 'copy'
    container 'crukcibioinformatics/crossmap'

    input:
    tuple val(gwas_id), file(gwas_vcf) from gwas_vcf_ch
    file gwas_lift_chain from gwas_lift_chain_ch.collect()
    file hg38_ref_genome from hg38_ref_genome_ch.collect()

    output:
    tuple val(gwas_id), file("${gwas_vcf.simpleName}.GRCh38.vcf") into tabix_index_ch

    script:
    """
    CrossMap.py vcf $gwas_lift_chain $gwas_vcf $hg38_ref_genome ${gwas_vcf.simpleName}.GRCh38.vcf
    """
}

process tabix_index_gwas{
    tag "${gwas_id}"
    publishDir "${params.outdir}/GRCh38_conv/", mode: 'copy'
    container = 'eqtlcatalogue/qtlmap:latest'

    input:
    tuple val(gwas_id), file(vcf_GRCh38) from tabix_index_ch

    output:
    tuple val(gwas_id), file("${vcf_GRCh38.simpleName}.GRCh38.sorted.vcf.gz"), file("${vcf_GRCh38.simpleName}.GRCh38.sorted.vcf.gz.tbi") into gwas_summstats_GRCh38

    script:
    """
    bcftools sort -o ${vcf_GRCh38.simpleName}.GRCh38.sorted.vcf.gz -O z $vcf_GRCh38
    tabix -p vcf ${vcf_GRCh38.simpleName}.GRCh38.sorted.vcf.gz
    """
}

process run_coloc{
    tag "${gwas_id}_${qtl_subset}"
    // publishDir "${params.outdir}/coloc_results_batch/", mode: 'copy'
    container 'kerimoff/coloc_main:latest'

    input:
    each batch_index from 1..params.n_batches
    tuple val(qtl_subset), file(eqtl_ss), file(eqtl_ss_index), val(gwas_id), file(gwas_sumstats), file(gwas_sumstats_index) from eqtl_summ_stats_ch.combine(gwas_summstats_GRCh38)
    file lead_pairs from gene_variant_list_ch.collect()

    output:
    set val(gwas_id), val("${gwas_id}_${qtl_subset}"), file("${gwas_id}_${qtl_subset}_${batch_index}_${params.n_batches}.tsv") into batch_files_merge_coloc_results

    script:
    """
    Rscript $baseDir/bin/run_coloc_batch.R \\
        --gwas_sumstats $gwas_sumstats \\
        --gwas_id $gwas_id \\
        --qtl_sumstats $eqtl_ss \\
        --qtl_subset $qtl_subset \\
        --lead_pairs $lead_pairs \\
        --window_coloc ${params.cis_window} \\
        --chunk '$batch_index ${params.n_batches}' \\
        --output_prefix ${gwas_id}_${qtl_subset}_${batch_index}_${params.n_batches}.tsv \\
        --outdir .
        
    """
}

process merge_coloc_results{
    publishDir "${params.outdir}/coloc_results_merged/${gwas_id}", mode: 'copy'
    container 'kerimoff/coloc_main:latest'

    input:
    tuple gwas_id, gwas_qtl_subset, file(gwas_qtl_subset_coloc_results_batch_files) from batch_files_merge_coloc_results.groupTuple(by: [1, 0], sort: true)

    output:
    tuple gwas_qtl_subset, file("${gwas_qtl_subset}.txt.gz")

    script:
    """
    awk 'NR == 1 || FNR > 1{print}' ${gwas_qtl_subset_coloc_results_batch_files.join(' ')} | bgzip -c > ${gwas_qtl_subset}.txt.gz
    """
}

workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops ... something went wrong" )
}
