#Summary tables that may be useful for reference
Amillepora_euk.emapper.annotations.tsv	--	annotations used in this study
all_deseqResults.tsv	--	single table of all DESeq results used in this study compiled
all_high_stress_go_mwu_results.tsv	--	gene ontology enrichment results for all type A stress samples together
all_zoox_type_counts.tsv	--	fold coverages for symbiont types. This was not used in publication but kept for reference
bioproject_stress_correlation_matrix.Rdata	--	correlation matrix for BioProjects' log2 fold change due to stress treatment
contextual_annots.tsv	--	table of contextual annotations labeling candidate genes for the general coral stress response based on consistent up- or downreulgation among the BioProjects
final_raw_gene_counts_matrix.tsv.zip	--	raw gene counts matrix (counted using featureCounts). This is the starting point for all the analyses after processing of the RNAseq reads
gene_contibutions.tsv	--	table of gene's contributions to PCA, DAPC, lasso regression, and random forest models
module_go_mwu_results.tsv	--	gene ontology enrichment for each WGCNA module
normalized_counts_project_controlled.Rdata	--	R object with normalized counts for the stress studies after controlling for BioProject. These were used for analyses of expression level such as PCA, DAPC, and WGCNA
wgcna_moduleMembership.tsv	--	table of module membership for each WGCNA module