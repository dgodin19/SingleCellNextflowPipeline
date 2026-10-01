process GENERATE_REPORT {
	publishDir path: { "${params.outdir}/${sample_id}/figures" }, mode: 'copy', pattern: '*.png'
	publishDir path: { "${params.outdir}/${sample_id}/tables" }, mode: 'copy', pattern: '*.tsv'
	publishDir path: { "${params.outdir}/${sample_id}/tables" }, mode: 'copy', pattern: '*.csv'

	input:
	tuple val(sample_id), path(qc_metrics), path(features), path(input_h5ad), path(markers)

	output:
	tuple val(sample_id), path("${sample_id}_qc_distributions.png"), path("${sample_id}_umap_clusters.png"), path("${sample_id}_top_markers.png"), emit: figures
	tuple val(sample_id), path("${sample_id}_qc_metrics.csv"), path("${sample_id}_highly_variable_genes.tsv"), path("${sample_id}_markers.tsv"), path("${sample_id}_cluster_counts.tsv"), path("${sample_id}_summary.tsv"), emit: tables

	script:
	"""
	generate_report.py "${sample_id}" "${qc_metrics}" "${features}" "${input_h5ad}" "${markers}" \\
		--output-prefix "${sample_id}"
	"""
}
