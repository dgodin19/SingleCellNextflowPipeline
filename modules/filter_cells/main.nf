process FILTER_CELLS {
	input:
	tuple val(sample_id), path(input_h5ad)

	output:
	tuple val(sample_id), path("${sample_id}_filtered.h5ad"), emit: adata
	tuple val(sample_id), path("${sample_id}_filtered_qc_metrics.csv"), emit: metrics

	script:
	"""
	filter_cells.py "${sample_id}" "${input_h5ad}" "${sample_id}_filtered.h5ad" "${sample_id}_filtered_qc_metrics.csv" \\
		--min-genes ${params.filter_min_genes} \\
		--max-genes ${params.filter_max_genes} \\
		--max-pct-mt ${params.filter_max_pct_mt}
	"""
}
