process CALCULATE_QC {
	input:
	tuple val(sample_id), path(matrix_input)

	output:
	tuple val(sample_id), path("${sample_id}_qc_metrics.csv"), emit: metrics
	tuple val(sample_id), path("${sample_id}_qc.h5ad"), emit: adata

	script:
	"""
	calculate_qc.py "${sample_id}" "${matrix_input}" "${sample_id}_qc_metrics.csv" "${sample_id}_qc.h5ad"
	"""
}
