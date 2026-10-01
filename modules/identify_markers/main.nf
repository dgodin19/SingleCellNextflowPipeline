process IDENTIFY_MARKERS {
	input:
	tuple val(sample_id), path(input_h5ad)

	output:
	tuple val(sample_id), path("${sample_id}_markers.h5ad"), emit: adata
	tuple val(sample_id), path("${sample_id}_markers.tsv"), emit: markers

	script:
	"""
	identify_markers.py "${sample_id}" "${input_h5ad}" "${sample_id}_markers.h5ad" "${sample_id}_markers.tsv" \\
		--n-top-markers ${params.identify_markers_n_top_markers}
	"""
}
