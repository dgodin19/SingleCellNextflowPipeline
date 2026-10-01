process DIMENSIONALITY_REDUCTION {
	input:
	tuple val(sample_id), path(input_h5ad)

	output:
	tuple val(sample_id), path("${sample_id}_reduced.h5ad"), emit: adata

	script:
	"""
	dimensionality_reduction.py "${sample_id}" "${input_h5ad}" "${sample_id}_reduced.h5ad" \\
		--n-pcs ${params.dimensionality_reduction_n_pcs} \\
		--n-neighbors ${params.dimensionality_reduction_n_neighbors} \\
		--umap-min-dist ${params.dimensionality_reduction_umap_min_dist}
	"""
}
