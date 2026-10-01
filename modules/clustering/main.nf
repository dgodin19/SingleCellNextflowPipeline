process CLUSTERING {
	input:
	tuple val(sample_id), path(input_h5ad)

	output:
	tuple val(sample_id), path("${sample_id}_clustered.h5ad"), emit: adata
	tuple val(sample_id), path("${sample_id}_cluster_assignments.tsv"), emit: assignments

	script:
	"""
	clustering.py "${sample_id}" "${input_h5ad}" "${sample_id}_clustered.h5ad" "${sample_id}_cluster_assignments.tsv" \\
		--resolution ${params.clustering_resolution} \\
		--random-seed ${params.clustering_random_seed}
	"""
}
