process NORMALIZE {
	input:
	tuple val(sample_id), path(input_h5ad)

	output:
	tuple val(sample_id), path("${sample_id}_normalized.h5ad"), emit: adata

	script:
	"""
	normalize.py "${sample_id}" "${input_h5ad}" "${sample_id}_normalized.h5ad" \\
		--target-sum ${params.normalization_target_sum}
	"""
}
