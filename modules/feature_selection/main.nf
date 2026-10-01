process FEATURE_SELECTION {
	input:
	tuple val(sample_id), path(input_h5ad)

	output:
	tuple val(sample_id), path("${sample_id}_feature_selected.h5ad"), emit: adata
	tuple val(sample_id), path("${sample_id}_highly_variable_genes.tsv"), emit: features

	script:
	"""
	feature_selection.py "${sample_id}" "${input_h5ad}" "${sample_id}_feature_selected.h5ad" "${sample_id}_highly_variable_genes.tsv" \\
		--n-top-genes ${params.feature_selection_n_top_genes}
	"""
}