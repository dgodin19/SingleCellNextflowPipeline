process CELLTYPIST {
	input:
	tuple val(sample_id), path(input_h5ad)

	output:
	tuple val(sample_id), path("${sample_id}_celltyped.h5ad"), emit: adata

	script:
	"""
	annotate_cell_types.py "${sample_id}" "${input_h5ad}" "${sample_id}_celltyped.h5ad" --model "${params.celltypist_model}"
	"""
}