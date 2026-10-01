#!/usr/bin/env nextflow

process VALIDATE_INPUT {
	input:
	tuple val(sample_id), path(matrix_input)

	output:
	tuple val(sample_id), path(matrix_input), emit: validated_matrix
	path 'matrix_validation.txt', emit: report

	script:
	"""
	validate_10x_matrix.py "${matrix_input}" matrix_validation.txt
	"""
}
