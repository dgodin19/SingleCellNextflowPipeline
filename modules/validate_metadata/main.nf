#!/usr/bin/env nextflow

process VALIDATE_METADATA {
	input:
	path metadata_csv, stageAs: 'samplesheet.csv'

	output:
	path 'validated_metadata.csv', emit: metadata

	script:
	"""
	validate_input.py samplesheet.csv validated_metadata.csv
	"""
}
