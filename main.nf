nextflow.enable.dsl=2

include { VALIDATE_METADATA } from './modules/validate_metadata'
include { VALIDATE_INPUT } from './modules/validate_input'
include { CALCULATE_QC } from './modules/calculate_qc'
include { FILTER_CELLS } from './modules/filter_cells'
include { NORMALIZE } from './modules/normalize'
include { CELLTYPIST } from './modules/celltypist'
include { FEATURE_SELECTION } from './modules/feature_selection'
include { DIMENSIONALITY_REDUCTION } from './modules/dimensionality_reduction'
include { CLUSTERING } from './modules/clustering'
include { IDENTIFY_MARKERS } from './modules/identify_markers'
include { GENERATE_REPORT } from './modules/generate_report'

workflow {
    VALIDATE_METADATA(params.metadata)
    matrix_inputs = VALIDATE_METADATA.out.metadata
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id, file("${projectDir}/${row.matrix_path}")) }
    VALIDATE_INPUT(matrix_inputs)
    CALCULATE_QC(VALIDATE_INPUT.out.validated_matrix)
    FILTER_CELLS(CALCULATE_QC.out.adata)
    NORMALIZE(FILTER_CELLS.out.adata)
    CELLTYPIST(NORMALIZE.out.adata)
    FEATURE_SELECTION(CELLTYPIST.out.adata)
    DIMENSIONALITY_REDUCTION(FEATURE_SELECTION.out.adata)
    CLUSTERING(DIMENSIONALITY_REDUCTION.out.adata)
    IDENTIFY_MARKERS(CLUSTERING.out.adata)
    report_inputs = CALCULATE_QC.out.metrics
        .join(FEATURE_SELECTION.out.features)
        .join(IDENTIFY_MARKERS.out.adata)
        .join(IDENTIFY_MARKERS.out.markers)
    GENERATE_REPORT(report_inputs)
}
