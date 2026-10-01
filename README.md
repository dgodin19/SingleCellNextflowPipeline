# Single-Cell Nextflow Pipeline

A lightweight Nextflow pipeline for processing 10x Genomics single-cell RNA-seq data, from input validation and quality control through clustering, cell-type annotation, marker-gene identification, and reporting.

## Directory Structure

```text
.
├── Dockerfile
├── assets
│   └── sample_metadata.csv
├── bin
│   ├── annotate_cell_types.py
│   ├── calculate_qc.py
│   ├── clustering.py
│   ├── dimensionality_reduction.py
│   ├── feature_selection.py
│   ├── filter_cells.py
│   ├── generate_report.py
│   ├── identify_markers.py
│   ├── normalize.py
│   ├── validate_10x_matrix.py
│   └── validate_input.py
├── data
│   ├── pbmc3k_filtered_gene_bc_matrices.tar.gz
│   └── sc5p_v2_hs_PBMC_10k_filtered_feature_bc_matrix.tar.gz
├── main.nf
├── modules
│   ├── calculate_qc
│   │   └── main.nf
│   ├── celltypist
│   │   └── main.nf
│   ├── clustering
│   │   └── main.nf
│   ├── dimensionality_reduction
│   │   └── main.nf
│   ├── feature_selection
│   │   └── main.nf
│   ├── filter_cells
│   │   └── main.nf
│   ├── generate_report
│   │   └── main.nf
│   ├── identify_markers
│   │   └── main.nf
│   ├── normalize
│   │   └── main.nf
│   ├── validate_input
│   │   └── main.nf
│   └── validate_metadata
│       └── main.nf
└── nextflow.config
```

## Workflow

```text
Input Validation
      ↓
Quality Control
      ↓
Cell Filtering
      ↓
Normalization
      ↓
Cell-Type Annotation
      ↓
Feature Selection
      ↓
Dimensionality Reduction
      ↓
Clustering
      ↓
Visualization
      ↓
Differential Marker Identification
      ↓
Reporting
```

Python dependencies are isolated in a Docker container and specified in the `Dockerfile`.

## Input Data

Example datasets that were used in the `data` directory:

* **10x Genomics Human PBMC 10k** — 10,000 PBMCs from a healthy donor. [10x Genomics PBMC 10k](https://www.10xgenomics.com/datasets/human-pbmc-from-a-healthy-donor-10-k-cells-v-2-2-standard-5-0-0)
* **10x Genomics PBMC 3k** — approximately 3,000 PBMCs from a healthy donor. [10x Genomics PBMC 3k dataset](https://www.10xgenomics.com/datasets/3-k-pbm-cs-from-a-healthy-donor-1-standard-1-1-0)

The pipeline uses a CSV metadata file to specify the input samples and their corresponding matrix paths. An example metadata file is provided at:

```text
assets/sample_metadata.csv
```

The metadata file contains the sample identifier, path to the 10x expression matrix, and experimental condition.

## Outputs

The pipeline will output the UMAP, marker genes, and the QC distribution. 

## Running the Pipeline

Build the Docker image:

```bash
docker build -t sc-nextflow:1.0 .
```

Run the pipeline with the Docker profile:

```bash
nextflow run main.nf -profile docker
```




