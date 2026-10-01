#!/usr/bin/env python3

import argparse
from pathlib import Path

import celltypist
from celltypist import models
import scanpy as sc


def annotate_cells(
    sample_id: str,
    input_path: Path,
    output_path: Path,
    model_name: str,
) -> None:
    adata = sc.read_h5ad(input_path)
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError('CellTypist requires at least one cell and one gene.')

    models.download_models(model=model_name)
    predictions = celltypist.annotate(adata, model=model_name)
    labels = predictions.predicted_labels['predicted_labels'].reindex(adata.obs_names)
    if labels.isna().any():
        raise ValueError('CellTypist did not return a prediction for every cell.')

    adata.obs['cell_type'] = labels.astype(str)
    adata.uns['celltypist'] = {
        'model': model_name,
        'label_key': 'cell_type',
    }
    adata.write_h5ad(output_path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Annotate cells with a CellTypist reference model.'
    )
    parser.add_argument('sample_id')
    parser.add_argument('input_h5ad', type=Path)
    parser.add_argument('output_h5ad', type=Path)
    parser.add_argument('--model', default='Immune_All_High.pkl')
    args = parser.parse_args()

    annotate_cells(args.sample_id, args.input_h5ad, args.output_h5ad, args.model)
    print(f'{args.sample_id}: annotated cells with {args.model}.')


if __name__ == '__main__':
    main()