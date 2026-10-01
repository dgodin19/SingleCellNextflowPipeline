#!/usr/bin/env python3

import argparse
from pathlib import Path
import scanpy as sc


def select_features(
    sample_id: str,
    input_path: Path,
    output_path: Path,
    features_path: Path,
    n_top_genes: int,
) -> int:
   

    if n_top_genes < 1:
        raise ValueError('n_top_genes must be at least 1.')

    adata = sc.read_h5ad(input_path)
    if adata.n_obs < 2 or adata.n_vars == 0:
        raise ValueError('Feature selection requires at least two cells and one gene.')
    if n_top_genes > adata.n_vars:
        raise ValueError(
            f'n_top_genes ({n_top_genes}) exceeds the number of genes ({adata.n_vars}).'
        )

    sc.pp.highly_variable_genes(
        adata,
        flavor='seurat',
        n_top_genes=n_top_genes,
        subset=False,
        inplace=True,
    )

    adata.obs['sample_id'] = sample_id
    adata.uns['feature_selection'] = {
        'method': 'highly_variable_genes',
        'flavor': 'seurat',
        'n_top_genes': n_top_genes,
    }

    selected = adata.var.loc[adata.var['highly_variable']].copy()
    selected.index.name = 'gene'
    selected = selected.reset_index()
    selected.insert(0, 'sample_id', sample_id)
    selected.to_csv(features_path, sep='\t', index=False)
    adata.write_h5ad(output_path)
    return selected.shape[0]


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Select highly variable genes from normalized single-cell data.'
    )
    parser.add_argument('sample_id')
    parser.add_argument('input_h5ad', type=Path)
    parser.add_argument('output_h5ad', type=Path)
    parser.add_argument('features_tsv', type=Path)
    parser.add_argument('--n-top-genes', type=int, default=2000)
    args = parser.parse_args()

    selected_count = select_features(
        args.sample_id,
        args.input_h5ad,
        args.output_h5ad,
        args.features_tsv,
        args.n_top_genes,
    )
    print(f'{args.sample_id}: selected {selected_count} highly variable genes.')


if __name__ == '__main__':
    main()