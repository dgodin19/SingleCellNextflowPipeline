#!/usr/bin/env python3

import argparse
import math
from pathlib import Path

import scanpy as sc


def reduce_dimensions(
    sample_id: str,
    input_path: Path,
    output_path: Path,
    n_pcs: int,
    n_neighbors: int,
    umap_min_dist: float,
) -> None:
    if n_pcs < 2:
        raise ValueError('n_pcs must be at least 2.')
    if n_neighbors < 2:
        raise ValueError('n_neighbors must be at least 2.')
    if not math.isfinite(umap_min_dist) or not 0 <= umap_min_dist <= 1:
        raise ValueError('umap_min_dist must be a finite value between 0 and 1.')

    adata = sc.read_h5ad(input_path)
    if adata.n_obs < 3:
        raise ValueError('Dimensionality reduction requires at least three cells.')
    if 'highly_variable' not in adata.var:
        raise ValueError('Input data must contain highly_variable gene annotations.')

    n_highly_variable = int(adata.var['highly_variable'].sum())
    if n_highly_variable < n_pcs + 1:
        raise ValueError(
            f'n_pcs ({n_pcs}) must be less than the number of highly variable '
            f'genes ({n_highly_variable}).'
        )
    if n_pcs >= adata.n_obs:
        raise ValueError(
            f'n_pcs ({n_pcs}) must be less than the number of cells ({adata.n_obs}).'
        )
    if n_neighbors >= adata.n_obs:
        raise ValueError(
            f'n_neighbors ({n_neighbors}) must be less than the number of '
            f'cells ({adata.n_obs}).'
        )

    adata.obs['sample_id'] = sample_id
    sc.pp.pca(adata, n_comps=n_pcs, use_highly_variable=True, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(adata, min_dist=umap_min_dist)
    adata.uns['dimensionality_reduction'] = {
        'pca': 'highly_variable_genes',
        'n_pcs': n_pcs,
        'n_neighbors': n_neighbors,
        'umap_min_dist': umap_min_dist,
    }
    adata.write_h5ad(output_path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Compute PCA, a neighbor graph, and UMAP for single-cell data.'
    )
    parser.add_argument('sample_id')
    parser.add_argument('input_h5ad', type=Path)
    parser.add_argument('output_h5ad', type=Path)
    parser.add_argument('--n-pcs', type=int, default=30)
    parser.add_argument('--n-neighbors', type=int, default=15)
    parser.add_argument('--umap-min-dist', type=float, default=0.5)
    args = parser.parse_args()

    reduce_dimensions(
        args.sample_id,
        args.input_h5ad,
        args.output_h5ad,
        args.n_pcs,
        args.n_neighbors,
        args.umap_min_dist,
    )
    print(f'{args.sample_id}: computed PCA, neighbors, and UMAP.')


if __name__ == '__main__':
    main()