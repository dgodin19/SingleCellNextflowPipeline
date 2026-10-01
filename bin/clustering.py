#!/usr/bin/env python3

import argparse
import math
from pathlib import Path

import scanpy as sc


def cluster_cells(
    sample_id: str,
    input_path: Path,
    output_path: Path,
    assignments_path: Path,
    resolution: float,
    random_seed: int,
) -> int:
    if not math.isfinite(resolution) or resolution <= 0:
        raise ValueError('resolution must be a finite number greater than zero.')
    if not 0 <= random_seed <= 2**32 - 1:
        raise ValueError('random_seed must be between 0 and 4294967295.')

    adata = sc.read_h5ad(input_path)
    if adata.n_obs < 2:
        raise ValueError('Clustering requires at least two cells.')
    if 'neighbors' not in adata.uns or 'connectivities' not in adata.obsp:
        raise ValueError(
            'Input data must contain a Scanpy neighbor graph from dimensionality reduction.'
        )

    adata.obs['sample_id'] = sample_id
    sc.tl.leiden(
        adata,
        resolution=resolution,
        random_state=random_seed,
        key_added='leiden',
        flavor='leidenalg',
    )
    adata.uns['clustering'] = {
        'method': 'leiden',
        'resolution': resolution,
        'random_seed': random_seed,
        'cluster_key': 'leiden',
    }

    assignments = adata.obs[['sample_id', 'leiden']].copy()
    assignments.insert(1, 'cell_id', adata.obs_names.to_numpy())
    assignments.to_csv(assignments_path, sep='\t', index=False)
    adata.write_h5ad(output_path)
    return int(adata.obs['leiden'].nunique())


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Cluster single-cell data with the Leiden algorithm.'
    )
    parser.add_argument('sample_id')
    parser.add_argument('input_h5ad', type=Path)
    parser.add_argument('output_h5ad', type=Path)
    parser.add_argument('assignments_tsv', type=Path)
    parser.add_argument('--resolution', type=float, default=1.0)
    parser.add_argument('--random-seed', type=int, default=0)
    args = parser.parse_args()

    cluster_count = cluster_cells(
        args.sample_id,
        args.input_h5ad,
        args.output_h5ad,
        args.assignments_tsv,
        args.resolution,
        args.random_seed,
    )
    print(f'{args.sample_id}: assigned cells to {cluster_count} Leiden clusters.')


if __name__ == '__main__':
    main()