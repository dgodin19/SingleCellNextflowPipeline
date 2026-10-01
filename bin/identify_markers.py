#!/usr/bin/env python3

import argparse
from pathlib import Path

import scanpy as sc


def identify_markers(
    sample_id: str,
    input_path: Path,
    output_path: Path,
    markers_path: Path,
    n_top_markers: int,
) -> int:
    if n_top_markers < 1:
        raise ValueError('n_top_markers must be at least 1.')

    adata = sc.read_h5ad(input_path)
    if adata.n_obs < 2 or adata.n_vars == 0:
        raise ValueError('Marker identification requires at least two cells and one gene.')
    if 'leiden' not in adata.obs:
        raise ValueError('Input data must contain Leiden cluster labels in adata.obs["leiden"].')
    if adata.obs['leiden'].isna().any():
        raise ValueError('Leiden cluster labels must not contain missing values.')
    if adata.obs['leiden'].nunique() < 2:
        raise ValueError('Marker identification requires at least two clusters.')

    n_genes = min(n_top_markers, adata.n_vars)
    adata.obs['sample_id'] = sample_id
    sc.tl.rank_genes_groups(
        adata,
        groupby='leiden',
        method='wilcoxon',
        n_genes=n_genes,
        use_raw=False,
        key_added='rank_genes_groups',
    )

    markers = sc.get.rank_genes_groups_df(
        adata,
        group=None,
        key='rank_genes_groups',
    )
    markers.insert(0, 'sample_id', sample_id)
    markers.insert(
        2,
        'rank',
        markers.groupby('group', sort=False).cumcount() + 1,
    )
    markers = markers.rename(columns={'group': 'cluster', 'names': 'gene'})
    markers.to_csv(markers_path, sep='\t', index=False)

    adata.uns['marker_identification'] = {
        'method': 'wilcoxon',
        'groupby': 'leiden',
        'n_top_markers': n_genes,
        'expression_source': 'X',
    }
    adata.write_h5ad(output_path)
    return len(markers)


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Rank marker genes for each Leiden cluster with the Wilcoxon test.'
    )
    parser.add_argument('sample_id')
    parser.add_argument('input_h5ad', type=Path)
    parser.add_argument('output_h5ad', type=Path)
    parser.add_argument('markers_tsv', type=Path)
    parser.add_argument('--n-top-markers', type=int, default=100)
    args = parser.parse_args()

    marker_count = identify_markers(
        args.sample_id,
        args.input_h5ad,
        args.output_h5ad,
        args.markers_tsv,
        args.n_top_markers,
    )
    print(f'{args.sample_id}: ranked {marker_count} cluster marker entries.')


if __name__ == '__main__':
    main()