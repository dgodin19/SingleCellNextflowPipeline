#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc


def passing_cell_mask(
	observations: pd.DataFrame,
	min_genes: int,
	max_genes: int,
	max_pct_mt: float,
) -> pd.Series:
	required_metrics = {'n_genes_by_counts', 'pct_counts_mt'}
	missing_metrics = required_metrics.difference(observations.columns)
	if missing_metrics:
		raise ValueError(
			'QC AnnData is missing required metrics: '
			+ ', '.join(sorted(missing_metrics))
		)

	detected_genes = pd.to_numeric(observations['n_genes_by_counts'], errors='coerce')
	mitochondrial_fraction = pd.to_numeric(observations['pct_counts_mt'], errors='coerce')
	return (
		detected_genes.ge(min_genes)
		& detected_genes.le(max_genes)
		& mitochondrial_fraction.le(max_pct_mt)
	)


def filter_cells(
	sample_id: str,
	input_path: Path,
	output_path: Path,
	metrics_path: Path,
	min_genes: int,
	max_genes: int,
	max_pct_mt: float,
) -> tuple[int, int]:
	

	if min_genes < 0 or max_genes < min_genes:
		raise ValueError('Gene thresholds must satisfy 0 <= min_genes <= max_genes.')
	if not 0 <= max_pct_mt <= 100:
		raise ValueError('max_pct_mt must be between 0 and 100.')

	adata = sc.read_h5ad(input_path)
	original_cell_count = adata.n_obs
	keep_cells = passing_cell_mask(adata.obs, min_genes, max_genes, max_pct_mt)
	filtered = adata[keep_cells.to_numpy()].copy()
	if filtered.n_obs == 0:
		raise ValueError(
			f'All {original_cell_count} cells failed QC filters for sample {sample_id}.'
		)

	filtered.obs['sample_id'] = sample_id
	metrics = filtered.obs[['total_counts', 'n_genes_by_counts', 'pct_counts_mt']].copy()
	metrics.insert(0, 'barcode', filtered.obs_names.astype(str))
	metrics.insert(0, 'sample_id', sample_id)
	metrics.to_csv(metrics_path, index=False)
	filtered.write_h5ad(output_path)
	return original_cell_count, filtered.n_obs


def main() -> None:
	parser = argparse.ArgumentParser(description='Filter cells using Scanpy QC metrics.')
	parser.add_argument('sample_id')
	parser.add_argument('input_h5ad', type=Path)
	parser.add_argument('output_h5ad', type=Path)
	parser.add_argument('metrics_csv', type=Path)
	parser.add_argument('--min-genes', type=int, default=200)
	parser.add_argument('--max-genes', type=int, default=2500)
	parser.add_argument('--max-pct-mt', type=float, default=5.0)
	args = parser.parse_args()

	original_cell_count, retained_cell_count = filter_cells(
		args.sample_id,
		args.input_h5ad,
		args.output_h5ad,
		args.metrics_csv,
		args.min_genes,
		args.max_genes,
		args.max_pct_mt,
	)
	print(
		f'{args.sample_id}: retained {retained_cell_count} of '
		f'{original_cell_count} cells.'
	)


if __name__ == '__main__':
	main()
