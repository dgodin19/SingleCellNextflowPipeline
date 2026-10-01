#!/usr/bin/env python3

import argparse
import math
from pathlib import Path


def normalize_adata(
	sample_id: str,
	input_path: Path,
	output_path: Path,
	target_sum: float,
) -> None:
	import scanpy as sc

	if not math.isfinite(target_sum) or target_sum <= 0:
		raise ValueError('target_sum must be a finite number greater than zero.')

	adata = sc.read_h5ad(input_path)
	if adata.n_obs == 0 or adata.n_vars == 0:
		raise ValueError('Cannot normalize an AnnData object with no cells or genes.')

	if 'counts' not in adata.layers:
		adata.layers['counts'] = adata.X.copy()

	adata.obs['sample_id'] = sample_id
	sc.pp.normalize_total(adata, target_sum=target_sum, inplace=True)
	sc.pp.log1p(adata)
	adata.uns['normalization'] = {
		'method': 'total_counts_then_log1p',
		'target_sum': target_sum,
		'counts_layer': 'counts',
	}
	adata.write_h5ad(output_path)


def main() -> None:
	parser = argparse.ArgumentParser(
		description='Normalize filtered single-cell counts with Scanpy.'
	)
	parser.add_argument('sample_id')
	parser.add_argument('input_h5ad', type=Path)
	parser.add_argument('output_h5ad', type=Path)
	parser.add_argument('--target-sum', type=float, default=10_000)
	args = parser.parse_args()

	normalize_adata(
		args.sample_id,
		args.input_h5ad,
		args.output_h5ad,
		args.target_sum,
	)
	print(
		f'{args.sample_id}: normalized to {args.target_sum:g} counts per cell '
		'and applied log1p.'
	)


if __name__ == '__main__':
	main()
