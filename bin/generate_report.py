#!/usr/bin/env python3

import argparse
from pathlib import Path
import shutil

import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


def require_columns(data: pd.DataFrame, columns: set[str], source: Path) -> None:
	missing = columns.difference(data.columns)
	if missing:
		raise ValueError(f'{source} is missing required columns: {", ".join(sorted(missing))}.')


def plot_qc(sample_id: str, metrics: pd.DataFrame, output_path: Path) -> None:
	columns = (
		('total_counts', 'Total counts per cell'),
		('n_genes_by_counts', 'Genes detected per cell'),
		('pct_counts_mt', 'Mitochondrial counts (%)'),
	)
	figure, axes = plt.subplots(1, 3, figsize=(13, 4))
	for axis, (column, title) in zip(axes, columns):
		values = pd.to_numeric(metrics[column], errors='coerce').dropna()
		values = values[np.isfinite(values)]
		if values.empty:
			raise ValueError(f'QC metric {column} contains no finite values.')
		axis.hist(values, bins=40, color='#287c78', edgecolor='white', linewidth=0.35)
		if column != 'pct_counts_mt' and (values > 0).all():
			axis.set_xscale('log')
		axis.set_title(title)
		axis.set_ylabel('Cells')
		axis.grid(axis='y', alpha=0.2)
	figure.suptitle(f'{sample_id}: quality-control distributions')
	figure.tight_layout()
	figure.savefig(output_path, dpi=160, bbox_inches='tight')
	plt.close(figure)


def plot_umap(sample_id: str, adata: sc.AnnData, output_path: Path) -> None:
	coordinates = np.asarray(adata.obsm['X_umap'])
	if coordinates.ndim != 2 or coordinates.shape != (adata.n_obs, 2):
		raise ValueError('X_umap must contain two coordinates for every cell.')
	if not np.isfinite(coordinates).all():
		raise ValueError('X_umap contains non-finite coordinates.')

	group_key = 'cell_type' if 'cell_type' in adata.obs else 'leiden'
	groups = adata.obs[group_key].astype(str).to_numpy()
	categories = sorted(pd.unique(groups))
	palette = plt.get_cmap('tab20', max(len(categories), 1))
	figure, axis = plt.subplots(figsize=(8, 6))
	for index, group in enumerate(categories):
		selected = groups == group
		axis.scatter(
			coordinates[selected, 0],
			coordinates[selected, 1],
			s=8,
			alpha=0.8,
			color=palette(index),
			label=group if group_key == 'cell_type' else f'Cluster {group}',
			linewidths=0,
			rasterized=True,
		)
		center = np.median(coordinates[selected], axis=0)
		axis.annotate(
			group if group_key == 'cell_type' else f'Cluster {group}',
			xy=center,
			ha='center',
			va='center',
			fontsize=8,
			fontweight='bold',
			bbox={'boxstyle': 'round,pad=0.2', 'facecolor': 'white', 'alpha': 0.8, 'edgecolor': 'none'},
		)
	group_label = 'cell type' if group_key == 'cell_type' else 'Leiden cluster'
	axis.set(title=f'{sample_id}: UMAP by {group_label}', xlabel='UMAP 1', ylabel='UMAP 2')
	axis.legend(title=group_label.title(), bbox_to_anchor=(1.02, 1), loc='upper left', frameon=False)
	axis.set_aspect('equal', adjustable='datalim')
	figure.tight_layout()
	figure.savefig(output_path, dpi=160, bbox_inches='tight')
	plt.close(figure)


def plot_markers(
	sample_id: str,
	adata: sc.AnnData,
	markers: pd.DataFrame,
	output_path: Path,
) -> None:
	if 'cell_type' not in adata.obs or adata.obs['cell_type'].isna().any():
		raise ValueError('Cell-type labels are required to plot marker panels.')

	cell_assignments = pd.DataFrame({
		'cluster': adata.obs['leiden'].astype(str),
		'cell_type': adata.obs['cell_type'].astype(str),
	})
	dominant_cell_types = (
		cell_assignments.groupby(['cluster', 'cell_type'], observed=True)
		.size()
		.rename('cell_count')
		.reset_index()
		.sort_values(
			['cluster', 'cell_count', 'cell_type'],
			ascending=[True, False, True],
		)
		.drop_duplicates('cluster')
		.set_index('cluster')['cell_type']
		.to_dict()
	)
	clusters = list(pd.unique(markers['cluster'].astype(str)))
	columns = min(4, max(len(clusters), 1))
	rows = (len(clusters) + columns - 1) // columns
	figure, axes = plt.subplots(rows, columns, figsize=(4.2 * columns, 2.8 * rows), squeeze=False)
	for axis, cluster in zip(axes.flat, clusters):
		selected = markers.loc[markers['cluster'].astype(str) == cluster].head(5).copy()
		adjusted_p = pd.to_numeric(selected['pvals_adj'], errors='coerce').fillna(1.0)
		significance = -np.log10(np.maximum(adjusted_p.to_numpy(), 1e-50))
		selected = selected.assign(significance=significance)
		axis.barh(selected['gene'].astype(str), selected['significance'], color='#d9824b')
		axis.invert_yaxis()
		axis.set_title(dominant_cell_types.get(cluster, 'Unknown cell type'))
		axis.set_xlabel('-log10 adjusted p-value')
		axis.grid(axis='x', alpha=0.2)
	for axis in list(axes.flat)[len(clusters):]:
		axis.set_visible(False)
	figure.suptitle(f'{sample_id}: top ranked markers by cell type')
	figure.tight_layout()
	figure.savefig(output_path, dpi=160, bbox_inches='tight')
	plt.close(figure)


def copy_if_needed(source: Path, destination: Path) -> None:
	if source.resolve() != destination.resolve():
		shutil.copyfile(source, destination)


def generate_report(
	sample_id: str,
	qc_metrics_path: Path,
	features_path: Path,
	input_h5ad: Path,
	markers_path: Path,
	output_prefix: str,
) -> None:
	qc_metrics = pd.read_csv(qc_metrics_path)
	features = pd.read_csv(features_path, sep='\t')
	markers = pd.read_csv(markers_path, sep='\t')
	require_columns(
		qc_metrics,
		{'sample_id', 'barcode', 'total_counts', 'n_genes_by_counts', 'pct_counts_mt'},
		qc_metrics_path,
	)
	require_columns(features, {'gene'}, features_path)
	require_columns(
		markers,
		{'cluster', 'rank', 'gene', 'scores', 'logfoldchanges', 'pvals_adj'},
		markers_path,
	)
	if qc_metrics.empty or markers.empty:
		raise ValueError('Report requires non-empty QC metrics and marker tables.')

	adata = sc.read_h5ad(input_h5ad)
	if 'leiden' not in adata.obs or adata.obs['leiden'].isna().any():
		raise ValueError('Clustered AnnData must contain non-missing Leiden labels.')
	if 'cell_type' not in adata.obs or adata.obs['cell_type'].isna().any():
		raise ValueError('Clustered AnnData must contain non-missing cell-type labels.')
	if 'X_umap' not in adata.obsm:
		raise ValueError('Clustered AnnData must contain a UMAP embedding in obsm["X_umap"].')

	prefix = Path(output_prefix)
	plot_qc(
		sample_id,
		qc_metrics,
		prefix.with_name(f'{prefix.name}_qc_distributions.png'),
	)
	plot_umap(
		sample_id,
		adata,
		prefix.with_name(f'{prefix.name}_umap_clusters.png'),
	)
	plot_markers(
		sample_id,
		adata,
		markers,
		prefix.with_name(f'{prefix.name}_top_markers.png'),
	)

	cell_type_counts = (
		adata.obs['cell_type'].astype(str)
		.value_counts()
		.sort_index()
		.rename_axis('cell_type')
		.reset_index(name='n_cells')
	)
	summary = pd.DataFrame([{
		'sample_id': sample_id,
		'n_cells': adata.n_obs,
		'n_genes': adata.n_vars,
		'n_clusters': adata.obs['leiden'].nunique(),
		'n_highly_variable_genes': len(features),
		'median_total_counts': pd.to_numeric(qc_metrics['total_counts']).median(),
		'median_genes_per_cell': pd.to_numeric(qc_metrics['n_genes_by_counts']).median(),
		'median_pct_counts_mt': pd.to_numeric(qc_metrics['pct_counts_mt']).median(),
	}])

	copy_if_needed(qc_metrics_path, Path(f'{output_prefix}_qc_metrics.csv'))
	copy_if_needed(features_path, Path(f'{output_prefix}_highly_variable_genes.tsv'))
	copy_if_needed(markers_path, Path(f'{output_prefix}_markers.tsv'))
	cell_type_counts.to_csv(f'{output_prefix}_cluster_counts.tsv', sep='\t', index=False)
	summary.to_csv(f'{output_prefix}_summary.tsv', sep='\t', index=False)


def main() -> None:
	parser = argparse.ArgumentParser(
		description='Generate per-sample single-cell figures and tables.'
	)
	parser.add_argument('sample_id')
	parser.add_argument('qc_metrics', type=Path)
	parser.add_argument('features', type=Path)
	parser.add_argument('clustered_h5ad', type=Path)
	parser.add_argument('markers', type=Path)
	parser.add_argument('--output-prefix', required=True)
	args = parser.parse_args()

	generate_report(
		args.sample_id,
		args.qc_metrics,
		args.features,
		args.clustered_h5ad,
		args.markers,
		args.output_prefix,
	)
	print(f'{args.sample_id}: wrote analysis figures and tables.')


if __name__ == '__main__':
	main()
