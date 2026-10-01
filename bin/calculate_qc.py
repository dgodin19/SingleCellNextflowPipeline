#!/usr/bin/env python3

import argparse
import shutil
import tarfile
import tempfile
from collections import defaultdict
from pathlib import Path, PurePosixPath
import scanpy as sc


def logical_name(name: str) -> str:
	return name[:-3] if name.endswith('.gz') else name


def find_matrix_files(matrix_root: Path) -> tuple[Path, Path, Path]:
	files_by_parent = defaultdict(dict)
	for path in matrix_root.rglob('*'):
		if path.is_file():
			files_by_parent[path.parent][logical_name(path.name)] = path

	layouts = []
	for parent, files in files_by_parent.items():
		feature_file = files.get('features.tsv') or files.get('genes.tsv')
		if feature_file and 'barcodes.tsv' in files and 'matrix.mtx' in files:
			layouts.append((parent, files['barcodes.tsv'], feature_file, files['matrix.mtx']))

	if len(layouts) != 1:
		raise ValueError(
			f'Expected one complete 10x matrix under {matrix_root}, found {len(layouts)}.'
		)
	return layouts[0][1:]


def stage_archive_matrix(archive_path: Path, destination: Path) -> Path:
	files_by_parent = defaultdict(dict)
	members = {}
	with tarfile.open(archive_path, mode='r:*') as archive:
		for member in archive.getmembers():
			if not member.isfile():
				continue
			member_path = PurePosixPath(member.name)
			files_by_parent[member_path.parent][logical_name(member_path.name)] = member
			members[member.name] = member

		layouts = []
		for parent, files in files_by_parent.items():
			feature_member = files.get('features.tsv') or files.get('genes.tsv')
			if feature_member and 'barcodes.tsv' in files and 'matrix.mtx' in files:
				layouts.append((parent, files['barcodes.tsv'], feature_member, files['matrix.mtx']))

		if len(layouts) != 1:
			raise ValueError(
				f'Expected one complete 10x matrix in {archive_path}, found {len(layouts)}.'
			)

		matrix_directory = destination / 'matrix'
		matrix_directory.mkdir()
		for _, *matrix_members in layouts:
			for member in matrix_members:
				source = archive.extractfile(member)
				if source is None:
					raise ValueError(f'Could not read {member.name} from {archive_path}.')
				target = matrix_directory / PurePosixPath(member.name).name
				with source, target.open('wb') as output:
					shutil.copyfileobj(source, output)

	return matrix_directory


def read_10x_matrix(matrix_input: Path, temporary_directory: Path):
	

	if matrix_input.is_dir():
		_, _, _ = find_matrix_files(matrix_input)
		matrix_directory = matrix_input
	elif matrix_input.is_file():
		matrix_directory = stage_archive_matrix(matrix_input, temporary_directory)
	else:
		raise FileNotFoundError(f'Matrix input does not exist: {matrix_input}')

	return sc.read_10x_mtx(
		matrix_directory,
		var_names='gene_symbols',
		make_unique=True,
		cache=False,
	)


def calculate_qc(sample_id: str, matrix_input: Path, metrics_output: Path, adata_output: Path) -> None:
	

	with tempfile.TemporaryDirectory(prefix='calculate_qc_') as temporary_directory:
		temporary_path = Path(temporary_directory)
		adata = read_10x_matrix(matrix_input, temporary_path)

		adata.var['mt'] = adata.var_names.astype(str).str.upper().str.startswith('MT-')
		sc.pp.calculate_qc_metrics(
			adata,
			qc_vars=['mt'],
			percent_top=None,
			log1p=False,
			inplace=True,
		)

		metrics = adata.obs[['total_counts', 'n_genes_by_counts', 'pct_counts_mt']].copy()
		metrics.insert(0, 'barcode', adata.obs_names.astype(str))
		metrics.insert(0, 'sample_id', sample_id)
		metrics.to_csv(metrics_output, index=False)

		adata.obs['sample_id'] = sample_id
		adata.write_h5ad(adata_output)


def main() -> None:
	parser = argparse.ArgumentParser(description='Calculate per-cell 10x QC metrics with Scanpy.')
	parser.add_argument('sample_id')
	parser.add_argument('matrix_input', type=Path)
	parser.add_argument('metrics_output', type=Path)
	parser.add_argument('adata_output', type=Path)
	args = parser.parse_args()

	calculate_qc(args.sample_id, args.matrix_input, args.metrics_output, args.adata_output)


if __name__ == '__main__':
	main()
