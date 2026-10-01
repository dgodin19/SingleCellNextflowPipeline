#!/usr/bin/env python3

import argparse
import csv
import gzip
import io
import tarfile
from contextlib import contextmanager
from decimal import Decimal, InvalidOperation
from pathlib import Path, PurePosixPath


def logical_name(name: str) -> str:
    return name[:-3] if name.endswith('.gz') else name


def collect_directory_files(root: Path):
    files = {}
    for path in root.rglob('*'):
        if path.is_file():
            relative_path = path.relative_to(root).as_posix()
            parent = PurePosixPath(relative_path).parent.as_posix()
            files.setdefault(parent, {}).setdefault(logical_name(path.name), []).append(path)
    return files


def collect_archive_files(archive: tarfile.TarFile):
    files = {}
    members = {}
    for member in archive.getmembers():
        if member.isfile():
            member_path = PurePosixPath(member.name)
            parent = member_path.parent.as_posix()
            name = logical_name(member_path.name)
            files.setdefault(parent, {}).setdefault(name, []).append(member.name)
            members[member.name] = member
    return files, members


@contextmanager
def open_text(source, archive=None):
    if archive is None:
        binary_handle = source.open('rb')
    else:
        binary_handle = archive.extractfile(source)
        if binary_handle is None:
            raise ValueError(f'Could not read archive member {source.name}.')

    if source.name.endswith('.gz'):
        binary_handle = gzip.GzipFile(fileobj=binary_handle)

    text_handle = io.TextIOWrapper(binary_handle, encoding='utf-8-sig', newline='')
    try:
        yield text_handle
    finally:
        text_handle.close()


def find_layout(files):
    layouts = []
    for parent, names in files.items():
        barcodes = names.get('barcodes.tsv', [])
        matrices = names.get('matrix.mtx', [])
        features = names.get('features.tsv', [])
        genes = names.get('genes.tsv', [])
        feature_files = features + genes
        if len(barcodes) == 1 and len(matrices) == 1 and len(feature_files) == 1:
            layouts.append((parent, barcodes[0], feature_files[0], matrices[0]))

    if not layouts:
        raise ValueError(
            'No complete 10x matrix found. barcodes.tsv, features.tsv or genes.tsv, '
            'and matrix.mtx must be together in one directory.'
        )
    if len(layouts) > 1:
        raise ValueError('Multiple complete 10x matrices found; provide one matrix directory/archive.')
    return layouts[0]


def count_barcodes(source, archive=None) -> int:
    count = 0
    seen = set()
    with open_text(source, archive) as handle:
        reader = csv.reader(handle, delimiter='\t')
        for row_number, row in enumerate(reader, start=1):
            if len(row) != 1 or not row[0].strip():
                raise ValueError(f'barcodes.tsv row {row_number} must contain one non-empty barcode.')
            barcode = row[0].strip()
            if barcode in seen:
                raise ValueError(f'barcodes.tsv contains duplicate barcode {barcode!r}.')
            seen.add(barcode)
            count += 1
    if count == 0:
        raise ValueError('barcodes.tsv is empty.')
    return count


def count_features(source, archive=None) -> int:
    count = 0
    with open_text(source, archive) as handle:
        reader = csv.reader(handle, delimiter='\t')
        for row_number, row in enumerate(reader, start=1):
            if len(row) not in (2, 3) or any(not value.strip() for value in row):
                raise ValueError(
                    f'features.tsv/genes.tsv row {row_number} must have two or three non-empty fields.'
                )
            count += 1
    if count == 0:
        raise ValueError('features.tsv/genes.tsv is empty.')
    return count


def validate_matrix(source, feature_count: int, barcode_count: int, archive=None) -> int:
    with open_text(source, archive) as handle:
        header = handle.readline().strip().lower().split()
        if (
            len(header) != 5
            or header[0] != '%%matrixmarket'
            or header[1:3] != ['matrix', 'coordinate']
            or header[3] not in ('integer', 'real')
            or header[4] != 'general'
        ):
            raise ValueError(
                'matrix.mtx must use Matrix Market coordinate integer or real general format.'
            )
        value_format = header[3]

        dimensions_line = None
        for line in handle:
            stripped = line.strip()
            if stripped and not stripped.startswith('%'):
                dimensions_line = stripped
                break

        if dimensions_line is None:
            raise ValueError('matrix.mtx is missing its dimensions line.')
        dimensions = dimensions_line.split()
        if len(dimensions) != 3:
            raise ValueError('matrix.mtx dimensions line must contain rows, columns, and entries.')
        try:
            matrix_rows, matrix_columns, expected_entries = map(int, dimensions)
        except ValueError as error:
            raise ValueError('matrix.mtx dimensions must be integers.') from error

        if matrix_rows != feature_count or matrix_columns != barcode_count:
            raise ValueError(
                'matrix.mtx dimensions do not match the feature and barcode counts '
                f'({matrix_rows}x{matrix_columns} vs {feature_count}x{barcode_count}).'
            )
        if matrix_rows < 1 or matrix_columns < 1 or expected_entries < 1:
            raise ValueError('matrix.mtx dimensions and entry count must be positive.')

        actual_entries = 0
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            fields = stripped.split()
            if len(fields) != 3:
                raise ValueError(f'matrix.mtx entry {line_number} must have three fields.')
            try:
                row_index, column_index = map(int, fields[:2])
                value = int(fields[2]) if value_format == 'integer' else Decimal(fields[2])
            except (InvalidOperation, ValueError) as error:
                raise ValueError(f'matrix.mtx entry {line_number} has an invalid value.') from error
            if not 1 <= row_index <= matrix_rows or not 1 <= column_index <= matrix_columns:
                raise ValueError(f'matrix.mtx entry {line_number} has an out-of-range index.')
            if isinstance(value, Decimal) and not value.is_finite():
                raise ValueError(f'matrix.mtx entry {line_number} has a non-finite value.')
            if value < 0:
                raise ValueError(f'matrix.mtx entry {line_number} has a negative count.')
            actual_entries += 1

        if actual_entries != expected_entries:
            raise ValueError(
                f'matrix.mtx declares {expected_entries} entries but contains {actual_entries}.'
            )
    return actual_entries


def validate_10x_matrix(matrix_path: Path):
    if matrix_path.is_dir():
        files = collect_directory_files(matrix_path)
        parent, barcode_source, feature_source, matrix_source = find_layout(files)
        directory = matrix_path / parent
        barcode_path = directory / barcode_source
        feature_path = directory / feature_source
        matrix_path = directory / matrix_source
        archive = None
        barcode_source, feature_source, matrix_source = barcode_path, feature_path, matrix_path
        return validate_sources(barcode_source, feature_source, matrix_source, archive)

    if not matrix_path.is_file():
        raise ValueError(f'Matrix input does not exist: {matrix_path}')

    archive = tarfile.open(matrix_path, mode='r:*')
    try:
        files, members = collect_archive_files(archive)
        _, barcode_name, feature_name, matrix_name = find_layout(files)
        barcode_member = members[barcode_name]
        feature_member = members[feature_name]
        matrix_member = members[matrix_name]
        return validate_sources(barcode_member, feature_member, matrix_member, archive)
    finally:
        archive.close()


def validate_sources(barcode_source, feature_source, matrix_source, archive=None):
    barcode_count = count_barcodes(barcode_source, archive)
    feature_count = count_features(feature_source, archive)
    entry_count = validate_matrix(matrix_source, feature_count, barcode_count, archive)
    return feature_count, barcode_count, entry_count


def main() -> None:
    parser = argparse.ArgumentParser(description='Validate a 10x Genomics matrix.')
    parser.add_argument('matrix_input', type=Path)
    parser.add_argument('report', type=Path)
    args = parser.parse_args()

    try:
        feature_count, barcode_count, entry_count = validate_10x_matrix(args.matrix_input)
    except (OSError, tarfile.TarError, UnicodeError, ValueError) as error:
        raise SystemExit(f'Invalid 10x matrix: {error}') from error

    report = (
        f'Valid 10x matrix: {args.matrix_input}\n'
        f'Features: {feature_count}\n'
        f'Barcodes: {barcode_count}\n'
        f'Matrix entries: {entry_count}\n'
    )
    args.report.write_text(report, encoding='utf-8')
    print(report, end='')


if __name__ == '__main__':
    main()