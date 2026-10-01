#!/usr/bin/env python3

import argparse
import csv
import shutil
from pathlib import Path, PurePosixPath, PureWindowsPath


def validate_metadata(metadata_path: Path) -> int:
    required_columns = {'sample_id', 'matrix_path', 'condition'}
    sample_ids = set()
    row_count = 0

    with metadata_path.open(newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise SystemExit('Input metadata CSV is empty.')

        if len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise SystemExit('Input metadata CSV contains duplicate column names.')

        missing_columns = required_columns.difference(reader.fieldnames)
        if missing_columns:
            raise SystemExit(
                'Input metadata CSV is missing required columns: '
                + ', '.join(sorted(missing_columns))
            )

        for row_number, row in enumerate(reader, start=2):
            if None in row:
                raise SystemExit(f'Row {row_number} has more values than the header.')

            values = {column: (row.get(column) or '').strip() for column in required_columns}
            empty_columns = [column for column, value in values.items() if not value]
            if empty_columns:
                raise SystemExit(
                    f'Row {row_number} has empty required fields: '
                    + ', '.join(sorted(empty_columns))
                )

            sample_id = values['sample_id']
            if sample_id in sample_ids:
                raise SystemExit(f'Row {row_number} has duplicate sample_id: {sample_id}')
            sample_ids.add(sample_id)

            matrix_path = values['matrix_path']
            matrix_path_parts = PurePosixPath(matrix_path)
            if (
                matrix_path_parts.is_absolute()
                or PureWindowsPath(matrix_path).is_absolute()
                or not matrix_path_parts.parts
                or '..' in matrix_path_parts.parts
                or '\\' in matrix_path
            ):
                raise SystemExit(
                    f'Row {row_number} matrix_path must be a safe, relative POSIX path.'
                )

            row_count += 1

    if row_count == 0:
        raise SystemExit('Input metadata CSV must contain at least one sample.')

    return row_count


def main() -> None:
    parser = argparse.ArgumentParser(description='Validate sample metadata CSV.')
    parser.add_argument('metadata_csv', type=Path)
    parser.add_argument('output_csv', type=Path)
    args = parser.parse_args()

    row_count = validate_metadata(args.metadata_csv)
    shutil.copyfile(args.metadata_csv, args.output_csv)
    print(f'Validated metadata for {row_count} sample(s).')


if __name__ == '__main__':
    main()