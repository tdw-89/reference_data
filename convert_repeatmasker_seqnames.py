#!/usr/bin/env python3
"""Replace RepeatMasker sequence names using an NCBI sequence report TSV.

Example:
    python3 convert_repeatmasker_seqnames.py sequence_report.tsv repeats/input.out

By default, the output is written beside the RepeatMasker file with
".refseq" inserted before its final extension. Use --output to choose a path.
"""

import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path


RECORD_PREFIX = re.compile(
    r"^[ \t]*[0-9]+[ \t]+\S+[ \t]+\S+[ \t]+\S+[ \t]+(?P<sequence>\S+)"
)


def normalize_header(value):
    return " ".join(value.strip().casefold().split())


def load_mapping(mapping_path):
    """Read sequence-name to RefSeq-accession pairs from an NCBI TSV."""
    mapping = {}

    with mapping_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        headers = {
            normalize_header(name): name
            for name in (reader.fieldnames or [])
            if name is not None
        }
        required = ("sequence name", "refseq seq accession")
        missing = [name for name in required if name not in headers]
        if missing:
            raise ValueError(
                "Mapping TSV is missing required column(s): " + ", ".join(missing)
            )

        sequence_column = headers["sequence name"]
        accession_column = headers["refseq seq accession"]

        for row_number, row in enumerate(reader, start=2):
            sequence_name = (row.get(sequence_column) or "").strip()
            accession = (row.get(accession_column) or "").strip()

            # NCBI sequence reports can include rows without a RefSeq accession.
            if not sequence_name or not accession:
                continue

            previous = mapping.get(sequence_name)
            if previous is not None and previous != accession:
                raise ValueError(
                    f"Conflicting RefSeq accessions for {sequence_name!r} "
                    f"in mapping TSV near line {row_number}: "
                    f"{previous!r} and {accession!r}"
                )
            mapping[sequence_name] = accession

    if not mapping:
        raise ValueError("No sequence-name to RefSeq-accession pairs were found")
    return mapping


def default_output_path(input_path):
    return input_path.with_name(f"{input_path.stem}.refseq{input_path.suffix}")


def convert_file(mapping_path, input_path, output_path):
    mapping = load_mapping(mapping_path)
    input_resolved = input_path.resolve()
    output_resolved = output_path.resolve()
    mapping_resolved = mapping_path.resolve()

    if output_resolved in (input_resolved, mapping_resolved):
        raise ValueError("Output path must differ from both input files")

    total_records = 0
    renamed_records = 0
    unmapped_names = Counter()

    # Stream the large RepeatMasker file and preserve its whitespace and line endings.
    with input_path.open("r", encoding="utf-8", newline="") as source, output_path.open(
        "w", encoding="utf-8", newline=""
    ) as destination:
        for line in source:
            match = RECORD_PREFIX.match(line)
            if match is None:
                destination.write(line)
                continue

            total_records += 1
            sequence_name = match.group("sequence")
            accession = mapping.get(sequence_name)
            if accession is None:
                unmapped_names[sequence_name] += 1
                destination.write(line)
                continue

            if accession != sequence_name:
                renamed_records += 1
                start, end = match.span("sequence")
                line = line[:start] + accession + line[end:]
            destination.write(line)

    return total_records, renamed_records, unmapped_names


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Replace the sequence-name field in RepeatMasker .out records "
            "with RefSeq accessions from an NCBI sequence_report.tsv."
        )
    )
    parser.add_argument("mapping_tsv", type=Path, help="NCBI sequence report TSV")
    parser.add_argument("repeatmasker_out", type=Path, help="RepeatMasker .out file")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help=(
            "output file path (default: beside the input with '.refseq' "
            "inserted before the final extension)"
        ),
    )
    args = parser.parse_args(argv)
    output_path = args.output or default_output_path(args.repeatmasker_out)

    try:
        total, renamed, unmapped = convert_file(
            args.mapping_tsv, args.repeatmasker_out, output_path
        )
    except (OSError, ValueError, csv.Error) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1

    print(f"Renamed {renamed:,} of {total:,} RepeatMasker records.")
    print(f"Output: {output_path}")
    if unmapped:
        names = sorted(unmapped)
        sample = ", ".join(names[:10])
        remainder = len(names) - min(len(names), 10)
        if remainder:
            sample += f", and {remainder:,} more"
        print(
            f"Left {len(unmapped):,} unmapped sequence name(s) unchanged "
            f"({sum(unmapped.values()):,} records): {sample}",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
