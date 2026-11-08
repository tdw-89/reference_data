#!/usr/bin/env python3
"""Replace FASTA defline sequence names using an NCBI sequence report TSV.

Example:
    python3 convert_fasta_seqnames.py sequence_report.tsv repeats/input.fasta

By default, the output is written beside the FASTA file with ".refseq"
inserted before its final extension. Use --output to choose a path.
"""

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

from convert_repeatmasker_seqnames import default_output_path, load_mapping


def line_ending(line):
    if line.endswith("\r\n"):
        return "\r\n"
    if line.endswith("\n") or line.endswith("\r"):
        return line[-1:]
    return ""


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

    # Stream large FASTA files and copy every sequence line unchanged.
    with input_path.open("r", encoding="utf-8", newline="") as source, output_path.open(
        "w", encoding="utf-8", newline=""
    ) as destination:
        for line in source:
            if not line.startswith(">"):
                destination.write(line)
                continue

            total_records += 1
            ending = line_ending(line)
            defline = line[1:]
            if ending:
                defline = defline[: -len(ending)]
            sequence_name = defline.strip()
            if not sequence_name:
                raise ValueError(f"Empty FASTA defline at record {total_records}")

            accession = mapping.get(sequence_name)
            if accession is None:
                unmapped_names[sequence_name] += 1
                destination.write(line)
            else:
                if accession != sequence_name:
                    renamed_records += 1
                destination.write(f">{accession}{ending}")

    return total_records, renamed_records, unmapped_names


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Replace FASTA defline sequence names with RefSeq accessions "
            "from an NCBI sequence_report.tsv."
        )
    )
    parser.add_argument("mapping_tsv", type=Path, help="NCBI sequence report TSV")
    parser.add_argument("fasta", type=Path, help="input FASTA file")
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
    output_path = args.output or default_output_path(args.fasta)

    try:
        total, renamed, unmapped = convert_file(args.mapping_tsv, args.fasta, output_path)
    except (OSError, ValueError, csv.Error) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1

    print(f"Renamed {renamed:,} of {total:,} FASTA deflines.")
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
