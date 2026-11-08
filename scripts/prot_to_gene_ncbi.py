#!/usr/bin/env python
"""Map NCBI CDS (protein) IDs to gene IDs from a GFF3 file, via shared GeneID Dbxref."""
import gzip
import sys
from collections import defaultdict
from urllib.parse import unquote

import polars as pl


def parse_attrs(field):
    attrs = {}
    for item in field.split(";"):
        key, sep, value = item.partition("=")
        if sep:
            attrs[key] = [unquote(v) for v in value.split(",")]
    return attrs


def main(input_path, output_path):
    opener = gzip.open if input_path.endswith(".gz") else open
    ids_by_gene = defaultdict(list)

    with opener(input_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            attrs = parse_attrs(cols[8])
            if "Dbxref" not in attrs or "ID" not in attrs:
                continue
            gene_id = next((x for x in attrs["Dbxref"] if x.startswith("GeneID:")), None)
            if gene_id is not None:
                ids_by_gene[gene_id].append(attrs["ID"][0])

    prot_to_gene = {}
    for ids in ids_by_gene.values():
        gene = next((i for i in ids if i.startswith("gene-")), None)
        if gene is None:
            continue
        for i in ids:
            if i.startswith("cds-"):
                prot_to_gene[i.removeprefix("cds-")] = gene.removeprefix("gene-")

    df = pl.DataFrame(
        {"gene": list(prot_to_gene.values()), "protein": list(prot_to_gene.keys())},
        schema={"gene": pl.String, "protein": pl.String},
    )
    df.write_csv(output_path)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit("Usage: prot_to_gene_ncbi.py <input_gff3[.gz]> <output_csv>")
    main(sys.argv[1], sys.argv[2])
