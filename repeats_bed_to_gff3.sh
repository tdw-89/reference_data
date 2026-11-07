#!/bin/zsh
# Convert a RepeatMasker BED (as written by each species' bed_conversion_command.sh,
# then lifted) into a GFF3 that BioinfoTools2's `add_features!` can load.
#
# Usage: ./repeats_bed_to_gff3.sh <in.bed> <out.gff3>
#
# Input columns: chrom, start (0-based), end (exclusive), <subfamily>#<class>/<family>,
# score, strand. Each row becomes one `repeat_region` with
#   ID=<subfamily>#<class>/<family>:<chrom>:<start>-<end>   (GFF 1-based closed)
# IDs must be unique, because BioinfoTools2 keys features by ID. Lifting can leave
# two rows with the same name and coordinates, so a repeated ID gets a `.<n>` suffix.
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: $0 <in.bed> <out.gff3>" >&2
    exit 1
fi

awk 'BEGIN { FS = OFS = "\t"; print "##gff-version 3" }
{
    start = $2 + 1
    id = $4 ":" $1 ":" start "-" $3
    if (id in seen) { id = id "." (++seen[id]) } else { seen[id] = 0 }
    print $1, "RepeatMasker", "repeat_region", start, $3, $5, $6, ".", "ID=" id
}' "$1" > "$2"
