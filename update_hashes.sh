#!/bin/sh
# Rewrites hash_list.txt as "path<TAB>md5". Re-runnable: only column 1 is read.
# Usage: ./update_hashes.sh [list_file]   (run from the directory the paths are relative to)
set -eu

list=${1:-hash_list.txt}
tmp=$(mktemp)
trap 'rm -f "$tmp" "$tmp.list"' EXIT

if command -v md5sum >/dev/null 2>&1; then
  md5cmd="md5sum"
else
  md5cmd="md5 -r"   # macOS: prints "hash path"
fi

awk 'NF {print $1}' "$list" > "$tmp.list"

# 4 files hashed in parallel; output order is restored by the awk join below
xargs -n 1 -P 4 $md5cmd < "$tmp.list" > "$tmp"

awk -v OFS='\t' '
  NR == FNR { h[$2] = $1; next }
  NF        { print $1, h[$1]; next }
            { print }
' "$tmp" "$list" > "$list.new"
mv "$list.new" "$list"
