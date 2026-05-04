#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash select_wells_sequences_25.sh /path/to/wells_sequences [output_prefix]
#
# Produces:
#   <output_prefix>.faa
#   <output_prefix>.metadata.tsv

WELLS_DIR="${1:-wells_sequences}"
OUT_PREFIX="${2:-selected_wells_major_refs_25}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python3 "${SCRIPT_DIR}/select_wells_references.py" \
  --wells-dir "${WELLS_DIR}" \
  --out-faa "${OUT_PREFIX}.faa" \
  --out-metadata "${OUT_PREFIX}.metadata.tsv" \
  --per-family 25 \
  --preset major \
  --source-fasta auto
#DMSORcompwoutMop.fasta
