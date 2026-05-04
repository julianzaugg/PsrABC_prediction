#!/usr/bin/env python3
"""
select_wells_references.py

Build a pipeline-compatible custom reference FASTA from the Wells et al. Dryad
MopB/DMSOR FASTA files.

This version is intended for a broader, manuscript-quality reference set. By
default it selects 25 representatives per major Wells family/clade from the
available Wells FASTA files.

Output FASTA headers follow the same family-first style as the built-in
pipeline references, e.g.:

  >PsrAPhsASrrA_Wells_001__original_id
  >NapA_Wells_001__original_id
  >FdhG_Wells_001__original_id

The first underscore-delimited field is the family/clade anchor used by the
pipeline's --custom-tree-refs metadata auto-parser.

Typical use:

  python3 select_wells_references.py \
    --wells-dir /path/to/wells_fastas \
    --out-faa data/wells_major_refs_25.faa \
    --out-metadata data/wells_major_refs_25.metadata.tsv \
    --per-family 25

Then run the pipeline with:

  bash 00_pipeline.sh ... --custom-tree-refs data/wells_major_refs_25.faa
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


# Family/clade definitions. Each clade is matched by exact sequence ID against
# the listed Wells files. The list deliberately includes raw, wcodes, and final
# variants where available so the same script works with either genomic source
# FASTAs (e.g. DMSORcompwoutMop.fasta) or metagenomic source FASTAs
# (e.g. MopBMAGs_50.fasta).
TARGET_FAMILY_FILES = OrderedDict({
    # Psr/Phs/Srr/Sre/Soe/Ttr/Arr groups relevant to the original PsrA question.
    "PsrAPhsASrrA": [
        "PsrAPhsASrrA.fasta",
        "PsrAPhsASrrAwcodes.fasta",
        "PsrAPhsAfinal.fasta",
    ],
    "TtrASrdA": [
        "TtrASrdAfinal.fasta",
        "TtrASrdAarchArrA.fasta",
        "TtrASrdAarchArrAwcodes.fasta",
    ],
    "bSreASoeA": [
        "bSreASoeAfinal.fasta",
        "bSreA.fasta",
        "bSreAwcodes.fasta",
    ],
    "aSreA": [
        "aSreA.fasta",
        "aSreAwcodes.fasta",
    ],
    "ArrAArxA": [
        "BactArrAArxA.fasta",
        "BactArrAArxAwcodes.fasta",
        "ArxAArrAfinal.fasta",
        "ArxAArrAcodes.fasta",
    ],

    # Major alternative MopB/DMSOR families that should act as explicit
    # non-Psr tree anchors.
    "NapA": [
        "NapA.fasta",
        "NapAwcodes.fasta",
        "NapAfinal.fasta",
    ],
    "NarG": [
        "NarG.fasta",
        "NarGwcodes.fasta",
        "NarGfinal.fasta",
    ],
    "NasCNarB": [
        "NasC.fasta",
        "NasCwcodes.fasta",
        "NasCNarBfinal.fasta",
    ],
    "DmsA": [
        "DmsAfinal.fasta",
    ],
    "DorATorA": [
        "DorATorAfinal.fasta",
    ],
    "BisC": [
        "BisCfinal.fasta",
    ],
    "AioAIdrA": [
        "AioA.fasta",
        "AioAwcodes.fasta",
        "AioAIdrAfinal.fasta",
        "IdrA.fasta",
    ],
    "ActB": [
        "ActB.fasta",
        "ActBwcodes.fasta",
        "ActBfinal.fasta",
    ],
    "AH": [
        "AH.fasta",
        "AHwcodes.fasta",
        "AHfinal.fasta",
    ],

    # Formate dehydrogenase / related outgroups and other broad branches.
    "FdhG": [
        "FdhG.fasta",
        "FdhGfullnames.fasta",
        "FdhGfinal.fasta",
    ],
    "FdhsFsdA": [
        "Fdhs.fasta",
        "Fdhsfinal.fasta",
        "FdhsFsdA.fasta",
        "FdhsFsdAfullnames.fasta",
    ],
    "FhcB": [
        "FhcB.fasta",
        "FhcBwcodes.fasta",
        "FhcBfinal.fasta",
    ],
    "FwdBFmdB": [
        "FwdB.fasta",
        "Bact_FwdB.fasta",
        "FwdB_FmdB.fasta",
        "FwdBfinal.fasta",
    ],
    "Nqo3": [
        "Nqo3.fasta",
        "Nqo3wcodes.fasta",
        "Nqo3final.fasta",
    ],
    "RhLPgtL": [
        "RhLPgtLfinal.fasta",
    ],
    "AspDMSOR": [
        "Aspfinal.fasta",
        "variousAspDMSORs.fasta",
        "variousAspDMSORswcodes.fasta",
    ],
    "SerDMSOR": [
        "variousSerDMSORs.fasta",
    ],

    # Aldehyde:ferredoxin oxidoreductase files are outside canonical MopB but
    # are present in the Wells download. Include only if requested with
    # --preset all or --include-family AOR.
    "AOR": [
        "AORswcodes.fasta",
    ],
})

DEFAULT_MAJOR_FAMILIES = [
    "PsrAPhsASrrA",
    "TtrASrdA",
    "bSreASoeA",
    "aSreA",
    "ArrAArxA",
    "NapA",
    "NarG",
    "NasCNarB",
    "DmsA",
    "DorATorA",
    "BisC",
    "AioAIdrA",
    "ActB",
    "AH",
    "FdhG",
    "FdhsFsdA",
    "FhcB",
    "FwdBFmdB",
    "Nqo3",
    "RhLPgtL",
    "AspDMSOR",
    "SerDMSOR",
]

KEY_PSR_FAMILIES = [
    "PsrAPhsASrrA",
    "TtrASrdA",
    "bSreASoeA",
    "aSreA",
    "ArrAArxA",
    "NapA",
    "FdhG",
    "DmsA",
]

AGGREGATE_BASENAMES = {
    "MopBMAGs.fasta",
    "MopBMAGs_50.fasta",
    "MopBMAGs_50aligned.fasta",
    "MopBMAGs_50aligned.faa",
    "allDMSORswoutMopaligned.fasta",
    "DMSORcompwoutMop.fasta",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Select Wells et al. family representatives for pipeline custom tree references."
    )
    p.add_argument("--wells-dir",
                   help="Directory containing Wells/Dryad FASTA files")
    p.add_argument("--out-faa",
                   help="Output FASTA path for selected references")
    p.add_argument("--out-metadata",
                   help="Output metadata TSV path")
    p.add_argument("--per-family", type=int, default=25,
                   help="Number of sequences to select per family [default: 25]")
    p.add_argument("--source-fasta", default="auto",
                   help=("Source FASTA within --wells-dir. Use 'auto' to prefer "
                         "DMSORcompwoutMop.fasta, then MopBMAGs_50.fasta [default: auto]"))
    p.add_argument("--preset", choices=["major", "key", "all"], default="major",
                   help="Family set to select [default: major]")
    p.add_argument("--include-family", action="append", default=None,
                   help=("Select only this family/clade; may be repeated. Valid values: "
                         + ", ".join(TARGET_FAMILY_FILES.keys())))
    p.add_argument("--exclude-family", action="append", default=[],
                   help="Exclude this family/clade; may be repeated")
    p.add_argument("--list-families", action="store_true",
                   help="Print available family/clade names and exit")
    p.add_argument("--strict", action="store_true",
                   help="Exit non-zero if fewer than --per-family sequences are available for any selected family")
    p.add_argument("--allow-duplicates", action="store_true",
                   help="Allow the same source sequence ID to be selected for multiple families")
    return p.parse_args()


def fasta_iter(path: Path) -> Iterable[Tuple[str, str, str]]:
    header = None
    seq_lines: List[str] = []
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    sid = header.split()[0]
                    yield sid, header, "".join(seq_lines)
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if header is not None:
            sid = header.split()[0]
            yield sid, header, "".join(seq_lines)


def ids_in_fasta(path: Path) -> set[str]:
    return {sid for sid, _header, _seq in fasta_iter(path)}


def wrap_seq(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def safe_label_part(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.|:-]+", "_", text).strip("_")


def resolve_source_fasta(wells_dir: Path, requested: str) -> Path:
    if requested != "auto":
        source_path = wells_dir / requested
        if not source_path.exists():
            sys.exit(f"[ERROR] Source FASTA not found: {source_path}")
        return source_path

    for fname in ["DMSORcompwoutMop.fasta", "MopBMAGs_50.fasta", "MopBMAGs.fasta"]:
        candidate = wells_dir / fname
        if candidate.exists():
            print(f"[*] --source-fasta auto selected: {fname}")
            return candidate
    sys.exit("[ERROR] Could not auto-detect source FASTA. Expected DMSORcompwoutMop.fasta or MopBMAGs_50.fasta")


def choose_families(args: argparse.Namespace) -> List[str]:
    if args.include_family:
        families = args.include_family
    elif args.preset == "key":
        families = KEY_PSR_FAMILIES
    elif args.preset == "all":
        families = list(TARGET_FAMILY_FILES.keys())
    else:
        families = DEFAULT_MAJOR_FAMILIES

    bad = [f for f in families + args.exclude_family if f not in TARGET_FAMILY_FILES]
    if bad:
        sys.exit("[ERROR] Unknown family/clade name(s): " + ", ".join(sorted(set(bad))))

    excluded = set(args.exclude_family)
    return [f for f in families if f not in excluded]


def build_family_id_map(wells_dir: Path, families: List[str]) -> Dict[str, set[str]]:
    fam_to_ids: Dict[str, set[str]] = {}
    for family in families:
        ids: set[str] = set()
        seen_files: List[str] = []
        missing_files: List[str] = []
        for fname in TARGET_FAMILY_FILES[family]:
            path = wells_dir / fname
            if not path.exists():
                missing_files.append(fname)
                continue
            if path.name in AGGREGATE_BASENAMES:
                continue
            seen_files.append(path.name)
            ids.update(ids_in_fasta(path))
        fam_to_ids[family] = ids
        if seen_files:
            print(f"[*] {family}: indexed {len(ids)} unique IDs from {', '.join(seen_files)}")
        else:
            print(f"[WARN] {family}: none of the expected family FASTA files were found", file=sys.stderr)
    return fam_to_ids


def select_evenly(records: List[Tuple[str, str, str]], n: int) -> List[Tuple[str, str, str]]:
    """Deterministically select up to n records spread across a sorted list."""
    if len(records) <= n:
        return records
    if n <= 1:
        return [records[0]]
    idxs: List[int] = []
    for i in range(n):
        idxs.append(round(i * (len(records) - 1) / (n - 1)))
    # Preserve order, de-duplicate rounded indices, then fill remaining positions.
    seen = set()
    unique: List[int] = []
    for idx in idxs:
        if idx not in seen:
            unique.append(idx)
            seen.add(idx)
    j = 0
    while len(unique) < n and j < len(records):
        if j not in seen:
            unique.append(j)
            seen.add(j)
        j += 1
    return [records[i] for i in sorted(unique[:n])]


def main() -> int:
    args = parse_args()

    if args.list_families:
        for family, files in TARGET_FAMILY_FILES.items():
            default_tag = ""
            if family in DEFAULT_MAJOR_FAMILIES:
                default_tag = " [major]"
            if family in KEY_PSR_FAMILIES:
                default_tag += " [key]"
            print(f"{family}{default_tag}\t" + ", ".join(files))
        return 0

    missing_required = [name for name in ["wells_dir", "out_faa", "out_metadata"] if getattr(args, name) in (None, "")]
    if missing_required:
        sys.exit("[ERROR] Missing required argument(s): " + ", ".join("--" + x.replace("_", "-") for x in missing_required))

    wells_dir = Path(args.wells_dir).expanduser().resolve()
    if not wells_dir.is_dir():
        sys.exit(f"[ERROR] Wells directory not found: {wells_dir}")
    if args.per_family < 1:
        sys.exit("[ERROR] --per-family must be >= 1")

    source_path = resolve_source_fasta(wells_dir, args.source_fasta)
    families = choose_families(args)
    print("[*] Selected family/clade set: " + ", ".join(families))

    fam_to_ids = build_family_id_map(wells_dir, families)
    source_records = {sid: (sid, header, seq) for sid, header, seq in fasta_iter(source_path)}
    print(f"[*] Source pool: {len(source_records)} sequences from {source_path.name}")

    selected_rows: List[dict[str, str]] = []
    selected_fasta_chunks: List[str] = []
    used_source_ids: set[str] = set()

    for family in families:
        family_ids = fam_to_ids.get(family, set())
        matches = [source_records[sid] for sid in sorted(source_records) if sid in family_ids]
        chosen = select_evenly(matches, args.per_family)
        print(f"[*] {family}: {len(matches)} source matches; selected {len(chosen)}")

        if args.strict and len(chosen) < args.per_family:
            sys.exit(f"[ERROR] {family}: only {len(chosen)} matches available; wanted {args.per_family}")

        family_counter = 0
        for sid, header, seq in chosen:
            if sid in used_source_ids and not args.allow_duplicates:
                print(f"[WARN] Skipping duplicated source ID across families: {sid}", file=sys.stderr)
                continue
            used_source_ids.add(sid)
            family_counter += 1

            accession = safe_label_part(sid)
            label = f"{family}_Wells_{family_counter:03d}"
            fasta_header = (
                f">{label}__{accession} "
                f"family={family} original_id={accession} original_header={safe_label_part(header)}"
            )
            selected_fasta_chunks.append(f"{fasta_header}\n{wrap_seq(seq)}\n")
            selected_rows.append({
                "label": label,
                "accession": accession,
                "organism": "Wells_et_al_reference",
                "protein": family,
                "clade": family,
                "source": f"Wells_Dryad_{source_path.name}",
                "note": f"Selected from {source_path.name}; original_header={header}",
                "status": "OK",
            })

    out_faa = Path(args.out_faa)
    out_meta = Path(args.out_metadata)
    out_faa.parent.mkdir(parents=True, exist_ok=True)
    out_meta.parent.mkdir(parents=True, exist_ok=True)

    with out_faa.open("w") as fh:
        for chunk in selected_fasta_chunks:
            fh.write(chunk)

    with out_meta.open("w") as fh:
        fh.write("label\taccession\torganism\tprotein\tclade\tsource\tnote\tstatus\n")
        for row in selected_rows:
            fh.write("\t".join(row[k] for k in [
                "label", "accession", "organism", "protein", "clade", "source", "note", "status"
            ]) + "\n")

    print(f"[DONE] Wrote {len(selected_rows)} selected Wells references")
    print(f"  FASTA   : {out_faa}")
    print(f"  Metadata: {out_meta}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
