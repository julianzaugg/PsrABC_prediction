#!/usr/bin/env bash
# =============================================================================
# PsrABC Identification and Phylogenetic Pipeline
# =============================================================================
# PURPOSE:
#   Identifies and classifies PsrABC (polysulfide reductase) operons from
#   Pyrodigal-annotated genome bins, distinguishing true PsrA from functional
#   homologues (PhsA, TtrA, SoeA, SreA) that share the Mo-bisPGD domain.
#
#   THREE INPUT MODES:
#   A) --ids supplied         : classify a specific list of pre-annotated candidates
#   B) --ids omitted          : scan ALL bins with PF00384 (Mo-bisPGD) only
#   C) --ids omitted + --hmss2: dual-gate scan — PF00384 UNION PsrAPhsASreA.hmm
#                               (RECOMMENDED — catches divergent PsrA missed by PF00384)
#
#   In Mode C, Step 0c runs PsrAPhsASreA.hmm across all proteins before the
#   neighbourhood extraction. The final candidate set is the UNION of PF00384
#   and HMSS2 hits. Each candidate is tagged with its discovery_source
#   (PF00384_only | HMSS2_only | both) which flows through to the final table.
#   Candidates discovered only by HMSS2 (PsrAPhsASreA.hmm) enter the same
#   scoring pipeline as PF00384 candidates — no score penalty is applied —
#   but their source is visible in classification_table.tsv for manual review.
#
#   RE-RUN CONTROL (sentinel-based):
#   Each step writes logs/stepNN.done on success. On re-run, steps with an
#   existing sentinel are skipped. Manual steps (DeepTMHMM, SignalP) are
#   considered done if their output files exist, regardless of sentinel.
#
#   --redo-from N  : clear sentinels for steps N..9 and re-run from step N
#   --redo-step N  : clear sentinel for step N only (and downstream N..9)
#
#   IQ-TREE MODES:
#   Default     : -bb 1000 -alrt 1000 (UFBoot + SH-aLRT, full accuracy)
#   --fast-tree : -fast (NNI only, no bootstraps — quick topology check)
#   Fast and full trees use different prefixes so they do not overwrite each
#   other. Step 9 summary always uses whichever tree was most recently built.
#
# BIN DIRECTORY STRUCTURE (expected):
#   BIN_DIR/BinName/BinName.faa  BinName.fna  BinName.gff
#
# PROTEIN ID FORMAT:
#   Pyrodigal names proteins after CONTIGS, not bins. A protein->bin index
#   is built at Step 0a. The bin name does NOT appear in protein IDs.
#
# PREREQUISITES:
#   conda install -c bioconda hmmer seqkit mafft trimal iqtree biopython requests
#
# USAGE:
#   # Full run, dual-gate discovery (Mode C — recommended)
#   bash 00_pipeline.sh \
#       --bindir    /path/to/bins \
#       --pfam      /path/to/Pfam-A.hmm \
#       --hmss2     /path/to/HMSS2/Hidden_Markov_Models/Inorganic_Sulfur_Metabolism \
#       --deeptmhmm /path/to/DeepTMHMM-Academic-License-v1.0 \
#       --outdir    results/psr_analysis
#
#   # PF00384-only scan (Mode B)
#   bash 00_pipeline.sh \
#       --bindir /path/to/bins --pfam /path/to/Pfam-A.hmm --outdir results/psr_analysis
#
#   # Re-run from step 7 (e.g. after adding references)
#   bash 00_pipeline.sh [same args] --redo-from 7
#
#   # Re-run only step 8 with fast mode
#   bash 00_pipeline.sh [same args] --redo-step 8 --fast-tree
#
# =============================================================================
set -euo pipefail

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; CYAN='\033[0;36m'; NC='\033[0m'
log()  { echo -e "${GREEN}[$(date '+%H:%M:%S')] $*${NC}"; }
warn() { echo -e "${YELLOW}[WARN] $*${NC}"; }
err()  { echo -e "${RED}[ERROR] $*${NC}" >&2; exit 1; }
skip() { echo -e "${CYAN}[$(date '+%H:%M:%S')] SKIP  $*${NC}"; }

# ---- Arguments ---------------------------------------------------------------
QUERY_IDS=""
BIN_DIR=""
PFAM_DB=""
OUTDIR="psr_analysis_$(date +%Y%m%d)"
THREADS=8
NEIGHBOURHOOD=10
PF00384_EVALUE="1e-5"
DEEPTMHMM_DIR=""
REDO_FROM=""
REDO_STEP=""
FAST_TREE=0
HMSS2_DIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ids)           QUERY_IDS="$2";        shift 2 ;;
        --bindir)        BIN_DIR="$2";           shift 2 ;;
        --pfam)          PFAM_DB="$2";           shift 2 ;;
        --outdir)        OUTDIR="$2";            shift 2 ;;
        --threads)       THREADS="$2";           shift 2 ;;
        --neighbourhood) NEIGHBOURHOOD="$2";     shift 2 ;;
        --evalue)        PF00384_EVALUE="$2";    shift 2 ;;
        --deeptmhmm)     DEEPTMHMM_DIR="$2";     shift 2 ;;
        --hmss2)         HMSS2_DIR="$2";         shift 2 ;;
        --redo-from)     REDO_FROM="$2";         shift 2 ;;
        --redo-step)     REDO_STEP="$2";         shift 2 ;;
        --fast-tree)     FAST_TREE=1;            shift ;;
        *) err "Unknown argument: $1" ;;
    esac
done

[[ -z "$BIN_DIR" ]] && err "Must supply --bindir"
[[ -z "$PFAM_DB" ]] && err "Must supply --pfam (path to Pfam-A.hmm)"
[[ -d "$BIN_DIR" ]] || err "Bin directory not found: $BIN_DIR"
[[ -f "$PFAM_DB" ]] || err "Pfam HMM not found: $PFAM_DB"
[[ -n "$QUERY_IDS" && ! -f "$QUERY_IDS" ]] && err "IDs file not found: $QUERY_IDS"

mkdir -p "$OUTDIR"/{00_scan,01_sequences,02_neighbourhoods,03_hmmer,04_topology,\
05_signalp,06_references,07_alignment,08_tree,09_summary,logs}

LOG="$OUTDIR/logs/pipeline.log"
exec > >(tee -a "$LOG") 2>&1

# =============================================================================
# Sentinel system
# =============================================================================
SENTINEL_DIR="$OUTDIR/logs"

_sentinel_name() {
    # If the argument is a plain integer, zero-pad it; otherwise use as-is.
    # e.g. 3 → step03.done   3b → step3b.done   0c → step0c.done
    local s="$1"
    if [[ "$s" =~ ^[0-9]+$ ]]; then
        printf "step%02d.done" "$s"
    else
        printf "step%s.done" "$s"
    fi
}
step_done() {
    [[ -f "$SENTINEL_DIR/$(_sentinel_name "$1")" ]]
}
step_mark() {
    touch "$SENTINEL_DIR/$(_sentinel_name "$1")"
    log "  ✓ Step $1 complete"
}
step_clear() {
    local from="$1" cleared=()
    for n in $(seq "$from" 9); do
        local f="$SENTINEL_DIR/step$(printf '%02d' "$n").done"
        [[ -f "$f" ]] && { rm "$f"; cleared+=("$n"); }
    done
    # Clear alphanumeric step sentinels when redoing from their parent step or earlier
    [[ "$from" -le 0 && -f "$SENTINEL_DIR/step0c.done" ]] && \
        rm "$SENTINEL_DIR/step0c.done" && cleared+=("0c")
    [[ "$from" -le 3 && -f "$SENTINEL_DIR/step3b.done" ]] && \
        rm "$SENTINEL_DIR/step3b.done" && cleared+=("3b")
    # Also clear fast-tree sentinel if redoing step 8 or later
    [[ "$from" -le 8 && -f "$SENTINEL_DIR/step8f.done" ]] && \
        rm "$SENTINEL_DIR/step8f.done" && cleared+=("8f")
    [[ ${#cleared[@]} -gt 0 ]] && log "  Cleared sentinels: ${cleared[*]}"
}

# Apply redo flags before anything else runs
if [[ -n "$REDO_FROM" ]]; then
    log "--- --redo-from $REDO_FROM : clearing steps $REDO_FROM..9 ---"
    step_clear "$REDO_FROM"
elif [[ -n "$REDO_STEP" ]]; then
    log "--- --redo-step $REDO_STEP : clearing steps $REDO_STEP..9 ---"
    step_clear "$REDO_STEP"
fi

# Determine discovery mode label for logging
if [[ -n "$QUERY_IDS" ]]; then
    MODE_LABEL="A (supplied IDs)"
elif [[ -n "$HMSS2_DIR" ]]; then
    MODE_LABEL="C (PF00384 + HMSS2 dual-gate)"
else
    MODE_LABEL="B (PF00384 scan only)"
fi

log "=== PsrABC Identification Pipeline ==="
log "Mode         : $MODE_LABEL"
log "Bin dir      : $BIN_DIR"
log "Pfam HMM     : $PFAM_DB"
log "Output dir   : $OUTDIR"
log "Threads      : $THREADS"
log "Neighbourhood: ±${NEIGHBOURHOOD} genes"
log "DeepTMHMM    : ${DEEPTMHMM_DIR:-not supplied}"
log "HMSS2 HMMs   : ${HMSS2_DIR:-not supplied}"
log "IQ-TREE mode : $([ "$FAST_TREE" -eq 1 ] && echo 'FAST (-fast)' || echo 'FULL (-bb 1000 -alrt 1000)')"

# =============================================================================
# Pfam fetch helper
# =============================================================================
pfam_fetch() {
    local acc="$1" out="$2" versioned
    versioned=$(grep -m1 "^ACC\s*${acc}\." "$PFAM_DB" | awk '{print $2}')
    [[ -z "$versioned" ]] && { warn "  $acc not found in $PFAM_DB"; return 1; }
    hmmfetch "$PFAM_DB" "$versioned" > "$out" && [[ -s "$out" ]] && return 0
    warn "  hmmfetch failed for $versioned"; rm -f "$out"; return 1
}

# =============================================================================
# STEP 0a: Build protein->bin index (always runs — fast, prerequisite for all)
# =============================================================================
log "--- STEP 0a: Building protein->bin index ---"

PROTEIN_INDEX="$OUTDIR/00_scan/protein_to_bin_index.tsv"
ALL_FAA_LIST="$OUTDIR/00_scan/all_faa_files.txt"

find "$BIN_DIR" -mindepth 2 -maxdepth 2 -name "*.faa" | sort > "$ALL_FAA_LIST"
N_FAA=$(wc -l < "$ALL_FAA_LIST")
log "  Found $N_FAA FAA files"
[[ "$N_FAA" -eq 0 ]] && err "No FAA files found — expected: BIN_DIR/BinName/BinName.faa"

python3 - <<'PYEOF' "$ALL_FAA_LIST" "$PROTEIN_INDEX"
import sys, os
faa_list_path, index_path = sys.argv[1], sys.argv[2]
with open(faa_list_path) as fh:
    faa_paths = [l.strip() for l in fh if l.strip()]
written = 0
with open(index_path, "w") as out:
    out.write("protein_id\tbin_name\tfaa_path\tgff_path\n")
    for faa_path in faa_paths:
        d        = os.path.dirname(faa_path)
        bin_name = os.path.basename(d)
        stem     = os.path.splitext(os.path.basename(faa_path))[0]
        gff      = os.path.join(d, stem + ".gff")
        if not os.path.exists(gff):
            print(f"  [WARN] GFF missing: {gff}", file=sys.stderr)
            gff = ""
        with open(faa_path) as faa:
            for line in faa:
                if line.startswith(">"):
                    out.write(f"{line[1:].split()[0]}\t{bin_name}\t{faa_path}\t{gff}\n")
                    written += 1
print(f"  Indexed {written} proteins from {len(faa_paths)} bins")
PYEOF

log "  Index built → $PROTEIN_INDEX"

ALL_PROTEINS_FAA="$OUTDIR/00_scan/all_proteins_combined.faa"
if [[ ! -s "$ALL_PROTEINS_FAA" ]]; then
    log "  Concatenating all FAA files..."
    xargs cat < "$ALL_FAA_LIST" > "$ALL_PROTEINS_FAA"
    log "  $(grep -c "^>" "$ALL_PROTEINS_FAA") total proteins"
else
    log "  Combined FAA exists — skipping concatenation"
fi

# =============================================================================
# STEP 0b (Mode B/C): PF00384 scan across all proteins
# =============================================================================
QUERY_IDS_RESOLVED="$QUERY_IDS"

PF00384_IDS="$OUTDIR/00_scan/PF00384_candidates.txt"
HMSS2_IDS="$OUTDIR/00_scan/HMSS2_prescreen_candidates.txt"
MERGED_IDS="$OUTDIR/00_scan/merged_candidates.txt"
DISCOVERY_SOURCE="$OUTDIR/00_scan/discovery_source.tsv"

if [[ -z "$QUERY_IDS" ]]; then
    if step_done 0; then
        skip "0b — PF00384 scan"
        QUERY_IDS_RESOLVED="$MERGED_IDS"
    else
        log "--- STEP 0b: PF00384 scan across all proteins ---"
        PF00384_HMM_SCAN="$OUTDIR/00_scan/PF00384.hmm"
        pfam_fetch PF00384 "$PF00384_HMM_SCAN" || err "Could not extract PF00384"
        hmmpress "$PF00384_HMM_SCAN" 2>/dev/null || true

        hmmsearch --cpu "$THREADS" \
            --tblout "$OUTDIR/00_scan/PF00384_all_hits.tbl" \
            -E "$PF00384_EVALUE" \
            "$PF00384_HMM_SCAN" "$ALL_PROTEINS_FAA" \
            > "$OUTDIR/00_scan/PF00384_all_hmmsearch.out"

        awk '!/^#/ {print $1}' "$OUTDIR/00_scan/PF00384_all_hits.tbl" \
            | sort -u > "$PF00384_IDS"
        N_PF=$(wc -l < "$PF00384_IDS")
        log "  $N_PF Mo-bisPGD candidates (PF00384) → $PF00384_IDS"
        [[ "$N_PF" -eq 0 ]] && err "No PF00384 hits found. Check --evalue or Pfam DB."

        step_mark 0
        # QUERY_IDS_RESOLVED set after Step 0c below (Mode C) or here (Mode B)
    fi
else
    skip "0b — Mode A (--ids supplied)"
fi

# =============================================================================
# STEP 0c (Mode C only): HMSS2 PsrAPhsASreA pre-screen across all proteins
#
# PURPOSE: Cast a second, independent net using a clade-specific HMM to catch
# divergent PsrA sequences that fall below PF00384's detection threshold.
# The final candidate set entering Step 1 is PF00384 UNION HMSS2 hits.
#
# Candidates are tagged by discovery_source:
#   both        — hit both PF00384 and PsrAPhsASreA.hmm (highest confidence)
#   PF00384_only — hit PF00384 only
#   HMSS2_only  — hit PsrAPhsASreA.hmm only (most important to review manually)
#
# This source tag flows through to classification_table.tsv as a column.
# HMSS2_only candidates are highlighted in the HTML output.
#
# NOTE: Step 0c runs ONLY in Mode C (--hmss2 supplied AND --ids omitted).
#       In Mode A or B, MERGED_IDS == PF00384_IDS and all sources == PF00384_only.
# =============================================================================
if [[ -z "$QUERY_IDS" ]]; then
    if [[ -n "$HMSS2_DIR" && -d "$HMSS2_DIR" ]]; then
        if step_done 0c; then
            skip "0c — HMSS2 PsrAPhsASreA pre-screen"
            QUERY_IDS_RESOLVED="$MERGED_IDS"
        else
            log "--- STEP 0c: HMSS2 PsrAPhsASreA pre-screen (Mode C dual-gate) ---"
            PSRA_HMM="$HMSS2_DIR/PsrAPhsASreA.hmm"
            [[ -f "$PSRA_HMM" ]] || err "PsrAPhsASreA.hmm not found in $HMSS2_DIR"

            hmmsearch --cpu "$THREADS" \
                --tblout "$OUTDIR/00_scan/HMSS2_PsrAPhsASreA_all_hits.tbl" \
                -E "$PF00384_EVALUE" \
                "$PSRA_HMM" "$ALL_PROTEINS_FAA" \
                > "$OUTDIR/00_scan/HMSS2_PsrAPhsASreA_hmmsearch.out"

            awk '!/^#/ {print $1}' "$OUTDIR/00_scan/HMSS2_PsrAPhsASreA_all_hits.tbl" \
                | sort -u > "$HMSS2_IDS"
            N_HMSS2=$(wc -l < "$HMSS2_IDS")
            log "  $N_HMSS2 candidates from PsrAPhsASreA.hmm → $HMSS2_IDS"

            # Merge PF00384 + HMSS2 hits; build discovery_source table
            python3 - <<'PYEOF' "$PF00384_IDS" "$HMSS2_IDS" "$MERGED_IDS" "$DISCOVERY_SOURCE"
import sys
pf_path, hmss2_path, merged_path, src_path = sys.argv[1:]

with open(pf_path)    as fh: pf_ids    = set(l.strip() for l in fh if l.strip())
with open(hmss2_path) as fh: hmss2_ids = set(l.strip() for l in fh if l.strip())

all_ids = pf_ids | hmss2_ids
n_both       = len(pf_ids & hmss2_ids)
n_pf_only    = len(pf_ids - hmss2_ids)
n_hmss2_only = len(hmss2_ids - pf_ids)

with open(merged_path, "w") as fh:
    for pid in sorted(all_ids):
        fh.write(pid + "\n")

with open(src_path, "w") as fh:
    fh.write("protein_id\tdiscovery_source\n")
    for pid in sorted(all_ids):
        in_pf    = pid in pf_ids
        in_hmss2 = pid in hmss2_ids
        if in_pf and in_hmss2:
            src = "both"
        elif in_pf:
            src = "PF00384_only"
        else:
            src = "HMSS2_only"
        fh.write(f"{pid}\t{src}\n")

print(f"  PF00384 hits     : {len(pf_ids)}")
print(f"  HMSS2 hits       : {len(hmss2_ids)}")
print(f"  Both             : {n_both}")
print(f"  PF00384_only     : {n_pf_only}")
print(f"  HMSS2_only       : {n_hmss2_only}  ← new candidates not in PF00384 scan")
print(f"  Total (union)    : {len(all_ids)}")
PYEOF

            N_MERGED=$(wc -l < "$MERGED_IDS")
            log "  $N_MERGED total candidates after union → $MERGED_IDS"
            log "  Discovery source table → $DISCOVERY_SOURCE"
            step_mark 0c
            QUERY_IDS_RESOLVED="$MERGED_IDS"
        fi
    else
        # Mode B: no HMSS2, PF00384 only — write trivial discovery_source and merged
        if [[ ! -f "$MERGED_IDS" ]]; then
            cp "$PF00384_IDS" "$MERGED_IDS"
            python3 - <<'PYEOF' "$PF00384_IDS" "$DISCOVERY_SOURCE"
import sys
pf_path, src_path = sys.argv[1:]
with open(pf_path) as fh:
    ids = [l.strip() for l in fh if l.strip()]
with open(src_path, "w") as fh:
    fh.write("protein_id\tdiscovery_source\n")
    for pid in ids:
        fh.write(f"{pid}\tPF00384_only\n")
PYEOF
            log "  Mode B: discovery_source set to PF00384_only for all candidates"
        fi
        QUERY_IDS_RESOLVED="$MERGED_IDS"
        [[ -n "$HMSS2_DIR" ]] && warn "HMSS2 dir not found: $HMSS2_DIR — running as Mode B"
    fi
else
    # Mode A: --ids supplied; write trivial discovery_source from supplied IDs
    if [[ ! -f "$DISCOVERY_SOURCE" ]]; then
        python3 - <<'PYEOF' "$QUERY_IDS" "$DISCOVERY_SOURCE"
import sys
ids_path, src_path = sys.argv[1:]
with open(ids_path) as fh:
    ids = [l.strip() for l in fh if l.strip() and not l.startswith("#")]
with open(src_path, "w") as fh:
    fh.write("protein_id\tdiscovery_source\n")
    for pid in ids:
        fh.write(f"{pid}\tsupplied\n")
PYEOF
    fi
fi

# =============================================================================
# STEP 1: Extract candidate sequences
# =============================================================================
CANDIDATES_FAA="$OUTDIR/01_sequences/candidates_psrA.faa"

if step_done 1; then
    skip "1 — sequence extraction"
else
    log "--- STEP 1: Extracting candidate sequences ---"
    seqkit grep --pattern-file "$QUERY_IDS_RESOLVED" "$ALL_PROTEINS_FAA" -w 0 \
        > "$CANDIDATES_FAA" 2> "$OUTDIR/logs/seqkit_extract.log"
    N_CAND=$(grep -c "^>" "$CANDIDATES_FAA" 2>/dev/null || true); N_CAND=${N_CAND:-0}
    N_QUERY=$(grep -c "^[^#]" "$QUERY_IDS_RESOLVED" 2>/dev/null || true); N_QUERY=${N_QUERY:-0}
    log "  Extracted $N_CAND / $N_QUERY sequences"
    [[ -s "$OUTDIR/logs/seqkit_extract.log" ]] && \
        warn "  seqkit: $(cat "$OUTDIR/logs/seqkit_extract.log")"
    [[ "$N_CAND" -eq 0 ]] && err "No sequences extracted. Check IDs and --bindir."
    step_mark 1
fi

# =============================================================================
# STEP 2: Genomic neighbourhood extraction
# =============================================================================
if step_done 2; then
    skip "2 — neighbourhood extraction"
else
    log "--- STEP 2: Genomic neighbourhoods (±${NEIGHBOURHOOD} genes) ---"
    python3 scripts/01_extract_neighbourhood.py \
        --ids    "$QUERY_IDS_RESOLVED" \
        --index  "$PROTEIN_INDEX" \
        --window "$NEIGHBOURHOOD" \
        --outdir "$OUTDIR/02_neighbourhoods"
    step_mark 2
fi

# =============================================================================
# STEP 3: HMMER searches
# =============================================================================
if step_done 3; then
    skip "3 — HMMER neighbourhood search"
else
    log "--- STEP 3: HMMER searches (PF03916/PF12800/PF13247/PF00384) ---"

    for acc in PF03916 PF14589 PF00384 PF12800 PF13247; do
        out="$OUTDIR/03_hmmer/${acc}.hmm"
        pfam_fetch "$acc" "$out" && hmmpress "$out" 2>/dev/null || true
    done

    ALL_NEIGHBOURS="$OUTDIR/03_hmmer/all_neighbours.faa"
    cat "$OUTDIR/02_neighbourhoods"/*_neighbourhood.faa > "$ALL_NEIGHBOURS" 2>/dev/null || \
        warn "No neighbourhood files found"

    for acc in PF03916 PF14589 PF12800 PF13247; do
        hmm="$OUTDIR/03_hmmer/${acc}.hmm"
        [[ -s "$hmm" ]] && hmmsearch --cpu "$THREADS" \
            --tblout "$OUTDIR/03_hmmer/${acc}_hits.tbl" -E 1e-5 \
            "$hmm" "$ALL_NEIGHBOURS" \
            > "$OUTDIR/03_hmmer/${acc}_hmmsearch.out"
    done

    hmmsearch --cpu "$THREADS" \
        --tblout "$OUTDIR/03_hmmer/PF00384_candidates_hits.tbl" -E 1e-5 \
        "$OUTDIR/03_hmmer/PF00384.hmm" "$CANDIDATES_FAA" \
        > "$OUTDIR/03_hmmer/PF00384_candidates.out"

    python3 scripts/02_parse_hmmer.py \
        --psrA_hits  "$OUTDIR/03_hmmer/PF00384_candidates_hits.tbl" \
        --nrfd_hits  "$OUTDIR/03_hmmer/PF03916_hits.tbl" \
        --nrfd2_hits "$OUTDIR/03_hmmer/PF14589_hits.tbl" \
        --psrB_hits  "$OUTDIR/03_hmmer/PF12800_hits.tbl" \
        --psrB_hits2 "$OUTDIR/03_hmmer/PF13247_hits.tbl" \
        --neighbourhood_faa "$ALL_NEIGHBOURS" \
        --psrA_faa   "$CANDIDATES_FAA" \
        --ids        "$QUERY_IDS_RESOLVED" \
        --outdir     "$OUTDIR/03_hmmer"
    step_mark 3
fi

# =============================================================================
# STEP 3b: HMSS2 HMM searches (annotation + neighbourhood subunit scoring)
# Profiles searched:
#   vs candidates.faa  : PsrAPhsASreA, SoeA, TtrA
#   vs neighbours.faa  : PsrBPhsBSreB, PsrCPhsCSreC, SoeB, SoeC, TtrB, TtrC
#
# NOTE: SoeA.hmm hits contribute to scoring in 06_build_summary.py (−3 penalty
# applied to candidates with a SoeA.hmm hit). All other HMSS2 hits are
# annotation columns only — not used in scoring.
# =============================================================================
HMSS2_OUT="$OUTDIR/03_hmmer/hmss2"

if [[ -n "$HMSS2_DIR" && -d "$HMSS2_DIR" ]]; then
    if step_done 3b; then
        skip "3b — HMSS2 HMM searches"
    else
        log "--- STEP 3b: HMSS2 HMM searches ---"
        mkdir -p "$HMSS2_OUT"

        # Profiles to run against candidate (catalytic subunit) FAA
        declare -a CAND_PROFILES=(PsrAPhsASreA SoeA TtrA)
        # Profiles to run against neighbourhood FAA
        declare -a NEIGH_PROFILES=(PsrBPhsBSreB PsrCPhsCSreC SoeB SoeC TtrB TtrC)

        for prof in "${CAND_PROFILES[@]}"; do
            hmm="$HMSS2_DIR/${prof}.hmm"
            if [[ -f "$hmm" ]]; then
                hmmsearch --cpu "$THREADS" \
                    --tblout "$HMSS2_OUT/${prof}_hits.tbl" -E 1e-5 \
                    "$hmm" "$CANDIDATES_FAA" \
                    > "$HMSS2_OUT/${prof}_hmmsearch.out" 2>&1
                log "  ${prof}: $(grep -vc '^#' "$HMSS2_OUT/${prof}_hits.tbl" 2>/dev/null || echo 0) hits"
            else
                warn "  ${prof}.hmm not found in $HMSS2_DIR — skipping"
            fi
        done

        for prof in "${NEIGH_PROFILES[@]}"; do
            hmm="$HMSS2_DIR/${prof}.hmm"
            if [[ -f "$hmm" ]]; then
                hmmsearch --cpu "$THREADS" \
                    --tblout "$HMSS2_OUT/${prof}_hits.tbl" -E 1e-5 \
                    "$hmm" "$OUTDIR/03_hmmer/all_neighbours.faa" \
                    > "$HMSS2_OUT/${prof}_hmmsearch.out" 2>&1
                log "  ${prof}: $(grep -vc '^#' "$HMSS2_OUT/${prof}_hits.tbl" 2>/dev/null || echo 0) hits"
            else
                warn "  ${prof}.hmm not found in $HMSS2_DIR — skipping"
            fi
        done

        python3 scripts/02b_parse_hmss2.py \
            --hmss2_dir    "$HMSS2_OUT" \
            --candidates   "$CANDIDATES_FAA" \
            --neighbours   "$OUTDIR/03_hmmer/all_neighbours.faa" \
            --ids          "$QUERY_IDS_RESOLVED" \
            --outdir       "$HMSS2_OUT"

        step_mark 3b
    fi
else
    [[ -n "$HMSS2_DIR" ]] && warn "HMSS2 dir not found: $HMSS2_DIR — skipping Step 3b"
    log "  (--hmss2 not supplied — Step 3b skipped)"
fi

# =============================================================================
# STEP 4: DeepTMHMM
# =============================================================================
NRFD_FAA="$OUTDIR/03_hmmer/nrfd_candidates.faa"
TOPO_OUT="$OUTDIR/04_topology"
DEEPTMHMM_OUTDIR="$(realpath "$TOPO_OUT")/deeptmhmm_out"
NRFD_FAA_ABS="$(realpath "$NRFD_FAA" 2>/dev/null || echo "$NRFD_FAA")"

if step_done 4 || [[ -d "$DEEPTMHMM_OUTDIR" ]]; then
    skip "4 — DeepTMHMM (output exists or sentinel set)"
    [[ ! -f "$SENTINEL_DIR/step04.done" ]] && step_mark 4
else
    log "--- STEP 4: DeepTMHMM TM topology ---"
    if [[ -n "$DEEPTMHMM_DIR" && -f "$DEEPTMHMM_DIR/predict.py" ]]; then
        log "  Running from $DEEPTMHMM_DIR"
        (
            cd "$DEEPTMHMM_DIR"
            python predict.py --fasta "$NRFD_FAA_ABS" --output-dir "$DEEPTMHMM_OUTDIR"
        ) > "$TOPO_OUT/deeptmhmm.log" 2>&1 \
            && { log "  DeepTMHMM complete"; step_mark 4; } \
            || { warn "  DeepTMHMM failed — check $TOPO_OUT/deeptmhmm.log"
                 touch "$TOPO_OUT/.manual_step_required"; }
    else
        warn "  DeepTMHMM not available."
        warn "  MANUAL: cd /path/to/DeepTMHMM-Academic-License-v1.0"
        warn "          python predict.py --fasta $NRFD_FAA_ABS --output-dir $DEEPTMHMM_OUTDIR"
        warn "  Then:   bash $0 [args] --redo-from 4"
        touch "$TOPO_OUT/.manual_step_required"
    fi
fi

python3 scripts/03_parse_topology.py \
    --tmhmm_results "$DEEPTMHMM_OUTDIR" \
    --nrfd_hits     "$OUTDIR/03_hmmer/PF03916_hits.tbl" \
    --ids           "$QUERY_IDS_RESOLVED" \
    --outdir        "$TOPO_OUT"

# =============================================================================
# STEP 5: SignalP 6
# =============================================================================
SIGNALP_RESULT="$OUTDIR/05_signalp/prediction_results.txt"

if step_done 5 || [[ -f "$SIGNALP_RESULT" ]]; then
    skip "5 — SignalP 6 (output exists or sentinel set)"
    [[ ! -f "$SENTINEL_DIR/step05.done" ]] && step_mark 5
else
    log "--- STEP 5: SignalP 6 TAT prediction ---"
    if command -v signalp6 &>/dev/null; then
        signalp6 \
            --fastafile  "$CANDIDATES_FAA" \
            --organism   other \
            --output_dir "$OUTDIR/05_signalp" \
            --format     txt \
            --mode       slow-sequential \
            > "$OUTDIR/05_signalp/signalp6.log" 2>&1 \
            && { log "  SignalP 6 complete"; step_mark 5; } \
            || warn "  SignalP 6 failed — check $OUTDIR/05_signalp/signalp6.log"
    else
        warn "  signalp6 not found."
        warn "  MANUAL: https://services.healthtech.dtu.dk/services/SignalP-6.0/"
        warn "  Upload: $CANDIDATES_FAA  |  Organism: other"
        warn "  Save prediction_results.txt to: $OUTDIR/05_signalp/"
        warn "  Then: bash $0 [args] --redo-from 5"
        touch "$OUTDIR/05_signalp/.manual_step_required"
    fi
fi

python3 scripts/04_parse_signalp.py \
    --signalp_dir "$OUTDIR/05_signalp" \
    --candidates  "$CANDIDATES_FAA" \
    --outdir      "$OUTDIR/05_signalp"

# =============================================================================
# STEP 6: Download reference sequences
# =============================================================================
REFS_FAA="$OUTDIR/06_references/references_all.faa"

if step_done 6 || [[ -f "$REFS_FAA" ]]; then
    skip "6 — reference sequences (already downloaded)"
    [[ ! -f "$SENTINEL_DIR/step06.done" ]] && step_mark 6
else
    log "--- STEP 6: Downloading reference sequences ---"
    python3 scripts/05_fetch_references.py --outdir "$OUTDIR/06_references"
    step_mark 6
fi

# =============================================================================
# STEP 7: MAFFT + TrimAl
# =============================================================================
COMBINED="$OUTDIR/07_alignment/combined_for_alignment.faa"
ALIGNED="$OUTDIR/07_alignment/combined_aligned.faa"
TRIMMED="$OUTDIR/07_alignment/combined_trimmed.faa"

if step_done 7; then
    skip "7 — MAFFT + TrimAl"
else
    log "--- STEP 7: MAFFT alignment + TrimAl ---"
    cat "$CANDIDATES_FAA" "$REFS_FAA" > "$COMBINED"
    N_SEQS=$(grep -c "^>" "$COMBINED" || true)
    log "  $N_SEQS sequences to align"
    [[ "$N_SEQS" -eq 0 ]] && err "combined_for_alignment.faa is empty"

    if [[ "$N_SEQS" -lt 200 ]]; then
        MAFFT_FLAGS=(--anysymbol --localpair --maxiterate 1000)
        log "  MAFFT: L-INS-i"
    else
        MAFFT_FLAGS=(--anysymbol --auto)
        log "  MAFFT: --auto"
    fi

    mafft "${MAFFT_FLAGS[@]}" --thread "$THREADS" --reorder "$COMBINED" \
        > "$ALIGNED" 2> "$OUTDIR/logs/mafft.log"

    N_ALIGNED=$(grep -c "^>" "$ALIGNED" 2>/dev/null || true); N_ALIGNED=${N_ALIGNED:-0}
    if [[ "$N_ALIGNED" -eq 0 ]]; then
        cat "$OUTDIR/logs/mafft.log" >&2; err "MAFFT produced empty output"
    fi
    log "  Aligned: $N_ALIGNED sequences"

    trimal -in "$ALIGNED" -out "$TRIMMED" -automated1 \
        -htmlout "$OUTDIR/07_alignment/trimal_report.html" \
        2> "$OUTDIR/logs/trimal.log"
    log "  TrimAl: $(grep -c "^>" "$TRIMMED" || true) sequences retained"
    step_mark 7
fi

# =============================================================================
# STEP 8: IQ-TREE
# Fast and full trees use separate prefixes and separate sentinels so both
# can coexist in the same output directory.
# =============================================================================
if [[ "$FAST_TREE" -eq 1 ]]; then
    TREE_PREFIX="$OUTDIR/08_tree/psr_phylogeny_fast"
    TREE_SENTINEL_FILE="$SENTINEL_DIR/step8f.done"
    TREE_LOG="$OUTDIR/logs/iqtree_fast.log"
    IQTREE_FLAGS=(-fast)
    TREE_MODE_LABEL="FAST"
else
    TREE_PREFIX="$OUTDIR/08_tree/psr_phylogeny"
    TREE_SENTINEL_FILE="$SENTINEL_DIR/step08.done"
    TREE_LOG="$OUTDIR/logs/iqtree_full.log"
    IQTREE_FLAGS=(-bb 1000 -alrt 1000)
    TREE_MODE_LABEL="FULL"
fi

if [[ -f "$TREE_SENTINEL_FILE" ]]; then
    skip "8 — IQ-TREE $TREE_MODE_LABEL"
else
    log "--- STEP 8: IQ-TREE $TREE_MODE_LABEL ---"
    IQTREE_BIN=""
    for bin in iqtree iqtree3 iqtree2; do
        command -v "$bin" &>/dev/null && { IQTREE_BIN="$bin"; break; }
    done
    [[ -z "$IQTREE_BIN" ]] && err "No IQ-TREE found (tried: iqtree iqtree3 iqtree2)"
    log "  Binary : $IQTREE_BIN"
    log "  Flags  : ${IQTREE_FLAGS[*]}"

    $IQTREE_BIN \
        -s "$TRIMMED" \
        --prefix "$TREE_PREFIX" \
        -m TEST \
        "${IQTREE_FLAGS[@]}" \
        -T "$THREADS" \
        --redo \
        > "$TREE_LOG" 2>&1

    log "  Tree: ${TREE_PREFIX}.treefile"
    touch "$TREE_SENTINEL_FILE"
    log "  ✓ Step 8 ($TREE_MODE_LABEL) complete"
fi

# =============================================================================
# STEP 9: Classification summary
# =============================================================================
if step_done 9; then
    skip "9 — classification summary"
else
    log "--- STEP 9: Classification summary ---"
    python3 scripts/06_build_summary.py \
        --ids              "$QUERY_IDS_RESOLVED" \
        --psrA_faa         "$CANDIDATES_FAA" \
        --nrfd_hits        "$OUTDIR/03_hmmer/PF03916_hits.tbl" \
        --topology_dir     "$OUTDIR/04_topology" \
        --signalp_dir      "$OUTDIR/05_signalp" \
        --treefile         "${TREE_PREFIX}.treefile" \
        --references       "$OUTDIR/06_references/reference_metadata.tsv" \
        --protein_index    "$PROTEIN_INDEX" \
        --discovery_source "$DISCOVERY_SOURCE" \
        --hmss2_dir        "$HMSS2_OUT" \
        --outdir           "$OUTDIR/09_summary"
    step_mark 9
fi

# =============================================================================
# DONE
# =============================================================================
log "=== Pipeline complete ==="
log ""
log "Key outputs:"
log "  Candidates      : $CANDIDATES_FAA"
log "  Discovery source: $DISCOVERY_SOURCE"
log "  Tree            : ${TREE_PREFIX}.treefile"
log "  Summary         : $OUTDIR/09_summary/classification_table.tsv"
log ""
log "Sentinel status:"
for n in 0 1 2 3 4 5 6 7 8 9; do
    f="$SENTINEL_DIR/$(_sentinel_name "$n")"
    [[ -f "$f" ]] \
        && echo -e "  ${GREEN}✓ step $n${NC}" \
        || echo -e "  ${YELLOW}✗ step $n${NC}"
done
[[ -f "$SENTINEL_DIR/step0c.done" ]] && echo -e "  ${GREEN}✓ step 0c (HMSS2 pre-screen)${NC}" \
    || { [[ -n "$HMSS2_DIR" ]] && echo -e "  ${YELLOW}✗ step 0c (HMSS2 pre-screen)${NC}"; }
[[ -f "$SENTINEL_DIR/step3b.done" ]] && echo -e "  ${GREEN}✓ step 3b (HMSS2 annotation)${NC}" \
    || { [[ -n "$HMSS2_DIR" ]] && echo -e "  ${YELLOW}✗ step 3b (HMSS2 annotation)${NC}"; }
[[ -f "$SENTINEL_DIR/step8f.done" ]] && echo -e "  ${GREEN}✓ step 8 fast${NC}"
echo ""
[[ -f "$OUTDIR/04_topology/.manual_step_required" ]] && \
    warn "  [!] DeepTMHMM not run — see Step 4 instructions above"
[[ -f "$OUTDIR/05_signalp/.manual_step_required" ]] && \
    warn "  [!] SignalP not run — see Step 5 instructions above"
