# CLAUDE.md — PsrABC Identification Pipeline

## Project Overview

A bioinformatics pipeline to identify and classify **PsrABC (polysulfide reductase)** operons from Pyrodigal-annotated metagenomic genome bins, and distinguish them from functional homologues sharing the Mo-bisPGD catalytic domain:

| Target | Clade | Key distinguishing features |
|--------|-------|-----------------------------|
| **PsrA** ← find these | PsrA | TAT signal + PsrC (8TM, NrfD) in neighbourhood |
| PhsA | PhsA | TAT + PhsC (5TM + b-haem) |
| TtrA | TtrA | TAT + TtrC (9TM) |
| SoeA | SoeA | **No TAT**, cytoplasmic, no PsrC neighbour |
| NapA/NarG/FdhG | various | outgroups / non-sulfur |

The pipeline operates in three discovery modes:

- **Mode A**: supplied list of pre-annotated psrA protein IDs
- **Mode B**: full PF00384 scan across all bins (catches annotation gaps)
- **Mode C** (recommended): dual-gate — PF00384 **UNION** PsrAPhsASreA.hmm (HMSS2). Catches divergent PsrA sequences that fall below PF00384's detection threshold. Requires `--hmss2`.

In Mode C, each candidate is tagged with a `discovery_source` (PF00384_only | HMSS2_only | both) which flows to the final table and HTML output. `HMSS2_only` candidates are highlighted in orange for manual review.

---

## Architecture

```
00_pipeline.sh                    # Main orchestrator (bash)
scripts/
  01_extract_neighbourhood.py     # Extract ±10 gene windows from GFF/FAA
  02_parse_hmmer.py               # Parse HMMER hits (PsrC, PsrB, Mo-bisPGD)
  02b_parse_hmss2.py              # Parse HMSS2 HMM hits — SoeA.hmm scored; rest annotation only
  03_parse_topology.py            # Classify TM helix count from DeepTMHMM output
  04_parse_signalp.py             # Parse SignalP 6 TAT predictions
  05_fetch_references.py          # Download reference sequences from UniProt/NCBI
  06_build_summary.py             # Integrate all evidence → classification table
```

### Pipeline steps

| Step | Name | Notes |
|------|------|-------|
| 0a | Build protein→bin index | Always runs; fast |
| 0b | PF00384 scan (Mode B/C) | Skipped in Mode A |
| 0c | HMSS2 pre-screen (Mode C only) | PsrAPhsASreA.hmm vs all proteins; union with PF00384 |
| 1  | Extract candidate sequences | seqkit from combined FAA |
| 2  | Genomic neighbourhood extraction | ±10 genes, contig-boundary aware |
| 3  | HMMER neighbourhood searches | PF03916, PF14589, PF12800, PF13247 |
| 3b | HMSS2 HMM annotation searches | Optional; requires --hmss2 |
| 4  | DeepTMHMM | May require manual run |
| 5  | SignalP 6 | May require manual run |
| 6  | Download reference sequences | UniProt/NCBI |
| 7  | MAFFT + TrimAl | --anysymbol required for FdhH selenocysteine |
| 8  | IQ-TREE | Full (-bb 1000 -alrt 1000) or --fast-tree |
| 9  | Classification summary | 06_build_summary.py |

### Key architectural decisions

1. **Protein→bin index** (`00_scan/protein_to_bin_index.tsv`): Pyrodigal names proteins after contigs (e.g. `NODE_3193_length_16941_cov_1.235637_2`), NOT after bins. The bin name is only the containing directory. The index maps every protein_id → (bin_name, faa_path, gff_path) and is built once at startup by scanning all FAA files.

2. **Dual-gate discovery (Mode C)**: Step 0c runs `PsrAPhsASreA.hmm` (from HMSS2) across all proteins in parallel with the PF00384 scan. The final candidate set is the union of both. `discovery_source.tsv` records the origin of each candidate. This is the primary mechanism for recovering divergent PsrA sequences missed by PF00384.

3. **Sentinel system** (`logs/stepNN.done`): Re-runs skip completed steps. Manual steps (DeepTMHMM, SignalP) are recognised as done if their output directories/files exist. Use `--redo-from N` or `--redo-step N` to force re-runs. Step names may be alphanumeric (e.g. `0c`, `3b`) — the `_sentinel_name()` helper handles zero-padding integers vs. passing through alphanumeric strings unchanged. `step_clear` handles `0c` and `3b` explicitly.

4. **Classification scoring** (in `06_build_summary.py → score_classification()`): Evidence is scored and summed. Tree score is deliberately capped at +2 to prevent phylogeny overriding biochemical evidence. **SoeA.hmm** hit (from HMSS2 Step 3b) applies a −3 score penalty — this is the **only** HMSS2 result that enters scoring. All other HMSS2 hits are annotation columns only.

5. **NrfD deduplication**: When multiple NrfD neighbours exist, the best one is selected by: PF14589 presence first (high-specificity), then topology classification confidence.

6. **HMSS2 integration (Steps 0c and 3b)**: Step 0c uses `PsrAPhsASreA.hmm` as a discovery gate. Step 3b runs the full HMSS2 profile set for annotation. `SoeA.hmm` is the only Step 3b profile that feeds back into scoring because it provides positive evidence for SoeA identity (converting ambiguous low-evidence cases from `AMBIGUOUS` to `LIKELY_SoeA_or_divergent`).

---

## Current State

### Implemented and working
- [x] Mode B full PF00384 scan across bins
- [x] **Mode C dual-gate discovery** (Step 0c) — PF00384 ∪ PsrAPhsASreA.hmm, `discovery_source` column in output
- [x] Protein→bin index building
- [x] Genomic neighbourhood extraction (GFF-based, contig-boundary aware)
- [x] HMMER searches: PF00384, PF03916, PF14589, PF12800, PF13247
- [x] DeepTMHMM integration (requires `cd` to install dir; output dir must not pre-exist)
- [x] SignalP 6 integration (organism=`other`, mode=`slow-sequential`)
- [x] Reference sequence download (18 sequences; two accessions corrected — see table below)
- [x] MAFFT alignment (`--anysymbol` for selenocysteine in FdhH)
- [x] TrimAl trimming
- [x] IQ-TREE (fast and full modes, separate prefixes/sentinels)
- [x] Classification summary with bin_name column, HTML output
- [x] **Step 3b: HMSS2 annotation HMM searches** (`--hmss2` flag, `02b_parse_hmss2.py`)
- [x] **SoeA.hmm scoring** — −3 penalty applied in `score_classification()`; converts ambiguous candidates to `LIKELY_SoeA_or_divergent`
- [x] **`discovery_source` column** — PF00384_only / HMSS2_only / both / supplied; orange HTML highlight for HMSS2_only
- [x] Sentinel system for alphanumeric step names (`_sentinel_name()` helper; `step_clear` handles `0c` and `3b`)
- [x] NrfD label clarified: `NrfD_unclassified` renamed `NrfD_PF03916(topology_ND)`
- [x] Reference accessions corrected: SreA Q9HGX4→Q8NKK1, ArrA_Chrysiogenes AAD05290.1→Q5Y818

---

## Tech Stack

| Tool | Version/Notes | Purpose |
|------|--------------|---------|
| Python | 3.x | All helper scripts |
| HMMER | `hmmsearch`, `hmmfetch`, `hmmpress` | Profile searches |
| seqkit | bioconda | Sequence extraction |
| MAFFT | bioconda, `--anysymbol --localpair` | Multiple sequence alignment |
| TrimAl | bioconda, `-automated1` | Alignment trimming |
| IQ-TREE | `iqtree` (not `iqtree2`) | ML phylogeny |
| DeepTMHMM | Academic licence, `predict.py` | TM topology prediction |
| SignalP 6 | `signalp6`, `--organism other` | TAT signal prediction |
| Pfam-A.hmm | Full database, pressed | Pfam HMM source |
| HMSS2 | Local HMM directory, `--hmss2` flag | Sulfur metabolism–specific HMMs |
| ete3 | optional, conda-forge | Tree clade assignment |
| biopython + requests | pip | Reference downloads |

### Pfam profiles used

| Profile | What it detects | Specificity |
|---------|----------------|-------------|
| PF00384 | Mo-bisPGD catalytic domain | All Mo-bisPGD enzymes (broad) |
| PF03916 | NrfD-like membrane subunit | PsrC, TtrC, PhsC (broad) |
| **PF14589** | NrfD_2 — polysulfide reductase subunit | PsrC-specific (HIGH) |
| PF12800 | NrfC-like 4Fe-4S | PsrB electron transfer subunit |
| PF13247 | 4Fe-4S dicluster | PsrB, broad backup |

### HMSS2 profiles used

Profiles run in two contexts:

**Step 0c (discovery gate — Mode C only)**

| Profile | Run against | Role |
|---------|-------------|------|
| PsrAPhsASreA.hmm | all_proteins_combined.faa | Second discovery gate; union with PF00384 hits |

**Step 3b (annotation + one scored profile)**

| Profile | Run against | Role |
|---------|-------------|------|
| PsrAPhsASreA.hmm | candidates.faa | Annotation column |
| **SoeA.hmm** | candidates.faa | **SCORED** — −3 penalty; positive SoeA evidence |
| TtrA.hmm | candidates.faa | Annotation column |
| PsrBPhsBSreB.hmm | neighbours.faa | Annotation column |
| PsrCPhsCSreC.hmm | neighbours.faa | Annotation column (cross-validation; Pfam outperforms) |
| SoeB.hmm | neighbours.faa | Annotation column |
| SoeC.hmm | neighbours.faa | Annotation column |
| TtrB.hmm | neighbours.faa | Annotation column |
| TtrC.hmm | neighbours.faa | Annotation column |

### Reference sequences

> ⚠️ **Two accessions were corrected.** Re-run `--redo-from 6` after verifying the corrected accessions manually before using results downstream.

| Label | Accession | Clade | Notes |
|-------|-----------|-------|-------|
| PsrA_Wolinella_succinogenes | **P31075** (NOT P31077 = PsrC!) | PsrA | |
| PsrA_Thermus_thermophilus | Q72LA4 | PsrA | PDB 2VPZ |
| PhsA_Salmonella_enterica_LT2 | P37600 | PhsA | |
| TtrA_Salmonella_enterica_LT2 | Q9Z4S6 | TtrA | |
| TtrA_Shewanella_ANA3 | WP_011715816.1 | TtrA | NCBI |
| SreA_Acidianus_ambivalens | **Q8NKK1** ← corrected | SreA | Verify gene=`sreA`, organism=*Acidianus ambivalens* at UniProt before use |
| SoeA_Allochromatium_vinosum | D3RNN8 | SoeA | |
| ArrA_Shewanella_ANA3 | Q7WTU0 | ArrA | |
| ArrA_Chrysiogenes_arsenatis | **Q5Y818** ← corrected | ArrA | Verify against Krafft & Macy 1998 before use |
| AioA_Alcaligenes_faecalis | Q7SIF4 | AioA | |
| NapA_Shewanella_oneidensis | Q8EIJ1 | NapA | |
| TorA_Ecoli | P33225 | TorA | |
| DmsA_Ecoli | P18775 | DmsA | |
| NarG_Ecoli | P09152 | NarG | |
| SerA_Thauera_selenatis | Q9S1H0 | SerA | |
| PcrA_Dechloromonas_aromatica | Q47CW6 | PcrA | |
| FdhG_Ecoli (outgroup) | P24183 | outgroup | |
| FdhH_Ecoli (outgroup) | P07658 | outgroup | contains selenocysteine U — requires `--anysymbol` in MAFFT |

---

## Next Steps

1. **First Mode C run**: add `--hmss2 /path/to/HMSS2/Hidden_Markov_Models/Inorganic_Sulfur_Metabolism` and run from scratch (or `--redo-from 0`). Inspect `09_summary/classification_table.html` — orange-bordered rows are `HMSS2_only` candidates; verify these against tree placement and neighbourhood evidence.

2. **Verify corrected reference accessions** before relying on tree-based clade assignments:
   - SreA: https://www.uniprot.org/uniprotkb/Q8NKK1/entry — confirm gene=`sreA`, organism=*Acidianus ambivalens*
   - ArrA_Chrysiogenes: https://www.uniprot.org/uniprotkb/Q5Y818/entry — confirm arsenate reductase from *Chrysiogenes arsenatis*
   - Then run `--redo-from 6` to re-download and rebuild the tree.

3. **Install ete3** (`conda install -c etetoolkit ete3`) for automatic tree clade assignments in Step 9.

4. **Validate SoeA.hmm scoring**: after first Mode C run, review whether the −3 SoeA.hmm penalty is correctly shifting `AMBIGUOUS` candidates to `LIKELY_SoeA_or_divergent` without overcalling. Examine `HMSS2_SOEA_*` columns in the output table.

5. **Consider additional PsrA reference sequences**: current tree has only 2 true PsrA anchors (*W. succinogenes* P31075, *T. thermophilus* Q72LA4). Adding more would strengthen clade assignments. When adding references, use `--redo-from 7` — pipeline does not auto-detect reference changes.

---

## Critical Context & Edge Cases

### ID format (do not break this)
- Protein IDs are **contig-based**: `NODE_3193_length_16941_cov_1.235637_2`
- The **last field only** (`_2`) is the gene number — strip only this to get the contig name
- The bin name is the **directory** containing the FAA/GFF, e.g. `SF2967_J5348.rosella_refine.165_0`
- Bin name **never appears** in the protein ID — always use the protein→bin index

### Sentinel system
- `step_done` / `step_mark` route through `_sentinel_name()` which zero-pads plain integers but passes alphanumeric names through unchanged
- Sentinel files: integers → `step03.done`, alphanumeric → `step0c.done`, `step3b.done`
- `step_clear` explicitly handles `0c` (cleared when redoing from step ≤0) and `3b` (cleared when redoing from step ≤3)
- **Do not use `printf '%02d'`** anywhere for step names — always use `_sentinel_name()`

### HMSS2 discovery gate (Step 0c) design
- Step 0c runs **only in Mode C** (--hmss2 supplied AND --ids omitted)
- It runs `PsrAPhsASreA.hmm` at the same E-value as PF00384 (`--evalue`, default 1e-5)
- `discovery_source.tsv` is written to `00_scan/` and passed to Step 9 via `--discovery_source`
- In Mode A or B, `discovery_source.tsv` is still written (values: `supplied` or `PF00384_only`) so the column is always present in the output table

### SoeA.hmm scoring design
- `SoeA.hmm` is the **only** HMSS2 profile that enters `score_classification()`
- Penalty is −3, applied before phylogeny scoring
- Rationale: converts `AMBIGUOUS` (score ~1–4, no TAT, no NrfD) to `LIKELY_SoeA_or_divergent` (score drops below the SoeA classification threshold)
- Does NOT override strong positive evidence — a TRUE_PsrA (score ≥8) with a SoeA.hmm hit would score ≥5 and remain at least LIKELY_PsrA
- All other HMSS2 profiles are annotation columns; they appear in output but cannot affect classification

### Pfam HMM accessions
- `hmmfetch Pfam-A.hmm PF00384` **fails** — must use versioned form `PF00384.28`
- The `pfam_fetch()` bash function handles this via `grep -m1 "^ACC\s*${acc}\."`

### NrfD evidence label meanings (in classification_table)
- `NrfD_PF14589(+)` — NrfD neighbour hit the high-specificity PsrC profile; strong PsrC evidence
- `NrfD_PF03916(topology_ND)` — broad HMM hit only; DeepTMHMM not yet run or result ambiguous
- When topology IS resolved, label comes from DeepTMHMM classification (e.g. `PsrC`, `TtrC`, `PhsC`)

### DeepTMHMM quirks
- Must `cd` to the install directory before running — model files are loaded relative to CWD
- Output directory must **not exist** before running — remove it before re-running
- Uses absolute paths for both `--fasta` and `--output-dir`

### SignalP 6 quirks
- Correct flag: `--organism other` (not `gram-`, not `gram+`)
- Correct flag: `--mode slow-sequential` (not `slow`)
- Output file: `prediction_results.txt`
- ID in output = full FASTA header string; parse only the first whitespace-delimited token

### MAFFT requires `--anysymbol`
- FdhH (P07658) contains selenocysteine (U) — MAFFT errors without this flag

### IQ-TREE binary name
- On this system: `iqtree` (not `iqtree2`)
- Pipeline tries `iqtree`, `iqtree3`, `iqtree2` in order

### IQ-TREE support string format
- Output with `-bb 1000 -alrt 1000` produces node labels like `95.6/100`
- ete3 must use `format=1` (or try formats 0–3 in sequence) to parse these

### Scoring thresholds
- `TRUE_PsrA` ≥ 8 points
- `LIKELY_PsrA` ≥ 5 points
- Tree evidence capped at **+2** (intentional — prevents phylogeny overriding biochemical evidence)
- **SoeA classification requires**: no TAT AND no NrfD in neighbourhood AND no PsrC topology — all three needed to prevent false SoeA calls from incomplete genomes

### PsrB evidence design
- PsrB is **supporting evidence only** — its absence never penalises classification
- Rationale: CISM four-cluster-protein subunits show high sequence divergence across subfamilies (Rothery et al. 2008, BBA 1778:1897); non-canonical Fe-S binding motifs reduce HMM sensitivity further (Driscoll et al. 2018, Chem. Sci.)
- `operon_completeness.tsv` reports `ABC_complete`, `AC_only`, `AB_only`, `A_only`

### Output files for manual review
After a full run, the most useful files are:
- `09_summary/classification_table.html` — colour-coded; HMSS2_only rows in orange; HMSS2 annotation columns in grey on the right
- `09_summary/classification_table.tsv` — for R/Python analysis; includes `discovery_source` column
- `00_scan/discovery_source.tsv` — full breakdown of which gate found each candidate
- `08_tree/psr_phylogeny.treefile` — load in iTOL or FigTree
- `03_hmmer/operon_completeness.tsv` — operon structure per candidate
- `03_hmmer/nrfd_hits.tsv` — includes `PF14589_high_specificity` column
- `03_hmmer/hmss2/hmss2_operon.tsv` — HMSS2 annotation per candidate (if `--hmss2` was used)
