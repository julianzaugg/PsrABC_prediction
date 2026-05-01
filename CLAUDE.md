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
  05_fetch_references.py          # Download curated reference sequences from UniProt/NCBI
  06_build_summary.py             # Integrate all evidence → classification table
```

### Pipeline steps

| Step | Name | Notes |
|------|------|-------|
| 0a | Build protein→bin index | Always runs; fast |
| 0b | PF00384 scan (Mode B/C) | Skipped in Mode A |
| 0c | HMSS2 PsrAPhsASreA pre-screen (Mode C only) | Union with PF00384 hits |
| 1  | Extract candidate sequences | seqkit from combined FAA |
| 2  | Genomic neighbourhood extraction | ±10 genes, contig-boundary aware |
| 3  | HMMER neighbourhood searches | PF03916, PF14589, PF12800, PF13247 |
| 3b | HMSS2 annotation searches | Optional; requires --hmss2 |
| 4  | DeepTMHMM | May require manual run |
| 5  | SignalP 6 | May require manual run |
| 6  | Download curated reference sequences | UniProt/NCBI via `05_fetch_references.py`; also builds `tree_references_all.faa` by appending `--custom-tree-refs` |
| 7  | MAFFT + TrimAl | Aligns candidates against `tree_references_all.faa` |
| 8  | IQ-TREE | Full (-bb 1000 -alrt 1000) or --fast-tree |
| 9  | Classification summary | `06_build_summary.py` |

### Key architectural decisions

1. **Protein→bin index** (`00_scan/protein_to_bin_index.tsv`): Pyrodigal names proteins after contigs (e.g. `NODE_3193_length_16941_cov_1.235637_2`), NOT after bins. The bin name is only the containing directory. The index maps every protein_id → (bin_name, faa_path, gff_path) and is built once at startup by scanning all FAA files.

2. **Dual-gate discovery (Mode C)**: Step 0c runs `PsrAPhsASreA.hmm` (from HMSS2) across all proteins in parallel with the PF00384 scan. The final candidate set is the union of both. `discovery_source.tsv` records the origin of each candidate.

3. **Custom tree references (`--custom-tree-refs`)**: Any FASTA supplied via this flag is appended to the standard 18-sequence reference set, and minimal metadata is auto-generated from the FASTA header. This is the mechanism for adding Wells et al. sequences. The combined file `06_references/tree_references_all.faa` is used for alignment and passed to `06_build_summary.py` via `--references`. See naming convention below.

4. **Sentinel system** (`logs/stepNN.done`): Re-runs skip completed steps. Manual steps (DeepTMHMM, SignalP) are recognised as done if their output directories/files exist. Use `--redo-from N` or `--redo-step N` to force re-runs. `--redo-from` takes integers only — do not pass alphanumeric step names (0c, 3b) to `--redo-from`. Step names may be alphanumeric — the `_sentinel_name()` helper handles zero-padding integers vs. passing through alphanumeric strings unchanged.

5. **Classification scoring** (in `06_build_summary.py → score_classification()`): Evidence is scored and summed. Tree score capped at +2. SoeA.hmm hit applies a conditional penalty (suppressed/−1/−3 depending on operon evidence). All other HMSS2 hits are annotation columns only. See scoring table below.

6. **NrfD deduplication**: When multiple NrfD neighbours exist, the best one is selected by: PF14589 presence first (high-specificity), then topology classification confidence.

7. **Tree clade keyword matching**: `score_classification()` scores tree placements by substring matching against `tree_clade`. Wells clade names in the data are: `PsrAPhsASrrA` (+2), `bSreASoeA` (−1), `ArrAArxA` (−2), `TtrASrdA` (−1). Outgroup clades (`FdhH_outgroup`, `NapA`, `NarG`, etc.) match no keyword and receive 0. The `elif` chain ensures only the first matching keyword applies, so `PsrAPhsASrrA` correctly scores +2 (PsrA match) not +1 (PhsA match).

---

## Current State

### Implemented and working
- [x] Mode C dual-gate discovery (Step 0c) — PF00384 ∪ PsrAPhsASreA.hmm, `discovery_source` column in output
- [x] Mode B full PF00384 scan (254 bins, 419,785 proteins → 410+ Mo-bisPGD candidates)
- [x] Protein→bin index building
- [x] Genomic neighbourhood extraction (GFF-based, contig-boundary aware)
- [x] HMMER searches: PF00384, PF03916, PF14589, PF12800, PF13247
- [x] DeepTMHMM integration (requires `cd` to install dir; output dir must not pre-exist)
- [x] SignalP 6 integration (organism=`other`, mode=`slow-sequential`)
- [x] Reference sequence download (18 curated sequences, all accessions verified in live run)
- [x] **Custom tree references** (`--custom-tree-refs`) — appended to standard references, auto-metadata from FASTA header; Wells et al. 2023 sequences integrated
- [x] MAFFT alignment (`--anysymbol` for selenocysteine in FdhH)
- [x] TrimAl trimming
- [x] IQ-TREE (fast and full modes, separate prefixes/sentinels)
- [x] Classification summary with `discovery_source` column, HTML output (HMSS2_only rows in orange)
- [x] Step 3b: HMSS2 annotation HMM searches (`--hmss2` flag, `02b_parse_hmss2.py`)
- [x] Conditional SoeA.hmm scoring (suppressed/−1/−3 based on PsrC operon evidence)
- [x] `parse_treefile_clades()` updated for custom references — longest-label-first matching, reads clade from metadata
- [x] Sentinel system for alphanumeric step names; `step_clear` handles `0c` and `3b`
- [x] NrfD label clarified: `NrfD_unclassified` renamed `NrfD_PF03916(topology_ND)`
- [x] Reference accessions verified in live run: SreA Q8NKK1, ArrA_Chrysiogenes Q5Y818

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

### Classification scoring

| Evidence | Score |
|----------|-------|
| Mo-bisPGD present (PF00384) | +2 |
| Mo-bisPGD absent | −2 |
| TAT signal peptide | +2 |
| PsrC topology (8TM) | +3 |
| NrfD PF14589 hit (no resolved topology) | +2 |
| NrfD PF03916 only (no resolved topology) | +1 |
| TtrC or SoeC topology | −1 |
| PsrB in neighbourhood | +1 |
| Tree clade: PsrA / PsrAPhsASrrA | +2 |
| Tree clade: PhsA | +1 |
| Tree clade: bSreASoeA / SoeA / TtrASrdA / TtrA | −1 |
| Tree clade: ArrAArxA / ArrA / ArxA | −2 |
| **SoeA.hmm hit, PsrC/PF14589 present** | 0 (suppressed) |
| **SoeA.hmm hit, NrfD+TAT present** | −1 (reduced) |
| **SoeA.hmm hit, no operon evidence** | −3 (full penalty) |

Classification thresholds: TRUE_PsrA ≥ 8 | LIKELY_PsrA ≥ 5 | PsrA_or_PhsA ≥ 2 + TAT + Mo | LIKELY_SoeA_or_divergent: no TAT + no NrfD + no PsrC topology | NOT_MoBisPGD_enzyme: no Mo-bisPGD | AMBIGUOUS: all other

### HMSS2 profiles used

**Step 0c (discovery gate — Mode C only)**

| Profile | Run against | Role |
|---------|-------------|------|
| PsrAPhsASreA.hmm | all_proteins_combined.faa | Second discovery gate; union with PF00384 hits |

**Step 3b (annotation + one scored profile)**

| Profile | Run against | Role |
|---------|-------------|------|
| PsrAPhsASreA.hmm | candidates.faa | Annotation column |
| **SoeA.hmm** | candidates.faa | **SCORED** — conditional penalty (see scoring table) |
| TtrA.hmm | candidates.faa | Annotation column |
| PsrBPhsBSreB.hmm | neighbours.faa | Annotation column |
| PsrCPhsCSreC.hmm | neighbours.faa | Annotation column (cross-validation; Pfam outperforms) |
| SoeB.hmm | neighbours.faa | Annotation column |
| SoeC.hmm | neighbours.faa | Annotation column |
| TtrB.hmm | neighbours.faa | Annotation column |
| TtrC.hmm | neighbours.faa | Annotation column |

### Reference sequences

18-sequence curated set downloaded automatically by `05_fetch_references.py`. All accessions verified in live run.

| Label | Accession | Clade |
|-------|-----------|-------|
| PsrA_Wolinella_succinogenes | P31075 (NOT P31077 = PsrC!) | PsrA |
| PsrA_Thermus_thermophilus | Q72LA4 | PsrA |
| PhsA_Salmonella_enterica_LT2 | P37600 | PhsA |
| TtrA_Salmonella_enterica_LT2 | Q9Z4S6 | TtrA |
| TtrA_Shewanella_ANA3 | WP_011715816.1 | TtrA |
| SreA_Acidianus_ambivalens | Q8NKK1 | SreA |
| SoeA_Allochromatium_vinosum | D3RNN8 | SoeA |
| ArrA_Shewanella_ANA3 | Q7WTU0 | ArrA |
| ArrA_Chrysiogenes_arsenatis | Q5Y818 | ArrA |
| AioA_Alcaligenes_faecalis | Q7SIF4 | AioA |
| NapA_Shewanella_oneidensis | Q8EIJ1 | NapA |
| TorA_Ecoli | P33225 | TorA |
| DmsA_Ecoli | P18775 | DmsA |
| NarG_Ecoli | P09152 | NarG |
| SerA_Thauera_selenatis | Q9S1H0 | SerA |
| PcrA_Dechloromonas_aromatica | Q47CW6 | PcrA |
| FdhG_Ecoli (outgroup) | P24183 | FdhG_outgroup |
| FdhH_Ecoli (outgroup) | P07658 | FdhH_outgroup — contains selenocysteine U |

### Custom tree references (Wells et al.)

Supplied via `--custom-tree-refs data/wells_selected_references.faa`. Wells et al. 2023 (*Microbiology Spectrum*, doi:10.1128/spectrum.04145-22; data at doi:10.5061/dryad.18931zd29) sequences are curated once and kept as a static file alongside the pipeline. The following clade names are in active use:

| Clade name in header | Tree score | Notes |
|---------------------|------------|-------|
| `PsrAPhsASrrA` | +2 | Contains `PsrA` keyword |
| `bSreASoeA` | −1 | Contains `Soe` keyword |
| `ArrAArxA` | −2 | Contains `Arr` keyword |
| `TtrASrdA` | −1 | Contains `Ttr` keyword |

Header naming convention: `>CladeName_Descriptive_organism_name__accession_or_id`

The auto-metadata awk script in Step 6 extracts the clade name as the first `_`-delimited field of the label (everything before the first underscore). This is used as both `protein` and `clade` in `reference_metadata_with_custom.tsv`, which `parse_treefile_clades()` reads via `meta.get("clade")`.

---

## Next Steps

1. **Resolve Unknown_10TM topology** for the ArrAArxA-placed candidates that have TAT + NrfD + PsrB (NODE_223 and 8 similar). These are the former best TRUE_PsrA candidates. If DeepTMHMM resolves them to PsrC (8TM), they will reach LIKELY_PsrA even with ArrA tree placement — the tree score (−2) is insufficient to outweigh Mo+2/TAT+2/PsrC+3/PsrB+1 = 8. Currently classified `PsrA_or_PhsA`.

2. **Consider raising ArrAArxA penalty to −3 or adding a hard rule**: currently ArrA placement yields only −2, so TAT + NrfD + PsrB candidates can still reach `PsrA_or_PhsA`. These are genuinely ambiguous, but if the environment has active arsenate reduction (which merits checking), these could be true arsenate reductases misplaced near ArrA. The `PsrA_or_PhsA` label is a reasonable flag for manual review.

3. **Install ete3** (`conda install -c etetoolkit ete3`) for automatic tree clade assignments in Step 9.

4. **Build a custom PsrB HMM** from P31076 (*W. succinogenes* PsrB) + BLAST hits — would outperform generic PF12800/PF13247 for PsrB detection, currently the weakest link.

---

## Critical Context & Edge Cases

### ID format (do not break this)
- Protein IDs are **contig-based**: `NODE_3193_length_16941_cov_1.235637_2`
- The **last field only** (`_2`) is the gene number
- Bin name is the **directory** containing the FAA/GFF; never appears in the protein ID
- Always use the protein→bin index for bin lookup

### Custom tree reference naming convention
- Format: `>CladeName_Descriptive_organism__accession_or_id`
- The clade name must be the **first underscore-delimited field** — this is what the awk metadata extractor and the scoring keyword matching both depend on
- Clade name should contain one of the scoring keywords: `PsrA`, `PhsA`, `Soe`, `Ttr`, `Arr`, `Arx` — or be an outgroup with no keyword (receives 0 tree score)
- Do not use `PsrA` as a prefix for non-PsrA clades — `PsrAPhsASrrA` is fine because the scoring correctly applies +2 (first match wins in the elif chain)
- `--custom-tree-refs` is reconstructed on every run after Step 6; re-running Step 7 is enough to incorporate changes to the FASTA file

### Sentinel system
- `_sentinel_name()` zero-pads plain integers, passes alphanumeric through unchanged
- `step_clear` handles `0c` (cleared from step ≤0) and `3b` (cleared from step ≤3)
- `--redo-from` accepts integers only — do NOT pass alphanumeric step names (0c, 3b) directly

### SoeA.hmm scoring design
- The MopB superfamily has ~15% cross-family sequence identity (Wells et al. 2023) — SoeA.hmm cross-hits in PsrA candidates are expected, particularly from bSreA/SoeA which is phylogenetically proximal
- Penalty is **suppressed** when `tm_canonical == PsrC` OR `NrfD_has_PF14589 == YES`
- Penalty is **−1** when NrfD present AND TAT present (moderate operon evidence)
- Penalty is **−3** otherwise
- Evidence string records which branch: `SoeA.hmm(+,suppressed_by_PsrC)`, `SoeA.hmm(+)` with −1, or `SoeA.hmm(+)` with −3
- A complete PsrABC operon (Mo+2, TAT+2, PsrC+3, PsrB+1, PsrA_tree+2 = 10) is unaffected by SoeA.hmm

### Tree clade keyword matching
- `score_classification()` uses substring matching: `"PsrA" in tree_clade`, `"Soe" in tree_clade`, etc.
- The `elif` chain means only the **first matching keyword** applies
- `PsrAPhsASrrA` matches `PsrA` first → +2 (not +1 for PhsA)
- `bSreASoeA` matches `Soe` → −1; does NOT match `PsrA` (correct)
- Outgroup clades (FdhH_outgroup, NapA, NarG) match nothing → 0 (correct)
- `parse_treefile_clades()` sorts reference labels longest-first to prevent short labels (e.g. `PsrA_Wolinella`) from matching leaves that should match longer labels (e.g. `PsrAPhsASrrA_Wells_01`)

### Pfam HMM accessions
- `hmmfetch Pfam-A.hmm PF00384` **fails** — must use versioned form `PF00384.28`
- The `pfam_fetch()` bash function handles this via `grep -m1 "^ACC\s*${acc}\."`

### NrfD evidence label meanings
- `NrfD_PF14589(+)` — high-specificity PsrC profile hit; strong PsrC evidence
- `NrfD_PF03916(topology_ND)` — broad HMM hit only; DeepTMHMM not yet run or result ambiguous
- When topology is resolved: label from DeepTMHMM (e.g. `PsrC`, `TtrC`, `PhsC`)

### DeepTMHMM quirks
- Must `cd` to the install directory before running — model files loaded relative to CWD
- Output directory must **not exist** before running — remove it before re-running
- Uses absolute paths for both `--fasta` and `--output-dir`

### SignalP 6 quirks
- Correct flag: `--organism other` (not `gram-`, not `gram+`)
- Correct flag: `--mode slow-sequential` (not `slow`)
- Output file: `prediction_results.txt`
- ID in output = full FASTA header string; parse only the first whitespace-delimited token

### MAFFT requires `--anysymbol`
- FdhH (P07658) contains selenocysteine (U) — MAFFT errors without this flag
- Any Wells sequences from FdhG/FdhH-adjacent families may also contain selenocysteine

### IQ-TREE
- Binary name on this system: `iqtree` (not `iqtree2`); pipeline tries `iqtree`, `iqtree3`, `iqtree2`
- Support string format with `-bb 1000 -alrt 1000`: node labels like `95.6/100`; ete3 must use `format=1`

### PsrB evidence design
- PsrB is **supporting evidence only** — absence never penalises classification
- `operon_completeness.tsv` reports `ABC_complete`, `AC_only`, `AB_only`, `A_only`

### Output files for manual review
- `09_summary/classification_table.html` — colour-coded; HMSS2_only rows in orange; HMSS2 columns in grey
- `09_summary/classification_table.tsv` — includes `discovery_source` and `tree_clade` columns
- `00_scan/discovery_source.tsv` — per-candidate gate origin
- `06_references/tree_references_all.faa` — exact reference set used in alignment
- `06_references/reference_metadata_with_custom.tsv` — metadata including auto-generated Wells entries
- `08_tree/psr_phylogeny.treefile` — load in iTOL or FigTree
- `03_hmmer/operon_completeness.tsv` — operon structure per candidate
- `03_hmmer/nrfd_hits.tsv` — includes `PF14589_high_specificity` column
- `03_hmmer/hmss2/hmss2_operon.tsv` — HMSS2 annotation per candidate
