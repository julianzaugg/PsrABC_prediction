# CLAUDE.md — PsrABC Identification Pipeline

## Project Overview

A bioinformatics pipeline to identify and classify **PsrABC-like Mo-bisPGD reductase operons** from Pyrodigal-annotated metagenomic genome bins. The pipeline now uses a **tree-gated classification framework**: phylogenetic placement against curated and Wells et al. DMSOR/MopB references defines the primary enzyme-family group, while local operon evidence refines confidence within compatible groups.

The pipeline distinguishes PsrA/PhsA/SrrA candidates from related Mo-bisPGD enzymes that can share A-B-C operon architecture:

| Target/group | Tree group | Key distinguishing features |
|--------------|------------|-----------------------------|
| **PsrA / PhsA / SrrA** ← main target | `PsrAPhsASrrA`, `PsrA`, `PhsA`, `SrrA` | Requires compatible tree placement; TAT + B-like neighbour + PsrC/NrfD-like membrane subunit increase confidence |
| ArrA / ArxA | `ArrAArxA`, `ArrA`, `ArxA` | ArrABC-like operon architecture; classified separately as `LIKELY_ArrA` |
| TtrA / SrdA | `TtrASrdA`, `TtrA`, `SrdA` | Ttr/Srd-like tree placement; classified as `LIKELY_TtrA_or_SrdA` |
| bSreA / SoeA / aSreA | `bSreASoeA`, `SoeA`, `SreA`, `aSreA` | Soe/Sre-like tree placement; classified as `LIKELY_SoeA_or_divergent` |
| Other DMSOR/MopB families | e.g. `FdhG`, `NapA`, `ActB`, `DmsA`, `NarG`, `NasCNarB` | Classified as `LIKELY_nonPsr_MopB`; operon-like cases are flagged as `LIKELY_nonPsr_MopB_operon_like` |

The pipeline operates in three discovery modes:

- **Mode A**: supplied list of pre-annotated psrA protein IDs
- **Mode B**: full PF00384 scan across all bins (catches annotation gaps)
- **Mode C** (recommended): dual-gate — PF00384 **UNION** PsrAPhsASreA.hmm (HMSS2). Catches divergent PsrA-like sequences that fall below PF00384's detection threshold. Requires `--hmss2`.

In Mode C, each candidate is tagged with a normalised `discovery_source` / discovery gate (`PF00384_gate` | `HMSS2_gate` | `both_gates`) which flows to the final table and HTML output.

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

1. **Protein→bin index** (`00_scan/protein_to_bin_index.tsv`): Pyrodigal names proteins after contigs (e.g. `NODE_3193_length_16941_cov_1.235637_2`), NOT after bins. The bin name is only the containing directory. The index maps every `protein_id` → (`bin_name`, `faa_path`, `gff_path`) and is built once at startup by scanning all FAA files.

2. **Dual-gate discovery (Mode C)**: Step 0c runs `PsrAPhsASreA.hmm` (from HMSS2) across all proteins in parallel with the PF00384 scan. The final candidate set is the union of both. `discovery_source.tsv` records the origin of each candidate. Final tables normalise these labels to `PF00384_gate`, `HMSS2_gate`, and `both_gates`.

3. **Custom tree references (`--custom-tree-refs`)**: Any FASTA supplied via this flag is appended to the standard curated reference set, and minimal metadata is auto-generated from the FASTA header. This is the mechanism for adding Wells et al. DMSOR/MopB references. The combined file `06_references/tree_references_all.faa` is used for alignment and passed to `06_build_summary.py` via `--references`. See naming convention below.

4. **Sentinel system** (`logs/stepNN.done`): Re-runs skip completed steps. Manual steps (DeepTMHMM, SignalP) are recognised as done if their output directories/files exist. Use `--redo-from N` or `--redo-step N` to force re-runs. `--redo-from` takes integers only — do not pass alphanumeric step names (0c, 3b) to `--redo-from`. Step names may be alphanumeric — the `_sentinel_name()` helper handles zero-padding integers vs. passing through alphanumeric strings unchanged.

5. **Tree-gated classification** (in `06_build_summary.py → score_classification()`): Tree placement defines the primary functional group. PsrA/PhsA/SrrA classifications require compatible placement in the Psr/Phs/Srr tree group. Local evidence (TAT signal, PsrB-like neighbour, PsrC/NrfD-like neighbour, PF14589, topology) refines confidence within compatible groups.

6. **NrfD deduplication**: When multiple NrfD neighbours exist, the best one is selected by: PF14589 presence first (high-specificity), then topology classification confidence.

7. **Nearest and second-nearest tree diagnostics**: `parse_treefile_clades()` assigns the nearest reference family by patristic distance and also reports the second-nearest reference family, distance margin, distance ratio, and assignment confidence. These are diagnostic outputs and do not currently alter classification.

8. **Non-Psr operon-like cases**: Candidates that place in known non-Psr MopB/DMSOR families but retain A+B+C-like local context are classified as `LIKELY_nonPsr_MopB_operon_like` rather than being forced into Psr categories.

---

## Current State

### Implemented and working
- [x] Mode C dual-gate discovery (Step 0c) — PF00384 ∪ PsrAPhsASreA.hmm, normalised discovery-gate labels in output
- [x] Mode B full PF00384 scan
- [x] Protein→bin index building
- [x] Genomic neighbourhood extraction (GFF-based, contig-boundary aware)
- [x] HMMER searches: PF00384, PF03916, PF14589, PF12800, PF13247
- [x] DeepTMHMM integration (requires `cd` to install dir; output dir must not pre-exist)
- [x] SignalP 6 integration (organism=`other`, mode=`slow-sequential`); `TATLIPO` categorical calls are retained as TAT-positive but zero/unavailable probability values are treated as missing
- [x] Reference sequence download (curated UniProt/NCBI references)
- [x] **Custom tree references** (`--custom-tree-refs`) — appended to standard references, auto-metadata from FASTA header; Wells et al. DMSOR/MopB references integrated
- [x] MAFFT alignment (`--anysymbol --auto`)
- [x] TrimAl trimming (`-gt 0.8 -st 0.001 -cons 60`)
- [x] IQ-TREE (fast and full modes, separate prefixes/sentinels)
- [x] Tree-gated classification summary with tree group, nearest/second-nearest reference diagnostics, tree confidence, and non-Psr operon-like flags
- [x] Step 3b: HMSS2 annotation HMM searches (`--hmss2` flag, `02b_parse_hmss2.py`)
- [x] Conditional SoeA.hmm evidence handling (suppressed/reduced/full penalty behaviour retained as evidence, but high-level family assignment is tree-gated)
- [x] `parse_treefile_clades()` updated for custom references — longest-label-first matching, reads clade from metadata, reports nearest and second-nearest reference families
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
| MAFFT | bioconda, `--anysymbol --auto` | Multiple sequence alignment |
| TrimAl | bioconda, `-gt 0.8 -st 0.001 -cons 60` | Alignment trimming |
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

### Classification framework (tree-gated)

Classification is now tree-gated rather than purely score-threshold based. The tree assignment defines the broad enzyme-family group, and local operon evidence refines confidence within that group.

#### Tree groups used by `06_build_summary.py`

| Tree group | Clade names recognised | Classification effect |
|------------|------------------------|-----------------------|
| `psr_phs_srr` | `PsrAPhsASrrA`, `PsrA`, `PhsA`, `SrrA` | Only this group can produce `TRUE_PsrA`, `LIKELY_PsrA`, or `PsrA_or_PhsA` |
| `arr_arx` | `ArrAArxA`, `ArrA`, `ArxA` | `LIKELY_ArrA` |
| `ttr_srd` | `TtrASrdA`, `TtrA`, `SrdA` | `LIKELY_TtrA_or_SrdA` |
| `soe_sre` | `bSreASoeA`, `aSreA`, `SoeA`, `SreA` | `LIKELY_SoeA_or_divergent` |
| `known_nonpsr_mopb` | `ActB`, `AH`, `AioAIdrA`, `AspDMSOR`, `BisC`, `DmsA`, `DorATorA`, `FdhG`, `FdhsFsdA`, `FhcB`, `FwdBFmdB`, `NapA`, `NarG`, `NasCNarB`, `Nqo3`, `RhLPgtL`, etc. | `LIKELY_nonPsr_MopB` or `LIKELY_nonPsr_MopB_operon_like` |
| unresolved/other | missing or unrecognised clade | conservative fallback to `AMBIGUOUS` unless evidence supports another non-Psr class |

#### Local Psr-like evidence retained for scoring/interpretation

| Evidence | Meaning |
|----------|---------|
| Mo-bisPGD present (PF00384) | Catalytic A-subunit domain confirmation; required for positive enzyme classification |
| TAT signal peptide | Supports periplasmic Psr/Phs/Ttr-like enzymes; `TATLIPO` categorical calls are TAT-positive but probability is treated as missing if unavailable |
| PsrC topology (8TM) | Strong PsrC-like membrane-subunit evidence |
| NrfD PF14589 hit | Strong polysulfide-reductase membrane-anchor evidence when topology is unresolved |
| NrfD PF03916 only | Broad NrfD-family membrane-anchor evidence |
| PsrB in neighbourhood | Fe-S electron-transfer subunit evidence |
| SoeA.hmm hit | Cross-reactive HMSS2 evidence; suppressed/reduced when PsrC/PF14589 or TAT+NrfD evidence contradicts SoeA identity |

#### Final classification labels

| Label | Meaning |
|-------|---------|
| `TRUE_PsrA` | Psr/Phs/Srr-compatible tree placement plus strongest PsrABC evidence |
| `LIKELY_PsrA` | Psr/Phs/Srr-compatible tree placement plus moderate/strong PsrABC evidence |
| `PsrA_or_PhsA` | Psr/Phs/Srr-compatible placement but subunit evidence insufficient to distinguish PsrA from PhsA/SrrA confidently |
| `LIKELY_ArrA` | ArrA/ArxA tree placement |
| `LIKELY_TtrA_or_SrdA` | TtrA/SrdA tree placement |
| `LIKELY_SoeA_or_divergent` | bSreA/SoeA/aSreA tree placement or Soe-like/divergent evidence |
| `LIKELY_nonPsr_MopB` | Known non-Psr DMSOR/MopB tree placement |
| `LIKELY_nonPsr_MopB_operon_like` | Known non-Psr tree placement plus A+B+C-like local operon architecture |
| `NOT_MoBisPGD_enzyme` | No PF00384 Mo-bisPGD confirmation |
| `AMBIGUOUS` | Conflicting or insufficient evidence |

#### Diagnostic output columns

`classification_table.tsv` now includes:

- `tree_group`
- `tree_distance`
- `tree_second_clade`
- `tree_second_distance`
- `tree_distance_margin`
- `tree_distance_ratio`
- `tree_assignment_confidence`
- `PsrABC_like_operon`
- `tree_incompatible_with_Psr`
- `ABC_like_nonPsr_operon`
- `PsrC_like_context`

The nearest/second-nearest tree diagnostics are currently reported for review and are not used as hard classification thresholds.

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

### Custom tree references (Wells et al. and other curated references)

Custom references are supplied with:

```bash
--custom-tree-refs data/selected_wells_major_refs_25.faa
```

Headers must follow the general style:

```text
>CladeName_Descriptive_label__accession_or_original_id
```

Examples:

```text
>PsrAPhsASrrA_Wells_001__GCA_000724175.1_42_1_3082_4
>TtrASrdA_Wells_001__GCA_000961845.1_59_60550_63631_4
>bSreASoeA_Wells_001__GCA_000981645.1_4_16765_19966_4
>FdhG_Wells_001__GCA_001303415.1_35_14725_17812_1
```

The first underscore-delimited field is used as the reference `clade` in `reference_metadata_with_custom.tsv`. `06_build_summary.py` uses this clade value for nearest-reference tree grouping.

For Wells et al. references, the current helper script selects approximately 25 sequences per major DMSOR/MopB family from the Wells/Dryad FASTA files. The expanded reference set is intended to anchor both Psr/Phs/Srr and non-Psr MopB families so that A+B+C-like operon architecture does not by itself force a PsrA assignment.

## Next Steps

1. **Inspect high-conflict cases**: prioritise candidates with `ABC_like_nonPsr_operon = YES`, `PsrC_like_context = YES`, or `tree_incompatible_with_Psr = YES`. These are biologically interesting non-Psr tree placements with PsrABC-like local architecture.

2. **Use full IQ-TREE mode for final manuscript classifications** rather than `--fast-tree`, and retain the expanded Wells reference set.

3. **Review low-confidence tree assignments** using `tree_assignment_confidence`, `tree_distance_margin`, and `tree_second_clade`. Low-confidence calls should be inspected manually, especially when the nearest and second-nearest clades imply different biology.

4. **Install ete3** (`conda install -c etetoolkit ete3`) for automatic tree clade assignments in Step 9.

## Critical Context & Edge Cases

### ID format (do not break this)
- Protein IDs are **contig-based**: `NODE_3193_length_16941_cov_1.235637_2`
- The **last field only** (`_2`) is the gene number
- Bin name is the **directory** containing the FAA/GFF; never appears in the protein ID
- Always use the protein→bin index for bin lookup

### Custom tree reference naming convention

Use `>CladeName_Descriptive_label__accession_or_original_id`. The first underscore-delimited field becomes the clade name used for tree grouping. Do not use spaces in FASTA headers.

Good examples:

```text
>PsrAPhsASrrA_Wells_001__GCA_000724175.1_42_1_3082_4
>NapA_Wells_001__GCA_...
>FdhG_Wells_001__GCA_...
```

Avoid generic names such as `custom_reference_001`, because they cannot anchor clade assignment.

### Sentinel system
- `_sentinel_name()` zero-pads plain integers, passes alphanumeric through unchanged
- `step_clear` handles `0c` (cleared from step ≤0) and `3b` (cleared from step ≤3)
- `--redo-from` accepts integers only — do NOT pass alphanumeric step names (0c, 3b) directly

### SoeA.hmm evidence design

`SoeA.hmm` is retained as useful annotation/cross-validation evidence, but high-level classification is now tree-gated. A SoeA.hmm hit does not override a Psr-compatible tree placement with strong PsrC/PF14589 evidence, and it does not rescue a candidate from a known non-Psr tree group.

The evidence string records whether the SoeA.hmm signal was suppressed or reduced by contradictory operon evidence.

### Tree grouping and diagnostics

`parse_treefile_clades()` assigns candidates to the nearest reference family by patristic distance and records second-nearest family diagnostics. Classification uses explicit tree groups rather than loose substring scoring. Psr classifications require `tree_group = psr_phs_srr`.

Diagnostic columns:

- `tree_distance`
- `tree_second_clade`
- `tree_second_distance`
- `tree_distance_margin`
- `tree_distance_ratio`
- `tree_assignment_confidence`

These diagnostics are for interpretation/manual review and do not currently change classifications.

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

SignalP 6 `TATLIPO` categorical predictions are treated as TAT-positive. If the parsed numeric TAT probability is zero or unavailable for a `TATLIPO` call, the final summary reports the probability as missing (`NA`) rather than a true zero.

### MAFFT and TrimAl settings

MAFFT is run with `--anysymbol --auto`. The `--anysymbol` flag is required to accommodate non-standard residues such as selenocysteine in some references. TrimAl uses `-gt 0.8 -st 0.001 -cons 60`, matching the Wells-style criterion of omitting columns with >20% gaps or very low similarity while retaining at least 60% of columns.

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
