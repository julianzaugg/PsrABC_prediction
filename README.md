# PsrABC Identification Pipeline

A bioinformatics pipeline for identifying and classifying **PsrABC (polysulfide reductase)** operons from Pyrodigal-annotated metagenomic genome bins. Distinguishes true PsrA from functional homologues (PhsA, TtrA, SoeA, SreA, bSreA) that share the Mo-bisPGD catalytic domain.

---

## Overview

Mo-bisPGD enzymes are a large and phylogenetically diverse family. PsrA cannot be reliably identified from sequence similarity or a single HMM — confident classification requires integrating multiple independent evidence types and constraining interpretation by phylogenetic placement.

The pipeline now uses a **tree-gated classification framework**: candidate A subunits are placed in a reference phylogeny containing curated enzymes plus expanded Wells et al. DMSOR/MopB representatives, and only candidates placed in the PsrA/PhsA/SrrA tree group can be classified as `TRUE_PsrA`, `LIKELY_PsrA`, or `PsrA_or_PhsA`.

| Evidence | Tool | Interpretation |
|----------|------|----------------|
| Mo-bisPGD domain | PF00384 (HMMER) | Confirms candidate belongs to the broad Mo-bisPGD enzyme family |
| TAT signal peptide | SignalP 6 | Supports periplasmic Psr/Phs/Ttr-like enzymes; `TATLIPO` categorical calls are TAT-positive |
| NrfD/PsrC neighbour | PF03916/PF14589 + DeepTMHMM | Supports A+B+C operon architecture; 8TM/PF14589 strengthens PsrC-like context |
| PsrB-like neighbour | PF12800/PF13247/HMSS2 | Supports Fe-S electron-transfer subunit context |
| Phylogenetic clade | MAFFT/TrimAl/IQ-TREE/ete3 | Defines the broad functional group used for classification |

The pipeline handles fragmented metagenomic assemblies where annotation is incomplete and divergent sequences are common.

---

## Discovery Modes

### Mode A — supplied IDs
```bash
bash 00_pipeline.sh --ids my_candidates.txt --bindir /bins --pfam Pfam-A.hmm --outdir results/
```

### Mode B — PF00384 scan only
```bash
bash 00_pipeline.sh --bindir /bins --pfam Pfam-A.hmm --outdir results/
```

### Mode C — dual-gate discovery (recommended)
```bash
bash 00_pipeline.sh \
    --bindir    /path/to/bins \
    --pfam      /path/to/Pfam-A.hmm \
    --hmss2     /path/to/HMSS2/Hidden_Markov_Models/Inorganic_Sulfur_Metabolism \
    --deeptmhmm /path/to/DeepTMHMM-Academic-License-v1.0 \
    --outdir    results/psr_analysis
```

Runs two independent HMM gates across all proteins and takes their union:
1. **PF00384** — broad Mo-bisPGD domain profile
2. **PsrAPhsASreA.hmm** (HMSS2) — clade-specific profile with higher sensitivity for divergent PsrA

Each candidate is tagged with its `discovery_source`, normalised in the final table as `PF00384_gate`, `HMSS2_gate`, or `both_gates`.

### Adding custom tree reference sequences
```bash
bash 00_pipeline.sh [other args] \
    --custom-tree-refs data/wells_selected_references.faa \
    --redo-from 7
```

Any FASTA supplied via `--custom-tree-refs` is appended to the standard reference set before alignment. This is the recommended way to add Wells et al. 2023 sequences or other curated representatives. The file is incorporated on every run after Step 6 — re-running from step 7 is sufficient when updating the file.

---

## Pipeline Steps

```
Step 0a   Build protein→bin index          (always)
Step 0b   PF00384 scan                     (Mode B/C)
Step 0c   HMSS2 PsrAPhsASreA pre-screen    (Mode C only)
Step 1    Extract candidate sequences
Step 2    Genomic neighbourhood extraction (±10 genes, contig-aware)
Step 3    HMMER neighbourhood searches     (PF03916, PF14589, PF12800, PF13247)
Step 3b   HMSS2 annotation searches        (optional, requires --hmss2)
Step 4    DeepTMHMM TM topology            (manual if not installed)
Step 5    SignalP 6 TAT prediction         (manual if not installed)
Step 6    Download reference sequences + build tree_references_all.faa
Step 7    MAFFT --anysymbol --auto + TrimAl -gt 0.8 -st 0.001 -cons 60
Step 8    IQ-TREE phylogeny
Step 9    Tree-gated classification summary
```

---

## Tree-gated Classification Framework

Classification is not purely score-threshold based. The candidate's placement in the reference tree defines the broad enzyme-family group, and local operon evidence refines confidence within compatible groups.

### Tree groups

| Tree group | Clades | Classification effect |
|------------|--------|-----------------------|
| `psr_phs_srr` | `PsrAPhsASrrA`, `PsrA`, `PhsA`, `SrrA` | Required for `TRUE_PsrA`, `LIKELY_PsrA`, and `PsrA_or_PhsA` |
| `arr_arx` | `ArrAArxA`, `ArrA`, `ArxA` | `LIKELY_ArrA` |
| `ttr_srd` | `TtrASrdA`, `TtrA`, `SrdA` | `LIKELY_TtrA_or_SrdA` |
| `soe_sre` | `bSreASoeA`, `aSreA`, `SoeA`, `SreA` | `LIKELY_SoeA_or_divergent` |
| `known_nonpsr_mopb` | `ActB`, `AH`, `AioAIdrA`, `AspDMSOR`, `BisC`, `DmsA`, `DorATorA`, `FdhG`, `FdhsFsdA`, `FhcB`, `FwdBFmdB`, `NapA`, `NarG`, `NasCNarB`, `Nqo3`, `RhLPgtL`, etc. | `LIKELY_nonPsr_MopB` or `LIKELY_nonPsr_MopB_operon_like` |

### Final labels

| Label | Meaning |
|-------|---------|
| `TRUE_PsrA` | Psr/Phs/Srr-compatible tree placement plus strongest PsrABC evidence |
| `LIKELY_PsrA` | Psr/Phs/Srr-compatible tree placement plus moderate/strong PsrABC evidence |
| `PsrA_or_PhsA` | Psr/Phs/Srr-compatible tree placement but insufficient evidence to resolve PsrA vs PhsA/SrrA |
| `LIKELY_ArrA` | ArrA/ArxA tree placement |
| `LIKELY_TtrA_or_SrdA` | TtrA/SrdA tree placement |
| `LIKELY_SoeA_or_divergent` | bSreA/SoeA/aSreA tree placement or Soe-like/divergent evidence |
| `LIKELY_nonPsr_MopB` | Known non-Psr DMSOR/MopB tree placement |
| `LIKELY_nonPsr_MopB_operon_like` | Known non-Psr tree placement plus A+B+C-like operon architecture |
| `NOT_MoBisPGD_enzyme` | No PF00384 Mo-bisPGD confirmation |
| `AMBIGUOUS` | Conflicting or insufficient evidence |

### Diagnostic columns

`classification_table.tsv` includes nearest/second-nearest reference diagnostics:

- `tree_group`
- `tree_distance`
- `tree_second_clade`
- `tree_second_distance`
- `tree_distance_margin`
- `tree_distance_ratio`
- `tree_assignment_confidence`

It also reports conflict/operon-context flags:

- `PsrABC_like_operon`
- `tree_incompatible_with_Psr`
- `ABC_like_nonPsr_operon`
- `PsrC_like_context`

The tree-distance diagnostics are currently for interpretation and manual review; they do not alter classifications.

---

## Custom Tree References: Naming Convention

Headers must follow:

```text
>CladeName_Descriptive_label__accession_or_original_id
```

The **first underscore-delimited field** becomes the clade name in the reference metadata and is used for tree grouping. Examples using Wells et al. DMSOR/MopB clades:

```text
>PsrAPhsASrrA_Wells_001__GCA_000724175.1_42_1_3082_4
>bSreASoeA_Wells_001__GCA_000981645.1_4_16765_19966_4
>ArrAArxA_Wells_001__GCA_...
>TtrASrdA_Wells_001__GCA_...
>NapA_Wells_001__GCA_...
>FdhG_Wells_001__GCA_...
```

Avoid generic labels such as `custom_reference_001`, because they cannot act as clade anchors.

The Wells helper script can be used to select approximately 25 representatives per major DMSOR/MopB family and write headers in this format. These references are supplied to the pipeline with `--custom-tree-refs`.

---

## Prerequisites

```bash
conda install -c bioconda hmmer seqkit mafft trimal iqtree biopython requests
conda install -c etetoolkit ete3   # optional, for automatic tree clade assignment
pip install requests biopython
```

External tools (academic licences — obtain separately):
- **DeepTMHMM** — https://dtu.biolib.com/DeepTMHMM (Academic License v1.0)
- **SignalP 6** — https://services.healthtech.dtu.dk/services/SignalP-6.0/
- **HMSS2** — local HMM directory, path supplied via `--hmss2`
- **Pfam-A.hmm** — full pressed database, path supplied via `--pfam`

---

## Bin Directory Structure

```
BIN_DIR/
  BinName_1/
    BinName_1.faa
    BinName_1.gff
  ...
```

Pyrodigal names proteins after **contigs**, not bins. Protein IDs look like `NODE_3193_length_16941_cov_1.235637_2`. The bin name never appears in the protein ID — the pipeline builds a protein→bin index at Step 0a to handle this.

---

## Re-run Control

```bash
# Re-run from step 7 (e.g. after updating custom tree references)
bash 00_pipeline.sh [same args] --redo-from 7

# Re-run step 8 in fast mode only (topology check)
bash 00_pipeline.sh [same args] --redo-step 8 --fast-tree

# Re-run discovery from scratch
bash 00_pipeline.sh [same args] --redo-from 0
```

`--redo-from` accepts integers only. IQ-TREE fast and full runs coexist with separate prefixes and sentinels. The `--custom-tree-refs` FASTA is always re-incorporated after Step 6 — running `--redo-from 7` is sufficient to update the tree after editing the file.

---

## Key Output Files

| File | Description |
|------|-------------|
| `09_summary/classification_table.html` | Colour-coded HTML summary |
| `09_summary/classification_table.tsv` | Full tree-gated classification table; includes discovery gate, tree group, tree diagnostics, and operon-context flags |
| `00_scan/discovery_source.tsv` | Per-candidate gate origin |
| `06_references/tree_references_all.faa` | Exact reference FASTA used in alignment (standard + custom) |
| `06_references/reference_metadata_with_custom.tsv` | Metadata for all references including auto-generated custom/Wells entries |
| `03_hmmer/operon_completeness.tsv` | PsrABC-like subunit counts per candidate |
| `03_hmmer/nrfd_hits.tsv` | NrfD/PsrC hits with PF14589 specificity flag |
| `03_hmmer/hmss2/hmss2_operon.tsv` | HMSS2 annotation per candidate |
| `08_tree/psr_phylogeny.treefile` | IQ-TREE maximum-likelihood tree |
| `04_topology/topology_summary.tsv` | DeepTMHMM TM helix classification per NrfD candidate |
| `05_signalp/tat_summary.tsv` | SignalP 6 TAT predictions per candidate |

---

## Manual Steps

**DeepTMHMM:**
```bash
cd /path/to/DeepTMHMM-Academic-License-v1.0
python predict.py \
    --fasta /abs/path/to/results/03_hmmer/nrfd_candidates.faa \
    --output-dir /abs/path/to/results/04_topology/deeptmhmm_out
bash 00_pipeline.sh [same args] --redo-from 4
```

**SignalP 6:**
- Upload `01_sequences/candidates_psrA.faa` to https://services.healthtech.dtu.dk/services/SignalP-6.0/
- Organism: `other`; save `prediction_results.txt` to `05_signalp/`
- `bash 00_pipeline.sh [same args] --redo-from 5`
