# PsrABC Identification Pipeline

A bioinformatics pipeline for identifying and classifying **PsrABC (polysulfide reductase)** operons from Pyrodigal-annotated metagenomic genome bins. Distinguishes true PsrA from functional homologues (PhsA, TtrA, SoeA, SreA) that share the Mo-bisPGD catalytic domain.

---

## Overview

Mo-bisPGD enzymes are a large and phylogenetically diverse family. PsrA alone cannot be reliably identified from sequence similarity or a single HMM — confident classification requires integrating multiple independent evidence types:

| Evidence | Tool | PsrA | SoeA | TtrA | PhsA |
|----------|------|------|------|------|------|
| Mo-bisPGD domain | PF00384 (HMMER) | ✓ | ✓ | ✓ | ✓ |
| TAT signal peptide | SignalP 6 | ✓ | **✗** | ✓ | ✓ |
| NrfD/PsrC neighbour | PF03916/PF14589 | ✓ 8TM | – | ✓ 9TM | ✓ 5TM+haem |
| Phylogenetic clade | IQ-TREE | ✓ | – | – | – |

The pipeline handles fragmented metagenomic assemblies where annotation is incomplete and divergent sequences are common.

---

## Discovery Modes

### Mode A — supplied IDs
```bash
bash 00_pipeline.sh --ids my_candidates.txt --bindir /bins --pfam Pfam-A.hmm --outdir results/
```
Use when you already have a curated list of candidate protein IDs (e.g. from a previous annotation run).

### Mode B — PF00384 scan
```bash
bash 00_pipeline.sh --bindir /bins --pfam Pfam-A.hmm --outdir results/
```
Scans all proteins in all bins with PF00384 (Mo-bisPGD). Catches unannotated candidates but may miss highly divergent PsrA sequences where PF00384 sensitivity is insufficient.

### Mode C — dual-gate discovery (recommended)
```bash
bash 00_pipeline.sh \
    --bindir    /path/to/bins \
    --pfam      /path/to/Pfam-A.hmm \
    --hmss2     /path/to/HMSS2/Hidden_Markov_Models/Inorganic_Sulfur_Metabolism \
    --deeptmhmm /path/to/DeepTMHMM-Academic-License-v1.0 \
    --outdir    results/psr_analysis
```

Runs two independent HMM gates across all proteins:
1. **PF00384** — broad Mo-bisPGD domain profile
2. **PsrAPhsASreA.hmm** (HMSS2) — clade-specific profile with higher sensitivity for divergent PsrA

The final candidate set entering downstream analysis is the **union** of both gates. Each candidate is tagged with its `discovery_source`:

| Source | Meaning |
|--------|---------|
| `both` | Hit both PF00384 and PsrAPhsASreA.hmm — highest initial confidence |
| `PF00384_only` | Standard Mo-bisPGD hit |
| `HMSS2_only` | Missed by PF00384; found only by the clade-specific HMM — **review carefully** |

`HMSS2_only` candidates are highlighted with an orange border in the HTML output. They enter the same scoring pipeline as all other candidates — no score penalty is applied — but their source is always visible for manual review.

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
Step 6    Download reference sequences
Step 7    MAFFT + TrimAl alignment
Step 8    IQ-TREE phylogeny
Step 9    Classification summary
```

Steps are sentinel-controlled — re-runs skip completed steps. Manual steps (DeepTMHMM, SignalP 6) are recognised as done if their output files exist.

---

## Classification Scoring

Each candidate is scored by summing evidence weights:

| Evidence | Score |
|----------|-------|
| Mo-bisPGD present (PF00384) | +2 |
| Mo-bisPGD absent | −2 |
| TAT signal peptide | +2 |
| PsrC topology (8TM, PF14589) | +3 (+1 extra for PF14589 vs PF03916) |
| TtrC or SoeC topology | −1 |
| PsrB in neighbourhood | +1 |
| Tree clade: PsrA | +2 (capped) |
| Tree clade: PhsA | +1 |
| Tree clade: SoeA/TtrA | −1 |
| **SoeA.hmm hit** (HMSS2) | **−3** |

Classification thresholds:

| Label | Threshold |
|-------|-----------|
| `TRUE_PsrA` | ≥ 8 |
| `LIKELY_PsrA` | ≥ 5 |
| `PsrA_or_PhsA` | ≥ 2, TAT present, Mo confirmed |
| `LIKELY_SoeA_or_divergent` | No TAT, no NrfD, no PsrC topology |
| `AMBIGUOUS` | All other cases |

**Tree evidence is capped at +2** to prevent phylogeny overriding biochemical evidence. **SoeA.hmm** is the only HMSS2 profile that contributes to scoring; all other HMSS2 results appear as annotation columns only.

SoeA classification requires **all three** negative conditions (no TAT, no NrfD neighbour, no PsrC topology) to avoid misclassifying incomplete-genome PsrA as SoeA.

PsrB absence **never penalises** classification — 4Fe-4S cluster proteins are routinely missed by HMM-based annotation in metagenomic data.

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
  BinName_2/
    BinName_2.faa
    BinName_2.gff
  ...
```

Pyrodigal names proteins after **contigs**, not bins. Protein IDs look like `NODE_3193_length_16941_cov_1.235637_2`. The bin name never appears in the protein ID — the pipeline builds a protein→bin index at Step 0a to handle this.

---

## Re-run Control

```bash
# Re-run from step 7 (e.g. after adding reference sequences)
bash 00_pipeline.sh [same args] --redo-from 7

# Re-run step 8 in fast mode only
bash 00_pipeline.sh [same args] --redo-step 8 --fast-tree

# Re-run discovery from scratch (e.g. after updating HMSS2 HMMs)
bash 00_pipeline.sh [same args] --redo-from 0
```

IQ-TREE fast (`--fast-tree`) and full runs coexist with separate prefixes (`psr_phylogeny_fast.*` vs `psr_phylogeny.*`) and separate sentinels (`step8f.done` vs `step08.done`).

---

## Key Output Files

| File | Description |
|------|-------------|
| `09_summary/classification_table.html` | Colour-coded HTML; HMSS2_only rows highlighted orange; annotation columns in grey |
| `09_summary/classification_table.tsv` | Full table for R/Python; includes `discovery_source` column |
| `00_scan/discovery_source.tsv` | Per-candidate gate origin (PF00384_only / HMSS2_only / both) |
| `03_hmmer/operon_completeness.tsv` | PsrABC subunit counts per candidate |
| `03_hmmer/nrfd_hits.tsv` | NrfD/PsrC hits with PF14589 specificity flag |
| `03_hmmer/hmss2/hmss2_operon.tsv` | HMSS2 annotation per candidate |
| `08_tree/psr_phylogeny.treefile` | IQ-TREE maximum-likelihood tree |
| `04_topology/topology_summary.tsv` | DeepTMHMM TM helix classification per NrfD candidate |
| `05_signalp/tat_summary.tsv` | SignalP 6 TAT predictions per candidate |

---

## Manual Steps

If DeepTMHMM or SignalP 6 are not installed locally, the pipeline will pause and print instructions:

**DeepTMHMM:**
```bash
cd /path/to/DeepTMHMM-Academic-License-v1.0
python predict.py \
    --fasta /abs/path/to/results/03_hmmer/nrfd_candidates.faa \
    --output-dir /abs/path/to/results/04_topology/deeptmhmm_out
# Then re-run:
bash 00_pipeline.sh [same args] --redo-from 4
```

**SignalP 6:**
- Upload `01_sequences/candidates_psrA.faa` to https://services.healthtech.dtu.dk/services/SignalP-6.0/
- Organism: `other`
- Save `prediction_results.txt` to `05_signalp/`
- Then: `bash 00_pipeline.sh [same args] --redo-from 5`

---

## Notes on HMSS2 Integration

HMSS2 (Hidden Markov models for Sulfur/Sulfate Metabolism) provides purpose-built profiles for sulfur metabolism enzymes. The pipeline uses HMSS2 in two distinct ways:

**Step 0c (discovery):** `PsrAPhsASreA.hmm` is run as an independent discovery gate alongside PF00384. This recovers divergent PsrA sequences that have diverged enough from the consensus Mo-bisPGD fold to fall below PF00384's detection threshold, but are still recognisable by a clade-specific profile trained on PsrA/PhsA/SreA sequences. These `HMSS2_only` candidates are the primary motivation for Mode C.

**Step 3b (annotation):** The full HMSS2 profile set is run against candidates and neighbourhood proteins for cross-validation. `SoeA.hmm` is the only profile whose hits enter scoring (−3 penalty), converting positive SoeA evidence from pure absence-of-evidence reasoning to a supported classification. All other Step 3b profiles appear as annotation columns in the output table for manual inspection.
