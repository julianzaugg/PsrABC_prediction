# PsrABC Identification Pipeline

A bioinformatics pipeline for identifying and classifying **PsrABC (polysulfide reductase)** operons from Pyrodigal-annotated metagenomic genome bins. Distinguishes true PsrA from functional homologues (PhsA, TtrA, SoeA, SreA, bSreA) that share the Mo-bisPGD catalytic domain.

---

## Overview

Mo-bisPGD enzymes are a large and phylogenetically diverse family. PsrA cannot be reliably identified from sequence similarity or a single HMM — confident classification requires integrating multiple independent evidence types:

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

Each candidate is tagged with its `discovery_source` (both | PF00384_only | HMSS2_only). `HMSS2_only` candidates are highlighted in orange in the HTML output for manual review.

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
Step 7    MAFFT + TrimAl alignment
Step 8    IQ-TREE phylogeny
Step 9    Classification summary
```

---

## Classification Scoring

Each candidate is scored by summing evidence weights:

| Evidence | Score |
|----------|-------|
| Mo-bisPGD present (PF00384) | +2 |
| Mo-bisPGD absent | −2 |
| TAT signal peptide | +2 |
| PsrC topology (8TM) | +3 |
| NrfD PF14589 hit (no topology) | +2 |
| NrfD PF03916 hit only (no topology) | +1 |
| TtrC or SoeC topology | −1 |
| PsrB in neighbourhood | +1 |
| Tree: PsrA / PsrAPhsASrrA clade | +2 (capped) |
| Tree: PhsA clade | +1 |
| Tree: bSreASoeA / SoeA / TtrASrdA / TtrA clade | −1 |
| Tree: ArrAArxA / ArrA / ArxA clade | −2 |
| SoeA.hmm hit, PsrC/PF14589 present | 0 (suppressed) |
| SoeA.hmm hit, NrfD + TAT present | −1 (reduced) |
| SoeA.hmm hit, no operon evidence | −3 (full penalty) |

Classification thresholds:

| Label | Condition |
|-------|-----------|
| `TRUE_PsrA` | Score ≥ 8 |
| `LIKELY_PsrA` | Score ≥ 5 |
| `PsrA_or_PhsA` | Score ≥ 2, TAT present, Mo confirmed |
| `LIKELY_SoeA_or_divergent` | No TAT, no NrfD, no PsrC topology |
| `NOT_MoBisPGD_enzyme` | No Mo-bisPGD domain |
| `AMBIGUOUS` | All other cases |

Tree evidence is capped at +2 to prevent phylogeny overriding biochemical evidence. SoeA classification requires all three negative conditions simultaneously to avoid false calls from incomplete genomes. PsrB absence never penalises classification.

---

## Custom Tree References: Naming Convention

Headers must follow: `>CladeName_Descriptive_organism__accession_or_id`

The **first underscore-delimited field** becomes the clade name in the reference metadata, which is then used for tree scoring keyword matching. Examples using Wells et al. 2023 clades:

```
>PsrAPhsASrrA_Desulfovibrio_vulgaris__WP_010940123
>bSreASoeA_Aquifex_aeolicus__O67847
>ArrAArxA_Geobacter_sulfurreducens__Q74BP2
>TtrASrdA_Salmonella_enterica__Q9Z4S8
```

Scoring keyword matching (substring, elif chain — first match wins):

| Clade name contains | Tree score |
|--------------------|------------|
| `PsrA` | +2 |
| `PhsA` | +1 |
| `Soe` | −1 |
| `Ttr` | −1 |
| `Arr` or `Arx` | −2 |
| none of the above | 0 |

Outgroup labels (FdhH_outgroup, NapA, NarG, TorA, etc.) match no keyword and receive 0 score, which is correct — they anchor the tree without biasing classification.

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
| `09_summary/classification_table.html` | Colour-coded HTML; HMSS2_only rows orange; annotation columns grey |
| `09_summary/classification_table.tsv` | Full table; includes `discovery_source` and `tree_clade` columns |
| `00_scan/discovery_source.tsv` | Per-candidate gate origin |
| `06_references/tree_references_all.faa` | Exact reference FASTA used in alignment (standard + custom) |
| `06_references/reference_metadata_with_custom.tsv` | Metadata for all references including auto-generated Wells entries |
| `03_hmmer/operon_completeness.tsv` | PsrABC subunit counts per candidate |
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
