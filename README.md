# PsrABC Identification Pipeline

A bioinformatics pipeline for identifying and classifying **PsrABC (polysulfide reductase)** operons from Pyrodigal-annotated metagenomic genome bins. Distinguishes true PsrA from functional homologues (PhsA, TtrA, SoeA, SreA, bSreA) that share the Mo-bisPGD catalytic domain.

---

## Overview

Mo-bisPGD enzymes are a large and phylogenetically diverse family. PsrA cannot be reliably identified from sequence similarity or a single HMM — confident classification requires integrating multiple independent evidence types and constraining interpretation by phylogenetic placement.

The pipeline uses a **tree-gated classification framework**: candidate A subunits are placed in a reference phylogeny containing curated enzymes plus expanded Wells et al. DMSOR/MopB representatives, and only candidates placed in the PsrA/PhsA/SrrA tree group can be classified as `TRUE_PsrA`, `LIKELY_PsrA`, or `PsrA_or_PhsA`.

| Evidence | Tool | Interpretation |
|----------|------|----------------|
| Mo-bisPGD domain | PF00384 (HMMER) | Confirms candidate belongs to the broad Mo-bisPGD enzyme family |
| TAT signal peptide | SignalP 6 (Teufel et al. 2022) | Supports periplasmic Psr/Phs/Ttr-like enzymes; `TATLIPO` categorical calls are treated as TAT-positive; zero numeric TAT probability associated with a `TATLIPO` call is masked as NA |
| NrfD/PsrC neighbour | PF03916/PF14589 + DeepTMHMM (Hallgren et al. 2022) | Supports A+B+C operon architecture; 8 TM helices (Jormakka et al. 2008) and/or PF14589 hits strengthen PsrC-like context |
| PsrB-like neighbour | PF12800/PF13247/HMSS2 | Supports Fe-S electron-transfer subunit context (Rothery et al. 2008) |
| Phylogenetic clade | MAFFT (Katoh & Standley 2013) / TrimAl (Capella-Gutiérrez et al. 2009) / IQ-TREE (Minh et al. 2020) / ete3 | Defines the broad functional group used for classification |

The pipeline handles fragmented metagenomic assemblies where annotation is incomplete and divergent sequences are common.

---

## Citations

If you use this pipeline, please cite the following tools and reference datasets.

### Tools

- **Pyrodigal** — Larralde, M. Pyrodigal: Python bindings and interface to Prodigal, an efficient method for gene prediction in prokaryotes. *Journal of Open Source Software* 7, 4296 (2022). https://doi.org/10.21105/joss.04296
- **HMMER** — Eddy, S.R. Accelerated profile HMM searches. *PLOS Computational Biology* 7, e1002195 (2011). http://hmmer.org
- **Pfam** — Mistry, J. et al. Pfam: The protein families database in 2021. *Nucleic Acids Research* 49, D412–D419 (2021). https://doi.org/10.1093/nar/gkaa913
- **MAFFT** — Katoh, K. & Standley, D.M. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Molecular Biology and Evolution* 30, 772–780 (2013). https://doi.org/10.1093/molbev/mst010
- **TrimAl** — Capella-Gutiérrez, S., Silla-Martínez, J.M. & Gabaldón, T. trimAl: a tool for automated alignment trimming in large-scale phylogenetic analysis. *Bioinformatics* 25, 1972–1973 (2009). https://doi.org/10.1093/bioinformatics/btp348
- **IQ-TREE 2** — Minh, B.Q. et al. IQ-TREE 2: New models and methods for phylogenetic inference. *Molecular Biology and Evolution* 37, 1530–1534 (2020). https://doi.org/10.1093/molbev/msaa015
- **SignalP 6.0** — Teufel, F. et al. SignalP 6.0 predicts all five types of signal peptides using protein language models. *Nature Biotechnology* 40, 1023–1025 (2022). https://doi.org/10.1038/s41587-021-01156-3
- **DeepTMHMM** — Hallgren, J. et al. DeepTMHMM predicts alpha and beta transmembrane proteins using deep neural networks. *bioRxiv* (2022). https://doi.org/10.1101/2022.04.08.487609
- **seqkit** — Shen, W. et al. SeqKit: A cross-platform and ultrafast toolkit for FASTA/Q file manipulation. *PLOS ONE* 11, e0163962 (2016). https://doi.org/10.1371/journal.pone.0163962
- **ete3** — Huerta-Cepas, J., Serra, F. & Bork, P. ETE 3: Reconstruction, analysis, and visualization of phylogenomic data. *Molecular Biology and Evolution* 33, 1635–1638 (2016). https://doi.org/10.1093/molbev/msw046

### Reference sequences and classification framework

- **Little et al. (2024)** — Little, C.J. et al. Ecological and evolutionary drivers of sulfur oxidation in the ocean. *Nature Microbiology* 9, 411–424 (2024). https://doi.org/10.1038/s41564-023-01560-2
  *(PhsA, TtrA, SoeA, ArrA, NapA, DmsA, TorA, NarG, SerA, PcrA reference sequences and clade definitions)*
- **Degré et al. (2026)** — Degré, A.G. et al. Genomic and metabolic diversity of Shewanella tetrathionate reductase operons. *Environmental Microbiology* (2026). https://doi.org/10.1111/1462-2920.70258
  *(TtrA Shewanella ANA-3 reference; WP_011715816.1)*
- **Jormakka et al. (2008)** — Jormakka, M. et al. Molecular mechanism of energy conservation in polysulfide respiration. *Nature Structural & Molecular Biology* 15, 730–737 (2008). https://doi.org/10.1038/nsmb.1434
  *(PDB 2VPZ; PsrC = 8 TM helices; PsrA TAT signal)*
- **Rothery et al. (2008)** — Rothery, R.A., Workun, G.J. & Weiner, J.H. The prokaryotic complex iron–sulfur molybdoenzyme family. *Biochimica et Biophysica Acta — Biomembranes* 1778, 1897–1929 (2008). https://doi.org/10.1016/j.bbamem.2007.09.002
  *(CISM/MopB family overview; PsrB iron–sulfur subunit divergence; TtrC = 9 TM helices)*
- **Stoffels et al. (2012)** — Stoffels, L. et al. Metagenomic and functional analysis of thiosulfate reductase during anaerobic growth of *Salmonella* Typhimurium. *PLOS Genetics* 8, e1002682 (2012). https://doi.org/10.1371/journal.pgen.1002682
  *(PhsC = 5 TM helices)*

### Expanded DMSOR/MopB reference tree

- **Wells et al. (2020)** — Wells, M. et al. Methane, arsenic, selenium and the origins of the DMSO reductase family. *Scientific Reports* 10, 10946 (2020). https://doi.org/10.1038/s41598-020-67892-9
- **Wells et al. (2023)** — Wells, M., Kim, M., Akob, D.M., Basu, P. & Stolz, J.F. Impact of the Dimethyl Sulfoxide Reductase Superfamily on the Evolution of Biogeochemical Cycles. *Microbiology Spectrum* 11, e04145-22 (2023). https://doi.org/10.1128/spectrum.04145-22
  *(Expanded DMSORcompwoutMop dataset; ~25 sequences per DMSOR/MopB family used to anchor non-Psr clades in the reference tree)*

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

Each candidate is tagged with its `discovery_source`, normalised in the final table as `PF00384_gate`, `HMSS2_gate`, or `both_gates`. Candidates entering through HMSS2 only (`HMSS2_gate`) are highlighted with an orange border in the HTML output and should be reviewed manually before acceptance.

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

Classification is not purely score-threshold based. The candidate's placement in the reference tree defines the broad enzyme-family group, and local operon evidence refines confidence within compatible groups. This design prevents non-Psr Mo-bisPGD enzymes with A-B-C operon architecture (e.g. ArrABC) from being misclassified as PsrA on the basis of operon evidence alone.

### Score components

Before the tree gate is applied, a score is computed from local evidence. The score influences confidence tiers within a given tree group but cannot promote a candidate across tree group boundaries.

| Evidence | Score |
|----------|-------|
| PF00384 Mo-bisPGD confirmed | +2 |
| TAT signal peptide | +2 |
| NrfD neighbour with PF14589 (high-specificity PsrC) | +2 |
| NrfD neighbour with PF03916 only (topology not determined) | +1 |
| PsrC topology (8 TM helices, DeepTMHMM) | +3 |
| PhsC topology (5 TM helices) | +1 |
| PsrB-like neighbour (PF12800/PF13247) | +1 |
| TtrC topology (9 TM helices) | −1 |
| SoeC topology (4–6 TM helices) | −1 |
| Incompatible tree group (Arr/Arx, Soe/Sre, Ttr/Srd, known non-Psr) | −2 |
| SoeA.hmm hit (HMSS2) with no contradicting PsrC/TAT evidence | −3 |
| SoeA.hmm hit (HMSS2) with TAT + NrfD present | −1 |
| SoeA.hmm hit (HMSS2) suppressed by PsrC or PF14589 evidence | 0 |

The SoeA.hmm score penalty is the only HMSS2 result that contributes to scoring. All other HMSS2 annotation columns are informational only.

### NrfD best-hit selection

When multiple NrfD-family proteins are found in a candidate's neighbourhood, the representative used for topology classification is selected by:

1. Hits with PF14589 (high-specificity NrfD_2) take priority
2. Hits with a classified (non-ambiguous) topology are preferred over unclassified hits
3. Hits with a PsrC classification are further preferred
4. The first hit otherwise

### Transmembrane topology and haem binding

Membrane subunit class is assigned solely from DeepTMHMM TM helix count (Jormakka et al. 2008; Rothery et al. 2008; Stoffels et al. 2012):

| TM helices | Classification |
|------------|----------------|
| 8 | PsrC |
| 9 | TtrC |
| 5 | PhsC (ambiguity with SoeC resolved by tree gate) |
| 4 or 6 | SoeC_like |
| ≥10 | Unclassified (consistent with other reductase families, e.g. ArrC) |

Note: PsrC and PhsC both coordinate b-type haems via conserved transmembrane histidine residues. b-type haem binding produces no detectable sequence motif (no covalent CXXCH attachment), so haem-motif searching is not used. Classification is based solely on TM helix count.

### TAT signal handling

SignalP 6.0 `TATLIPO` categorical predictions are treated as TAT-positive. However, where a `TATLIPO` call is associated with a numeric TAT probability of 0.0 or missing, this value is masked as `NA` in the output table to avoid misleading quantitative comparisons. The `has_TAT_signal` column is still reported as `YES` in these cases.

### Tree groups

| Tree group | Clades | Classification effect |
|------------|--------|-----------------------|
| `psr_phs_srr` | `PsrAPhsASrrA`, `PsrA`, `PhsA`, `SrrA` | Required for `TRUE_PsrA`, `LIKELY_PsrA`, and `PsrA_or_PhsA` |
| `arr_arx` | `ArrAArxA`, `ArrA`, `ArxA` | `LIKELY_ArrA`; hard override when strong operon evidence is present (ArrABC operon architecture is identical to PsrABC — TAT + NrfD/PsrB does not discriminate between the two) |
| `ttr_srd` | `TtrASrdA`, `TtrA`, `SrdA` | `LIKELY_TtrA_or_SrdA` |
| `soe_sre` | `bSreASoeA`, `aSreA`, `SoeA`, `SreA` | `LIKELY_SoeA_or_divergent` |
| `known_nonpsr_mopb` | `ActB`, `AH`, `AioAIdrA`, `AspDMSOR`, `BisC`, `DmsA`, `DorATorA`, `FdhG`, `FdhsFsdA`, `FhcB`, `FwdBFmdB`, `NapA`, `NarG`, `NasCNarB`, `Nqo3`, `RhLPgtL`, etc. | `LIKELY_nonPsr_MopB` or `LIKELY_nonPsr_MopB_operon_like` |

### Final labels

| Label | Meaning |
|-------|---------|
| `TRUE_PsrA` | Psr/Phs/Srr-compatible tree placement plus strongest PsrABC evidence (score ≥ 8, strong operon) |
| `LIKELY_PsrA` | Psr/Phs/Srr-compatible tree placement plus moderate/strong PsrABC evidence (score ≥ 5, moderate operon) |
| `PsrA_or_PhsA` | Psr/Phs/Srr-compatible or unresolved tree placement, moderate operon evidence but insufficient to resolve PsrA vs PhsA/SrrA |
| `LIKELY_ArrA` | ArrA/ArxA tree placement; hard override when TAT + NrfD or PsrB is present |
| `LIKELY_TtrA_or_SrdA` | TtrA/SrdA tree placement |
| `LIKELY_SoeA_or_divergent` | bSreA/SoeA/aSreA tree placement, or unresolved placement with no TAT/NrfD/PsrC evidence |
| `LIKELY_nonPsr_MopB` | Known non-Psr DMSOR/MopB tree placement |
| `LIKELY_nonPsr_MopB_operon_like` | Known non-Psr tree placement with A+B+C-like operon architecture (Mo-bisPGD A subunit, B-like Fe-S neighbour, NrfD-family membrane subunit); flagged for manual interpretation |
| `NOT_MoBisPGD_enzyme` | No PF00384 Mo-bisPGD confirmation |
| `AMBIGUOUS` | Conflicting or insufficient evidence |

### Diagnostic columns

`classification_table.tsv` includes nearest/second-nearest reference diagnostics:

| Column | Description |
|--------|-------------|
| `tree_group` | Broad enzyme-family group assigned from nearest reference clade |
| `tree_distance` | Patristic distance to nearest reference family |
| `tree_second_clade` | Second-nearest labelled reference family |
| `tree_second_distance` | Patristic distance to second-nearest reference family |
| `tree_distance_margin` | `tree_second_distance − tree_distance`; larger values indicate a more confident placement |
| `tree_distance_ratio` | `tree_distance / tree_second_distance`; lower values indicate a more confident placement |
| `tree_assignment_confidence` | `HIGH` (margin ≥ 0.20 or ratio ≤ 0.75), `MEDIUM` (margin ≥ 0.05 or ratio ≤ 0.90), `LOW` otherwise; heuristic, for flagging manual review only — not used as a hard classification threshold |

Operon context flags (also in `classification_table.tsv`):

| Column | Description |
|--------|-------------|
| `PsrABC_like_operon` | Strict flag: TAT + PsrB + (PsrC topology or PF14589 or NrfD) all present |
| `tree_incompatible_with_Psr` | `PsrABC_like_operon=YES` but tree group is Arr/Arx, Soe/Sre, Ttr/Srd, or known non-Psr MopB |
| `ABC_like_nonPsr_operon` | A+B+C-like local architecture (Mo-bisPGD + NrfD + PsrB) combined with known non-Psr MopB tree placement |
| `PsrC_like_context` | PsrC topology (8 TM) or PF14589 hit present |

These diagnostics are for interpretation and manual review; they do not directly alter classifications.

### MAFFT `--anysymbol` flag

The `--anysymbol` flag is passed to MAFFT to accommodate non-standard residues (e.g. selenocysteine, encoded as U) present in some reference sequences including formate dehydrogenase outgroups.

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

Avoid generic labels such as `custom_reference_001`, because they cannot act as clade anchors. Longer, more specific labels are matched before shorter ones during tree leaf annotation to avoid partial-prefix ambiguity.

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

`--redo-from` accepts integers only. IQ-TREE fast and full runs coexist with separate prefixes and sentinels (`psr_phylogeny_fast.*` vs `psr_phylogeny.*`). The `--custom-tree-refs` FASTA is always re-incorporated after Step 6 — running `--redo-from 7` is sufficient to update the tree after editing the file.

---

## Key Output Files

| File | Description |
|------|-------------|
| `09_summary/classification_table.html` | Colour-coded HTML summary; `HMSS2_gate` candidates highlighted with orange border |
| `09_summary/classification_table.tsv` | Full tree-gated classification table; includes discovery gate, tree group, tree diagnostics, and operon-context flags |
| `00_scan/discovery_source.tsv` | Per-candidate gate origin (`PF00384_gate`, `HMSS2_gate`, `both_gates`, `supplied`) |
| `06_references/tree_references_all.faa` | Exact reference FASTA used in alignment (standard + custom) |
| `06_references/reference_metadata_with_custom.tsv` | Metadata for all references including auto-generated custom/Wells entries |
| `03_hmmer/operon_completeness.tsv` | PsrABC-like subunit counts per candidate (`ABC_complete`, `AC_only`, `AB_only`, `A_only`) |
| `03_hmmer/nrfd_hits.tsv` | NrfD/PsrC hits with PF14589 specificity flag |
| `03_hmmer/hmss2/hmss2_operon.tsv` | HMSS2 annotation per candidate (annotation only, except `SoeA.hmm` which contributes to scoring) |
| `08_tree/psr_phylogeny.treefile` | IQ-TREE maximum-likelihood tree (full run) |
| `08_tree/psr_phylogeny_fast.treefile` | IQ-TREE fast tree (if `--fast-tree` used) |
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
