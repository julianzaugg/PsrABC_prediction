# PsrABC Identification & Phylogenetic Pipeline

## Overview

This pipeline identifies true **PsrABC** (polysulfide reductase) operons from
metagenomic bin FAA/GFF files, distinguishing them from functional homologues
(PhsABC, TtrABC, SoeABC, SreABC) that share the same Mo-bisPGD catalytic domain
but differ in subunit architecture and substrate specificity.

### Classification Logic

```
                     ┌─ PsrA annotation ─────────────────────────────────────┐
                     │                                                         │
            PF00384? ─────── NO ──────────────────────► NOT_MoBisPGD_enzyme   │
                     │                                                         │
                    YES                                                        │
                     │                                                         │
            TAT signal? ─── NO ───────────────────────► LIKELY_SoeA           │
                     │                                                         │
                    YES                                                        │
                     │                                                         │
  NrfD (PF03916)     │                                                         │
  in neighbourhood? ─ NO ─────────────────────────────► PsrA_missing_PsrC    │
                     │                                   (genome incomplete?)  │
                    YES                                                        │
                     │                                                         │
    TM count of     8TM ─────────────────────────────► TRUE_PsrA  ✓          │
    NrfD hit?       9TM ─────────────────────────────► LIKELY_TtrA           │
                    5TM + haem ──────────────────────► LIKELY_PhsA           │
                    5TM no haem ─────────────────────► PhsA or SoeC         │
                     └───────────────────────────────────────────────────────┘
```

Confirmed by: **phylogenetic tree position** (nearest reference clade)

---

## Reference Sequences

> ⚠️ Accessions marked **corrected** were wrong in earlier versions. Verify at
> UniProt before re-running. If you have a previous run, re-run from Step 6
> (`--redo-from 6`) after verifying.

| Label | Organism | Protein | Accession | Clade | Note |
|-------|----------|---------|-----------|-------|------|
| PsrA_Wolinella_succinogenes | *W. succinogenes* DSM1740 | PsrA | **P31075** (UniProt) | PsrA | Krafft et al. 1992. **NOT P31077** (= PsrC subunit!) |
| PsrA_Thermus_thermophilus | *T. thermophilus* HB8 | PsrA | **Q72LA4** (UniProt) | PsrA | PDB 2VPZ; Jormakka et al. 2008 |
| PhsA_Salmonella_enterica_LT2 | *S. enterica* Typhimurium LT2 | PhsA | **P37600** (UniProt) | PhsA | Heinzinger et al. 1995 |
| TtrA_Salmonella_enterica_LT2 | *S. enterica* Typhimurium LT2 | TtrA | **Q9Z4S6** (UniProt) | TtrA | Hensel et al. 1999 |
| TtrA_Shewanella_ANA3 | *Shewanella* sp. ANA-3 | TtrA | **WP_011715816.1** (NCBI) | TtrA | Degré et al. 2026 |
| SreA_Acidianus_ambivalens | *A. ambivalens* | SreA | **Q8NKK1** (UniProt) ← corrected | SreA | Laska et al. 2003. Verify gene = `sreA` at UniProt |
| SoeA_Allochromatium_vinosum | *A. vinosum* DSM180 | SoeA | **D3RNN8** (UniProt) | SoeA | Dahl et al. 2013; no TAT signal |
| ArrA_Shewanella_ANA3 | *Shewanella* sp. ANA-3 | ArrA | **Q7WTU0** (UniProt) | ArrA | Little et al. 2024 |
| ArrA_Chrysiogenes_arsenatis | *C. arsenatis* | ArrA | **Q5Y818** (UniProt) ← corrected | ArrA | Krafft & Macy 1998. Verify at UniProt |
| AioA_Alcaligenes_faecalis | *Alcaligenes faecalis* | AioA | **Q7SIF4** (UniProt) | AioA | Little et al. 2024 |
| NapA_Shewanella_oneidensis | *S. oneidensis* MR-1 | NapA | **Q8EIJ1** (UniProt) | NapA | TAT-exported |
| TorA_Ecoli | *E. coli* | TorA | **P33225** (UniProt) | TorA | TAT-exported |
| DmsA_Ecoli | *E. coli* | DmsA | **P18775** (UniProt) | DmsA | TAT-exported |
| NarG_Ecoli | *E. coli* K-12 | NarG | **P09152** (UniProt) | NarG | No TAT signal |
| SerA_Thauera_selenatis | *T. selenatis* | SerA | **Q9S1H0** (UniProt) | SerA | TAT-exported |
| PcrA_Dechloromonas_aromatica | *D. aromatica* | PcrA | **Q47CW6** (UniProt) | PcrA | TAT-exported |
| FdhG_Ecoli_K12 (outgroup) | *E. coli* K-12 | FdhG | **P24183** (UniProt) | Outgroup | Formate DH-N alpha |
| FdhH_Ecoli_K12 (outgroup) | *E. coli* K-12 | FdhH | **P07658** (UniProt) | Outgroup | Contains selenocysteine — requires `--anysymbol` in MAFFT |

---

## Prerequisites

### Conda environment

```bash
conda create -n biotools -c bioconda -c conda-forge \
    hmmer seqkit mafft trimal iqtree biopython requests
conda activate biotools
```

### Pfam database (required for HMMER)

```bash
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
```

### Optional local tools (web alternatives available)

- **DeepTMHMM** (TM topology prediction)
  ```bash
  pip install biolib && biolib install DTU/DeepTMHMM
  # OR use web: https://dtu.biolib.com/DeepTMHMM
  ```

- **SignalP 6.0** (TAT signal prediction)
  - Academic download: https://services.healthtech.dtu.dk/services/SignalP-6.0/
  - OR use web server (same URL)

- **ete3** (tree-based clade assignment)
  ```bash
  conda install -c etetoolkit ete3
  ```

---

## Input Files

| File | Description |
|------|-------------|
| `my_psrA_ids.txt` | One Pyrodigal protein ID per line (e.g. `NODE_3193_length_16941_cov_1.235637_2`) |
| `/path/to/bins/` | Directory tree containing `BinName.faa` + `BinName.gff` files |
| `Pfam-A.hmm` | Pfam database (hmmpress'd) |
| `hmss2_hmms/` | (Optional) HMSS2 HMM directory for sulfur-specific annotation |

---

## Usage

```bash
# Basic run (all steps)
bash 00_pipeline.sh \
    --ids    my_psrA_ids.txt \
    --bindir /path/to/bins/ \
    --pfam   /path/to/Pfam-A.hmm \
    --outdir results/psr_$(date +%Y%m%d) \
    --threads 16

# With optional HMSS2 annotation (adds sulfur HMM results as extra annotation columns)
bash 00_pipeline.sh \
    --ids    my_psrA_ids.txt \
    --bindir /path/to/bins/ \
    --pfam   /path/to/Pfam-A.hmm \
    --outdir results/psr_$(date +%Y%m%d) \
    --threads 16 \
    --hmss2  /path/to/HMSS2/Inorganic_Sulfur_Metabolism/

# Re-run from a specific step (e.g. after fixing references or adding --hmss2)
bash 00_pipeline.sh [args] --redo-from 6

# Re-run a single step only
bash 00_pipeline.sh [args] --redo-step 9

# Fast tree (exploration only; use full mode for publication)
bash 00_pipeline.sh [args] --fast-tree
```

---

## Manual Steps (Web Tools)

### Step 4 — DeepTMHMM (TM topology)

1. Go to: https://dtu.biolib.com/DeepTMHMM
2. Upload: `results/.../03_hmmer/nrfd_candidates.faa`
3. Download results as **3-line format**
4. Save as: `results/.../04_topology/deeptmhmm_results.3line`
5. Re-run: `bash 00_pipeline.sh [args] --redo-from 4`

### Step 5 — SignalP 6.0 (TAT signal)

1. Go to: https://services.healthtech.dtu.dk/services/SignalP-6.0/
2. Upload: `results/.../01_sequences/candidates_psrA.faa`
3. Select organism: **Other**
4. Download `prediction_results.txt`
5. Save as: `results/.../05_signalp/prediction_results.txt`
6. Re-run: `bash 00_pipeline.sh [args] --redo-from 5`

---

## Output Directory Structure

```
results/psr_YYYYMMDD/
├── 01_sequences/
│   └── candidates_psrA.faa              # Query PsrA candidate sequences
├── 02_neighbourhoods/
│   ├── PROTID_neighbourhood.faa         # ±10 gene neighbourhoods per candidate
│   └── PROTID_neighbourhood.tsv
├── 03_hmmer/
│   ├── PF00384_candidates_hits.tbl      # Mo-bisPGD confirmation
│   ├── PF03916_hits.tbl                 # NrfD-like hits in neighbourhoods (broad)
│   ├── PF14589_hits.tbl                 # NrfD_2 hits (high-specificity PsrC)
│   ├── PF12800_hits.tbl                 # PsrB hits (NrfC-like)
│   ├── PF13247_hits.tbl                 # PsrB hits (4Fe-4S, broad backup)
│   ├── nrfd_candidates.faa              # NrfD sequences → send to DeepTMHMM
│   ├── psrB_candidates.faa              # PsrB candidate sequences
│   ├── psrA_mo_domain_check.tsv         # PF00384 result per candidate
│   ├── nrfd_hits.tsv                    # NrfD hits with PF14589 flag
│   ├── psrB_hits.tsv
│   ├── operon_completeness.tsv          # Per-candidate ABC operon summary
│   └── hmss2/                           # (only if --hmss2 supplied)
│       ├── hmss2_candidates.tsv         # HMSS2 hits on catalytic subunits
│       ├── hmss2_neighbours.tsv         # HMSS2 hits on neighbourhood proteins
│       └── hmss2_operon.tsv             # Per-candidate summary → fed into Step 9
├── 04_topology/
│   └── topology_summary.tsv            # TM count → PsrC/TtrC/PhsC classification
├── 05_signalp/
│   └── tat_summary.tsv                 # TAT signal predictions
├── 06_references/
│   ├── references_all.faa
│   └── reference_metadata.tsv
├── 07_alignment/
│   ├── combined_for_alignment.faa
│   ├── combined_aligned.faa            # MAFFT output
│   └── combined_trimmed.faa            # TrimAl output
├── 08_tree/
│   ├── psr_phylogeny.treefile          # ← MAIN TREE (load in iTOL/FigTree)
│   ├── psr_phylogeny.iqtree
│   └── psr_phylogeny.contree
├── 09_summary/
│   ├── classification_table.tsv        # ← MAIN OUTPUT
│   └── classification_table.html       # Colour-coded view
└── logs/
    ├── pipeline.log
    ├── step01.done, step02.done ...     # Sentinels; delete to force re-run
    └── step3b.done                      # HMSS2 step sentinel
```

---

## Classification Table Column Guide

The main output is `09_summary/classification_table.tsv` (also as `.html`).
The HTML view colour-codes rows by classification; HMSS2 annotation columns
appear grey-shaded on the right.

### Primary evidence columns

- **`prot_id`** — Pyrodigal protein ID (contig-based, e.g. `NODE_3193_length_16941_cov_1.235637_2`). The final number is the gene index on that contig. The bin name is never in this ID — see `bin_name`.

- **`bin_name`** — The MAG bin this candidate belongs to. Populated from the protein→bin index built at startup; `NA` if the index was not supplied to Step 9.

- **`has_PF00384`** — `YES`/`NO`. Hit to the Mo-bisPGD catalytic domain profile (E ≤ 1×10⁻⁵). This is a hard prerequisite — any `NO` here means the sequence is not a confirmed Mo-bisPGD enzyme and should be treated with scepticism regardless of other columns.

- **`PF00384_evalue`** — E-value of the best Mo-bisPGD hit. True hits are typically well below 1×10⁻²⁰.

- **`has_TAT_signal`** — `YES`/`NO`/`NOT_RUN`. TAT/Tat-SPI signal peptide predicted by SignalP 6.0. **The single most important discriminator**: PsrA, TtrA, and PhsA are periplasmic (TAT+); SoeA is cytoplasmic (TAT−). A `NO` here with a confirmed Mo-bisPGD domain strongly suggests SoeA or a cytoplasmic homologue.

- **`TAT_probability`** — Raw SignalP 6.0 probability for the TAT(Tat/SPI) class (0–1). Values >0.5 are called positive. `NA` if SignalP was not run.

- **`signalp_prediction`** — Full SignalP 6.0 prediction label: `OTHER`, `SP(Sec/SPI)`, `TAT(Tat/SPI)`, `TATLIPO(Tat/SPII)`, `LIPO(Sec/SPII)`, or `PILIN(Sec/SPIII)`. `NOT_RUN` if SignalP was not executed. A `SP(Sec/SPI)` result (Sec rather than TAT export) is unusual for PsrA and warrants further investigation.

### Neighbourhood / operon columns

- **`NrfD_in_neighbourhood`** — `YES`/`NO`. Whether any protein within ±10 genes on the same contig hit a NrfD-family membrane subunit profile (PF03916 or PF14589, E ≤ 1×10⁻⁵). Detects the PsrC-type membrane anchor.

- **`n_NrfD_neighbours`** — Count of NrfD-family hits in the neighbourhood. >1 may indicate a duplicate, a mixed operon, or contig assembly artefact.

- **`NrfD_ids`** — Semicolon-separated protein IDs of all NrfD hits found in the neighbourhood.

- **`NrfD_has_PF14589`** — `YES`/`NO`. Whether any NrfD neighbour hit the high-specificity NrfD_2 profile (PF14589). A `YES` is strong evidence for a true PsrC subunit specifically, as opposed to TtrC or PhsC.

- **`best_NrfD_id`** — The single NrfD hit selected for topology classification. Selection priority: PF14589 hit first → resolved topology second → first hit otherwise.

- **`membrane_subunit_class`** — DeepTMHMM-based classification of the best NrfD hit:

  | Value | Meaning | Implication |
  |-------|---------|-------------|
  | `PsrC_8TM` | 8 TM helices, no haem | Strong PsrA evidence |
  | `TtrC_9TM` | 9 TM helices | Suggests TtrA |
  | `PhsC_5TM_haem` | 5 TM + CXXCH motif | Suggests PhsA |
  | `PhsC_or_SoeC` | 5 TM, no haem | Ambiguous; check tree |
  | `SoeC` | 4–6 TM, variable | Suggests SoeA context |
  | `NrfD_PF03916(topology_ND)` | HMM hit confirmed, DeepTMHMM not yet run | **Run DeepTMHMM to resolve** |
  | `NrfD_PF14589(+)` | High-specificity hit, topology not yet run | High-confidence PsrC pending topology |
  | `topology_not_run` / `NOT_RUN` | DeepTMHMM not executed | Run Step 4 |
  | `no_NrfD_found` | No NrfD protein in neighbourhood | Missing PsrC; incomplete genome? |

- **`PsrB_in_neighbourhood`** — `YES`/`NO`. Hit to PF12800 (NrfC-like) or PF13247 (4Fe-4S dicluster) in the neighbourhood, detecting the PsrB electron transfer subunit. **Absence is not penalised** — 4Fe-4S subunits are highly sequence-divergent and frequently missed (Rothery et al. 2008).

- **`n_PsrB_neighbours`** — Count of PsrB candidate hits.

- **`PsrB_ids`** — Protein IDs of PsrB candidates found.

### Phylogeny column

- **`tree_clade`** — Nearest reference sequence in the ML phylogeny by branch-length distance (assigned by ete3). Values are reference labels (e.g. `PsrA_Wolinella_succinogenes`, `SoeA_Allochromatium_vinosum`). Returns `inspect_tree_manually` if ete3 is not installed or assignment failed. Tree evidence contributes at most +2 to the score — deliberately capped so biochemical evidence dominates.

### Classification columns

- **`classification`** — Final integrated call:

  | Value | Meaning | Score threshold |
  |-------|---------|-----------------|
  | `TRUE_PsrA` | Mo-bisPGD + TAT + PsrC (8TM) all confirmed | ≥ 8 |
  | `LIKELY_PsrA` | Most evidence consistent with PsrA; one line missing or unresolved | ≥ 5 |
  | `PsrA_or_PhsA` | Mo-bisPGD + TAT confirmed; membrane subunit class ambiguous | ≥ 2 + TAT + Mo |
  | `LIKELY_SoeA_or_divergent` | Mo-bisPGD present, **no TAT**, no NrfD neighbour (all three required) | — |
  | `NOT_MoBisPGD_enzyme` | Failed PF00384 — not a Mo-bisPGD enzyme | — |
  | `AMBIGUOUS` | Conflicting or incomplete evidence; inspect manually | — |

- **`confidence`** — `HIGH`, `MEDIUM`, or `LOW`. Reflects how cleanly the evidence aligned, not the raw score.

- **`evidence`** — Pipe-separated tokens showing exactly what contributed to the score. Read this when a classification is unexpected. Tokens and their score contributions:

  | Token | Score |
  |-------|-------|
  | `Mo-bisPGD(+)` | +2 |
  | `Mo-bisPGD(-)` | −2 |
  | `TAT(+)` | +2 |
  | `TAT(-)` | 0 |
  | `PsrC_8TM(+)` | +3 |
  | `TtrC_9TM` | −1 |
  | `PhsC_5TM_haem` | 0 |
  | `SoeC` | −1 |
  | `NrfD_PF14589(+)` | +2 |
  | `NrfD_PF03916(topology_ND)` | +1 |
  | `NrfD(-)` | −1 |
  | `PsrB(+)` | +1 |
  | `PsrB(not_found)` | 0 |
  | `Tree:LABEL` | +2 (PsrA), +1 (PhsA), −1 (SoeA/TtrA); capped at +2 total |

### HMSS2 annotation columns (grey in HTML — annotation only, not scored)

Only present when the pipeline was run with `--hmss2`. These columns record hits to purpose-built sulfur metabolism HMMs from the HMSS2 database. They do **not** influence `classification` or `confidence` — they are for cross-validation and manual inspection only.

- **`HMSS2_PsrAPhsASreA`** / **`_evalue`** — Hit to the HMSS2 PsrA/PhsA/SreA catalytic HMM. More specific than PF00384 (which detects all Mo-bisPGD enzymes). A `YES` supports a genuine polysulfide/thiosulfate/sulfur reductase identity.

- **`HMSS2_SoeA`** / **`_evalue`** — Hit to the HMSS2 SoeA HMM. The most diagnostically valuable HMSS2 column: a `YES` converts a `LIKELY_SoeA_or_divergent` call from purely negative evidence into a positive identification.

- **`HMSS2_TtrA`** / **`_evalue`** — Hit to the HMSS2 TtrA HMM. Helps distinguish TtrA from PsrA when topology data is unavailable.

- **`HMSS2_neigh_PsrBPhsBSreB`** / **`_evalue`** — Neighbourhood hit to the HMSS2 PsrB/PhsB/SreB HMM. More specific than PF12800/PF13247 for the electron transfer subunit.

- **`HMSS2_neigh_PsrCPhsCSreC`** / **`_evalue`** — Neighbourhood hit to the HMSS2 PsrC/PhsC/SreC HMM. Note: Pfam PF03916+PF14589 outperform this in practice — treat as supporting evidence only.

- **`HMSS2_neigh_SoeB`**, **`HMSS2_neigh_SoeC`**, **`HMSS2_neigh_TtrB`**, **`HMSS2_neigh_TtrC`** — Neighbourhood hits to respective HMSS2 subunit HMMs. Useful for confirming context: e.g. a `LIKELY_SoeA` candidate with a `SoeC` neighbourhood hit is a more confident SoeA call.

---

## Interpreting the Tree

Load `08_tree/psr_phylogeny.treefile` in **iTOL** (https://itol.embl.de),
**FigTree** (https://github.com/rambaut/figtree), or **R ggtree**.

Query sequences retain their Pyrodigal IDs.
Reference sequences are labelled `CLADE_ORGANISM__ACCESSION`.

**Key clades:**
- Clustering with `PsrA_Wolinella_succinogenes` / `PsrA_Thermus_thermophilus` → **TRUE PsrA**
- Clustering with `PhsA_Salmonella_enterica_LT2` → **PhsA-type** (thiosulfate reductase)
- Clustering with `SoeA_Allochromatium_vinosum` → **SoeA-type** (sulfite dehydrogenase, NO TAT)
- Clustering with `TtrA_Salmonella_enterica_LT2` → **TtrA-type** (tetrathionate reductase)

---

## Key Diagnostic Rules Summary

| Evidence | PsrA | PhsA | TtrA | SoeA |
|----------|------|------|------|------|
| PF00384 (Mo-bisPGD) | ✓ | ✓ | ✓ | ✓ |
| TAT signal | ✓ | ✓ | ✓ | **✗** |
| NrfD (PF03916) in operon | ✓ | ✓ | ✓ | ✓ |
| Membrane subunit TM count | **8** | **5** | **9** | variable |
| Haem in membrane subunit | **✗** | **✓** (b-haem) | ✗ | variable |
| PF14589 (NrfD_2, high-specificity) | **✓** | variable | variable | ✗ |

---

## Citation

If you use this pipeline, please cite:
- Little et al. (2024) Nat. Microbiol. doi:10.1038/s41564-023-01560-2
- Degré et al. (2026) Environ. Microbiol. doi:10.1111/1462-2920.70258
- Jormakka et al. (2008) Nat. Struct. Mol. Biol. 15:730 (PDB 2VPZ, PsrC topology)
- Rothery et al. (2008) BBA 1778:1897 (CISM family; PsrB subunit divergence)
- Teufel et al. (2022) Nat. Biotechnol. (SignalP 6.0)
- Hallgren et al. (2022) bioRxiv (DeepTMHMM)
- Minh et al. (2020) Mol. Biol. Evol. (IQ-TREE 2)
- Katoh & Standley (2013) Mol. Biol. Evol. (MAFFT)
- Capella-Gutiérrez et al. (2009) Bioinformatics (TrimAl)
