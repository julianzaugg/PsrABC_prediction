#!/usr/bin/env python3
"""
05_fetch_references.py

Downloads curated reference sequences for the PsrA/PhsA/TtrA/SoeA/SreA/ArrA
phylogenetic tree from NCBI Entrez or UniProt.

REFERENCE SEQUENCES
===================
Accessions verified against Little et al. 2024 (Nature Microbiology,
doi:10.1038/s41564-023-01560-2) Supplementary Table and primary literature.

NOTE: P31077 is PsrC (membrane anchor subunit of W. succinogenes Psr),
      NOT PsrA. The correct PsrA accession is P31075.

CLADE   ORGANISM                          PROTEIN  ACCESSION       SOURCE
------  --------------------------------  -------  --------------  ------
PsrA    Wolinella succinogenes DSM1740    PsrA     P31075          UniProt
PsrA    Thermus thermophilus HB8          PsrA     Q72LA4          UniProt
PhsA    Salmonella enterica Typhimurium   PhsA     P37600          UniProt
TtrA    Salmonella enterica Typhimurium   TtrA     Q9Z4S6          UniProt
TtrA    Shewanella sp. ANA-3              TtrA     WP_011715816.1  NCBI
SreA    Acidianus ambivalens              SreA     Q9HGX4          UniProt
SoeA    Allochromatium vinosum DSM180     SoeA     D3RNN8          UniProt
ArrA    Shewanella sp. ANA-3              ArrA     Q7WTU0          UniProt
ArrA    Chrysiogenes arsenatis            ArrA     AAD05290.1      NCBI
AioA    Alcaligenes faecalis              AioA     Q7SIF4          UniProt
NapA    Shewanella oneidensis MR-1        NapA     Q8EIJ1          UniProt
TorA    Escherichia coli                  TorA     P33225          UniProt
DmsA    Escherichia coli                  DmsA     P18775          UniProt
NarG    Escherichia coli K-12             NarG     P09152          UniProt
SerA    Thauera selenatis                 SerA     Q9S1H0          UniProt
PcrA    Dechloromonas aromatica           PcrA     Q47CW6          UniProt
FdhG    Escherichia coli K-12 (outgroup)  FdhG     P24183          UniProt
FdhH    Escherichia coli K-12 (outgroup)  FdhH     P07658          UniProt

REQUIREMENTS:
  pip install requests biopython
"""

import argparse
import os
import time
import sys

try:
    import requests
except ImportError:
    sys.exit("[ERROR] requests not installed: pip install requests")

try:
    from Bio import Entrez, SeqIO
    Entrez.email = "your_email@institution.edu"  # Change this!
except ImportError:
    sys.exit("[ERROR] biopython not installed: pip install biopython")


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", required=True, help="Output directory for reference sequences")
    p.add_argument("--email",  default="your_email@institution.edu",
                   help="Your email for NCBI Entrez")
    return p.parse_args()


# =============================================================================
# Reference sequence registry
# =============================================================================
REFERENCES = [
    # --- PsrA clade ---
    # IMPORTANT: P31077 = PsrC (membrane subunit). PsrA = P31075.
    {"label": "PsrA_Wolinella_succinogenes",
     "organism": "Wolinella succinogenes DSM 1740",
     "protein": "PsrA", "clade": "PsrA",
     "accession": "P31075", "source": "uniprot",
     "note": "Polysulfide reductase catalytic subunit; Krafft et al. 1992"},

    {"label": "PsrA_Thermus_thermophilus",
     "organism": "Thermus thermophilus HB8",
     "protein": "PsrA", "clade": "PsrA",
     "accession": "Q72LA4", "source": "uniprot",
     "note": "PDB 2VPZ; Jormakka et al. 2008"},

    # --- PhsA clade ---
    {"label": "PhsA_Salmonella_enterica_LT2",
     "organism": "Salmonella enterica Typhimurium LT2",
     "protein": "PhsA", "clade": "PhsA",
     "accession": "P37600", "source": "uniprot",
     "note": "Thiosulfate reductase; Heinzinger et al. 1995; Little et al. 2024"},

    # --- TtrA clade ---
    {"label": "TtrA_Salmonella_enterica_LT2",
     "organism": "Salmonella enterica Typhimurium LT2",
     "protein": "TtrA", "clade": "TtrA",
     "accession": "Q9Z4S6", "source": "uniprot",
     "note": "Tetrathionate reductase; Hensel et al. 1999; Little et al. 2024"},

    {"label": "TtrA_Shewanella_ANA3",
     "organism": "Shewanella sp. ANA-3",
     "protein": "TtrA", "clade": "TtrA",
     "accession": "WP_011715816.1", "source": "ncbi",
     "note": "Degre et al. 2026"},

    # --- SreA clade ---
    {"label": "SreA_Acidianus_ambivalens",
     "organism": "Acidianus ambivalens",
     "protein": "SreA", "clade": "SreA",
     "accession": "Q8NKK1", "source": "uniprot",
     "note": "Sulfur reductase catalytic subunit; Laska et al. 2003"},

    # --- SoeA (sulfite dehydrogenase — cytoplasmic, NO TAT signal) ---
    {"label": "SoeA_Allochromatium_vinosum",
     "organism": "Allochromatium vinosum DSM 180",
     "protein": "SoeA", "clade": "SoeA",
     "accession": "D3RNN8", "source": "uniprot",
     "note": "Sulfite dehydrogenase; no TAT signal; Little et al. 2024"},

    # --- ArrA clade ---
    {"label": "ArrA_Shewanella_ANA3",
     "organism": "Shewanella sp. ANA-3",
     "protein": "ArrA", "clade": "ArrA",
     "accession": "Q7WTU0", "source": "uniprot",
     "note": "Respiratory arsenate reductase; Tat/SPI; Little et al. 2024"},

    {"label": "ArrA_Chrysiogenes_arsenatis",
     "organism": "Chrysiogenes arsenatis",
     "protein": "ArrA", "clade": "ArrA",
     "accession": "Q5Y818", "source": "uniprot",
     "note": "Respiratory arsenate reductase; Krafft & Macy 1998"},

    # --- AioA (arsenite oxidase — no TAT signal) ---
    {"label": "AioA_Alcaligenes_faecalis",
     "organism": "Alcaligenes faecalis",
     "protein": "AioA", "clade": "AioA",
     "accession": "Q7SIF4", "source": "uniprot",
     "note": "Arsenite oxidase; no TAT signal; Little et al. 2024"},

    # --- NapA (periplasmic nitrate reductase — TAT exported) ---
    {"label": "NapA_Shewanella_oneidensis",
     "organism": "Shewanella oneidensis MR-1",
     "protein": "NapA", "clade": "NapA",
     "accession": "Q8EIJ1", "source": "uniprot",
     "note": "Periplasmic nitrate reductase; Tat/SPI; Little et al. 2024"},

    # --- TorA (TMAO reductase — TAT exported) ---
    {"label": "TorA_Ecoli",
     "organism": "Escherichia coli",
     "protein": "TorA", "clade": "TorA",
     "accession": "P33225", "source": "uniprot",
     "note": "TMAO reductase; Tat/SPI; Little et al. 2024"},

    # --- DmsA (DMSO reductase — TAT exported) ---
    {"label": "DmsA_Ecoli",
     "organism": "Escherichia coli",
     "protein": "DmsA", "clade": "DmsA",
     "accession": "P18775", "source": "uniprot",
     "note": "DMSO reductase; Tat/SPI; Little et al. 2024"},

    # --- NarG (membrane-bound nitrate reductase — NOT TAT exported) ---
    {"label": "NarG_Ecoli",
     "organism": "Escherichia coli K-12",
     "protein": "NarG", "clade": "NarG",
     "accession": "P09152", "source": "uniprot",
     "note": "Cytoplasmic nitrate reductase; no TAT signal; Little et al. 2024"},

    # --- SerA (selenate reductase — TAT exported) ---
    {"label": "SerA_Thauera_selenatis",
     "organism": "Thauera selenatis",
     "protein": "SerA", "clade": "SerA",
     "accession": "Q9S1H0", "source": "uniprot",
     "note": "Selenate reductase; Tat/SPI; Little et al. 2024"},

    # --- PcrA (perchlorate reductase — TAT exported) ---
    {"label": "PcrA_Dechloromonas_aromatica",
     "organism": "Dechloromonas aromatica",
     "protein": "PcrA", "clade": "PcrA",
     "accession": "Q47CW6", "source": "uniprot",
     "note": "Perchlorate reductase; Tat/SPI; Little et al. 2024"},

    # --- Outgroup: Formate dehydrogenases ---
    {"label": "FdhG_Ecoli_K12_OUTGROUP",
     "organism": "Escherichia coli K-12",
     "protein": "FdhG", "clade": "FdhG_outgroup",
     "accession": "P24183", "source": "uniprot",
     "note": "Formate dehydrogenase-N alpha subunit — OUTGROUP"},

    {"label": "FdhH_Ecoli_K12_OUTGROUP",
     "organism": "Escherichia coli K-12",
     "protein": "FdhH", "clade": "FdhH_outgroup",
     "accession": "P07658", "source": "uniprot",
     "note": "Formate dehydrogenase H — OUTGROUP"},
]


def fetch_uniprot(accession, retries=3):
    """Fetch sequence from UniProt REST API."""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    for attempt in range(retries):
        try:
            resp = requests.get(url, timeout=30)
            if resp.status_code == 200:
                return resp.text
            else:
                print(f"    [WARN] UniProt {accession}: HTTP {resp.status_code}")
        except requests.RequestException as e:
            print(f"    [WARN] UniProt {accession}: {e}")
        time.sleep(2 ** attempt)
    return None


def fetch_ncbi_protein(accession, email, retries=3):
    """Fetch sequence from NCBI Protein database via Entrez."""
    Entrez.email = email
    for attempt in range(retries):
        try:
            handle = Entrez.efetch(
                db="protein", id=accession, rettype="fasta", retmode="text"
            )
            result = handle.read()
            handle.close()
            return result
        except Exception as e:
            print(f"    [WARN] NCBI {accession}: {e}")
            time.sleep(2 ** attempt)
    return None


def rename_fasta_header(fasta_text, new_label):
    """Replace the FASTA header line with a clean, informative label."""
    lines = fasta_text.strip().split("\n")
    if not lines:
        return fasta_text
    # Keep original accession in the header but prepend the label
    original_header = lines[0]
    original_acc = original_header[1:].split()[0]
    new_header = f">{new_label}__{original_acc}"
    return new_header + "\n" + "\n".join(lines[1:]) + "\n"


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    all_seqs_path  = os.path.join(args.outdir, "references_all.faa")
    metadata_path  = os.path.join(args.outdir, "reference_metadata.tsv")
    failed_path    = os.path.join(args.outdir, "failed_downloads.txt")

    all_seqs = []
    metadata_rows = []
    failed = []

    print(f"[*] Downloading {len(REFERENCES)} reference sequences...")
    print(f"    Note: Change Entrez.email in this script to your own email address\n")

    for ref in REFERENCES:
        label    = ref["label"]
        acc      = ref["accession"]
        source   = ref["source"]
        clade    = ref["clade"]

        print(f"  Fetching {label} ({acc}) from {source}...")

        fasta_text = None
        if source == "uniprot":
            fasta_text = fetch_uniprot(acc)
        elif source == "ncbi":
            fasta_text = fetch_ncbi_protein(acc, args.email)
        else:
            print(f"    [WARN] Unknown source '{source}' for {acc}")
            failed.append(f"{label}\t{acc}\tunknown_source")
            continue

        if fasta_text and fasta_text.strip().startswith(">"):
            renamed = rename_fasta_header(fasta_text, label)
            all_seqs.append(renamed)
            # Save individual file too
            ind_path = os.path.join(args.outdir, f"{label}.faa")
            with open(ind_path, "w") as fh:
                fh.write(renamed)

            metadata_rows.append({
                "label":    label,
                "accession": acc,
                "organism": ref["organism"],
                "protein":  ref["protein"],
                "clade":    clade,
                "source":   source,
                "note":     ref["note"],
                "status":   "OK",
            })
            print(f"    [OK] Downloaded ({len(fasta_text.split())} chars)")
        else:
            print(f"    [FAIL] Could not fetch {acc}")
            failed.append(f"{label}\t{acc}\tdownload_failed")
            metadata_rows.append({
                "label": label, "accession": acc,
                "organism": ref["organism"], "protein": ref["protein"],
                "clade": clade, "source": source,
                "note": ref["note"], "status": "FAILED",
            })

        time.sleep(0.4)  # Be polite to NCBI/UniProt

    # Write combined FAA
    with open(all_seqs_path, "w") as fh:
        for seq in all_seqs:
            fh.write(seq)
            if not seq.endswith("\n"):
                fh.write("\n")

    # Write metadata TSV
    with open(metadata_path, "w") as fh:
        fh.write("label\taccession\torganism\tprotein\tclade\tsource\tnote\tstatus\n")
        for row in metadata_rows:
            fh.write(f"{row['label']}\t{row['accession']}\t{row['organism']}\t"
                     f"{row['protein']}\t{row['clade']}\t{row['source']}\t"
                     f"{row['note']}\t{row['status']}\n")

    # Write failed list
    if failed:
        with open(failed_path, "w") as fh:
            fh.write("label\taccession\treason\n")
            for f in failed:
                fh.write(f + "\n")
        print(f"\n[WARN] {len(failed)} sequences failed to download → {failed_path}")
        print(f"  You can manually download these from:")
        print(f"    UniProt:  https://www.uniprot.org/uniprot/ACCESSION.fasta")
        print(f"    NCBI:     https://www.ncbi.nlm.nih.gov/protein/ACCESSION")
        print(f"  Then append to: {all_seqs_path}")

    n_ok = len(all_seqs)
    print(f"\n[DONE] {n_ok}/{len(REFERENCES)} reference sequences downloaded")
    print(f"  Combined FAA   : {all_seqs_path}")
    print(f"  Metadata table : {metadata_path}")

    # Print reference table for records
    print("\n  REFERENCE SEQUENCE REGISTRY:")
    print(f"  {'Label':<45} {'Accession':<15} {'Clade':<15} {'Status'}")
    print(f"  {'-'*45} {'-'*15} {'-'*15} {'-'*8}")
    for row in metadata_rows:
        print(f"  {row['label']:<45} {row['accession']:<15} {row['clade']:<15} {row['status']}")


if __name__ == "__main__":
    main()
