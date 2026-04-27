#!/usr/bin/env python3
"""
04_parse_signalp.py

Parse SignalP 6.0 output to identify TAT signal peptides in PsrA candidates.

TAT signal = key diagnostic for PsrA/TtrA/PhsA vs SoeA:
  - PsrA, TtrA, PhsA : periplasmic, exported via TAT → SHOULD have TAT signal
  - SoeA              : cytoplasmic, NOT exported → NO TAT signal

SignalP 6.0 actual output format (prediction_results.txt):
  # SignalP-6.0   Organism: Other   Timestamp: ...
  # ID    Prediction    OTHER    SP(Sec/SPI)    LIPO(Sec/SPII)    TAT(Tat/SPI)    TATLIPO(Tat/SPII)    PILIN(Sec/SPIII)    CS Position
  NODE_223_...  # 5878 # 8892 # ... ID=6_8;...    OTHER    0.838    0.161    0.000    0.000    0.000    0.000

IMPORTANT: The sequence ID in SignalP output includes the full FASTA header
description (everything after ">"), not just the protein ID. We extract only
the first whitespace-delimited token to match against our protein IDs.

Column order (0-indexed after splitting on whitespace, after stripping comment header):
  col 0   : full header string (we take only the first word as protein_id)
  col -7  : Prediction label (OTHER / SP(Sec/SPI) / TAT(Tat/SPI) / etc.)
  col -6  : OTHER probability
  col -5  : SP(Sec/SPI) probability
  col -4  : LIPO(Sec/SPII) probability
  col -3  : TAT(Tat/SPI) probability  ← what we want
  col -2  : TATLIPO(Tat/SPII) probability
  col -1  : PILIN(Sec/SPIII) probability
  (CS Position may also be appended if predicted)

We parse from the right using negative indices to be robust to varying ID lengths.
"""

import argparse
import os
import glob


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--signalp_dir", required=True, help="Directory with SignalP 6 output")
    p.add_argument("--candidates",  required=True, help="Candidate PsrA FAA")
    p.add_argument("--outdir",      required=True)
    return p.parse_args()


def find_signalp_output(signalp_dir):
    """Find the SignalP 6 output file. Prioritises the standard filename."""
    exact = os.path.join(signalp_dir, "prediction_results.txt")
    if os.path.exists(exact):
        return exact
    for pat in ["summary.signalp5", "*.signalp6", "signalp6_output.txt"]:
        matches = glob.glob(os.path.join(signalp_dir, pat))
        if matches:
            return matches[0]
    return None


# Known SignalP 6 prediction labels — used to identify the Prediction column
_SP6_LABELS = {"OTHER", "SP(Sec/SPI)", "LIPO(Sec/SPII)",
               "TAT(Tat/SPI)", "TATLIPO(Tat/SPII)", "PILIN(Sec/SPIII)"}

def parse_signalp6_output(filepath):
    """
    Parse SignalP 6.0 prediction_results.txt.

    The ID field spans everything before the Prediction label, and may contain
    spaces because SignalP uses the full FASTA header as the ID. We extract
    the protein ID as the first whitespace-delimited word of that field.

    Columns after the Prediction label (always the last 6 float columns):
      OTHER  SP(Sec/SPI)  LIPO(Sec/SPII)  TAT(Tat/SPI)  TATLIPO  PILIN

    Returns dict: {protein_id: {"prediction": str, "TAT_prob": float, "has_TAT": bool}}
    """
    results = {}
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            
            # Use tab split to keep ID + metadata as a single block in parts[0]
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            
            # The short ID is the first word in the first tab-column
            # This matches NODE_4296..._7 in your FAA file
            full_id_string = parts[0]
            prot_id = full_id_string.split()[0] 
            
            # Prediction is always the second tab-column
            prediction = parts[1].strip()
            
            # Probability parsing - SignalP 6 order: OTHER, SP, LIPO, TAT, TATLIPO
            # These are indices 2, 3, 4, 5, 6 in the tab-split list
            tat_prob = 0.0
            try:
                if len(parts) > 5:
                    tat_prob = float(parts[5])
            except (ValueError, IndexError):
                tat_prob = 0.0

            results[prot_id] = {
                "prediction": prediction,
                "TAT_prob": tat_prob,
                "has_TAT": "TAT" in prediction.upper()
            }
    return results


def load_faa_ids(faa_path):
    ids = []
    with open(faa_path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    candidate_ids = load_faa_ids(args.candidates)
    print(f"[*] {len(candidate_ids)} candidate proteins to check for TAT signal")

    sp_file = find_signalp_output(args.signalp_dir)

    sp_results = {}
    if sp_file:
        print(f"[*] Parsing SignalP 6 output: {sp_file}")
        sp_results = parse_signalp6_output(sp_file)
        print(f"    {len(sp_results)} proteins parsed from SignalP output")
    else:
        print(f"[WARN] No SignalP 6 output found in {args.signalp_dir}")
        print(f"       MANUAL STEP: https://services.healthtech.dtu.dk/services/SignalP-6.0/")
        print(f"       Upload: {args.candidates}  |  Organism: other")
        print(f"       Save prediction_results.txt to: {args.signalp_dir}/")

    out_path = os.path.join(args.outdir, "tat_summary.tsv")
    n_tat = n_no_tat = n_unknown = 0

    with open(out_path, "w") as fh:
        fh.write("protein_id\tsignalp_prediction\tTAT_probability\thas_TAT\tinterpretation\n")
        for prot_id in candidate_ids:
            if prot_id in sp_results:
                info    = sp_results[prot_id]
                has_tat = info["has_TAT"]
                pred    = info["prediction"]
                prob    = f"{info['TAT_prob']:.4f}"
                if has_tat:
                    interp = "periplasmic_TAT_export—consistent_with_PsrA/TtrA/PhsA"
                    n_tat += 1
                elif "OTHER" in pred.upper():
                    interp = "no_signal_peptide—consistent_with_SoeA_or_cytoplasmic"
                    n_no_tat += 1
                elif "SP" in pred.upper():
                    interp = "Sec_signal_not_TAT—unusual_for_PsrA_check_further"
                    n_no_tat += 1
                else:
                    interp = f"no_TAT—prediction={pred}"
                    n_no_tat += 1
            else:
                pred    = "NOT_RUN"
                prob    = "NA"
                has_tat = False
                interp  = "SignalP_not_run—manual_step_required"
                n_unknown += 1

            fh.write(f"{prot_id}\t{pred}\t{prob}\t{'YES' if has_tat else 'NO'}\t{interp}\n")

    print(f"\n[DONE] TAT summary:")
    print(f"  TAT signal     : {n_tat}  → periplasmic (PsrA/TtrA/PhsA type)")
    print(f"  No TAT signal  : {n_no_tat}  → cytoplasmic or Sec-exported")
    print(f"  Not run        : {n_unknown}")
    print(f"  Output → {out_path}")


if __name__ == "__main__":
    main()
