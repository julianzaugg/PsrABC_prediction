#!/usr/bin/env python3
"""
03_parse_topology.py

Parses DeepTMHMM output (3-line format or TSV summary) and classifies
NrfD-like proteins by TM helix count into:

  PsrC  : 8 TM helices, no haem (NrfD subtype — polysulfide reductase anchor)
  TtrC  : 9 TM helices           (tetrathionate reductase anchor)
  PhsC  : 5 TM helices + haem    (thiosulfate reductase anchor — note: haem
                                   prediction requires separate check)
  SoeC  : variable, typically 4-6 TM helices (sulfite dehydrogenase anchor)
  Other : unexpected TM count

DeepTMHMM 3-line output format:
  >PROTEIN_ID  topology_type  [topology_string]
  AAAAAAA...   (sequence)
  MMMIIIOOOO   (topology: M=TM, I=inside, O=outside, S=signal)

DeepTMHMM TSV summary (if downloaded from web):
  protein_id  \t  topology_type  \t  n_TM  \t  topology_string

We also check for CXXCH haem-binding motifs in the sequence to help
distinguish PhsC from PsrC (PhsC has 2 b-type haems).

USAGE:
  python3 03_parse_topology.py \
      --tmhmm_results /path/to/deeptmhmm_output_dir/ \
      --nrfd_hits PF03916_hits.tbl \
      --ids original_ids.txt \
      --outdir 04_topology/
"""

import argparse
import os
import re
import glob


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--tmhmm_results", required=True,
                   help="Directory with DeepTMHMM output files")
    p.add_argument("--nrfd_hits",     required=True,
                   help="PF03916 HMMER tblout from step 3")
    p.add_argument("--ids",           required=True,
                   help="Original query psrA IDs file")
    p.add_argument("--outdir",        required=True)
    return p.parse_args()


def classify_by_tm_count(n_tm, has_haem_motif):
    """
    Rule-based classification based on published subunit topologies.
    Ref: Jormakka et al. 2008 (PsrC=8TM), Rothery et al. 2008 (TtrC=9TM),
         Stoffels et al. 2012 (PhsC=5TM+haem)
    """
    if n_tm == 8 and not has_haem_motif:
        return "PsrC", "HIGH"
    elif n_tm == 8 and has_haem_motif:
        return "PsrC_or_ambiguous", "MEDIUM"
    elif n_tm == 9:
        return "TtrC", "HIGH"
    elif n_tm == 5 and has_haem_motif:
        return "PhsC", "HIGH"
    elif n_tm == 5 and not has_haem_motif:
        return "PhsC_or_SoeC", "MEDIUM"
    elif n_tm in (4, 6):
        return "SoeC_like", "LOW"
    elif n_tm <= 3:
        return "Low_TM_not_PsrC", "LOW"
    else:
        return f"Unknown_{n_tm}TM", "LOW"


def count_tm_helices(topology_string):
    """Count TM segments from DeepTMHMM topology string (M=transmembrane)."""
    # Each continuous run of 'M' characters is one TM helix
    helices = re.findall(r"M+", topology_string.upper())
    return len(helices)


def has_haem_binding_motif(sequence):
    """Check for CXXCH (c-type haem) or CXXCHxxM (b-type haem) motifs."""
    # c-type haem: CXXCH
    c_type = bool(re.search(r"C.{2}CH", sequence.upper()))
    return c_type


def load_faa(faa_path):
    seqs = {}
    current_id = None
    current_seq = []
    with open(faa_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id:
                    seqs[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id:
        seqs[current_id] = "".join(current_seq)
    return seqs


def parse_deeptmhmm_3line(filepath):
    """
    Parse DeepTMHMM 3-line format output.
    Returns dict: {protein_id: {"topology": str, "n_tm": int, "type": str}}
    """
    results = {}
    with open(filepath) as fh:
        lines = [l.rstrip() for l in fh if l.strip()]

    i = 0
    while i < len(lines):
        if lines[i].startswith(">"):
            header = lines[i]
            prot_id = header[1:].split()[0]
            # Header may contain topology type after protein ID
            topo_type_match = re.search(r"\|\s*(\S+)", header)
            topo_type = topo_type_match.group(1) if topo_type_match else "unknown"

            seq_line   = lines[i+1] if i+1 < len(lines) else ""
            label_line = lines[i+2] if i+2 < len(lines) else ""

            n_tm = count_tm_helices(label_line)
            results[prot_id] = {
                "topology_string": label_line,
                "topology_type":   topo_type,
                "n_tm":            n_tm,
                "sequence":        seq_line,
            }
            i += 3
        else:
            i += 1
    return results


def parse_deeptmhmm_tsv(filepath):
    """
    Parse DeepTMHMM TSV summary (alternative download format from web server).
    Expected columns: ID, topology_type, n_TM, topology_string
    """
    results = {}
    with open(filepath) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 3:
                continue
            prot_id   = parts[0]
            topo_type = parts[1] if len(parts) > 1 else "unknown"
            n_tm_str  = parts[2] if len(parts) > 2 else "0"
            topo_str  = parts[3] if len(parts) > 3 else ""
            try:
                n_tm = int(n_tm_str)
            except ValueError:
                n_tm = count_tm_helices(topo_str)
            results[prot_id] = {
                "topology_string": topo_str,
                "topology_type":   topo_type,
                "n_tm":            n_tm,
                "sequence":        "",
            }
    return results


def find_deeptmhmm_output(results_dir):
    """Auto-detect the output format from DeepTMHMM."""
    # Look for 3-line format files
    for pattern in ["*.3line", "*.pred", "predicted_topologies.3line",
                    "TMRs.gff3", "deeptmhmm_results.3line"]:
        matches = glob.glob(os.path.join(results_dir, pattern))
        if matches:
            return matches[0], "3line"
    # Look for TSV
    for pattern in ["*.tsv", "*.txt"]:
        matches = glob.glob(os.path.join(results_dir, pattern))
        for m in matches:
            # Check if it looks like a topology TSV
            with open(m) as f:
                first = f.readline()
            if "TM" in first.upper() or "topology" in first.lower():
                return m, "tsv"
    return None, None


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Load original IDs for cross-reference
    with open(args.ids) as fh:
        query_ids = [l.strip() for l in fh if l.strip() and not l.startswith("#")]

    # Find and parse DeepTMHMM output
    output_file, fmt = find_deeptmhmm_output(args.tmhmm_results)

    tm_results = {}
    if output_file and fmt == "3line":
        print(f"[*] Parsing DeepTMHMM 3-line format: {output_file}")
        tm_results = parse_deeptmhmm_3line(output_file)
    elif output_file and fmt == "tsv":
        print(f"[*] Parsing DeepTMHMM TSV format: {output_file}")
        tm_results = parse_deeptmhmm_tsv(output_file)
    else:
        print(f"[WARN] No DeepTMHMM output found in {args.tmhmm_results}")
        print(f"       Manual step required — see pipeline instructions")
        print(f"       Creating placeholder topology summary...")

    # Load NrfD candidate sequences for haem motif checking
    nrfd_faa_path = os.path.join(os.path.dirname(args.nrfd_hits), "nrfd_candidates.faa")
    nrfd_seqs = {}
    if os.path.exists(nrfd_faa_path):
        nrfd_seqs = load_faa(nrfd_faa_path)

    # Classify each NrfD candidate
    rows = []
    for prot_id, info in tm_results.items():
        seq = nrfd_seqs.get(prot_id, info.get("sequence", ""))
        haem = has_haem_binding_motif(seq) if seq else False
        n_tm = info["n_tm"]
        subunit_class, confidence = classify_by_tm_count(n_tm, haem)
        rows.append({
            "protein_id":    prot_id,
            "n_tm":          n_tm,
            "topology_type": info.get("topology_type", "NA"),
            "haem_motif":    "YES" if haem else "NO",
            "subunit_class": subunit_class,
            "confidence":    confidence,
        })

    # Also add placeholder rows for NrfD hits without topology results
    with open(args.nrfd_hits) as fh:
        for line in fh:
            if line.startswith("#"): continue
            nrfd_id = line.split()[0]
            if nrfd_id not in {r["protein_id"] for r in rows}:
                seq = nrfd_seqs.get(nrfd_id, "")
                haem = has_haem_binding_motif(seq) if seq else False
                rows.append({
                    "protein_id":    nrfd_id,
                    "n_tm":          "ND",
                    "topology_type": "NOT_RUN",
                    "haem_motif":    "YES" if haem else "NO",
                    "subunit_class": "UNCLASSIFIED_run_DeepTMHMM",
                    "confidence":    "NONE",
                })

    # Write topology summary table
    out_path = os.path.join(args.outdir, "topology_summary.tsv")
    with open(out_path, "w") as fh:
        fh.write("protein_id\tn_TM_helices\ttopology_type\thaem_CXXCH_motif\t"
                 "subunit_classification\tconfidence\n")
        for row in rows:
            fh.write(f"{row['protein_id']}\t{row['n_tm']}\t{row['topology_type']}\t"
                     f"{row['haem_motif']}\t{row['subunit_class']}\t{row['confidence']}\n")

    print(f"\n[DONE] Topology classification: {len(rows)} NrfD-like proteins classified")
    print(f"       Classification rules:")
    print(f"         8 TM, no haem → PsrC   (polysulfide reductase)")
    print(f"         9 TM          → TtrC   (tetrathionate reductase)")
    print(f"         5 TM + haem   → PhsC   (thiosulfate reductase)")
    print(f"         5 TM, no haem → PhsC or SoeC (ambiguous)")
    print(f"       Summary table → {out_path}")


if __name__ == "__main__":
    main()
