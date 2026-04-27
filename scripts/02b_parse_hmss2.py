#!/usr/bin/env python3
"""
02b_parse_hmss2.py

Parse HMSS2 HMM search results for PsrABC/TtrABC/SoeABC profiles.

These results are ANNOTATION ONLY — they are not used in the final
classification scoring. They appear as extra columns in the summary
table so that hits can be inspected manually and used to validate or
question the primary evidence (PF-based HMMs, DeepTMHMM, SignalP).

Profiles expected (all at E <= 1e-5):

  Run against candidates.faa (catalytic subunit):
    PsrAPhsASreA.hmm  — PsrA/PhsA/SreA clade catalytic subunit
    SoeA.hmm          — SoeA sulfite dehydrogenase catalytic subunit
    TtrA.hmm          — TtrA tetrathionate reductase catalytic subunit

  Run against neighbours.faa (neighbourhood proteins):
    PsrBPhsBSreB.hmm  — PsrB/PhsB/SreB electron transfer subunit
    PsrCPhsCSreC.hmm  — PsrC/PhsC/SreC membrane anchor subunit
    SoeB.hmm          — SoeB electron transfer subunit
    SoeC.hmm          — SoeC membrane anchor subunit
    TtrB.hmm          — TtrB electron transfer subunit
    TtrC.hmm          — TtrC membrane anchor subunit

Output:
  hmss2_candidates.tsv  — one row per candidate, catalytic-subunit hits
  hmss2_neighbours.tsv  — one row per neighbourhood hit, subunit identity
  hmss2_operon.tsv      — per-candidate summary of neighbour hits (joined
                          back to psrA via "from:" tag in FASTA header)
"""

import argparse
import os
import re
import glob


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--hmss2_dir",  required=True,
                   help="Directory containing HMSS2 tblout files (*_hits.tbl)")
    p.add_argument("--candidates", required=True, help="Candidate PsrA FAA")
    p.add_argument("--neighbours", required=True, help="Combined neighbourhood FAA")
    p.add_argument("--ids",        required=True, help="Query IDs file")
    p.add_argument("--outdir",     required=True)
    return p.parse_args()


EVALUE_CUTOFF = 1e-5

# Profiles run against the candidate (catalytic) FAA
CAND_PROFILES = ["PsrAPhsASreA", "SoeA", "TtrA"]

# Profiles run against the neighbourhood FAA
NEIGH_PROFILES = ["PsrBPhsBSreB", "PsrCPhsCSreC", "SoeB", "SoeC", "TtrB", "TtrC"]

# Human-readable labels for each profile
PROFILE_LABELS = {
    "PsrAPhsASreA": "PsrA/PhsA/SreA_catalytic",
    "SoeA":         "SoeA_catalytic",
    "TtrA":         "TtrA_catalytic",
    "PsrBPhsBSreB": "PsrB/PhsB/SreB_electron_transfer",
    "PsrCPhsCSreC": "PsrC/PhsC/SreC_membrane_anchor",
    "SoeB":         "SoeB_electron_transfer",
    "SoeC":         "SoeC_membrane_anchor",
    "TtrB":         "TtrB_electron_transfer",
    "TtrC":         "TtrC_membrane_anchor",
}


def parse_tblout(path):
    """
    Returns dict: {target_id: best_evalue}
    Only retains hits at or below EVALUE_CUTOFF.
    """
    hits = {}
    if not os.path.exists(path):
        return hits
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            if len(fields) < 5:
                continue
            target = fields[0]
            try:
                evalue = float(fields[4])
            except ValueError:
                continue
            if evalue <= EVALUE_CUTOFF:
                if target not in hits or evalue < hits[target]:
                    hits[target] = evalue
    return hits


def load_faa_ids(faa_path):
    """Returns list of protein IDs from a FASTA file."""
    ids = []
    with open(faa_path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def load_neighbour_to_psra(faa_path):
    """
    Parse the 'from:PSRA_ID' tag embedded by 01_extract_neighbourhood.py
    in neighbourhood FASTA headers.
    Returns dict: {neighbour_id: psrA_id}
    """
    mapping = {}
    with open(faa_path) as fh:
        for line in fh:
            if line.startswith(">"):
                neigh_id = line[1:].split()[0]
                m = re.search(r"from:(\S+)", line)
                psra_id = m.group(1) if m else "unknown"
                mapping[neigh_id] = psra_id
    return mapping


def best_evalue_str(ev):
    return f"{ev:.2e}" if ev is not None else "NA"


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    with open(args.ids) as fh:
        query_ids = [l.strip() for l in fh if l.strip() and not l.startswith("#")]

    candidate_ids  = load_faa_ids(args.candidates)
    neighbour_map  = load_neighbour_to_psra(args.neighbours)   # neigh_id → psrA_id

    # -------------------------------------------------------------------------
    # Load all tblout results
    # -------------------------------------------------------------------------
    cand_hits  = {}   # profile → {prot_id: evalue}
    neigh_hits = {}   # profile → {prot_id: evalue}

    for prof in CAND_PROFILES:
        tbl = os.path.join(args.hmss2_dir, f"{prof}_hits.tbl")
        cand_hits[prof] = parse_tblout(tbl)
        n = len(cand_hits[prof])
        if n:
            print(f"  [{prof}] {n} candidate hits")
        else:
            print(f"  [{prof}] no hits (or file missing)")

    for prof in NEIGH_PROFILES:
        tbl = os.path.join(args.hmss2_dir, f"{prof}_hits.tbl")
        neigh_hits[prof] = parse_tblout(tbl)
        n = len(neigh_hits[prof])
        if n:
            print(f"  [{prof}] {n} neighbourhood hits")
        else:
            print(f"  [{prof}] no hits (or file missing)")

    # -------------------------------------------------------------------------
    # Write 1: hmss2_candidates.tsv
    # One row per candidate. Columns: one YES/evalue per catalytic profile.
    # -------------------------------------------------------------------------
    cand_out = os.path.join(args.outdir, "hmss2_candidates.tsv")
    cand_cols = ["protein_id"] + [f"HMSS2_{p}" for p in CAND_PROFILES] + \
                [f"HMSS2_{p}_evalue" for p in CAND_PROFILES]

    with open(cand_out, "w") as fh:
        fh.write("\t".join(cand_cols) + "\n")
        for prot_id in candidate_ids:
            row = [prot_id]
            for prof in CAND_PROFILES:
                ev = cand_hits[prof].get(prot_id)
                row.append("YES" if ev is not None else "NO")
            for prof in CAND_PROFILES:
                ev = cand_hits[prof].get(prot_id)
                row.append(best_evalue_str(ev))
            fh.write("\t".join(row) + "\n")

    print(f"\n  Candidate HMSS2 hits  → {cand_out}")

    # -------------------------------------------------------------------------
    # Write 2: hmss2_neighbours.tsv
    # One row per neighbourhood protein that hit any HMSS2 profile.
    # -------------------------------------------------------------------------
    neigh_out = os.path.join(args.outdir, "hmss2_neighbours.tsv")
    all_neigh_hits = {}   # neigh_id → {prof: evalue}
    for prof, hits in neigh_hits.items():
        for nid, ev in hits.items():
            all_neigh_hits.setdefault(nid, {})[prof] = ev

    with open(neigh_out, "w") as fh:
        header = ["neighbour_id", "nearest_psrA"] + \
                 [f"HMSS2_{p}" for p in NEIGH_PROFILES] + \
                 [f"HMSS2_{p}_evalue" for p in NEIGH_PROFILES]
        fh.write("\t".join(header) + "\n")
        for nid in sorted(all_neigh_hits.keys()):
            psra = neighbour_map.get(nid, "unknown")
            row  = [nid, psra]
            for prof in NEIGH_PROFILES:
                ev = all_neigh_hits[nid].get(prof)
                row.append("YES" if ev is not None else "NO")
            for prof in NEIGH_PROFILES:
                ev = all_neigh_hits[nid].get(prof)
                row.append(best_evalue_str(ev))
            fh.write("\t".join(row) + "\n")

    print(f"  Neighbour HMSS2 hits  → {neigh_out}")

    # -------------------------------------------------------------------------
    # Write 3: hmss2_operon.tsv
    # One row per query candidate. For each neighbourhood profile, report
    # whether any neighbour of that candidate hit the profile, and the best
    # evalue seen. This is the table that feeds into 06_build_summary.py.
    # -------------------------------------------------------------------------
    # Build psrA → {prof: best_evalue} from neighbourhood hits
    psra_neigh_summary = {qid: {} for qid in query_ids}
    for prof, hits in neigh_hits.items():
        for nid, ev in hits.items():
            psra_id = neighbour_map.get(nid, "unknown")
            if psra_id in psra_neigh_summary:
                existing = psra_neigh_summary[psra_id].get(prof)
                if existing is None or ev < existing:
                    psra_neigh_summary[psra_id][prof] = ev

    operon_out = os.path.join(args.outdir, "hmss2_operon.tsv")
    operon_cols = ["psrA_id"] + \
                  [f"HMSS2_{p}" for p in CAND_PROFILES] + \
                  [f"HMSS2_{p}_evalue" for p in CAND_PROFILES] + \
                  [f"HMSS2_neigh_{p}" for p in NEIGH_PROFILES] + \
                  [f"HMSS2_neigh_{p}_evalue" for p in NEIGH_PROFILES]

    with open(operon_out, "w") as fh:
        fh.write("\t".join(operon_cols) + "\n")
        for qid in query_ids:
            row = [qid]
            # Catalytic hits on the candidate itself
            for prof in CAND_PROFILES:
                ev = cand_hits[prof].get(qid)
                row.append("YES" if ev is not None else "NO")
            for prof in CAND_PROFILES:
                ev = cand_hits[prof].get(qid)
                row.append(best_evalue_str(ev))
            # Neighbourhood subunit hits
            neigh_data = psra_neigh_summary.get(qid, {})
            for prof in NEIGH_PROFILES:
                ev = neigh_data.get(prof)
                row.append("YES" if ev is not None else "NO")
            for prof in NEIGH_PROFILES:
                ev = neigh_data.get(prof)
                row.append(best_evalue_str(ev))
            fh.write("\t".join(row) + "\n")

    print(f"  Per-candidate operon  → {operon_out}")

    # -------------------------------------------------------------------------
    # Console summary
    # -------------------------------------------------------------------------
    print(f"\n  HMSS2 ANNOTATION SUMMARY ({len(query_ids)} candidates):")
    print(f"  {'Profile':<30} {'Hits'}")
    print(f"  {'-'*30} {'-'*6}")
    for prof in CAND_PROFILES:
        n = sum(1 for qid in query_ids if cand_hits[prof].get(qid) is not None)
        print(f"  {PROFILE_LABELS[prof]:<30} {n}")
    for prof in NEIGH_PROFILES:
        n = sum(1 for qid in query_ids
                if psra_neigh_summary.get(qid, {}).get(prof) is not None)
        print(f"  {PROFILE_LABELS[prof]:<30} {n}  (candidates with ≥1 neighbour hit)")
    print(f"\n  NOTE: HMSS2 results are annotation only — not used in scoring.")
    print(f"        Cross-reference against primary evidence in classification_table.tsv")


if __name__ == "__main__":
    main()
