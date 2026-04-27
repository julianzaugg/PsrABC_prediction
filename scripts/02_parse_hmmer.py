#!/usr/bin/env python3
"""
02_parse_hmmer.py

Parse HMMER tblout files for:
  - PF00384 hits against PsrA candidates  (confirms Mo-bisPGD domain)
  - PF03916 hits in neighbourhood FAA     (identifies NrfD/PsrC-like subunits)

Annotates which neighbourhood proteins hit PF03916 and links them back
to the psrA candidate they were found near.
Extracts NrfD sequences to a combined FAA for topology prediction.
"""

import argparse
import os
import re
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--psrA_hits",       required=True, help="hmmsearch tblout: PF00384 vs candidates")
    p.add_argument("--nrfd_hits",       required=True, help="hmmsearch tblout: PF03916 vs neighbours (PsrC)")
    p.add_argument("--nrfd2_hits",      default=None,  help="hmmsearch tblout: PF14589 vs neighbours (NrfD_2, more specific PsrC)")
    p.add_argument("--psrB_hits",       default=None,  help="hmmsearch tblout: PF12800 vs neighbours (PsrB, NrfC-like)")
    p.add_argument("--psrB_hits2",      default=None,  help="hmmsearch tblout: PF13247 vs neighbours (PsrB, 4Fe-4S dicluster)")
    p.add_argument("--neighbourhood_faa", required=True, help="Combined neighbourhood FAA")
    p.add_argument("--psrA_faa",        required=True, help="Candidate PsrA FAA")
    p.add_argument("--ids",             required=True, help="Original query IDs file")
    p.add_argument("--outdir",          required=True, help="Output directory")
    return p.parse_args()


def parse_tblout(tblout_path, evalue_cutoff=1e-5):
    """
    Returns dict: {target_name: best_evalue}
    from an hmmsearch --tblout file.
    """
    hits = {}
    with open(tblout_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 5:
                continue
            target = fields[0]
            evalue = float(fields[4])  # full sequence E-value
            if evalue <= evalue_cutoff:
                if target not in hits or evalue < hits[target]:
                    hits[target] = evalue
    return hits


def load_faa(faa_path):
    """Returns dict: protein_id → (header_line, sequence)"""
    seqs = {}
    current_id = None
    current_header = None
    current_seq = []
    with open(faa_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id:
                    seqs[current_id] = (current_header, "".join(current_seq))
                current_header = line
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id:
        seqs[current_id] = (current_header, "".join(current_seq))
    return seqs


def write_faa(seqs_dict, path, id_set=None):
    """Write sequences to FASTA. If id_set given, only write those IDs."""
    with open(path, "w") as fh:
        for seq_id, (header, seq) in seqs_dict.items():
            if id_set and seq_id not in id_set:
                continue
            fh.write(header + "\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + "\n")


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Load original query IDs
    with open(args.ids) as fh:
        query_ids = set(l.strip() for l in fh if l.strip() and not l.startswith("#"))

    # --- Parse PF00384 hits (Mo-bisPGD confirmation) ---
    print("[*] Parsing PF00384 hits (Mo-bisPGD domain confirmation)...")
    mo_hits = parse_tblout(args.psrA_hits)
    print(f"    {len(mo_hits)} candidate PsrA sequences with PF00384 hit (E≤1e-5)")

    confirmed_psrA = set()
    for prot_id in query_ids:
        if prot_id in mo_hits:
            confirmed_psrA.add(prot_id)
        else:
            print(f"    [WARN] {prot_id}: NO PF00384 hit — may not be a true Mo-bisPGD enzyme")

    # --- Parse PF03916 hits (NrfD-like membrane subunit / PsrC, broad) ---
    print("\n[*] Parsing PF03916 hits (NrfD-like / PsrC-type, broad)...")
    nrfd_hits = parse_tblout(args.nrfd_hits)
    print(f"    {len(nrfd_hits)} neighbourhood proteins with PF03916 hit (E≤1e-5)")

    # --- Parse PF14589 hits (NrfD_2 — more specific for polysulfide reductase PsrC) ---
    # PF14589 is a narrower profile than PF03916 and more diagnostic for true PsrC.
    # Proteins hitting PF14589 are flagged as high-specificity PsrC candidates.
    nrfd2_hits = {}
    if args.nrfd2_hits and os.path.exists(args.nrfd2_hits) and os.path.getsize(args.nrfd2_hits) > 0:
        nrfd2_hits = parse_tblout(args.nrfd2_hits)
        print(f"    {len(nrfd2_hits)} neighbourhood proteins with PF14589 hit (E≤1e-5) — HIGH specificity NrfD_2")
    else:
        print(f"    [INFO] No PF14589 results — skipping high-specificity NrfD_2 search")

    # Merge: any protein hitting either PF03916 or PF14589 is a NrfD candidate
    # Track which profiles each hit matched for downstream reporting
    all_nrfd_hits = {}   # protein_id → {"evalue": float, "profiles": set}
    for pid, ev in nrfd_hits.items():
        all_nrfd_hits.setdefault(pid, {"evalue": ev, "profiles": set()})
        all_nrfd_hits[pid]["profiles"].add("PF03916")
        all_nrfd_hits[pid]["evalue"] = min(all_nrfd_hits[pid]["evalue"], ev)
    for pid, ev in nrfd2_hits.items():
        all_nrfd_hits.setdefault(pid, {"evalue": ev, "profiles": set()})
        all_nrfd_hits[pid]["profiles"].add("PF14589")
        all_nrfd_hits[pid]["evalue"] = min(all_nrfd_hits[pid]["evalue"], ev)
    print(f"    {len(all_nrfd_hits)} unique NrfD candidates (PF03916 + PF14589 combined)")

    # --- Parse PsrB hits (PF12800 and/or PF13247) ---
    # PsrB = electron transfer subunit, 4x[4Fe-4S] clusters
    # NOTE: PsrB is useful supporting evidence but NOT used to reject a candidate
    # if absent — 4Fe-4S proteins are easily missed by annotation and HMMs.
    # PF12800 (NrfC-like) is the more specific profile for PsrB.
    # PF13247 (4Fe-4S dicluster) is broader but catches more divergent PsrB.
    psrB_hits = {}
    for tbl_arg, label in [(args.psrB_hits, "PF12800"), (args.psrB_hits2, "PF13247")]:
        if tbl_arg and os.path.exists(tbl_arg) and os.path.getsize(tbl_arg) > 0:
            hits = parse_tblout(tbl_arg)
            print(f"    {len(hits)} neighbourhood proteins with {label} hit (E≤1e-5)")
            for h, ev in hits.items():
                if h not in psrB_hits or ev < psrB_hits[h]:
                    psrB_hits[h] = ev
        else:
            print(f"    [INFO] No {label} results file — skipping PsrB search via {label}")
    print(f"    {len(psrB_hits)} unique PsrB candidates (combined PF12800+PF13247)")

    # Load neighbourhood FAA
    neighbour_seqs = load_faa(args.neighbourhood_faa)

    # For each NrfD hit, trace back to the psrA it was found near
    nrfd_to_psrA = {}
    for seq_id, (header, seq) in neighbour_seqs.items():
        if seq_id not in all_nrfd_hits:
            continue
        m = re.search(r"from:(\S+)", header)
        psrA_id = m.group(1) if m else "unknown"
        nrfd_to_psrA[seq_id] = psrA_id

    # For each PsrB hit, trace back to the psrA it was found near
    psrB_to_psrA = {}
    for seq_id, (header, seq) in neighbour_seqs.items():
        if seq_id not in psrB_hits:
            continue
        m = re.search(r"from:(\S+)", header)
        psrA_id = m.group(1) if m else "unknown"
        psrB_to_psrA[seq_id] = psrA_id

    # Write NrfD candidate sequences to FAA for topology prediction
    nrfd_faa_path = os.path.join(args.outdir, "nrfd_candidates.faa")
    write_faa(neighbour_seqs, nrfd_faa_path, id_set=set(nrfd_to_psrA.keys()))
    print(f"\n    NrfD/PsrC candidates → {nrfd_faa_path}  (send to DeepTMHMM)")

    # Write PsrB candidates to FAA
    psrB_faa_path = os.path.join(args.outdir, "psrB_candidates.faa")
    write_faa(neighbour_seqs, psrB_faa_path, id_set=set(psrB_to_psrA.keys()))
    print(f"    PsrB candidates      → {psrB_faa_path}")

    # --- Write summary tables ---

    # PsrA confirmation table
    psrA_table = os.path.join(args.outdir, "psrA_mo_domain_check.tsv")
    with open(psrA_table, "w") as fh:
        fh.write("prot_id\tPF00384_hit\tPF00384_evalue\n")
        for prot_id in sorted(query_ids):
            hit = prot_id in mo_hits
            ev = f"{mo_hits[prot_id]:.2e}" if hit else "NA"
            fh.write(f"{prot_id}\t{'YES' if hit else 'NO'}\t{ev}\n")
    print(f"\n    PsrA Mo-domain check → {psrA_table}")

    # NrfD (PsrC) hits table — includes which profiles matched and PF14589 flag
    nrfd_table = os.path.join(args.outdir, "nrfd_hits.tsv")
    with open(nrfd_table, "w") as fh:
        fh.write("nrfd_protein_id\tbest_evalue\tprofiles_matched\tPF14589_high_specificity\tnearest_psrA\n")
        for nrfd_id in sorted(nrfd_to_psrA.keys()):
            info     = all_nrfd_hits.get(nrfd_id, {"evalue": float("nan"), "profiles": set()})
            ev       = f"{info['evalue']:.2e}"
            profiles = ";".join(sorted(info["profiles"]))
            pf14589  = "YES" if "PF14589" in info["profiles"] else "NO"
            fh.write(f"{nrfd_id}\t{ev}\t{profiles}\t{pf14589}\t{nrfd_to_psrA[nrfd_id]}\n")
    print(f"    NrfD/PsrC hits table  → {nrfd_table}")

    # PsrB hits table
    psrB_table = os.path.join(args.outdir, "psrB_hits.tsv")
    with open(psrB_table, "w") as fh:
        fh.write("psrB_protein_id\tbest_evalue\tnearest_psrA\n")
        for psrB_id in sorted(psrB_to_psrA.keys()):
            ev = psrB_hits.get(psrB_id, float("nan"))
            fh.write(f"{psrB_id}\t{ev:.2e}\t{psrB_to_psrA[psrB_id]}\n")
    print(f"    PsrB hits table       → {psrB_table}")

    # Operon completeness summary
    psrA_has_nrfd  = {p: [] for p in query_ids}
    psrA_has_psrB  = {p: [] for p in query_ids}
    for nrfd_id, psrA_id in nrfd_to_psrA.items():
        if psrA_id in psrA_has_nrfd:
            psrA_has_nrfd[psrA_id].append(nrfd_id)
    for psrB_id, psrA_id in psrB_to_psrA.items():
        if psrA_id in psrA_has_psrB:
            psrA_has_psrB[psrA_id].append(psrB_id)

    operon_table = os.path.join(args.outdir, "operon_completeness.tsv")
    with open(operon_table, "w") as fh:
        fh.write("psrA_id\thas_PF00384\tNrfD_PsrC_count\tNrfD_ids\tNrfD_has_PF14589\t"
                 "PsrB_count\tPsrB_ids\toperon_note\n")
        for prot_id in sorted(query_ids):
            has_mo  = prot_id in mo_hits
            nrfds   = psrA_has_nrfd.get(prot_id, [])
            psrBs   = psrA_has_psrB.get(prot_id, [])
            n_nrfd  = len(nrfds)
            n_psrB  = len(psrBs)
            nrfd_str = ";".join(nrfds) if nrfds else "none"
            psrB_str = ";".join(psrBs) if psrBs else "none"
            # Flag if any NrfD hit for this psrA has PF14589 (high-specificity)
            has_pf14589 = any(
                "PF14589" in all_nrfd_hits.get(n, {}).get("profiles", set())
                for n in nrfds
            )
            pf14589_str = "YES" if has_pf14589 else "NO"

            if has_mo and n_nrfd >= 1 and n_psrB >= 1:
                note = "ABC_complete"
            elif has_mo and n_nrfd >= 1 and n_psrB == 0:
                note = "AC_only—PsrB_not_found"
            elif has_mo and n_nrfd == 0 and n_psrB >= 1:
                note = "AB_only—PsrC_not_found"
            elif has_mo and n_nrfd == 0 and n_psrB == 0:
                note = "A_only—PsrBC_missing"
            elif not has_mo:
                note = "WARNING_no_Mo_domain"
            else:
                note = "PARTIAL"
            fh.write(f"{prot_id}\t{'YES' if has_mo else 'NO'}\t{n_nrfd}\t{nrfd_str}\t"
                     f"{pf14589_str}\t{n_psrB}\t{psrB_str}\t{note}\n")
    print(f"    Operon completeness   → {operon_table}")
    print(f"\n    NOTE: PsrB absence is LOW concern (4Fe-4S proteins frequently missed)")
    print(f"          PF14589=YES in NrfD hit = high-specificity PsrC evidence")
    print(f"\n[DONE] HMMER parsing complete")


if __name__ == "__main__":
    main()
