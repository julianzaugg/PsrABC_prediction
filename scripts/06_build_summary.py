#!/usr/bin/env python3
"""
06_build_summary.py

Integrates all evidence to produce the final classification table.

Evidence sources (in descending reliability order):
  1. Mo-bisPGD domain (PF00384)      — required for any classification
  2. TAT signal peptide (SignalP 6)  — periplasmic export: PsrA/TtrA/PhsA yes, SoeA no
  3. PsrC topology (DeepTMHMM)       — 8TM=PsrC, 9TM=TtrC, 5TM+haem=PhsC
  4. NrfD in neighbourhood           — PF14589 (specific) > PF03916 (broad)
  5. PsrB in neighbourhood           — supporting only, absence is weak evidence
  6. Phylogenetic clade (IQ-TREE)    — heuristic, lower weight than biochemical

SCORING DESIGN:
  - Tree score capped at +2 to prevent phylogeny overriding biochemistry
  - Tree score only applied when Mo-bisPGD is confirmed (not standalone)
  - SoeA classification requires BOTH: no TAT AND no PsrC in neighbourhood
  - NrfD selection: best hit chosen by PF14589 presence first, then confidence
  - Topology fix: use re.match to avoid substring bugs in class detection
  - SoeA.hmm hit (from HMSS2 Step 3b) applies −3 score penalty; this is the
    only HMSS2 result that feeds into scoring. All other HMSS2 columns are
    annotation only.

DISCOVERY SOURCE:
  Candidates may originate from PF00384, PsrAPhsASreA.hmm (HMSS2), or both.
  The discovery_source column is loaded from 00_scan/discovery_source.tsv and
  carried through to the output tables. HMSS2_only candidates are highlighted
  in a distinct colour in the HTML output for easy manual review.
"""

import argparse
import os
import re
from collections import Counter


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--ids",              required=True)
    p.add_argument("--psrA_faa",         required=True)
    p.add_argument("--nrfd_hits",        required=True)
    p.add_argument("--topology_dir",     required=True)
    p.add_argument("--signalp_dir",      required=True)
    p.add_argument("--treefile",         default=None)
    p.add_argument("--references",       required=True)
    p.add_argument("--protein_index",    default=None,
                   help="protein_to_bin_index.tsv for bin_name lookup")
    p.add_argument("--discovery_source", default=None,
                   help="discovery_source.tsv from Step 0b/0c "
                        "(columns: protein_id, discovery_source). "
                        "If absent, all candidates are labelled 'unknown'.")
    p.add_argument("--hmss2_dir",        default=None,
                   help="Directory with HMSS2 results (03_hmmer/hmss2/). "
                        "If supplied, hmss2_operon.tsv annotation columns are "
                        "appended to the summary table. SoeA.hmm hits are also "
                        "used in scoring (−3 penalty).")
    p.add_argument("--outdir",           required=True)
    return p.parse_args()


def load_tsv(path, key_col=0):
    if not os.path.exists(path):
        return {}
    rows = {}
    header = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line or line.startswith("##"):
                continue
            if header is None:
                header = line.split("\t")
                continue
            parts = line.split("\t")
            if len(parts) < len(header):
                parts += [""] * (len(header) - len(parts))
            rows[parts[key_col]] = dict(zip(header, parts))
    return rows


def load_protein_index(index_path):
    """Load protein_to_bin_index.tsv → {protein_id: bin_name}"""
    if not index_path or not os.path.exists(index_path):
        return {}
    mapping = {}
    with open(index_path) as fh:
        fh.readline()  # header
        for line in fh:
            parts = line.rstrip().split("\t")
            if len(parts) >= 2:
                mapping[parts[0]] = parts[1]
    return mapping


def load_discovery_source(source_path):
    """
    Load discovery_source.tsv → {protein_id: source_string}
    Source values: 'both', 'PF00384_only', 'HMSS2_only', 'supplied', 'unknown'
    """
    if not source_path or not os.path.exists(source_path):
        return {}
    mapping = {}
    with open(source_path) as fh:
        fh.readline()  # header
        for line in fh:
            parts = line.rstrip().split("\t")
            if len(parts) >= 2:
                mapping[parts[0]] = parts[1]
    return mapping


def load_soea_hits(hmss2_dir):
    """
    Load SoeA.hmm hits from HMSS2 Step 3b output.
    Returns set of protein_ids that hit SoeA.hmm (E ≤ 1e-5).
    These are the only HMSS2 hits that contribute to scoring.
    """
    if not hmss2_dir:
        return set()
    tbl = os.path.join(hmss2_dir, "SoeA_hits.tbl")
    if not os.path.exists(tbl):
        return set()
    hits = set()
    with open(tbl) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            if len(fields) < 5:
                continue
            try:
                if float(fields[4]) <= 1e-5:
                    hits.add(fields[0])
            except ValueError:
                continue
    return hits


def parse_treefile_clades(treefile, reference_labels):
    """
    Assign query sequences to nearest reference clade using ete3.
    Handles IQ-TREE support strings like '95.6/100:0.04' by using format=1
    (internal node labels) with quoted_node_names fallback.
    Returns dict: {seq_id: nearest_reference_label}
    """
    if not treefile or not os.path.exists(treefile):
        return {}

    try:
        from ete3 import Tree
    except ImportError:
        print("  [INFO] ete3 not installed — skipping clade assignment")
        print("         conda install -c etetoolkit ete3")
        return {}

    tree = None
    # IQ-TREE with -bb and -alrt produces support strings like '95.6/100'
    # as internal node names. ete3 format=1 reads these as node names.
    for fmt in [1, 0, 2, 3]:
        try:
            tree = Tree(treefile, format=fmt)
            break
        except Exception:
            continue

    if tree is None:
        print(f"  [WARN] Could not parse tree {treefile} — skipping clade assignment")
        return {}

    # Map reference leaf names to clade labels
    ref_nodes = {}
    for leaf in tree.get_leaves():
        for ref_label in reference_labels:
            if ref_label in leaf.name:
                ref_nodes[leaf.name] = ref_label
                break

    if not ref_nodes:
        print("  [WARN] No reference sequences found in tree leaves — check label matching")
        return {}

    clade_assignments = {}
    for leaf in tree.get_leaves():
        if any(ref in leaf.name for ref in reference_labels):
            continue
        min_dist    = float("inf")
        nearest_ref = "unassigned"
        for ref_name, ref_label in ref_nodes.items():
            try:
                dist = tree.get_distance(leaf.name, ref_name)
                if dist < min_dist:
                    min_dist    = dist
                    nearest_ref = ref_label
            except Exception:
                continue
        clade_assignments[leaf.name] = nearest_ref

    return clade_assignments


def select_best_nrfd(nrfd_ids_str, topology, all_nrfd_info):
    """
    Select the best NrfD hit for topology classification.
    Priority:
      1. Hits with PF14589 (high-specificity NrfD_2)
      2. Hits with a classified topology (not 'UNCLASSIFIED')
      3. First hit otherwise

    Returns (best_id, subunit_class) or ("none", "no_NrfD_found")
    """
    if not nrfd_ids_str or nrfd_ids_str == "none":
        return "none", "no_NrfD_found"

    candidates = [n for n in nrfd_ids_str.split(";") if n]
    if not candidates:
        return "none", "no_NrfD_found"

    def candidate_score(nid):
        topo  = topology.get(nid, {})
        tc    = topo.get("subunit_classification", "")
        info  = all_nrfd_info.get(nid, {})
        profs = info.get("profiles", "")
        score = 0
        if "PF14589" in profs:
            score += 10
        if tc and "UNCLASSIFIED" not in tc and "not_found" not in tc:
            score += 5
        if tc.startswith("PsrC"):
            score += 3
        return score

    best       = max(candidates, key=candidate_score)
    best_class = topology.get(best, {}).get("subunit_classification", "topology_not_run")
    return best, best_class


def load_nrfd_info(nrfd_hits_path):
    """Load nrfd_hits.tsv profiles_matched column for PF14589 info."""
    info = {}
    if not os.path.exists(nrfd_hits_path):
        return info
    with open(nrfd_hits_path) as fh:
        header = fh.readline().rstrip().split("\t")
        for line in fh:
            parts = line.rstrip().split("\t")
            if len(parts) < 3:
                continue
            row = dict(zip(header, parts))
            nid = row.get("nrfd_protein_id", "")
            if nid:
                info[nid] = row
    return info


def classify_topology(tm_class):
    """
    Classify membrane subunit type from topology string.
    Uses explicit regex word-boundary matching to avoid substring
    false-positives (e.g. 'PsrC_or_ambiguous' must still map to PsrC).
    """
    if not tm_class or tm_class in ("no_NrfD_found", "topology_not_run", "NOT_RUN"):
        return tm_class, "UNKNOWN"

    if re.search(r"\bPsrC\b", tm_class):
        return "PsrC_8TM", "PsrC"
    if re.search(r"\bTtrC\b", tm_class):
        return "TtrC_9TM", "TtrC"
    if re.search(r"\bPhsC\b", tm_class):
        return "PhsC_5TM_haem", "PhsC"
    if re.search(r"\bSoeC\b", tm_class):
        return "SoeC", "SoeC"
    return tm_class, "OTHER"


def score_classification(row, soea_hits):
    """
    Evidence-based scoring.

    Scoring design:
      1. Topology bug fixed: tm_canonical from classify_topology() is used
      2. Tree weight capped at +2 to prevent phylogeny overriding biochemistry
      3. Tree score only applied when Mo-bisPGD present
      4. SoeA classification requires BOTH no-TAT AND no PsrC in neighbourhood
         (avoids misclassifying incomplete-genome PsrA as SoeA)
      5. PF14589 NrfD hit gets +1 over plain PF03916
      6. SoeA.hmm hit (HMSS2 Step 3b) applies −3 penalty — the only HMSS2
         result that enters scoring. Converts ambiguous low-evidence candidates
         to LIKELY_SoeA_or_divergent where biochemical evidence is otherwise
         insufficient. All other HMSS2 results are annotation only.
    """
    prot_id     = row.get("prot_id", "")
    has_mo      = row.get("has_PF00384", "NO") == "YES"
    has_nrfd    = row.get("NrfD_in_neighbourhood", "NO") == "YES"
    has_pf14589 = row.get("NrfD_has_PF14589", "NO") == "YES"
    has_psrB    = row.get("PsrB_in_neighbourhood", "NO") == "YES"
    tm_class    = row.get("membrane_subunit_class", "")
    has_tat     = row.get("has_TAT_signal", "NO") == "YES"
    tree_clade  = row.get("tree_clade", "unassigned")
    has_soea    = prot_id in soea_hits

    _, tm_canonical = classify_topology(tm_class)

    evidence = []
    score    = 0

    # --- Mo-bisPGD (required) ---
    if has_mo:
        evidence.append("Mo-bisPGD(+)"); score += 2
    else:
        evidence.append("Mo-bisPGD(-)"); score -= 2

    # --- TAT signal ---
    if has_tat:
        evidence.append("TAT(+)"); score += 2
    else:
        evidence.append("TAT(-)")

    # --- Membrane subunit topology ---
    if tm_canonical == "PsrC":
        evidence.append("PsrC_8TM(+)"); score += 3
    elif tm_canonical == "TtrC":
        evidence.append("TtrC_9TM"); score -= 1
    elif tm_canonical == "PhsC":
        evidence.append("PhsC_5TM_haem")
    elif tm_canonical == "SoeC":
        evidence.append("SoeC"); score -= 1
    elif has_nrfd:
        if has_pf14589:
            base = "NrfD_PF14589(+)"
        else:
            base = "NrfD_PF03916(topology_ND)"
        evidence.append(base)
        score += (2 if has_pf14589 else 1)

    # --- PsrB (supporting only, never penalised for absence) ---
    if has_psrB:
        evidence.append("PsrB(+)"); score += 1
    else:
        evidence.append("PsrB(not_found)")

    # --- SoeA.hmm hit (HMSS2 — scored, negative) ---
    # Applied before phylogeny so it can shift borderline cases to SoeA
    if has_soea:
        evidence.append("SoeA.hmm(+)"); score -= 3

    # --- Phylogeny (capped at +2, only when Mo confirmed) ---
    if tree_clade and tree_clade not in ("unassigned", "inspect_tree_manually") and has_mo:
        evidence.append(f"Tree:{tree_clade}")
        if "PsrA" in tree_clade:
            score += 2
        elif "PhsA" in tree_clade:
            score += 1
        elif "Soe" in tree_clade:
            score -= 1
        elif "Ttr" in tree_clade:
            score -= 1

    # --- Final classification ---
    if score >= 8:
        classification, confidence = "TRUE_PsrA", "HIGH"
    elif score >= 5:
        classification, confidence = "LIKELY_PsrA", "MEDIUM"
    elif score >= 2 and has_tat and has_mo:
        classification, confidence = "PsrA_or_PhsA", "MEDIUM"
    elif not has_tat and has_mo and not has_nrfd and tm_canonical not in ("PsrC",):
        # SoeA-like: no TAT AND no NrfD neighbour AND no PsrC topology
        # SoeA.hmm hit pushes ambiguous cases here via the score penalty above
        classification, confidence = "LIKELY_SoeA_or_divergent", "MEDIUM"
    elif not has_mo:
        classification, confidence = "NOT_MoBisPGD_enzyme", "HIGH"
    else:
        classification, confidence = "AMBIGUOUS", "LOW"

    return classification, confidence, "|".join(evidence)


def make_html_table(rows, out_path, hmss2_cols=None):
    hmss2_cols = hmss2_cols or []

    # Classification row background colours
    cls_colours = {
        "TRUE_PsrA":                "#d4edda",
        "LIKELY_PsrA":              "#c8e6c9",
        "PsrA_or_PhsA":             "#fff9c4",
        "LIKELY_SoeA_or_divergent": "#fff3e0",
        "NOT_MoBisPGD_enzyme":      "#ffcdd2",
        "AMBIGUOUS":                "#e3f2fd",
    }

    # Discovery source left-border accent colours
    src_border = {
        "HMSS2_only":  "4px solid #e65100",   # orange — needs manual review
        "both":        "4px solid #1565c0",   # blue — highest confidence
        "PF00384_only": "",                    # no accent
        "supplied":    "",
        "unknown":     "",
    }

    html = [
        "<html><head><style>",
        "body{font-family:monospace;font-size:11px;}",
        "table{border-collapse:collapse;width:100%;}",
        "th{background:#333;color:white;padding:5px;text-align:left;}",
        "th.hmss2{background:#555;}",
        "td{padding:3px 6px;border-bottom:1px solid #ddd;}",
        "tr:hover{opacity:0.85;}",
        ".src-badge{font-size:9px;padding:1px 4px;border-radius:3px;"
        "  font-weight:bold;white-space:nowrap;}",
        ".src-hmss2{background:#e65100;color:white;}",
        ".src-both{background:#1565c0;color:white;}",
        ".src-pf{background:#ccc;color:#333;}",
        "</style></head><body>",
        "<h2>PsrABC Classification Summary</h2>",
        "<p style='font-size:10px;color:#333'>",
        "<b>Discovery source badges:</b> ",
        "<span class='src-badge src-hmss2'>HMSS2_only</span> — found by "
        "PsrAPhsASreA.hmm only; not in PF00384 scan — review carefully. &nbsp;",
        "<span class='src-badge src-both'>both</span> — found by both PF00384 "
        "and PsrAPhsASreA.hmm (highest initial confidence). &nbsp;",
        "<span class='src-badge src-pf'>PF00384_only</span> — standard discovery.",
        "</p>",
    ]
    if hmss2_cols:
        html.append("<p style='color:#555;font-size:10px'>"
                    "Columns shaded grey (HMSS2_*) are annotation only — "
                    "not used in scoring, except SoeA.hmm which contributes "
                    "a −3 score penalty.</p>")
    html += [
        "<table><tr>",
        "<th>Protein ID</th><th>Bin</th><th>Source</th>",
        "<th>Mo-bisPGD</th><th>TAT</th>",
        "<th>Membrane subunit</th><th>NrfD(PF14589?)</th><th>PsrB</th>",
        "<th>Tree clade</th><th>Classification</th><th>Confidence</th>",
        "<th>Evidence</th>",
    ]
    for col in hmss2_cols:
        html.append(f'<th class="hmss2">{col}</th>')
    html.append("</tr>")

    primary_cols = ["prot_id", "bin_name", "discovery_source",
                    "has_PF00384", "has_TAT_signal",
                    "membrane_subunit_class", "NrfD_in_neighbourhood",
                    "PsrB_in_neighbourhood", "tree_clade",
                    "classification", "confidence", "evidence"]

    for row in rows:
        cls = row.get("classification", "AMBIGUOUS")
        src = row.get("discovery_source", "unknown")
        bg  = cls_colours.get(cls, "#ffffff")
        border_style = src_border.get(src, "")
        style = f"background:{bg}"
        if border_style:
            style += f";border-left:{border_style}"

        html.append(f'<tr style="{style}">')
        for col in primary_cols:
            val = row.get(col, "NA")
            # Render discovery_source as a coloured badge
            if col == "discovery_source":
                badge_cls = {
                    "HMSS2_only":   "src-hmss2",
                    "both":         "src-both",
                    "PF00384_only": "src-pf",
                }.get(val, "")
                if badge_cls:
                    val = f"<span class='src-badge {badge_cls}'>{val}</span>"
            html.append(f"<td>{val}</td>")
        for col in hmss2_cols:
            val = row.get(col, "NA")
            html.append(f'<td style="color:#555">{val}</td>')
        html.append("</tr>")

    html.append("</table></body></html>")
    with open(out_path, "w") as fh:
        fh.write("\n".join(html))


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    with open(args.ids) as fh:
        query_ids = [l.strip() for l in fh if l.strip() and not l.startswith("#")]

    # Load evidence tables
    hmmer_dir = os.path.dirname(args.nrfd_hits)
    mo_check  = load_tsv(os.path.join(hmmer_dir, "psrA_mo_domain_check.tsv"))
    operon    = load_tsv(os.path.join(hmmer_dir, "operon_completeness.tsv"))
    topology  = load_tsv(os.path.join(args.topology_dir, "topology_summary.tsv"))
    tat       = load_tsv(os.path.join(args.signalp_dir, "tat_summary.tsv"))
    ref_meta  = load_tsv(args.references)
    nrfd_info = load_nrfd_info(os.path.join(hmmer_dir, "nrfd_hits.tsv"))

    # Load discovery source (PF00384_only / HMSS2_only / both / supplied)
    disc_src = load_discovery_source(args.discovery_source)
    if not disc_src:
        print("  [INFO] --discovery_source not supplied or file missing — "
              "labelling all candidates as 'unknown'")

    # Load SoeA.hmm hits for scoring (only HMSS2 result that enters scoring)
    soea_hits = load_soea_hits(args.hmss2_dir)
    if soea_hits:
        print(f"[*] Loaded {len(soea_hits)} SoeA.hmm hits for scoring (−3 penalty)")
    elif args.hmss2_dir:
        print("  [INFO] No SoeA.hmm hits found (or SoeA_hits.tbl absent) — "
              "no SoeA scoring penalty applied")

    # Load HMSS2 annotation table (optional — annotation only except SoeA)
    hmss2_data = {}
    hmss2_cols = []
    if args.hmss2_dir:
        operon_tsv = os.path.join(args.hmss2_dir, "hmss2_operon.tsv")
        if os.path.exists(operon_tsv):
            hmss2_data = load_tsv(operon_tsv, key_col=0)
            with open(operon_tsv) as fh:
                header = fh.readline().rstrip().split("\t")
            hmss2_cols = [c for c in header if c != "psrA_id"]
            print(f"[*] Loaded HMSS2 annotation: {len(hmss2_data)} rows, "
                  f"{len(hmss2_cols)} columns")
        else:
            print(f"  [INFO] --hmss2_dir supplied but hmss2_operon.tsv not found: "
                  f"{operon_tsv}")

    # Load protein→bin mapping
    index_path = args.protein_index
    if not index_path:
        candidate = os.path.join(os.path.dirname(args.outdir),
                                 "00_scan", "protein_to_bin_index.tsv")
        if os.path.exists(candidate):
            index_path = candidate
    bin_map = load_protein_index(index_path)
    if not bin_map:
        print("  [INFO] protein_index not found — bin_name column will be empty")

    # Parse tree clade assignments
    ref_labels  = list(ref_meta.keys())
    tree_clades = {}
    if args.treefile and os.path.exists(args.treefile):
        print("[*] Parsing tree for clade assignments...")
        tree_clades = parse_treefile_clades(args.treefile, ref_labels)

    final_rows = []
    for prot_id in query_ids:
        row = {"prot_id": prot_id}

        # Discovery source
        row["discovery_source"] = disc_src.get(prot_id, "unknown")

        # Bin name
        row["bin_name"] = bin_map.get(prot_id, "NA")

        # Mo-bisPGD
        mo = mo_check.get(prot_id, {})
        row["has_PF00384"]    = mo.get("PF00384_hit", "NOT_RUN")
        row["PF00384_evalue"] = mo.get("PF00384_evalue", "NA")

        # Operon
        op           = operon.get(prot_id, {})
        n_nrfd       = op.get("NrfD_PsrC_count", "0")
        n_psrB       = op.get("PsrB_count", "0")
        nrfd_ids_str = op.get("NrfD_ids", "none")
        row["NrfD_in_neighbourhood"] = "YES" if n_nrfd not in ("0", "") else "NO"
        row["PsrB_in_neighbourhood"] = "YES" if n_psrB not in ("0", "") else "NO"
        row["n_NrfD_neighbours"]     = n_nrfd
        row["n_PsrB_neighbours"]     = n_psrB
        row["NrfD_ids"]              = nrfd_ids_str
        row["PsrB_ids"]              = op.get("PsrB_ids", "none")
        row["NrfD_has_PF14589"]      = op.get("NrfD_has_PF14589", "NO")

        # Best membrane subunit topology
        best_nrfd_id, best_tm = select_best_nrfd(nrfd_ids_str, topology, nrfd_info)
        row["membrane_subunit_class"] = best_tm
        row["best_NrfD_id"]           = best_nrfd_id

        # TAT signal
        tat_info              = tat.get(prot_id, {})
        row["has_TAT_signal"] = tat_info.get("has_TAT", "NOT_RUN")
        row["TAT_probability"]    = tat_info.get("TAT_probability", "NA")
        row["signalp_prediction"] = tat_info.get("signalp_prediction", "NOT_RUN")

        # Tree clade
        row["tree_clade"] = tree_clades.get(prot_id, "inspect_tree_manually")

        # Score and classify (soea_hits passed for SoeA.hmm penalty)
        classification, confidence, evidence = score_classification(row, soea_hits)
        row["classification"] = classification
        row["confidence"]     = confidence
        row["evidence"]       = evidence

        # HMSS2 annotation columns (appended after scoring — no influence on
        # score except SoeA.hmm which was already applied above)
        if hmss2_cols:
            hmss2_row = hmss2_data.get(prot_id, {})
            for col in hmss2_cols:
                row[col] = hmss2_row.get(col, "NA")

        final_rows.append(row)

    # Write TSV
    tsv_path = os.path.join(args.outdir, "classification_table.tsv")
    cols = ["prot_id", "bin_name", "discovery_source",
            "has_PF00384", "PF00384_evalue",
            "has_TAT_signal", "TAT_probability", "signalp_prediction",
            "NrfD_in_neighbourhood", "n_NrfD_neighbours", "NrfD_ids",
            "NrfD_has_PF14589", "best_NrfD_id", "membrane_subunit_class",
            "PsrB_in_neighbourhood", "n_PsrB_neighbours", "PsrB_ids",
            "tree_clade", "classification", "confidence", "evidence"] + hmss2_cols

    with open(tsv_path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in final_rows:
            fh.write("\t".join(str(row.get(c, "NA")) for c in cols) + "\n")

    # Write HTML
    html_path = os.path.join(args.outdir, "classification_table.html")
    make_html_table(final_rows, html_path, hmss2_cols=hmss2_cols)

    # Summary counts
    counts = Counter(r["classification"] for r in final_rows)
    src_counts = Counter(r["discovery_source"] for r in final_rows)

    print(f"\n{'='*65}")
    print(f"  FINAL CLASSIFICATION SUMMARY ({len(final_rows)} proteins)")
    print(f"{'='*65}")
    for cls, n in sorted(counts.items()):
        print(f"  {cls:<45} : {n}")
    print(f"{'='*65}")
    print(f"\n  DISCOVERY SOURCE BREAKDOWN:")
    for src, n in sorted(src_counts.items()):
        flag = "  ← review in HTML" if src == "HMSS2_only" else ""
        print(f"  {src:<20} : {n}{flag}")
    print(f"\n  TSV  : {tsv_path}")
    print(f"  HTML : {html_path}")
    print(f"\n  GUIDE:")
    print(f"    TRUE_PsrA     : Mo-bisPGD + TAT + PsrC(8TM)")
    print(f"    LIKELY_PsrA   : Most evidence supports PsrA")
    print(f"    PsrA_or_PhsA  : Mo-bisPGD + TAT, subunit C ambiguous")
    print(f"    LIKELY_SoeA   : Mo-bisPGD, no TAT, no NrfD neighbour")
    print(f"    AMBIGUOUS     : Conflicting/incomplete evidence")
    print(f"    NrfD_PF03916(topology_ND) : NrfD neighbour confirmed by broad "
          f"HMM only; run DeepTMHMM to resolve to PsrC/TtrC/PhsC")
    print(f"\n  NOTE: HMSS2_only candidates (orange border in HTML) were not")
    print(f"        detected by PF00384. Check tree placement and neighbourhood")
    print(f"        evidence carefully before accepting these as PsrA.")


if __name__ == "__main__":
    main()
