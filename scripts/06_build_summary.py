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
  - Expanded reference-tree placement gates high-level family labels
  - Tree score only applied when Mo-bisPGD is confirmed (not standalone)
  - SoeA classification requires BOTH: no TAT AND no PsrC in neighbourhood
  - NrfD selection: best hit chosen by PF14589 presence first, then confidence
  - Topology fix: use re.match to avoid substring bugs in class detection
  - SoeA.hmm hit (from HMSS2 Step 3b) applies conditional −3/−1/0 penalty;
    this is the only HMSS2 result that feeds into scoring. All other HMSS2
    columns are annotation only.
  - ArrA/ArxA tree placement with operon evidence (TAT + NrfD or PsrB) triggers
    a hard LIKELY_ArrA override, bypassing score-based classification. The ArrA
    operon architecture (ArrABC) is analogous to PsrABC — TAT + NrfD/PsrB does
    NOT discriminate between the two. Without operon evidence, ArrA placement
    applies only the −2 score penalty and falls through normal scoring.

DISCOVERY SOURCE:
  Candidates may originate from PF00384, PsrAPhsASreA.hmm (HMSS2), or both.
  The discovery_source column is loaded from 00_scan/discovery_source.tsv and
  normalised to clearer gate labels (PF00384_gate, HMSS2_gate, both_gates).
  This avoids misleading combinations such as HMSS2_gate + has_PF00384=YES,
  where the candidate entered through HMSS2 but was later confirmed by PF00384.
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
    Source values: 'both', 'PF00384_gate', 'HMSS2_gate', 'supplied', 'unknown'
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




def normalise_discovery_source(src):
    """Return clearer discovery-gate labels for output tables.

    The original pipeline labels describe how a protein entered the candidate
    set. A row can therefore be HMSS2_gate but still later have PF00384=YES,
    which is technically correct but visually confusing. These labels make the
    column explicitly about the entry gate.
    """
    s = str(src or "unknown").strip()
    mapping = {
        "PF00384_only": "PF00384_gate",
        "PF00384_gate": "PF00384_gate",
        "PF00384": "PF00384_gate",
        "HMSS2_only": "HMSS2_gate",
        "HMSS2_gate": "HMSS2_gate",
        "HMSS2": "HMSS2_gate",
        "both": "both_gates",
        "PF00384_and_HMSS2": "both_gates",
        "PF00384+HMSS2": "both_gates",
        "supplied": "supplied",
        "unknown": "unknown",
        "": "unknown",
    }
    return mapping.get(s, s)


def clean_tat_probability(has_tat, probability, prediction):
    """Mask misleading zero probabilities for categorical TATLIPO calls.

    SignalP can report TATLIPO as a categorical prediction even when the parsed
    TAT probability field is 0.0 or missing/not comparable. For heatmaps and
    Excel summaries, treating those as true zero-probability TAT calls is
    misleading. Keep has_TAT_signal=YES but report TAT_probability as NA.
    """
    pred = str(prediction or "").upper()
    prob = str(probability or "NA")
    if str(has_tat).upper() == "YES" and pred == "TATLIPO":
        try:
            if float(prob) == 0.0:
                return "NA"
        except Exception:
            return "NA"
    return probability


def derive_context_flags(row):
    """Derive simple interpretability flags from final evidence fields."""
    has_mo      = row.get("has_PF00384", "NO") == "YES"
    has_tat     = row.get("has_TAT_signal", "NO") == "YES"
    has_nrfd    = row.get("NrfD_in_neighbourhood", "NO") == "YES"
    has_pf14589 = row.get("NrfD_has_PF14589", "NO") == "YES"
    has_psrB    = row.get("PsrB_in_neighbourhood", "NO") == "YES"
    tm_class    = row.get("membrane_subunit_class", "")
    tree_group  = tree_group_from_clade(row.get("tree_clade", "unassigned"))
    _, tm_canonical = classify_topology(tm_class)

    # Strict PsrABC-like flag: canonical Psr-like operon evidence including TAT.
    # This is used for the existing Psr/conflict annotation.
    psrABC_like = (
        has_mo and has_tat and has_psrB and
        (tm_canonical == "PsrC" or has_pf14589 or has_nrfd)
    )

    # Broader A+B+C-like local context, independent of TAT. This captures
    # biologically interesting non-Psr MopB/DMSOR placements that nevertheless
    # retain an adjacent B-like Fe-S subunit and C-like/NrfD membrane subunit.
    abc_like_nonpsr = has_mo and has_psrB and has_nrfd
    psrC_like_context = tm_canonical == "PsrC" or has_pf14589

    tree_incompat = psrABC_like and tree_group in {
        "arr_arx", "soe_sre", "ttr_srd", "known_nonpsr_mopb"
    }
    return {
        "tree_group": tree_group,
        "PsrABC_like_operon": "YES" if psrABC_like else "NO",
        "ABC_like_nonPsr_operon": "YES" if (abc_like_nonpsr and tree_group == "known_nonpsr_mopb") else "NO",
        "PsrC_like_context": "YES" if psrC_like_context else "NO",
        "tree_incompatible_with_Psr": "YES" if tree_incompat else "NO",
    }

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


def classify_tree_assignment_confidence(nearest_dist, second_dist):
    """Classify nearest-reference assignment confidence from tree distances.

    This is diagnostic only. It does not alter scoring/classification.
    Distances are tree branch-length distances from a candidate to the nearest
    labelled reference family and the second-nearest labelled reference family.
    """
    try:
        nearest = float(nearest_dist)
        second = float(second_dist)
    except Exception:
        return "NA"

    if second <= 0 or nearest < 0:
        return "NA"

    margin = second - nearest
    ratio = nearest / second

    # Conservative, intentionally heuristic thresholds. These should be used to
    # flag rows for manual inspection, not as hard taxonomic boundaries.
    if margin >= 0.20 or ratio <= 0.75:
        return "HIGH"
    if margin >= 0.05 or ratio <= 0.90:
        return "MEDIUM"
    return "LOW"


def parse_treefile_clades(treefile, reference_metadata):
    """
    Assign query sequences to the clade/family of their nearest reference leaf.

    reference_metadata is the dict loaded from reference_metadata.tsv, keyed by
    the reference FASTA label. Built-in references use labels such as:
        PsrA_Wolinella_succinogenes
    and FASTA leaves such as:
        PsrA_Wolinella_succinogenes__sp|P31075|PSRA_WOLSU

    Custom references should use the same family-first style:
        PsrA_My_reference__accession
        PsrAPhsASrrA_Wells_01__GCA_...
        bSreASoeA_Wells_01__GCA_...

    Returns dict keyed by query sequence ID. Each value contains nearest and
    second-nearest reference-family diagnostics:
        tree_clade, tree_distance, tree_second_clade, tree_second_distance,
        tree_distance_margin, tree_distance_ratio, tree_assignment_confidence

    IMPORTANT: this is still a nearest-reference heuristic, not proof of robust
    monophyletic clade membership. The second-nearest diagnostics help identify
    cases where a candidate is close to two labelled families.
    """
    empty = {}
    if not treefile or not os.path.exists(treefile):
        return empty

    try:
        from ete3 import Tree
    except ImportError:
        print("  [INFO] ete3 not installed — skipping clade assignment")
        print("         conda install -c etetoolkit ete3")
        return empty

    tree = None
    for fmt in [1, 0, 2, 3]:
        try:
            tree = Tree(treefile, format=fmt)
            break
        except Exception:
            continue

    if tree is None:
        print(f"  [WARN] Could not parse tree {treefile} — skipping clade assignment")
        return empty

    # Longest labels first avoids accidental partial matches when labels share
    # prefixes, e.g. PsrA_* and PsrAPhsASrrA_*.
    reference_labels = sorted(reference_metadata.keys(), key=len, reverse=True)

    ref_nodes = {}
    for leaf in tree.get_leaves():
        for ref_label in reference_labels:
            if leaf.name == ref_label or leaf.name.startswith(ref_label + "__") or ref_label in leaf.name:
                meta = reference_metadata.get(ref_label, {})
                ref_clade = meta.get("clade") or meta.get("protein") or ref_label
                ref_nodes[leaf.name] = ref_clade
                break

    if not ref_nodes:
        print("  [WARN] No reference sequences found in tree leaves — check label matching")
        return empty

    diagnostics = {}
    for leaf in tree.get_leaves():
        if leaf.name in ref_nodes:
            continue

        # For each reference family/clade, retain only the closest labelled
        # reference leaf. Then compare the closest and second-closest families.
        best_by_clade = {}
        for ref_name, ref_clade in ref_nodes.items():
            try:
                dist = float(tree.get_distance(leaf.name, ref_name))
            except Exception:
                continue
            if ref_clade not in best_by_clade or dist < best_by_clade[ref_clade]:
                best_by_clade[ref_clade] = dist

        if not best_by_clade:
            diagnostics[leaf.name] = {
                "tree_clade": "unassigned",
                "tree_distance": "NA",
                "tree_second_clade": "NA",
                "tree_second_distance": "NA",
                "tree_distance_margin": "NA",
                "tree_distance_ratio": "NA",
                "tree_assignment_confidence": "NA",
            }
            continue

        ordered = sorted(best_by_clade.items(), key=lambda kv: kv[1])
        nearest_clade, nearest_dist = ordered[0]
        if len(ordered) > 1:
            second_clade, second_dist = ordered[1]
            margin = second_dist - nearest_dist
            ratio = nearest_dist / second_dist if second_dist > 0 else float("nan")
            conf = classify_tree_assignment_confidence(nearest_dist, second_dist)
        else:
            second_clade, second_dist = "NA", "NA"
            margin, ratio, conf = "NA", "NA", "NA"

        def fmt(x):
            if isinstance(x, str):
                return x
            try:
                return f"{float(x):.6g}"
            except Exception:
                return "NA"

        diagnostics[leaf.name] = {
            "tree_clade": nearest_clade,
            "tree_distance": fmt(nearest_dist),
            "tree_second_clade": second_clade,
            "tree_second_distance": fmt(second_dist),
            "tree_distance_margin": fmt(margin),
            "tree_distance_ratio": fmt(ratio),
            "tree_assignment_confidence": conf,
        }

    return diagnostics


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


def tree_group_from_clade(tree_clade):
    """
    Collapse tree_clade labels into broad groups used by the classifier.

    There is still only one tree. This helper asks whether the candidate's
    nearest labelled reference clade is compatible with a Psr/Phs/Srr,
    Arr/Arx, Soe/Sre, Ttr/Srd, or known non-Psr MopB interpretation.
    """
    c = str(tree_clade or "").strip()
    if c.lower() in {"", "na", "nan", "none", "unassigned", "inspect_tree_manually"}:
        return "unresolved"

    # Exact labels from the expanded Wells-derived reference set and built-ins.
    psr_like = {"PsrAPhsASrrA", "PsrA", "PhsA", "SrrA"}
    arr_like = {"ArrAArxA", "ArrA", "ArxA"}
    soe_like = {"bSreASoeA", "aSreA", "SoeA", "SreA", "Sre"}
    ttr_like = {"TtrASrdA", "TtrA", "SrdA"}
    known_nonpsr = {
        "ActB", "AH", "AioA", "AioAIdrA", "IdrA", "AspDMSOR",
        "BisC", "DmsA", "DorATorA", "DorA", "TorA", "FdhG",
        "FdhsFsdA", "Fdhs", "FsdA", "FhcB", "FwdBFmdB", "FwdB", "FmdB",
        "NapA", "NarG", "NasCNarB", "NasC", "NarB", "Nqo3", "RhLPgtL",
        "PcrA",
    }

    if c in psr_like:
        return "psr_phs_srr"
    if c in arr_like:
        return "arr_arx"
    if c in soe_like:
        return "soe_sre"
    if c in ttr_like:
        return "ttr_srd"
    if c in known_nonpsr:
        return "known_nonpsr_mopb"

    # Conservative pattern fallbacks for user-provided custom labels.
    if re.search(r"(^|[_/|-])(PsrA|PhsA|SrrA)([_/|-]|$)", c):
        return "psr_phs_srr"
    if re.search(r"(^|[_/|-])(ArrA|ArxA)([_/|-]|$)", c):
        return "arr_arx"
    if re.search(r"(^|[_/|-])(SoeA|SreA|bSreA|aSreA)([_/|-]|$)", c):
        return "soe_sre"
    if re.search(r"(^|[_/|-])(TtrA|SrdA)([_/|-]|$)", c):
        return "ttr_srd"
    if re.search(r"(^|[_/|-])(ActB|AioA|IdrA|BisC|DmsA|DorA|TorA|FdhG|FdhH|FsdA|FhcB|FwdB|FmdB|NapA|NarG|NasC|NarB|Nqo3|PcrA)([_/|-]|$)", c):
        return "known_nonpsr_mopb"

    return "other"


def score_classification(row, soea_hits):
    """
    Tree-gated evidence classifier.

    The local operon/domain score is still calculated and reported through the
    evidence string, but high-level PsrA/PhsA/SrrA calls are constrained by the
    candidate's placement in the expanded reference tree.

    In practice, a candidate with Mo-bisPGD + TAT + nearby B/C-like proteins can
    look operon-compatible with PsrABC, ArrABC, or other DMSOR/MopB systems.
    Therefore, tree placement is used as a gate:
      * Psr/Phs/Srr tree placement permits TRUE_PsrA / LIKELY_PsrA.
      * Arr/Arx tree placement yields LIKELY_ArrA when operon evidence exists.
      * Soe/Sre tree placement yields LIKELY_SoeA_or_divergent.
      * Ttr/Srd tree placement yields LIKELY_TtrA_or_SrdA.
      * Known non-Psr MopB placements yield LIKELY_nonPsr_MopB.
      * Unresolved/other placements are capped at PsrA_or_PhsA or AMBIGUOUS.
    """
    prot_id     = row.get("prot_id", "")
    has_mo      = row.get("has_PF00384", "NO") == "YES"
    has_nrfd    = row.get("NrfD_in_neighbourhood", "NO") == "YES"
    has_pf14589 = row.get("NrfD_has_PF14589", "NO") == "YES"
    has_psrB    = row.get("PsrB_in_neighbourhood", "NO") == "YES"
    tm_class    = row.get("membrane_subunit_class", "")
    has_tat     = row.get("has_TAT_signal", "NO") == "YES"
    tree_clade  = row.get("tree_clade", "unassigned")
    tree_group  = tree_group_from_clade(tree_clade)
    has_soea    = prot_id in soea_hits

    _, tm_canonical = classify_topology(tm_class)

    evidence = []
    score    = 0

    # --- Mo-bisPGD domain ---
    if has_mo:
        evidence.append("Mo-bisPGD(+)"); score += 2
    else:
        evidence.append("Mo-bisPGD(-)"); score -= 2

    # --- TAT signal ---
    if has_tat:
        evidence.append("TAT(+)"); score += 2
    else:
        evidence.append("TAT(-)")

    # --- Membrane subunit topology / NrfD-family neighbour ---
    if tm_canonical == "PsrC":
        evidence.append("PsrC_8TM(+)"); score += 3
    elif tm_canonical == "TtrC":
        evidence.append("TtrC_9TM"); score -= 1
    elif tm_canonical == "PhsC":
        evidence.append("PhsC_5TM_haem"); score += 1
    elif tm_canonical == "SoeC":
        evidence.append("SoeC"); score -= 1
    elif has_nrfd:
        if has_pf14589:
            evidence.append("NrfD_PF14589(+)"); score += 2
        else:
            evidence.append("NrfD_PF03916(topology_ND)"); score += 1

    # --- PsrB-like/electron-transfer neighbour ---
    if has_psrB:
        evidence.append("PsrB(+)"); score += 1
    else:
        evidence.append("PsrB(not_found)")

    # --- HMSS2 SoeA evidence. Keep as a score modifier, but do not allow it
    #     to override the tree gate. ---
    if has_soea:
        if tm_canonical == "PsrC" or has_pf14589:
            evidence.append("SoeA.hmm(+,suppressed_by_PsrC)")
        elif has_nrfd and has_tat:
            evidence.append("SoeA.hmm(+)"); score -= 1
        else:
            evidence.append("SoeA.hmm(+)"); score -= 3

    # --- Tree placement. This is recorded as evidence and used below as a gate. ---
    if tree_clade and tree_clade not in ("unassigned", "inspect_tree_manually"):
        evidence.append(f"Tree:{tree_clade}")
        evidence.append(f"Tree_gate:{tree_group}")
        if has_mo:
            if tree_group == "psr_phs_srr":
                score += 2
            elif tree_group in {"arr_arx", "ttr_srd", "soe_sre", "known_nonpsr_mopb"}:
                # Penalise incompatible placements in the Psr-operon score, but
                # final class is handled by the tree gate below.
                score -= 2
    else:
        evidence.append("Tree:unresolved")
        evidence.append("Tree_gate:unresolved")

    has_operon_evidence = has_tat and (has_nrfd or has_psrB)
    strong_psr_operon   = has_tat and has_mo and (tm_canonical == "PsrC" or has_pf14589 or has_nrfd) and has_psrB
    moderate_psr_operon = has_tat and has_mo and (has_nrfd or has_psrB or tm_canonical in {"PsrC", "PhsC"})

    # Explicit conflict/context flags.
    # PsrABC_like_operon is strict and includes TAT, preserving the previous
    # Psr-oriented interpretation flag. ABC_like_nonPsr_operon is broader and
    # captures A+B+C-like operon context among candidates whose A subunit places
    # with known non-Psr MopB/DMSOR families.
    psrABC_like_operon = strong_psr_operon
    abc_like_nonpsr_operon = has_mo and has_nrfd and has_psrB and tree_group == "known_nonpsr_mopb"
    psrC_like_context = tm_canonical == "PsrC" or has_pf14589
    tree_incompatible_with_psr = psrABC_like_operon and tree_group in {
        "arr_arx", "soe_sre", "ttr_srd", "known_nonpsr_mopb"
    }
    if psrABC_like_operon:
        evidence.append("PsrABC_like_operon(+)")
    if abc_like_nonpsr_operon:
        evidence.append("ABC_like_nonPsr_operon(+)")
    if psrC_like_context:
        evidence.append("PsrC_like_context(+)")
    if tree_incompatible_with_psr:
        evidence.append("Tree_incompatible_with_Psr(+)")

    # --- Final tree-gated classification ---
    if not has_mo:
        classification, confidence = "NOT_MoBisPGD_enzyme", "HIGH"

    elif tree_group == "psr_phs_srr":
        if score >= 8 and strong_psr_operon:
            classification, confidence = "TRUE_PsrA", "HIGH"
        elif score >= 5 and moderate_psr_operon:
            classification, confidence = "LIKELY_PsrA", "MEDIUM"
        elif moderate_psr_operon:
            classification, confidence = "PsrA_or_PhsA", "MEDIUM"
        else:
            classification, confidence = "AMBIGUOUS", "LOW"

    elif tree_group == "arr_arx":
        if has_operon_evidence:
            classification, confidence = "LIKELY_ArrA", "LOW"
        else:
            classification, confidence = "AMBIGUOUS", "LOW"

    elif tree_group == "soe_sre":
        classification, confidence = "LIKELY_SoeA_or_divergent", "MEDIUM"

    elif tree_group == "ttr_srd":
        classification, confidence = "LIKELY_TtrA_or_SrdA", "MEDIUM"

    elif tree_group == "known_nonpsr_mopb":
        # Tree places this as a non-Psr MopB/DMSOR family. Preserve local
        # operon evidence, but block Psr labels. Split A+B+C-like contexts into
        # a separate manual-review class because these are biologically distinct
        # from isolated catalytic subunit hits.
        if abc_like_nonpsr_operon:
            classification, confidence = "LIKELY_nonPsr_MopB_operon_like", "MEDIUM"
        else:
            classification, confidence = "LIKELY_nonPsr_MopB", "MEDIUM"

    elif tree_group in {"unresolved", "other"}:
        # Conservative fallback when the tree is absent or not represented by a
        # known reference family. Do not allow TRUE/LIKELY_PsrA without a
        # Psr/Phs/Srr-compatible tree placement.
        if score >= 5 and moderate_psr_operon:
            classification, confidence = "PsrA_or_PhsA", "LOW"
        elif not has_tat and not has_nrfd and tm_canonical not in ("PsrC",):
            classification, confidence = "LIKELY_SoeA_or_divergent", "LOW"
        else:
            classification, confidence = "AMBIGUOUS", "LOW"

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
        "LIKELY_ArrA":              "#e8d5f5",   # purple — ArrA tree + operon evidence
        "LIKELY_TtrA_or_SrdA":      "#d1c4e9",   # lavender — Ttr/Srd tree placement
        "LIKELY_nonPsr_MopB_operon_like": "#bdbdbd",   # darker grey — non-Psr tree with A+B+C operon context
        "LIKELY_nonPsr_MopB":             "#e0e0e0",   # grey — non-Psr MopB/DMSOR tree placement
        "NOT_MoBisPGD_enzyme":      "#ffcdd2",
        "AMBIGUOUS":                "#e3f2fd",
    }

    # Discovery source left-border accent colours
    src_border = {
        "HMSS2_gate":  "4px solid #e65100",   # orange — needs manual review
        "both_gates":  "4px solid #1565c0",   # blue — highest confidence
        "PF00384_gate": "",                    # no accent
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
        "<span class='src-badge src-hmss2'>HMSS2_gate</span> — found by "
        "PsrAPhsASreA.hmm only; not in PF00384 scan — review carefully. &nbsp;",
        "<span class='src-badge src-both'>both_gates</span> — found by both PF00384 "
        "and PsrAPhsASreA.hmm (highest initial confidence). &nbsp;",
        "<span class='src-badge src-pf'>PF00384_gate</span> — standard discovery.",
        "</p>",
    ]
    if hmss2_cols:
        html.append("<p style='color:#555;font-size:10px'>"
                    "Columns shaded grey (HMSS2_*) are annotation only — "
                    "not used in scoring, except SoeA.hmm which contributes "
                    "a conditional score penalty. "
                    "<b style='color:#6a0dad'>LIKELY_ArrA</b> (purple): "
                    "ArrA/ArxA tree placement combined with TAT + NrfD/PsrB operon "
                    "evidence — the ArrABC operon is architecturally identical to "
                    "PsrABC; these require manual review and ecological context "
                    "before interpretation.</p>")
    html += [
        "<table><tr>",
        "<th>Protein ID</th><th>Bin</th><th>Source</th>",
        "<th>Mo-bisPGD</th><th>TAT</th>",
        "<th>Membrane subunit</th><th>NrfD(PF14589?)</th><th>PsrB</th>",
        "<th>Tree clade</th><th>Tree dist.</th><th>Second clade</th><th>Second dist.</th>",
        "<th>Tree margin</th><th>Tree ratio</th><th>Tree assignment confidence</th>",
        "<th>Tree group</th><th>PsrABC-like operon</th>",
        "<th>ABC-like non-Psr operon</th><th>PsrC-like context</th>",
        "<th>Tree incompatible with Psr</th><th>Classification</th><th>Confidence</th>",
        "<th>Evidence</th>",
    ]
    for col in hmss2_cols:
        html.append(f'<th class="hmss2">{col}</th>')
    html.append("</tr>")

    primary_cols = ["prot_id", "bin_name", "discovery_source",
                    "has_PF00384", "has_TAT_signal",
                    "membrane_subunit_class", "NrfD_in_neighbourhood",
                    "PsrB_in_neighbourhood", "tree_clade", "tree_distance",
                    "tree_second_clade", "tree_second_distance", "tree_distance_margin",
                    "tree_distance_ratio", "tree_assignment_confidence", "tree_group",
                    "PsrABC_like_operon", "tree_incompatible_with_Psr",
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
                    "HMSS2_gate":   "src-hmss2",
                    "both_gates":   "src-both",
                    "PF00384_gate": "src-pf",
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

    # Load discovery source (PF00384_gate / HMSS2_gate / both / supplied)
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

    # Parse tree clade assignments and nearest/second-nearest diagnostics
    tree_diagnostics = {}
    if args.treefile and os.path.exists(args.treefile):
        print("[*] Parsing tree for clade assignments and distance diagnostics...")
        tree_diagnostics = parse_treefile_clades(args.treefile, ref_meta)

    final_rows = []
    for prot_id in query_ids:
        row = {"prot_id": prot_id}

        # Discovery source
        row["discovery_source"] = normalise_discovery_source(disc_src.get(prot_id, "unknown"))

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
        row["signalp_prediction"] = tat_info.get("signalp_prediction", "NOT_RUN")
        row["TAT_probability"] = clean_tat_probability(
            row["has_TAT_signal"],
            tat_info.get("TAT_probability", "NA"),
            row["signalp_prediction"],
        )

        # Tree clade and nearest/second-nearest diagnostics
        tree_diag = tree_diagnostics.get(prot_id, {})
        row["tree_clade"] = tree_diag.get("tree_clade", "inspect_tree_manually")
        row["tree_distance"] = tree_diag.get("tree_distance", "NA")
        row["tree_second_clade"] = tree_diag.get("tree_second_clade", "NA")
        row["tree_second_distance"] = tree_diag.get("tree_second_distance", "NA")
        row["tree_distance_margin"] = tree_diag.get("tree_distance_margin", "NA")
        row["tree_distance_ratio"] = tree_diag.get("tree_distance_ratio", "NA")
        row["tree_assignment_confidence"] = tree_diag.get("tree_assignment_confidence", "NA")

        # Context flags used for interpretability. These are also recomputed
        # inside score_classification when constructing the evidence string,
        # but storing them as columns makes downstream review/filtering easier.
        row.update(derive_context_flags(row))

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
            "tree_clade", "tree_distance", "tree_second_clade",
            "tree_second_distance", "tree_distance_margin", "tree_distance_ratio",
            "tree_assignment_confidence", "tree_group", "PsrABC_like_operon",
            "ABC_like_nonPsr_operon", "PsrC_like_context",
            "tree_incompatible_with_Psr",
            "classification", "confidence", "evidence"] + hmss2_cols

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
        flag = "  ← review in HTML" if src == "HMSS2_gate" else ""
        print(f"  {src:<20} : {n}{flag}")
    print(f"\n  TSV  : {tsv_path}")
    print(f"  HTML : {html_path}")
    print(f"\n  GUIDE:")
    print(f"    TRUE_PsrA          : Psr/Phs/Srr tree placement + strong PsrABC evidence")
    print(f"    LIKELY_PsrA        : Psr/Phs/Srr tree placement + moderate/strong operon evidence")
    print(f"    PsrA_or_PhsA       : Psr-compatible or unresolved tree, but lower-confidence subunit assignment")
    print(f"    LIKELY_ArrA        : ArrA/ArxA tree placement + compatible operon evidence")
    print(f"    LIKELY_TtrA_or_SrdA: Ttr/Srd tree placement")
    print(f"    LIKELY_nonPsr_MopB_operon_like : Known non-Psr tree placement with A+B+C-like local operon context")
    print(f"    LIKELY_nonPsr_MopB : Known non-Psr MopB/DMSOR tree placement, e.g. FdhG/NapA/ActB")
    print(f"    LIKELY_SoeA        : Soe/Sre tree placement or no-TAT/no-NrfD divergent case")
    print(f"    AMBIGUOUS          : Conflicting/incomplete evidence")
    print(f"    NrfD_PF03916(topology_ND) : NrfD neighbour confirmed by broad "
          f"HMM only; run DeepTMHMM to resolve to PsrC/TtrC/PhsC")
    print(f"    tree_assignment_confidence : Diagnostic only. LOW means nearest and")
    print(f"        second-nearest labelled reference families are close in tree distance.")
    print(f"\n  NOTE: HMSS2_gate candidates (orange border in HTML) were not")
    print(f"        detected by PF00384. Check tree placement and neighbourhood")
    print(f"        evidence carefully before accepting these as PsrA.")


if __name__ == "__main__":
    main()
