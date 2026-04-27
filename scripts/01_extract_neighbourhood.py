#!/usr/bin/env python3
"""
01_extract_neighbourhood.py

Extracts genomic neighbourhoods (±WINDOW genes, same contig) for a list of
query protein IDs from Pyrodigal-annotated genome bins.

IMPORTANT — ID FORMAT:
    Pyrodigal names proteins after the CONTIG they come from, not the bin.
    For example, bin "SF2967_J5348.rosella_refine.165_0" may contain contig
    "NODE_3193_length_16941_cov_1.235637", whose proteins are named:
        NODE_3193_length_16941_cov_1.235637_1
        NODE_3193_length_16941_cov_1.235637_2  ...
    The bin name does NOT appear in the protein ID, so we cannot derive the
    FAA/GFF path by parsing the protein ID. Instead, we use the protein→bin
    index built by 00_pipeline.sh (Step 0a), which maps each protein_id to
    its bin_name, faa_path, and gff_path.

GFF format (Pyrodigal):
    NODE_3193_length_16941_cov_1.235637  pyrodigal_v2.3.0  CDS  1  186  ...
    ID field in attributes: ID=NODE_3193_length_16941_cov_1.235637_1;...

USAGE:
    python3 01_extract_neighbourhood.py \
        --ids   query_ids.txt \
        --index protein_to_bin_index.tsv \
        --window 10 \
        --outdir 02_neighbourhoods/
"""

import argparse
import os
import re


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--ids",    required=True,
                   help="File with query protein IDs, one per line")
    p.add_argument("--index",  required=True,
                   help="protein_to_bin_index.tsv from pipeline Step 0a "
                        "(columns: protein_id, bin_name, faa_path, gff_path)")
    p.add_argument("--window", type=int, default=10,
                   help="Number of genes either side of target to extract (default: 10)")
    p.add_argument("--outdir", required=True,
                   help="Output directory for neighbourhood files")
    return p.parse_args()


def load_index(index_path):
    """
    Load protein->bin index TSV.
    Returns dict: protein_id -> {"bin_name": str, "faa_path": str, "gff_path": str}
    """
    index = {}
    with open(index_path) as fh:
        fh.readline()  # skip header line
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            prot_id, bin_name, faa_path, gff_path = parts[0], parts[1], parts[2], parts[3]
            index[prot_id] = {
                "bin_name": bin_name,
                "faa_path": faa_path,
                "gff_path": gff_path,
            }
    return index


def load_gff(gff_path):
    """
    Parse a Pyrodigal GFF file.
    Returns:
        features  : list of dicts in genomic order, each with:
                    contig, gene_id, start, end, strand
        gene_index: dict mapping gene_id -> position in features list
    """
    features   = []
    gene_index = {}

    with open(gff_path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            contig = parts[0]
            start  = int(parts[3])
            end    = int(parts[4])
            strand = parts[6]
            attrs  = parts[8]

            m = re.search(r"ID=([^;]+)", attrs)
            if not m:
                continue
            gene_id = m.group(1)

            idx = len(features)
            features.append({
                "contig": contig,
                "gene_id": gene_id,
                "start":  start,
                "end":    end,
                "strand": strand,
            })
            gene_index[gene_id] = idx

    return features, gene_index


def load_faa(faa_path):
    """Returns dict: protein_id -> amino acid sequence string"""
    seqs        = {}
    current_id  = None
    current_seq = []
    with open(faa_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id:
                    seqs[current_id] = "".join(current_seq)
                current_id  = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id:
        seqs[current_id] = "".join(current_seq)
    return seqs


def safe_filename(prot_id):
    """Replace characters unsafe for filenames."""
    return re.sub(r"[^\w.\-]", "_", prot_id)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Load query IDs
    with open(args.ids) as fh:
        query_ids = [l.strip() for l in fh if l.strip() and not l.startswith("#")]
    print(f"[*] {len(query_ids)} query proteins to process")

    # Load protein->bin index (built by pipeline Step 0a)
    print(f"[*] Loading protein->bin index: {args.index}")
    index = load_index(args.index)
    print(f"    {len(index)} proteins in index")

    # Cache GFF and FAA per file path to avoid re-loading for proteins in same bin
    gff_cache = {}   # gff_path -> (features, gene_index)
    faa_cache = {}   # faa_path -> seq_dict

    summary_rows  = []
    n_ok          = 0
    n_no_index    = 0
    n_no_gff      = 0

    for prot_id in query_ids:

        # Step 1: look up which bin/file this protein comes from
        if prot_id not in index:
            print(f"  [WARN] {prot_id}: not found in protein index — "
                  f"does this ID exist in a FAA file under --bindir?")
            n_no_index += 1
            continue

        info     = index[prot_id]
        bin_name = info["bin_name"]
        faa_path = info["faa_path"]
        gff_path = info["gff_path"]

        if not gff_path or not os.path.exists(gff_path):
            print(f"  [WARN] {prot_id}: GFF not found at '{gff_path}'")
            n_no_gff += 1
            continue

        # Step 2: load GFF and FAA for this bin (cached)
        if gff_path not in gff_cache:
            gff_cache[gff_path] = load_gff(gff_path)
        features, gene_index = gff_cache[gff_path]

        if faa_path not in faa_cache:
            faa_cache[faa_path] = load_faa(faa_path)
        seqs = faa_cache[faa_path]

        # Step 3: locate target gene in the GFF
        if prot_id not in gene_index:
            print(f"  [WARN] {prot_id}: ID not found in GFF {gff_path}")
            print(f"         Possible ID mismatch between FAA header and GFF ID= field.")
            n_no_gff += 1
            continue

        target_idx    = gene_index[prot_id]
        target_contig = features[target_idx]["contig"]

        # Step 4: collect neighbours on the SAME contig only (never cross boundaries)
        neighbours = []
        for offset in range(-args.window, args.window + 1):
            i = target_idx + offset
            if i < 0 or i >= len(features):
                continue
            feat = features[i]
            if feat["contig"] != target_contig:
                continue
            neighbours.append((offset, feat))

        # Step 5: write neighbourhood files
        fname   = safe_filename(prot_id)
        out_faa = os.path.join(args.outdir, f"{fname}_neighbourhood.faa")
        out_tsv = os.path.join(args.outdir, f"{fname}_neighbourhood.tsv")

        with open(out_faa, "w") as ffa, open(out_tsv, "w") as ftsv:
            ftsv.write("offset\tgene_id\tcontig\tstart\tend\tstrand\tis_target\tbin_name\n")
            for offset, feat in neighbours:
                gid       = feat["gene_id"]
                is_target = (offset == 0)
                ftsv.write(
                    f"{offset:+d}\t{gid}\t{feat['contig']}\t"
                    f"{feat['start']}\t{feat['end']}\t{feat['strand']}\t"
                    f"{'YES' if is_target else 'no'}\t{bin_name}\n"
                )
                if gid in seqs:
                    tag = "[target]" if is_target else f"[offset={offset:+d}]"
                    ffa.write(f">{gid} {tag} from:{prot_id} bin:{bin_name}\n")
                    seq = seqs[gid]
                    for j in range(0, len(seq), 60):
                        ffa.write(seq[j:j + 60] + "\n")
                else:
                    # Present in GFF but not FAA: RNA gene, incomplete ORF, etc.
                    print(f"    [INFO] {gid}: in GFF but not FAA (RNA/incomplete ORF?)")

        print(f"  [OK] {prot_id}  bin={bin_name}  contig={target_contig}  "
              f"{len(neighbours)} neighbours written")
        summary_rows.append({
            "prot_id":      prot_id,
            "bin_name":     bin_name,
            "contig":       target_contig,
            "n_neighbours": len(neighbours),
            "faa_path":     faa_path,
            "gff_path":     gff_path,
        })
        n_ok += 1

    # Summary TSV
    summary_path = os.path.join(args.outdir, "neighbourhood_summary.tsv")
    with open(summary_path, "w") as fh:
        fh.write("prot_id\tbin_name\tcontig\tn_neighbours\tfaa_path\tgff_path\n")
        for row in summary_rows:
            fh.write(
                f"{row['prot_id']}\t{row['bin_name']}\t{row['contig']}\t"
                f"{row['n_neighbours']}\t{row['faa_path']}\t{row['gff_path']}\n"
            )

    print(f"\n[DONE] Neighbourhood extraction complete:")
    print(f"  Successful   : {n_ok}")
    print(f"  Not in index : {n_no_index}  (protein ID absent from all FAA files)")
    print(f"  GFF issues   : {n_no_gff}   (missing GFF or protein ID not in GFF)")
    print(f"  Summary      : {summary_path}")


if __name__ == "__main__":
    main()
