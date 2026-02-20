#!/usr/bin/env python3
"""
run_plma.py

Compute PLMA (Partial Local Multiple Alignment) for paralog families.
Fetches protein sequences, adds orthologs, runs paloma-D via Docker,
parses .agraph output, and generates plma.json for each pair.

Usage:
    python scripts/run_plma.py                   # Process missing families only
    python scripts/run_plma.py --force            # Regenerate all
    python scripts/run_plma.py SMARCA2_SMARCA4    # Specific pair
"""

import argparse
import json
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
import requests

# --- Paths ---
PROJECT_DIR = Path(__file__).resolve().parent.parent
INPUT_DIR = PROJECT_DIR / "input"
DATA_DIR = PROJECT_DIR / "data"
PAIRS_DIR = DATA_DIR / "pairs"
PLMA_WORK_DIR = DATA_DIR / "cache" / "plma"
FEATURES_CSV = INPUT_DIR / "ens111_human_allFeatures.csv"
PAIRS_CSV = INPUT_DIR / "pairs.csv"

# Sequence cache (shared with previous pipeline)
SEQ_DIR = Path("/Users/olivierdennler/Documents/data/mut_disease_paralogs/clinical_protein_fasta")

# OrthoInspector data
ORTHO_CSV = Path("/Users/olivierdennler/Documents/data/SLI_2023/ens111_humanGene_orthoInspectorOrthologs_1472Sp.csv")
TAXID_LIST = [10090, 6239, 4932]  # mouse, C. elegans, yeast

# Known sequences that crash paloma-D (discovered empirically)
PALOMA_BLACKLIST = {"Q9U7E0", "F8VPZ5", "Q93781", "Q86D18"}  # cause segfault/weight mismatch

# Docker
DOCKER_CONTAINER = "04b80bf721e5"

# paloma-D parameters
PALOMA_T = 10
PALOMA_M = 1
PALOMA_BIG_M = 20
PALOMA_Q = 2


def log(msg: str):
    print(f"[plma] {msg}", flush=True)


# ---------------------------------------------------------------------------
# 1. Family discovery
# ---------------------------------------------------------------------------

def get_families_for_pairs(pairs_df, features_df):
    """Group pairs into families using the family_id from features CSV.

    Returns dict: family_id -> {genes: set, uniprots: {gene: uniprot}, pair_ids: list}
    """
    families = {}
    for _, row in pairs_df.iterrows():
        gene_a, gene_b = row["A1"], row["A2"]
        pair_id = f"{gene_a}_{gene_b}"

        # Look up in features CSV (try both orientations)
        feat = features_df[
            ((features_df["A1"] == gene_a) & (features_df["A2"] == gene_b)) |
            ((features_df["A1"] == gene_b) & (features_df["A2"] == gene_a))
        ]
        if feat.empty:
            log(f"  WARNING: {pair_id} not found in features CSV, skipping")
            continue

        r = feat.iloc[0]
        fam_id = str(int(float(r["family_id"])))
        # Ensure correct gene->uniprot mapping
        if r["A1"] == gene_a:
            ua, ub = r["uniprot_A1"], r["uniprot_A2"]
        else:
            ua, ub = r["uniprot_A2"], r["uniprot_A1"]

        if fam_id not in families:
            families[fam_id] = {"genes": set(), "uniprots": {}, "pair_ids": []}
        families[fam_id]["genes"].update([gene_a, gene_b])
        families[fam_id]["uniprots"][gene_a] = ua
        families[fam_id]["uniprots"][gene_b] = ub
        families[fam_id]["pair_ids"].append(pair_id)

    return families


def expand_family_genes(family_id, features_df):
    """Get ALL genes in a family from features CSV (not just those in our pairs)."""
    fam_rows = features_df[features_df["family_id"] == int(family_id)]
    genes = set()
    uniprots = {}
    for _, r in fam_rows.iterrows():
        genes.add(r["A1"])
        genes.add(r["A2"])
        uniprots[r["A1"]] = r["uniprot_A1"]
        uniprots[r["A2"]] = r["uniprot_A2"]
    return genes, uniprots


# ---------------------------------------------------------------------------
# 2. Sequence fetching
# ---------------------------------------------------------------------------

def fetch_uniprot_fasta(uniprot_id: str) -> Path:
    """Get FASTA file for a UniProt ID, downloading if needed."""
    fasta_path = SEQ_DIR / f"{uniprot_id}.fasta"
    if fasta_path.exists() and fasta_path.stat().st_size > 0:
        return fasta_path

    log(f"    Downloading {uniprot_id} from UniProt...")
    url = f"https://rest.uniprot.org/uniprotkb/search?query={uniprot_id}&format=fasta"
    resp = requests.get(url, timeout=30)
    if resp.status_code == 200 and resp.text.startswith(">"):
        fasta_path.write_text(resp.text)
        return fasta_path

    # Try direct accession
    url2 = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    resp2 = requests.get(url2, timeout=30)
    if resp2.status_code == 200 and resp2.text.startswith(">"):
        fasta_path.write_text(resp2.text)
        return fasta_path

    log(f"    WARNING: Could not download {uniprot_id}")
    return None


def get_fasta_length(fasta_path: Path) -> int:
    """Get sequence length from FASTA file."""
    seq = ""
    for line in fasta_path.read_text().splitlines():
        if not line.startswith(">"):
            seq += line.strip()
    return len(seq)


def resolve_multi_uniprot(uniprot_ids_str: str) -> str:
    """If multiple UniProt IDs (pipe-separated), pick the one with existing FASTA or longest."""
    if "|" not in str(uniprot_ids_str):
        return str(uniprot_ids_str)

    ids = str(uniprot_ids_str).split("|")
    # Prefer existing files
    for uid in ids:
        fasta_path = SEQ_DIR / f"{uid.strip()}.fasta"
        if fasta_path.exists():
            return uid.strip()
    # Otherwise return first
    return ids[0].strip()


# ---------------------------------------------------------------------------
# 3. Ortholog fetching
# ---------------------------------------------------------------------------

def load_orthologs_for_genes(gene_uniprots: dict, ortho_df: pd.DataFrame) -> list:
    """Get ortholog UniProt IDs for all genes in a family."""
    ortho_ids = set()
    for gene, uniprot in gene_uniprots.items():
        filtered = ortho_df[
            (ortho_df["uniprot_id"] == uniprot) &
            (ortho_df["species"].isin(TAXID_LIST))
        ]
        for ortho_list_str in filtered["orthologs"]:
            ortho_list = str(ortho_list_str).strip("[]").replace("'", "").split(", ")
            for oid in ortho_list:
                oid = oid.strip()
                if oid and oid != "nan":
                    ortho_ids.add(oid)
    return list(ortho_ids)


# ---------------------------------------------------------------------------
# 4. Build multi-FASTA
# ---------------------------------------------------------------------------

def build_family_fasta(family_id, gene_uniprots, ortho_ids, work_dir):
    """Build the multi-FASTA file for a family (paralogs + orthologs)."""
    work_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = work_dir / f"family_{family_id}.fasta"

    seen_files = set()
    fasta_parts = []

    # Add paralog sequences
    for gene, uniprot in gene_uniprots.items():
        fp = fetch_uniprot_fasta(uniprot)
        if fp and str(fp) not in seen_files:
            seen_files.add(str(fp))
            fasta_parts.append(fp.read_text().rstrip())

    # Add ortholog sequences (skip blacklisted)
    for oid in ortho_ids:
        if oid in PALOMA_BLACKLIST:
            log(f"    Skipping blacklisted {oid}")
            continue
        fp = fetch_uniprot_fasta(oid)
        if fp and str(fp) not in seen_files:
            seen_files.add(str(fp))
            fasta_parts.append(fp.read_text().rstrip())

    fasta_path.write_text("\n".join(fasta_parts) + "\n")
    log(f"  Family FASTA: {len(seen_files)} sequences written to {fasta_path.name}")
    return fasta_path


# ---------------------------------------------------------------------------
# 5. Run paloma-D via Docker
# ---------------------------------------------------------------------------

def run_paloma_docker(fasta_path: Path, work_dir: Path) -> Path:
    """Run paloma-D inside Docker container, return .agraph path."""
    t, m, M, q = PALOMA_T, PALOMA_M, PALOMA_BIG_M, PALOMA_Q
    agraph_name = f"{fasta_path.stem}_t{t}m{m}M{M}_q{q}.agraph"
    agraph_local = work_dir / agraph_name

    # Docker workflow: cleanup old -> copy in -> exec paloma-D -> copy out -> cleanup
    cmd_in_docker = f". ~/.bashrc && paloma-D -i {fasta_path.name} -t {t} -m {m} -M {M} -q {q} --oplma && exit"
    cleanup_pattern = fasta_path.stem
    cmd_cleanup = f"rm -f {cleanup_pattern}*"

    cmd = (
        f"docker start {DOCKER_CONTAINER} && "
        f'docker exec {DOCKER_CONTAINER} /bin/bash -c "{cmd_cleanup}" ; '
        f"docker cp {fasta_path} {DOCKER_CONTAINER}:/ && "
        f"docker start {DOCKER_CONTAINER} && "
        f'docker exec {DOCKER_CONTAINER} /bin/bash -c "{cmd_in_docker}" && '
        f"docker cp {DOCKER_CONTAINER}:/{agraph_name} {work_dir}/ && "
        f"docker start {DOCKER_CONTAINER} && "
        f'docker exec {DOCKER_CONTAINER} /bin/bash -c "{cmd_cleanup}"'
    )

    # Remove any existing agraph from a previous attempt
    if agraph_local.exists():
        agraph_local.unlink()

    log(f"  Running paloma-D in Docker ({_count_fasta_seqs(fasta_path)} seqs)...")
    result = subprocess.run(cmd, shell=True, capture_output=True, timeout=600)

    if agraph_local.exists():
        log(f"  paloma-D produced {agraph_name} ({agraph_local.stat().st_size} bytes)")
        return agraph_local

    # Log error
    stderr = result.stderr.decode("utf-8", errors="replace")
    if stderr:
        log(f"  paloma-D stderr: {stderr[:500]}")
    log(f"  ERROR: paloma-D did not produce {agraph_name}")
    return None


def _count_fasta_seqs(fasta_path: Path) -> int:
    return sum(1 for line in fasta_path.read_text().splitlines() if line.startswith(">"))


def run_paloma_with_retry(fasta_path, work_dir, human_ids, gene_uniprots=None):
    """Run paloma-D, retry with human-only if it fails."""
    result = run_paloma_docker(fasta_path, work_dir)
    if result:
        return result

    if not gene_uniprots:
        return None

    # Retry: rebuild FASTA from individual human-only files
    log(f"  Retrying paloma-D with human-only sequences...")
    fasta_parts = []
    seen = set()
    for gene, uniprot in gene_uniprots.items():
        fp = SEQ_DIR / f"{uniprot}.fasta"
        if fp.exists() and str(fp) not in seen:
            seen.add(str(fp))
            fasta_parts.append(fp.read_text().rstrip())

    fasta_path.write_text("\n".join(fasta_parts) + "\n")
    log(f"  Rebuilt with {len(fasta_parts)} human-only sequences")

    result = run_paloma_docker(fasta_path, work_dir)
    return result


# ---------------------------------------------------------------------------
# 6. Parse .agraph file
# ---------------------------------------------------------------------------

def parse_agraph(agraph_path: Path):
    """Parse paloma-D .agraph output into block data.

    Returns:
        blocks: {block_id: {seq_num_str: {start, end, length, seq}}}
        seq_info: {seq_num_str: {id, desc, seq, length}}
    """
    text = agraph_path.read_text()
    lines = text.splitlines()

    # Parse sequences
    seq_info = {}
    i = 0
    while i < len(lines):
        line = lines[i].rstrip()
        if line == "Sequences:":
            i += 1
            while i < len(lines) and not lines[i].startswith("Aligned"):
                line = lines[i].rstrip()
                # Match "  N:  # sequence"
                m = re.match(r"^\s+(\d+):\s+# sequence", line)
                if m:
                    seq_num = m.group(1)
                    seq_data = {"num": seq_num}
                    i += 1
                    while i < len(lines):
                        sl = lines[i].rstrip()
                        if re.match(r"^\s+\d+:\s+# sequence", sl) or sl.startswith("Aligned"):
                            break
                        id_m = re.match(r'\s+id:\s+(.+)', sl)
                        desc_m = re.match(r'\s+descr:\s+"?\s*(.+?)"?\s*$', sl)
                        seq_m = re.match(r'\s+seq:\s+(\S+)', sl)
                        weight_m = re.match(r'\s+weight:', sl)
                        if id_m:
                            seq_data["id"] = id_m.group(1).strip()
                        elif desc_m:
                            seq_data["desc"] = desc_m.group(1).strip().rstrip('"')
                        elif seq_m:
                            seq_data["seq"] = seq_m.group(1)
                            seq_data["length"] = len(seq_m.group(1))
                        elif weight_m:
                            pass
                        else:
                            pass
                        i += 1
                    seq_info[seq_num] = seq_data
                    continue
                i += 1
            continue
        i += 1

    # Parse aligned sites -> blocks
    blocks = {}
    block_num = 0
    i = 0
    while i < len(lines):
        line = lines[i].rstrip()
        if line.strip() == "Aligned sites:":
            i += 1
            # Skip contig line "  N:  # alignment contig"
            while i < len(lines):
                line = lines[i].rstrip()
                contig_m = re.match(r"^\s+\d+:\s+# alignment contig", line)
                if contig_m:
                    i += 1
                    # Parse blocks within this contig
                    while i < len(lines):
                        line = lines[i].rstrip()
                        block_m = re.match(r"^\s+(\d+):\s+# block", line)
                        if block_m:
                            block_num += 1
                            block_id = f"B{block_num}"
                            block_start_pos = int(block_m.group(1))
                            block_members = {}  # seq_num -> list of positions
                            i += 1
                            while i < len(lines):
                                sl = lines[i].rstrip()
                                # Match lines like "      - {1: [pos], 3: [pos], ...}"
                                if re.match(r"^\s+- \{", sl):
                                    # Parse {seq: [pos], ...}
                                    pairs = re.findall(r"(\d+):\s*\[(\d+)\]", sl)
                                    for sn, pos in pairs:
                                        if sn not in block_members:
                                            block_members[sn] = []
                                        block_members[sn].append(int(pos))
                                    i += 1
                                else:
                                    break
                            # Convert position lists to block info
                            block_data = {}
                            for sn, positions in block_members.items():
                                start = min(positions)
                                end = max(positions)
                                length = len(positions)
                                # Extract AA sequence from full sequence
                                seq_str = ""
                                if sn in seq_info and "seq" in seq_info[sn]:
                                    full_seq = seq_info[sn]["seq"]
                                    for p in sorted(positions):
                                        if 1 <= p <= len(full_seq):
                                            seq_str += full_seq[p - 1]
                                block_data[sn] = {
                                    "start": start,
                                    "end": end,
                                    "length": length,
                                    "seq": seq_str,
                                }
                            if block_data:
                                blocks[block_id] = block_data
                            continue
                        elif re.match(r"^\s+\d+:\s+# alignment contig", line):
                            # New contig
                            break
                        else:
                            i += 1
                    continue
                else:
                    i += 1
            break
        i += 1

    return blocks, seq_info


# ---------------------------------------------------------------------------
# 7. Build plma.json (reuse extract_plma.py logic)
# ---------------------------------------------------------------------------

def categorize_blocks(blocks, gene_a_seq, gene_b_seq, human_seqs):
    """Categorize blocks by conservation pattern relative to the pair."""
    categorized = []
    for bloc_id, members in blocks.items():
        seq_nums = set(members.keys())
        humans_in = seq_nums & human_seqs
        a_in = gene_a_seq in humans_in
        b_in = gene_b_seq in humans_in
        other_humans = humans_in - {gene_a_seq, gene_b_seq}

        if a_in and b_in:
            category = "shared_with_family" if other_humans else "pair_exclusive"
        elif a_in:
            category = "a_with_family" if other_humans else "specific_a"
        elif b_in:
            category = "b_with_family" if other_humans else "specific_b"
        else:
            category = "family_only"

        categorized.append({
            "id": bloc_id,
            "category": category,
            "members": members,
            "n_seqs": len(seq_nums),
            "n_human": len(humans_in),
        })
    return categorized


def build_plma_json_from_agraph(agraph_path, family_id, gene_a, gene_b, uniprot_a, uniprot_b):
    """Build plma.json data from .agraph file."""
    blocks, seq_info = parse_agraph(agraph_path)

    if not blocks:
        log(f"  WARNING: No blocks found in {agraph_path.name}")
        return None

    # Build seq_descriptions matching extract_plma format
    seq_descriptions = {}
    for sn, info in seq_info.items():
        full_id = info.get("id", "")
        desc = info.get("desc", "")
        # Extract uniprot from "sp|XXXXX|NAME_SPECIES"
        uniprot = ""
        if "|" in full_id:
            parts = full_id.split("|")
            if len(parts) >= 2:
                uniprot = parts[1]
        # Extract gene name from description
        gene_name = ""
        gn_match = re.search(r"GN=(\S+)", desc)
        if gn_match:
            gene_name = gn_match.group(1)

        seq_descriptions[sn] = {
            "uniprot": uniprot,
            "gene": gene_name,
            "desc": f"{full_id} {desc}",
        }

    # Identify human sequences
    human_seqs = set()
    for sn, info in seq_descriptions.items():
        desc = info.get("desc", "")
        if "OX=9606" in desc or "_HUMAN" in desc:
            human_seqs.add(sn)
    n_orthologs = len(seq_descriptions) - len(human_seqs)
    log(f"  {len(human_seqs)} human + {n_orthologs} ortholog sequences, {len(blocks)} blocks")

    # Find seq numbers for gene A and gene B
    gene_a_seq = gene_b_seq = None
    for sn, info in seq_descriptions.items():
        if info["uniprot"] == uniprot_a or info["gene"] == gene_a:
            gene_a_seq = sn
        if info["uniprot"] == uniprot_b or info["gene"] == gene_b:
            gene_b_seq = sn

    if gene_a_seq is None:
        for sn, info in seq_descriptions.items():
            if gene_a.upper() in info.get("desc", "").upper():
                gene_a_seq = sn
                break
    if gene_b_seq is None:
        for sn, info in seq_descriptions.items():
            if gene_b.upper() in info.get("desc", "").upper():
                gene_b_seq = sn
                break

    if gene_a_seq is None or gene_b_seq is None:
        log(f"  ERROR: Could not map genes to sequences (A={gene_a_seq}, B={gene_b_seq})")
        log(f"  Looking for {gene_a}={uniprot_a}, {gene_b}={uniprot_b}")
        log(f"  Available: {[(sn, d['uniprot'], d['gene']) for sn, d in seq_descriptions.items()]}")
        return None

    # Build sequences list
    sequences = []
    for sn in sorted(seq_info.keys(), key=lambda x: int(x)):
        info = seq_descriptions.get(sn, {})
        sequences.append({
            "num": sn,
            "gene": info.get("gene", ""),
            "uniprot": info.get("uniprot", ""),
            "length": seq_info[sn].get("length", 0),
            "is_gene_a": sn == gene_a_seq,
            "is_gene_b": sn == gene_b_seq,
            "is_human": sn in human_seqs,
        })

    # Categorize blocks
    categorized = categorize_blocks(blocks, gene_a_seq, gene_b_seq, human_seqs)

    # Build block JSON
    blocks_json = []
    for block in categorized:
        b = {
            "id": block["id"],
            "category": block["category"],
            "n_seqs": block["n_seqs"],
            "n_human": block["n_human"],
            "positions": {},
        }
        for sn, info in block["members"].items():
            b["positions"][sn] = {
                "start": info["start"],
                "end": info["end"],
                "length": info["length"],
                "seq": info.get("seq", ""),
            }
        blocks_json.append(b)

    # Sort blocks (same topological sort as extract_plma.py)
    _sort_blocks(blocks_json, sequences)

    # Category summary
    summary = {
        "shared_with_family": 0, "pair_exclusive": 0,
        "specific_a": 0, "specific_b": 0,
        "a_with_family": 0, "b_with_family": 0,
        "family_only": 0,
    }
    for b in blocks_json:
        cat = b["category"]
        if cat in ("shared_with_family", "pair_exclusive"):
            la = b["positions"].get(gene_a_seq, {}).get("length", 0)
            lb = b["positions"].get(gene_b_seq, {}).get("length", 0)
            summary[cat] += (la + lb) // 2
        elif cat in ("specific_a", "a_with_family"):
            summary[cat] += b["positions"].get(gene_a_seq, {}).get("length", 0)
        elif cat in ("specific_b", "b_with_family"):
            summary[cat] += b["positions"].get(gene_b_seq, {}).get("length", 0)
        elif cat == "family_only":
            lengths = [p["length"] for p in b["positions"].values()]
            summary[cat] += max(lengths) if lengths else 0

    return {
        "family_id": str(family_id),
        "gene_a": gene_a,
        "gene_b": gene_b,
        "gene_a_seq": gene_a_seq,
        "gene_b_seq": gene_b_seq,
        "n_blocks": len(blocks_json),
        "n_sequences": len(sequences),
        "sequences": sequences,
        "blocks": blocks_json,
        "summary": summary,
    }


def _sort_blocks(blocks_json, sequences):
    """Sort blocks to minimize positional inversions (same logic as extract_plma.py)."""
    from collections import defaultdict
    import heapq

    n_blk = len(blocks_json)
    if n_blk <= 1:
        return

    sl_map = {s["num"]: s["length"] for s in sequences if s.get("length")}

    norm_pos = []
    for b in blocks_json:
        fracs = []
        for sn, p in b["positions"].items():
            sl = sl_map.get(sn, 0)
            if sl > 0 and p.get("start") is not None:
                fracs.append(p["start"] / sl)
        norm_pos.append(sum(fracs) / len(fracs) if fracs else 1.0)

    adj = defaultdict(set)
    in_deg = [0] * n_blk

    for i in range(n_blk):
        pi = blocks_json[i]["positions"]
        for j in range(i + 1, n_blk):
            pj = blocks_json[j]["positions"]
            shared = set(pi.keys()) & set(pj.keys())
            if not shared:
                continue
            v_ij = v_ji = 0
            for sn in shared:
                si, sj = pi[sn]["start"], pj[sn]["start"]
                if si < sj:
                    v_ij += 1
                elif sj < si:
                    v_ji += 1
            if v_ij > v_ji:
                adj[i].add(j)
                in_deg[j] += 1
            elif v_ji > v_ij:
                adj[j].add(i)
                in_deg[i] += 1

    heap = [(norm_pos[i], i) for i in range(n_blk) if in_deg[i] == 0]
    heapq.heapify(heap)
    order = []
    visited = set()
    while heap:
        _, node = heapq.heappop(heap)
        if node in visited:
            continue
        visited.add(node)
        order.append(node)
        for nb in adj[node]:
            in_deg[nb] -= 1
            if in_deg[nb] == 0:
                heapq.heappush(heap, (norm_pos[nb], nb))

    if len(order) < n_blk:
        remaining = [(norm_pos[i], i) for i in range(n_blk) if i not in visited]
        remaining.sort()
        while len(order) < n_blk:
            best = None
            for _, idx in remaining:
                if idx not in visited:
                    best = idx
                    break
            if best is None:
                break
            visited.add(best)
            order.append(best)
            for nb in adj[best]:
                in_deg[nb] -= 1
                if in_deg[nb] == 0 and nb not in visited:
                    heapq.heappush(heap, (norm_pos[nb], nb))
            while heap:
                _, node = heapq.heappop(heap)
                if node in visited:
                    continue
                visited.add(node)
                order.append(node)
                for nb in adj[node]:
                    in_deg[nb] -= 1
                    if in_deg[nb] == 0 and nb not in visited:
                        heapq.heappush(heap, (norm_pos[nb], nb))

    blocks_json[:] = [blocks_json[i] for i in order]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Compute PLMA for paralog families")
    parser.add_argument("pair_id", nargs="?", help="Specific pair to process")
    parser.add_argument("--force", action="store_true", help="Regenerate all (ignore existing)")
    args = parser.parse_args()

    # Load data
    log("Loading features CSV...")
    features_df = pd.read_csv(FEATURES_CSV, usecols=[
        "A1", "A2", "family_id", "uniprot_A1", "uniprot_A2",
    ], low_memory=False)

    log("Loading pairs CSV...")
    pairs_df = pd.read_csv(PAIRS_CSV)
    if args.pair_id:
        parts = args.pair_id.split("_", 1)
        pairs_df = pairs_df[
            ((pairs_df["A1"] == parts[0]) & (pairs_df["A2"] == parts[1])) |
            ((pairs_df["A1"] == parts[1]) & (pairs_df["A2"] == parts[0]))
        ]
        if pairs_df.empty:
            log(f"ERROR: Pair {args.pair_id} not found in pairs.csv")
            sys.exit(1)

    log("Loading OrthoInspector data...")
    ortho_df = pd.read_csv(ORTHO_CSV, low_memory=False)

    # Group into families
    families = get_families_for_pairs(pairs_df, features_df)
    log(f"Found {len(families)} unique families for {len(pairs_df)} pairs")

    # Ensure work directory
    PLMA_WORK_DIR.mkdir(parents=True, exist_ok=True)
    SEQ_DIR.mkdir(parents=True, exist_ok=True)

    success = 0
    errors = 0

    for fam_id, fam_info in families.items():
        log(f"\n{'='*60}")
        log(f"Family {fam_id}: {sorted(fam_info['genes'])}")
        log(f"  Pairs: {fam_info['pair_ids']}")

        # Check if all plma.json exist (skip unless --force)
        if not args.force:
            all_exist = all(
                (PAIRS_DIR / pid / "plma.json").exists()
                for pid in fam_info["pair_ids"]
            )
            if all_exist:
                log(f"  All plma.json exist, skipping (use --force to regenerate)")
                success += len(fam_info["pair_ids"])
                continue

        # Expand to full family from features CSV
        all_genes, all_uniprots = expand_family_genes(fam_id, features_df)
        # Merge with our pair uniprots (more reliable)
        all_uniprots.update(fam_info["uniprots"])
        log(f"  Full family: {len(all_genes)} genes")

        # 2-gene families: skip PLMA (handled in JS as "2-gene family" note)
        if len(all_genes) <= 2:
            log(f"  2-gene family — skipping PLMA computation")
            for pair_id in fam_info["pair_ids"]:
                out_dir = PAIRS_DIR / pair_id
                out_dir.mkdir(parents=True, exist_ok=True)
                plma_file = out_dir / "plma.json"
                # Write empty marker so we know it was processed
                plma_file.write_text(json.dumps({"family_id": fam_id, "n_blocks": 0, "note": "2-gene family"}))
                log(f"  Wrote {pair_id}/plma.json (2-gene family marker)")
                success += 1
            continue

        # Fetch ortholog IDs
        log(f"  Fetching orthologs for {len(all_uniprots)} genes...")
        ortho_ids = load_orthologs_for_genes(all_uniprots, ortho_df)
        log(f"  Found {len(ortho_ids)} ortholog sequences")

        # Build multi-FASTA
        fam_work = PLMA_WORK_DIR / f"family_{fam_id}"
        fasta_path = build_family_fasta(fam_id, all_uniprots, ortho_ids, fam_work)

        # Run paloma-D (with retry on failure)
        human_ids = set(all_uniprots.values())
        agraph_path = run_paloma_with_retry(fasta_path, fam_work, human_ids, all_uniprots)
        if agraph_path is None:
            log(f"  ERROR: paloma-D failed for family {fam_id}")
            errors += len(fam_info["pair_ids"])
            continue

        # Build plma.json for each pair in this family
        for pair_id in fam_info["pair_ids"]:
            gene_a, gene_b = pair_id.split("_", 1)
            uniprot_a = fam_info["uniprots"].get(gene_a, all_uniprots.get(gene_a, ""))
            uniprot_b = fam_info["uniprots"].get(gene_b, all_uniprots.get(gene_b, ""))

            log(f"  Building plma.json for {pair_id}...")
            plma_data = build_plma_json_from_agraph(
                agraph_path, fam_id, gene_a, gene_b, uniprot_a, uniprot_b
            )

            if plma_data is None:
                errors += 1
                continue

            out_dir = PAIRS_DIR / pair_id
            out_dir.mkdir(parents=True, exist_ok=True)
            out_file = out_dir / "plma.json"
            with open(out_file, "w") as f:
                json.dump(plma_data, f)

            log(f"  Wrote {pair_id}/plma.json: {plma_data['n_blocks']} blocks, {plma_data['n_sequences']} sequences")
            success += 1

    log(f"\nDone! Success: {success}, Errors: {errors}")


if __name__ == "__main__":
    main()
