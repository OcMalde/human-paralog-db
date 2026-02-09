#!/usr/bin/env python3
"""
extract_plma.py

Parse PLMA dot files and generate plma.json for each pair in the database.
Extracts alignment blocks, categorizes them by conservation pattern,
and outputs data for the alignment visualization.

Usage:
    python scripts/extract_plma.py
    python scripts/extract_plma.py SMARCA2_SMARCA4   # Specific pair
"""

import argparse
import functools
import json
import os
import re
import sys
from pathlib import Path

import pandas as pd

# Project paths
PROJECT_DIR = Path(__file__).resolve().parent.parent
INPUT_DIR = PROJECT_DIR / "input"
DATA_DIR = PROJECT_DIR / "data"
PAIRS_DIR = DATA_DIR / "pairs"
FEATURES_CSV = INPUT_DIR / "ens111_human_allFeatures.csv"
PLMA_DIR = Path("/Users/olivierdennler/Documents/data/SLI_2023/families_fasta")


def log(msg):
    print(f"[plma] {msg}", flush=True)


def extract_block_sequences(plma_dot_file):
    """Parse a PLMA dot file and extract block sequences.

    Returns:
        dict_bloc_infos: {block_id: {seq_num: {start, end, seq}}}
        seq_number_name: {seq_num_str: uniprot_id}
        seq_lengths: {seq_num_str: int}
        seq_descriptions: {seq_num_str: full_label}
    """
    with open(plma_dot_file, 'r') as f:
        lines = f.readlines()

    concat = functools.reduce(lambda a, b: a.strip() + b.strip(), lines)

    # Parse cluster_1: sequence names
    seq_number_name = {}
    seq_descriptions = {}
    parts = concat.split('subgraph cluster_')
    if len(parts) < 2:
        return {}, {}, {}, {}

    cluster1 = parts[1]
    nodes = cluster1.split(';\"')[1:]
    for node in nodes:
        if '[label = "' not in node:
            continue
        label_part = node.split('[label = "')[1].split('"]')[0]
        seq_num = label_part.split(": ")[0].strip()
        seq_desc = label_part.split(": ")[1] if ": " in label_part else ""
        # Extract uniprot ID from "sp|XXXXX|NAME_HUMAN ..."
        uniprot = ""
        if "|" in seq_desc:
            pipe_parts = seq_desc.split("|")
            if len(pipe_parts) >= 2:
                uniprot = pipe_parts[1].replace("\\", "")
        # Extract gene name from "GN=XXXX"
        gene_name = ""
        gn_match = re.search(r'GN=(\S+)', seq_desc)
        if gn_match:
            gene_name = gn_match.group(1)

        seq_number_name[seq_num] = uniprot
        seq_descriptions[seq_num] = {
            'uniprot': uniprot,
            'gene': gene_name,
            'desc': seq_desc.replace("\\", ""),
        }

    # Parse cluster_2: sequence lengths
    seq_lengths = {}
    if len(parts) >= 3:
        cluster2 = parts[2]
        # Pattern: "(seq_num, length, 0)"
        length_nodes = re.findall(r'\((\d+),\s*(\d+),\s*0\)', cluster2)
        for sn, sl in length_nodes:
            seq_lengths[sn] = int(sl)

    # Parse cluster_3+: alignment blocks
    dict_bloc_infos = {}
    block_clusters = parts[3:] if len(parts) > 3 else []

    for bloc_idx, sub_cluster in enumerate(block_clusters):
        bloc_id = f"B{bloc_idx + 1}"
        blocs_infos = {}

        # Get content before closing brace
        sub_data = sub_cluster.split('}')[0]
        sub_nodes = sub_data.split(';"')[1:]

        for sub_node in sub_nodes:
            if '[label = "' not in sub_node:
                continue
            # Parse "(seq_num, start, length)" [label = "SEQUENCE"]
            coord_part = sub_node.split('" [label = "')[0].strip('"')
            seq_part = sub_node.split('" [label = "')[1].split('"]')[0]

            # Extract coordinates
            coord_clean = coord_part.strip('()')
            coord_fields = [c.strip() for c in coord_clean.split(',')]
            if len(coord_fields) < 3:
                continue

            seq_num = coord_fields[0]
            start = int(coord_fields[1])
            length = int(coord_fields[2])
            end = start + length - 1

            blocs_infos[seq_num] = {
                'start': start,
                'end': end,
                'length': length,
                'seq': seq_part,
            }

        if blocs_infos:
            dict_bloc_infos[bloc_id] = blocs_infos

    return dict_bloc_infos, seq_number_name, seq_lengths, seq_descriptions


def categorize_blocks(blocks, gene_a_seq_num, gene_b_seq_num):
    """Categorize each block by conservation pattern relative to the pair.

    Categories:
    - shared_with_family: Both A and B present, plus other family members
    - pair_exclusive: Only A and B present (no other family members)
    - specific_a: Only A present (not B)
    - specific_b: Only B present (not A)
    - family_only: Neither A nor B present
    """
    categorized = []
    for bloc_id, members in blocks.items():
        seq_nums = set(members.keys())
        a_in = gene_a_seq_num in seq_nums
        b_in = gene_b_seq_num in seq_nums

        if a_in and b_in:
            if len(seq_nums) == 2:
                category = "pair_exclusive"
            else:
                category = "shared_with_family"
        elif a_in and not b_in:
            category = "specific_a"
        elif b_in and not a_in:
            category = "specific_b"
        else:
            category = "family_only"

        categorized.append({
            'id': bloc_id,
            'category': category,
            'members': members,
            'n_seqs': len(seq_nums),
        })

    return categorized


def build_plma_json(pair_id, gene_a, gene_b, uniprot_a, uniprot_b, family_id):
    """Build the plma.json data for a pair."""
    # Find PLMA dot file
    fam_id_str = str(int(float(family_id)))
    dot_file = PLMA_DIR / f"family_{fam_id_str}_t10m1M20.dot"
    if not dot_file.exists():
        log(f"  WARNING: PLMA dot file not found: {dot_file}")
        return None

    log(f"  Parsing {dot_file.name}...")
    blocks, seq_names, seq_lengths, seq_descs = extract_block_sequences(str(dot_file))

    if not blocks:
        log(f"  WARNING: No blocks found in {dot_file.name}")
        return None

    # Find which sequence numbers correspond to gene A and gene B
    gene_a_seq = None
    gene_b_seq = None

    for sn, info in seq_descs.items():
        if info['uniprot'] == uniprot_a or info['gene'] == gene_a:
            gene_a_seq = sn
        if info['uniprot'] == uniprot_b or info['gene'] == gene_b:
            gene_b_seq = sn

    if gene_a_seq is None:
        log(f"  WARNING: Could not find {gene_a} ({uniprot_a}) in PLMA sequences")
        # Try gene name matching from descriptions
        for sn, info in seq_descs.items():
            desc_upper = info.get('desc', '').upper()
            if gene_a.upper() in desc_upper:
                gene_a_seq = sn
                break
    if gene_b_seq is None:
        log(f"  WARNING: Could not find {gene_b} ({uniprot_b}) in PLMA sequences")
        for sn, info in seq_descs.items():
            desc_upper = info.get('desc', '').upper()
            if gene_b.upper() in desc_upper:
                gene_b_seq = sn
                break

    if gene_a_seq is None or gene_b_seq is None:
        log(f"  ERROR: Could not map genes to sequences (A={gene_a_seq}, B={gene_b_seq})")
        log(f"  Available sequences: {seq_descs}")
        return None

    # Build sequence list
    sequences = []
    for sn in sorted(seq_descs.keys(), key=lambda x: int(x)):
        info = seq_descs[sn]
        sequences.append({
            'num': sn,
            'gene': info['gene'],
            'uniprot': info['uniprot'],
            'length': seq_lengths.get(sn, 0),
            'is_gene_a': sn == gene_a_seq,
            'is_gene_b': sn == gene_b_seq,
        })

    # Categorize blocks
    categorized = categorize_blocks(blocks, gene_a_seq, gene_b_seq)

    # Build block data for JSON, including AA sequences
    blocks_json = []
    for block in categorized:
        b = {
            'id': block['id'],
            'category': block['category'],
            'n_seqs': block['n_seqs'],
            'positions': {},
        }
        for sn, info in block['members'].items():
            b['positions'][sn] = {
                'start': info['start'],
                'end': info['end'],
                'length': info['length'],
                'seq': info.get('seq', ''),
            }
        blocks_json.append(b)

    # Sort blocks by PLMA order (B1, B2, B3...) which is the alignment order
    def block_sort_key(b):
        # Extract numeric part from block id (e.g., "B12" -> 12)
        m = re.match(r'B(\d+)', b['id'])
        return int(m.group(1)) if m else 0

    blocks_json.sort(key=block_sort_key)

    # Compute category summary
    summary = {
        'shared_with_family': 0,
        'pair_exclusive': 0,
        'specific_a': 0,
        'specific_b': 0,
        'family_only': 0,
    }
    for b in blocks_json:
        cat = b['category']
        # Sum aa count for the pair's sequences
        if cat in ('shared_with_family', 'pair_exclusive'):
            # Use average of A and B lengths
            la = b['positions'].get(gene_a_seq, {}).get('length', 0)
            lb = b['positions'].get(gene_b_seq, {}).get('length', 0)
            summary[cat] += (la + lb) // 2
        elif cat == 'specific_a':
            summary[cat] += b['positions'].get(gene_a_seq, {}).get('length', 0)
        elif cat == 'specific_b':
            summary[cat] += b['positions'].get(gene_b_seq, {}).get('length', 0)
        elif cat == 'family_only':
            lengths = [p['length'] for p in b['positions'].values()]
            summary[cat] += max(lengths) if lengths else 0

    plma_data = {
        'family_id': fam_id_str,
        'gene_a': gene_a,
        'gene_b': gene_b,
        'gene_a_seq': gene_a_seq,
        'gene_b_seq': gene_b_seq,
        'n_blocks': len(blocks_json),
        'n_sequences': len(sequences),
        'sequences': sequences,
        'blocks': blocks_json,
        'summary': summary,
    }

    return plma_data


def main():
    parser = argparse.ArgumentParser(description='Extract PLMA data for pairs')
    parser.add_argument('pair_id', nargs='?', help='Specific pair to process')
    args = parser.parse_args()

    # Load features CSV for family_id and uniprot mappings
    log("Loading features CSV...")
    df = pd.read_csv(FEATURES_CSV, usecols=[
        'A1', 'A2', 'family_id', 'uniprot_A1', 'uniprot_A2',
    ], low_memory=False)
    log(f"  Loaded {len(df)} pairs")

    # Load pairs list
    pairs_csv = INPUT_DIR / "pairs.csv"
    pairs_df = pd.read_csv(pairs_csv)

    if args.pair_id:
        pair_ids = [args.pair_id]
    else:
        pair_ids = [f"{r['A1']}_{r['A2']}" for _, r in pairs_df.iterrows()]

    log(f"Processing {len(pair_ids)} pairs...")

    success = 0
    errors = 0

    for pair_id in pair_ids:
        parts = pair_id.split('_', 1)
        if len(parts) != 2:
            log(f"  Skipping invalid pair_id: {pair_id}")
            errors += 1
            continue

        gene_a, gene_b = parts

        # Find in features CSV
        row = df[(df['A1'] == gene_a) & (df['A2'] == gene_b)]
        if row.empty:
            row = df[(df['A1'] == gene_b) & (df['A2'] == gene_a)]
        if row.empty:
            log(f"  {pair_id}: NOT FOUND in features CSV")
            errors += 1
            continue

        r = row.iloc[0]
        family_id = r['family_id']
        # Make sure uniprot mapping is correct for gene_a and gene_b
        if r['A1'] == gene_a:
            uniprot_a = r['uniprot_A1']
            uniprot_b = r['uniprot_A2']
        else:
            uniprot_a = r['uniprot_A2']
            uniprot_b = r['uniprot_A1']

        log(f"[{pair_id}] family={int(family_id)}, {gene_a}={uniprot_a}, {gene_b}={uniprot_b}")

        plma_data = build_plma_json(pair_id, gene_a, gene_b, uniprot_a, uniprot_b, family_id)

        if plma_data is None:
            errors += 1
            continue

        # Write plma.json
        out_dir = PAIRS_DIR / pair_id
        out_dir.mkdir(parents=True, exist_ok=True)
        out_file = out_dir / "plma.json"
        with open(out_file, 'w') as f:
            json.dump(plma_data, f)

        log(f"  Wrote {out_file.name}: {plma_data['n_blocks']} blocks, {plma_data['n_sequences']} sequences")
        success += 1

    log(f"Done! Success: {success}, Errors: {errors}")


if __name__ == "__main__":
    main()
