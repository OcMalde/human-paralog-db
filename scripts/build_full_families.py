#!/usr/bin/env python3
"""
Build full family data from the features CSV.

This script reads the complete paralog pair dataset and builds:
1. Family clusters (connected components of gene pairs)
2. Sequence identity data for all pairs
3. A mapping of all genes in each family relevant to pairs in our DB

Output: data/full_families.json
"""

import csv
import json
import sqlite3
from collections import defaultdict
from pathlib import Path

# Use local config
import sys
sys.path.insert(0, str(Path(__file__).parent))
from lib.config import FEATURES_CSV, DB_PATH, BASE_DIR, log


def build_union_find(pairs):
    """Build union-find structure to find connected components (families)."""
    parent = {}

    def find(x):
        if x not in parent:
            parent[x] = x
        if parent[x] != x:
            parent[x] = find(parent[x])  # Path compression
        return parent[x]

    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    # Union all pairs
    for g1, g2 in pairs:
        union(g1, g2)

    # Build family groups
    families = defaultdict(set)
    for gene in parent:
        root = find(gene)
        families[root].add(gene)

    return families, find


def main():
    log("Building full families from features CSV...")

    # Check if features CSV exists
    if not FEATURES_CSV.exists():
        log(f"ERROR: Features CSV not found: {FEATURES_CSV}")
        return

    # Read all pairs from features CSV
    log(f"Reading {FEATURES_CSV}...")
    all_pairs = []
    pair_identities = {}  # (g1, g2) -> max_sequence_identity

    with open(FEATURES_CSV, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            g1 = row.get('A1', '').strip()
            g2 = row.get('A2', '').strip()
            if g1 and g2:
                all_pairs.append((g1, g2))
                # Store max sequence identity for the pair
                try:
                    identity = float(row.get('max_sequence_identity', 0))
                except (ValueError, TypeError):
                    identity = 0
                # Store in both directions for easy lookup
                pair_identities[(g1, g2)] = identity
                pair_identities[(g2, g1)] = identity

    log(f"Loaded {len(all_pairs)} pairs from CSV")

    # Build family clusters
    log("Building family clusters...")
    families, find_root = build_union_find(all_pairs)
    log(f"Found {len(families)} distinct families")

    # Get genes from our database that we need families for
    log(f"Reading pairs from database: {DB_PATH}")
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    db_pairs = conn.execute("SELECT pair_id, gene_a, gene_b FROM pairs").fetchall()
    conn.close()

    db_genes = set()
    for row in db_pairs:
        db_genes.add(row['gene_a'])
        db_genes.add(row['gene_b'])

    log(f"Found {len(db_genes)} unique genes in database: {db_genes}")

    # For each gene in our DB, find its full family
    relevant_families = {}  # family_root -> { genes, pair_identities }

    for gene in db_genes:
        root = find_root(gene)
        if root not in relevant_families:
            family_genes = families[root]

            # Build identity matrix for this family
            family_identities = {}
            for g in family_genes:
                gene_identities = {}
                for other_g in family_genes:
                    if g != other_g:
                        identity = pair_identities.get((g, other_g))
                        if identity is not None:
                            gene_identities[other_g] = identity
                if gene_identities:
                    family_identities[g] = gene_identities

            relevant_families[root] = {
                'genes': sorted(family_genes),
                'identities': family_identities
            }

    log(f"Extracted {len(relevant_families)} relevant families")

    # Build output: gene -> family data mapping
    output = {
        'families': {},  # gene -> family_id
        'family_data': {}  # family_id -> { genes, identities }
    }

    family_id = 0
    for root, family_info in relevant_families.items():
        fid = f"family_{family_id}"
        output['family_data'][fid] = family_info

        # Map all genes in this family to the family_id
        for gene in family_info['genes']:
            output['families'][gene] = fid

        log(f"  {fid}: {len(family_info['genes'])} genes")
        family_id += 1

    # Write output
    output_path = BASE_DIR / "data" / "full_families.json"
    with open(output_path, 'w') as f:
        json.dump(output, f, separators=(',', ':'))

    log(f"Wrote {output_path}")
    log("Done!")


if __name__ == "__main__":
    main()