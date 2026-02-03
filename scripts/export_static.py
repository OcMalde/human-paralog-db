#!/usr/bin/env python3
"""
export_static.py

Exports database contents to static JSON files for GitHub Pages / CDN hosting.

Usage:
    python scripts/export_static.py                    # Export all pairs
    python scripts/export_static.py SMARCA2_SMARCA4    # Export specific pair(s)
    python scripts/export_static.py --clean            # Clean and re-export

Output structure:
    data/
    ├── index.json              # List of all pairs with metadata
    ├── family_index.json       # Gene -> pair IDs mapping
    ├── meta/
    │   └── boxplot_stats.json  # Global statistics for boxplots
    └── pairs/
        └── GENE1_GENE2/
            ├── report.json     # Full report data
            ├── summary.json    # Summary data
            └── pdb.json        # PDB variants (base64)
"""

import argparse
import base64
import gzip
import json
import shutil
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional

# ============= Configuration =============
BASE_DIR = Path(__file__).parent.parent
DB_PATH = BASE_DIR / "data" / "pairs.db"
STRUCTURES_DIR = BASE_DIR / "data" / "structures"
OUTPUT_DIR = BASE_DIR / "data" / "pairs"
META_DIR = BASE_DIR / "data" / "meta"


def log(msg: str):
    print(f"[export] {msg}")


def get_db() -> sqlite3.Connection:
    """Get database connection."""
    if not DB_PATH.exists():
        raise FileNotFoundError(f"Database not found at {DB_PATH}")
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    return conn


def export_pair(conn: sqlite3.Connection, pair_id: str, output_dir: Path) -> bool:
    """Export a single pair to JSON files."""

    # Get pair info
    pair = conn.execute(
        "SELECT * FROM pairs WHERE pair_id = ?", (pair_id,)
    ).fetchone()

    if not pair:
        log(f"  Pair {pair_id} not found in database")
        return False

    # Get report data
    report = conn.execute(
        "SELECT data_json, summary_json FROM report_data WHERE pair_id = ?", (pair_id,)
    ).fetchone()

    if not report or not report['data_json']:
        log(f"  No report data for {pair_id}")
        return False

    # Create output directory
    pair_dir = output_dir / pair_id
    pair_dir.mkdir(parents=True, exist_ok=True)

    # Parse and write report.json
    report_data = json.loads(report['data_json'])
    with open(pair_dir / "report.json", 'w') as f:
        json.dump(report_data, f, separators=(',', ':'))

    # Write summary.json
    if report['summary_json']:
        summary_data = json.loads(report['summary_json'])
        with open(pair_dir / "summary.json", 'w') as f:
            json.dump(summary_data, f, separators=(',', ':'))

    # Generate PDB variants
    pdb_variants = generate_pdb_variants(conn, pair_id, pair['acc_a'], pair['acc_b'])
    if pdb_variants:
        with open(pair_dir / "pdb.json", 'w') as f:
            json.dump(pdb_variants, f, separators=(',', ':'))

    return True


def generate_pdb_variants(conn: sqlite3.Connection, pair_id: str, acc_a: str, acc_b: str) -> Dict[str, str]:
    """Generate PDB variant files (base64 encoded)."""

    # Get protein structures
    prot_a = conn.execute("SELECT pdb_b64 FROM proteins WHERE acc = ?", (acc_a,)).fetchone()
    prot_b = conn.execute("SELECT pdb_b64 FROM proteins WHERE acc = ?", (acc_b,)).fetchone()

    if not prot_a or not prot_b:
        return {}

    pdb_a = base64.b64decode(prot_a['pdb_b64']).decode('utf-8', errors='ignore')
    pdb_b = base64.b64decode(prot_b['pdb_b64']).decode('utf-8', errors='ignore')

    # Get alignment transformation
    aln = conn.execute(
        "SELECT u_matrix, t_vector FROM alignments WHERE pair_id = ?", (pair_id,)
    ).fetchone()

    U = json.loads(aln['u_matrix']) if aln and aln['u_matrix'] else None
    T = json.loads(aln['t_vector']) if aln and aln['t_vector'] else None

    variants = {}

    # Full combined PDB
    combined = create_aligned_pdb(pdb_a, pdb_b, U, T, True, True)
    variants['pdb64_full'] = base64.b64encode(combined.encode()).decode()

    # Chain A only
    chain_a = create_aligned_pdb(pdb_a, pdb_b, U, T, True, False)
    variants['pdb64_a'] = base64.b64encode(chain_a.encode()).decode()

    # Chain B only
    chain_b = create_aligned_pdb(pdb_a, pdb_b, U, T, False, True)
    variants['pdb64_b'] = base64.b64encode(chain_b.encode()).decode()

    return variants


def create_aligned_pdb(pdb_a: str, pdb_b: str, U: List[float], T: List[float],
                       include_chain_a: bool = True, include_chain_b: bool = True) -> str:
    """Create combined PDB with chain B transformed."""
    lines = []

    # Chain A - unchanged
    if include_chain_a:
        for line in pdb_a.splitlines():
            if line.startswith(("ATOM  ", "HETATM")):
                lines.append(line[:21] + "A" + line[22:])

    # Chain B - transformed
    if include_chain_b:
        if U and T and len(U) == 9 and len(T) == 3:
            for line in pdb_b.splitlines():
                if not line.startswith(("ATOM  ", "HETATM")):
                    continue
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    X = U[0]*x + U[1]*y + U[2]*z + T[0]
                    Y = U[3]*x + U[4]*y + U[5]*z + T[1]
                    Z = U[6]*x + U[7]*y + U[8]*z + T[2]
                    new_line = line[:21] + "B" + line[22:30] + f"{X:8.3f}{Y:8.3f}{Z:8.3f}" + line[54:]
                    lines.append(new_line)
                except:
                    lines.append(line[:21] + "B" + line[22:])
        else:
            for line in pdb_b.splitlines():
                if line.startswith(("ATOM  ", "HETATM")):
                    lines.append(line[:21] + "B" + line[22:])

    lines.append("END")
    return "\n".join(lines)


def create_index_file(conn: sqlite3.Connection, output_path: Path):
    """Create index.json with all pairs and metadata."""

    pairs = conn.execute("""
        SELECT p.pair_id, p.gene_a, p.gene_b, p.acc_a, p.acc_b,
               a.tm_score, a.fident
        FROM pairs p
        LEFT JOIN alignments a ON p.pair_id = a.pair_id
        ORDER BY p.pair_id
    """).fetchall()

    index_data = []
    for p in pairs:
        index_data.append({
            "id": p['pair_id'],
            "geneA": p['gene_a'],
            "geneB": p['gene_b'],
            "accA": p['acc_a'],
            "accB": p['acc_b'],
            "tm": float(p['tm_score']) if p['tm_score'] else None,
            "fident": float(p['fident']) if p['fident'] else None,
        })

    with open(output_path, 'w') as f:
        json.dump(index_data, f, indent=2)

    log(f"Created index.json with {len(index_data)} pairs")


def create_family_index(conn: sqlite3.Connection, output_path: Path):
    """Create family_index.json mapping gene -> pair IDs."""

    pairs = conn.execute("SELECT pair_id, gene_a, gene_b FROM pairs").fetchall()

    family_index = defaultdict(list)
    for p in pairs:
        family_index[p['gene_a']].append(p['pair_id'])
        family_index[p['gene_b']].append(p['pair_id'])

    # Sort and deduplicate
    family_index = {gene: sorted(set(pair_ids)) for gene, pair_ids in family_index.items()}

    with open(output_path, 'w') as f:
        json.dump(family_index, f, separators=(',', ':'))

    log(f"Created family_index.json with {len(family_index)} genes")


def create_search_index(conn: sqlite3.Connection, output_path: Path):
    """Create search_index.json for client-side search."""

    pairs = conn.execute("""
        SELECT p.pair_id, p.gene_a, p.gene_b, p.acc_a, p.acc_b
        FROM pairs p
    """).fetchall()

    # Simple search index: list of searchable terms per pair
    search_data = []
    for p in pairs:
        search_data.append({
            "id": p['pair_id'],
            "terms": [
                p['gene_a'].lower(),
                p['gene_b'].lower(),
                p['acc_a'].lower() if p['acc_a'] else '',
                p['acc_b'].lower() if p['acc_b'] else '',
                p['pair_id'].lower(),
            ]
        })

    with open(output_path, 'w') as f:
        json.dump(search_data, f, separators=(',', ':'))

    log(f"Created search_index.json with {len(search_data)} entries")


def main():
    parser = argparse.ArgumentParser(description='Export database to static JSON files')
    parser.add_argument('pair_ids', nargs='*', help='Specific pair IDs to export')
    parser.add_argument('--clean', action='store_true', help='Remove existing exports first')
    args = parser.parse_args()

    if not DB_PATH.exists():
        log(f"ERROR: Database not found at {DB_PATH}")
        log("Run populate_db.py and generate_report_data.py first")
        sys.exit(1)

    conn = get_db()

    # Clean if requested
    if args.clean and OUTPUT_DIR.exists():
        log(f"Cleaning {OUTPUT_DIR}...")
        shutil.rmtree(OUTPUT_DIR)

    # Create output directories
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    META_DIR.mkdir(parents=True, exist_ok=True)

    # Determine pairs to export
    if args.pair_ids:
        pair_ids = args.pair_ids
    else:
        pair_ids = [r[0] for r in conn.execute("SELECT pair_id FROM pairs").fetchall()]

    log(f"Exporting {len(pair_ids)} pairs...")

    success = 0
    errors = 0

    for i, pair_id in enumerate(pair_ids):
        log(f"[{i+1}/{len(pair_ids)}] {pair_id}")
        try:
            if export_pair(conn, pair_id, OUTPUT_DIR):
                success += 1
            else:
                errors += 1
        except Exception as e:
            log(f"  ERROR: {e}")
            errors += 1

    # Create index files
    log("\nCreating index files...")
    create_index_file(conn, BASE_DIR / "data" / "index.json")
    create_family_index(conn, BASE_DIR / "data" / "family_index.json")
    create_search_index(conn, BASE_DIR / "data" / "search_index.json")

    conn.close()

    log(f"\n{'='*50}")
    log(f"Export complete!")
    log(f"  Success: {success}")
    log(f"  Errors: {errors}")
    log(f"  Output: {OUTPUT_DIR}")

    # Calculate total size
    total_size = sum(f.stat().st_size for f in OUTPUT_DIR.rglob('*') if f.is_file())
    log(f"  Total size: {total_size / 1024 / 1024:.2f} MB")


if __name__ == "__main__":
    main()