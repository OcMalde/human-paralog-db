#!/usr/bin/env python3
"""
populate_db.py

Populates the pairs.db database from a pairs input file.
Automatically maps gene names to UniProt IDs.
Fetches AlphaMissense PDB files and runs foldseek for alignments.

Usage:
    python scripts/populate_db.py                           # Use default input/pairs.csv
    python scripts/populate_db.py --pairs input/pairs.csv   # Specific pairs file
    python scripts/populate_db.py --all                     # All pairs from features CSV
    python scripts/populate_db.py SMARCA2_SMARCA4           # Specific pair(s)

Requirements:
    - foldseek installed and in PATH (or set FOLDSEEK env var)
    - Internet connection for fetching AM PDB files
    - Input data files in input/ directory
"""

import argparse
import base64
import json
import os
import re
import sqlite3
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Optional, Tuple, List, Dict

import pandas as pd
import requests

# ============= Configuration =============
BASE_DIR = Path(__file__).parent.parent
DB_PATH = Path(os.environ.get("PARALOG_DB", str(BASE_DIR / "data" / "pairs.db")))
STRUCTURES_DIR = BASE_DIR / "data" / "structures"
CACHE_DIR = BASE_DIR / "data" / "cache"

# Input data files (should be provided/linked)
INPUT_DIR = BASE_DIR / "input"
PAIRS_CSV = INPUT_DIR / "pairs.csv"
FEATURES_CSV = INPUT_DIR / "ens111_human_allFeatures.csv"
GENE_MAPPING_CSV = INPUT_DIR / "all_genes_ids.csv"

# Foldseek binary
FOLDSEEK = os.environ.get("FOLDSEEK", "foldseek")

# AlphaMissense PDB URL
AM_PDB_URL = "https://alphamissense.hegelab.org/pdb/AF-{acc}-F1-AM_v4.pdb"

# Annotation API URLs
TED_API_URL = "https://ted.cathdb.info/api/v1/uniprot/summary/{acc}?skip=0&limit=100"
EBI_API_URL = "https://www.ebi.ac.uk/proteins/api/features/{acc}"


def log(msg: str):
    print(f"[populate] {msg}")


# ============= Gene to UniProt Mapping =============

_gene_to_uniprot_cache: Optional[Dict[str, str]] = None


def load_gene_to_uniprot_mapping() -> Dict[str, str]:
    """Load gene name to canonical UniProt ID mapping."""
    global _gene_to_uniprot_cache
    if _gene_to_uniprot_cache is not None:
        return _gene_to_uniprot_cache

    _gene_to_uniprot_cache = {}

    if not GENE_MAPPING_CSV.exists():
        log(f"WARNING: Gene mapping file not found at {GENE_MAPPING_CSV}")
        return _gene_to_uniprot_cache

    log(f"Loading gene-to-UniProt mapping from {GENE_MAPPING_CSV}...")
    df = pd.read_csv(GENE_MAPPING_CSV)

    def get_canonical_uniprot(acc_str):
        """Extract canonical UniProt ID from string like '[A0A..., P51531]'."""
        if pd.isna(acc_str):
            return None
        # Extract all accessions using regex
        accs = re.findall(r'[A-Z][A-Z0-9]{5}', str(acc_str))
        # Prefer Swiss-Prot IDs (6 chars, start with P, Q, O, not A0A)
        for acc in accs:
            if acc[0] in 'PQO' and not acc.startswith('A0A'):
                return acc
        return accs[0] if accs else None

    for _, row in df.iterrows():
        gene_name = row.get('gene_name')
        if pd.notna(gene_name):
            uniprot = get_canonical_uniprot(row.get('UniProt_Acc'))
            if uniprot:
                _gene_to_uniprot_cache[str(gene_name).strip()] = uniprot

    log(f"  Loaded {len(_gene_to_uniprot_cache)} gene-to-UniProt mappings")
    return _gene_to_uniprot_cache


def get_uniprot_for_gene(gene_name: str) -> Optional[str]:
    """Get UniProt accession for a gene name."""
    mapping = load_gene_to_uniprot_mapping()
    return mapping.get(gene_name)


# ============= Database Setup =============

def create_tables(conn):
    """Create database tables."""
    conn.executescript("""
        CREATE TABLE IF NOT EXISTS pairs (
            pair_id TEXT PRIMARY KEY,
            gene_a TEXT,
            gene_b TEXT,
            acc_a TEXT,
            acc_b TEXT
        );

        CREATE TABLE IF NOT EXISTS alignments (
            pair_id TEXT PRIMARY KEY,
            qaln TEXT,
            taln TEXT,
            tm_score REAL,
            fident REAL,
            lddt REAL,
            u_matrix TEXT,
            t_vector TEXT,
            FOREIGN KEY (pair_id) REFERENCES pairs(pair_id)
        );

        CREATE TABLE IF NOT EXISTS proteins (
            acc TEXT PRIMARY KEY,
            length INTEGER,
            bfactors_json TEXT,
            pdb_b64 TEXT
        );

        CREATE TABLE IF NOT EXISTS annotations (
            acc TEXT,
            source TEXT,
            data_json TEXT,
            PRIMARY KEY (acc, source)
        );

        CREATE TABLE IF NOT EXISTS domain_alignments (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            pair_id TEXT,
            domain_a_name TEXT,
            domain_a_range TEXT,
            domain_b_name TEXT,
            domain_b_range TEXT,
            tm_score REAL,
            fident REAL,
            dam_percent REAL,
            pdb_b64 TEXT,
            FOREIGN KEY (pair_id) REFERENCES pairs(pair_id)
        );

        CREATE TABLE IF NOT EXISTS report_data (
            pair_id TEXT PRIMARY KEY,
            data_json TEXT,
            summary_json TEXT,
            generated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        );

        CREATE INDEX IF NOT EXISTS idx_pairs_gene_a ON pairs(gene_a);
        CREATE INDEX IF NOT EXISTS idx_pairs_gene_b ON pairs(gene_b);
        CREATE INDEX IF NOT EXISTS idx_annotations_acc ON annotations(acc);
    """)
    conn.commit()


# ============= Structure Fetching =============

def fetch_am_pdb(acc: str) -> Optional[Tuple[str, List[float], int]]:
    """
    Fetch AlphaMissense PDB from hegelab.
    Returns (pdb_content, bfactors, length) or None.
    """
    # Check cache first
    STRUCTURES_DIR.mkdir(parents=True, exist_ok=True)
    cache_path = STRUCTURES_DIR / f"AF-{acc}-F1-AM.pdb"

    if cache_path.exists():
        pdb_content = cache_path.read_text()
        log(f"    Loaded from cache: {cache_path.name}")
    else:
        url = AM_PDB_URL.format(acc=acc)
        try:
            log(f"    Fetching: {url}")
            resp = requests.get(url, timeout=30)
            if not resp.ok:
                log(f"    HTTP {resp.status_code} for {acc}")
                return None
            pdb_content = resp.text
            # Cache it
            cache_path.write_text(pdb_content)
            log(f"    Cached: {cache_path.name}")
        except Exception as e:
            log(f"    Failed to fetch PDB for {acc}: {e}")
            return None

    # Extract B-factors from CA atoms
    bfactors = []
    seen_residues = set()

    for line in pdb_content.split('\n'):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            try:
                res_num = int(line[22:26].strip())
                if res_num not in seen_residues:
                    seen_residues.add(res_num)
                    bf = float(line[60:66].strip())
                    bfactors.append(bf)
            except:
                pass

    return pdb_content, bfactors, len(bfactors)


# ============= Annotation Fetching =============

def fetch_ted_domains(acc: str) -> Optional[dict]:
    """Fetch TED domain annotations."""
    url = TED_API_URL.format(acc=acc)
    try:
        resp = requests.get(url, timeout=30, headers={"Accept": "application/json"})
        if resp.ok:
            return resp.json()
    except Exception as e:
        log(f"    TED fetch failed for {acc}: {e}")
    return None


def fetch_ebi_features(acc: str) -> Optional[dict]:
    """Fetch EBI/UniProt feature annotations."""
    url = EBI_API_URL.format(acc=acc)
    try:
        resp = requests.get(url, timeout=30, headers={"Accept": "application/json"})
        if resp.ok:
            return resp.json()
    except Exception as e:
        log(f"    EBI fetch failed for {acc}: {e}")
    return None


def store_annotation(conn, acc: str, source: str, data: Optional[dict]):
    """Store annotation in database."""
    if data is None:
        return
    conn.execute("""
        INSERT OR REPLACE INTO annotations (acc, source, data_json)
        VALUES (?, ?, ?)
    """, (acc, source, json.dumps(data)))


# ============= Foldseek Alignment =============

def run_foldseek(pdb1_path: Path, pdb2_path: Path) -> Optional[Tuple[str, str, float, float, List[float], List[float]]]:
    """
    Run foldseek easy-search to align two structures.
    Returns (qaln, taln, tm_score, fident, U_matrix, T_vector) or None.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        out_tsv = Path(tmpdir) / "result.tsv"
        tmp_db = Path(tmpdir) / "tmp"

        fields = "query,target,fident,alntmscore,qaln,taln,u,t"

        args = [
            FOLDSEEK, "easy-search",
            str(pdb1_path), str(pdb2_path),
            str(out_tsv), str(tmp_db),
            "-a",
            "--max-seqs", "1",
            "--format-output", fields,
        ]

        try:
            result = subprocess.run(args, capture_output=True, timeout=120)
            if result.returncode != 0:
                log(f"    Foldseek error: {result.stderr.decode()[:200]}")
                return None

            if not out_tsv.exists() or out_tsv.stat().st_size == 0:
                log(f"    Foldseek produced no output")
                return None

            # Parse output
            line = out_tsv.read_text().strip().split('\n')[0]
            cols = line.split('\t')

            if len(cols) < 6:
                return None

            fident = float(cols[2]) if cols[2] else None
            tm_score = float(cols[3]) if cols[3] else None
            qaln = cols[4]
            taln = cols[5]

            # Parse rotation matrix and translation vector
            U = [float(x) for x in cols[6].split(",")] if len(cols) > 6 and cols[6] else None
            T = [float(x) for x in cols[7].split(",")] if len(cols) > 7 and cols[7] else None

            return qaln, taln, tm_score, fident, U, T

        except subprocess.TimeoutExpired:
            log("    Foldseek timeout")
            return None
        except Exception as e:
            log(f"    Foldseek error: {e}")
            return None


# ============= Main Population =============

def load_pairs_from_csv(pairs_csv: Path) -> List[Tuple[str, str]]:
    """Load pairs from input CSV file."""
    df = pd.read_csv(pairs_csv)

    # Find gene columns
    a1_col = next((c for c in ['A1', 'gene_a', 'Gene1'] if c in df.columns), None)
    a2_col = next((c for c in ['A2', 'gene_b', 'Gene2'] if c in df.columns), None)

    if not a1_col or not a2_col:
        raise ValueError(f"Could not find gene columns in {pairs_csv}. Found: {list(df.columns)}")

    pairs = []
    for _, row in df.iterrows():
        g1 = str(row[a1_col]).strip()
        g2 = str(row[a2_col]).strip()
        if g1 and g2 and g1 != 'nan' and g2 != 'nan':
            pairs.append((g1, g2))

    return pairs


def populate_pair(conn, gene_a: str, gene_b: str, proteins_cache: dict) -> bool:
    """Populate database for a single pair. Returns True if successful."""

    # Get UniProt accessions
    acc_a = get_uniprot_for_gene(gene_a)
    acc_b = get_uniprot_for_gene(gene_b)

    if not acc_a:
        log(f"  WARNING: No UniProt ID found for {gene_a}")
        return False
    if not acc_b:
        log(f"  WARNING: No UniProt ID found for {gene_b}")
        return False

    pair_id = f"{gene_a}_{gene_b}"
    log(f"  {pair_id}: {gene_a} ({acc_a}) vs {gene_b} ({acc_b})")

    # Insert pair
    conn.execute("""
        INSERT OR REPLACE INTO pairs (pair_id, gene_a, gene_b, acc_a, acc_b)
        VALUES (?, ?, ?, ?, ?)
    """, (pair_id, gene_a, gene_b, acc_a, acc_b))

    # Fetch/cache proteins and structures
    pdb_paths = {}
    for acc, gene in [(acc_a, gene_a), (acc_b, gene_b)]:
        if acc not in proteins_cache:
            log(f"  Fetching structure for {acc} ({gene})...")
            result = fetch_am_pdb(acc)
            if result:
                pdb_content, bfactors, length = result
                pdb_path = STRUCTURES_DIR / f"AF-{acc}-F1-AM.pdb"
                proteins_cache[acc] = (pdb_path, bfactors, length)

                # Store in database
                conn.execute("""
                    INSERT OR REPLACE INTO proteins (acc, length, bfactors_json, pdb_b64)
                    VALUES (?, ?, ?, ?)
                """, (acc, length, json.dumps(bfactors), base64.b64encode(pdb_content.encode()).decode()))

                # Fetch and store annotations
                log(f"    Fetching annotations for {acc}...")
                ted_data = fetch_ted_domains(acc)
                ebi_data = fetch_ebi_features(acc)
                store_annotation(conn, acc, "ted", ted_data)
                store_annotation(conn, acc, "ebi", ebi_data)
            else:
                proteins_cache[acc] = None
                log(f"    Failed to fetch {acc}")

        if proteins_cache.get(acc):
            pdb_paths[acc] = proteins_cache[acc][0]

    # Run foldseek alignment if both structures available
    if acc_a in pdb_paths and acc_b in pdb_paths:
        log(f"  Running foldseek alignment...")
        aln_result = run_foldseek(pdb_paths[acc_a], pdb_paths[acc_b])

        if aln_result:
            qaln, taln, tm_score, fident, U, T = aln_result
            conn.execute("""
                INSERT OR REPLACE INTO alignments (pair_id, qaln, taln, tm_score, fident, u_matrix, t_vector)
                VALUES (?, ?, ?, ?, ?, ?, ?)
            """, (pair_id, qaln, taln, tm_score, fident,
                  json.dumps(U) if U else None,
                  json.dumps(T) if T else None))
            tm_str = f"{tm_score:.3f}" if tm_score else "N/A"
            fid_str = f"{fident:.3f}" if fident else "N/A"
            log(f"    Alignment: {len(qaln)} cols, TM={tm_str}, fident={fid_str}")
            return True
        else:
            log(f"    Alignment failed")
            return False
    else:
        log(f"    Skipping alignment (missing structures)")
        return False


def main():
    parser = argparse.ArgumentParser(description='Populate pairs database')
    parser.add_argument('--pairs', type=Path, default=PAIRS_CSV, help='Input pairs CSV file')
    parser.add_argument('--all', action='store_true', help='Process all pairs from features CSV')
    parser.add_argument('pair_ids', nargs='*', help='Specific pair IDs (GENE1_GENE2)')
    args = parser.parse_args()

    # Check for foldseek
    try:
        result = subprocess.run([FOLDSEEK, "--help"], capture_output=True, timeout=5)
        log(f"Found foldseek: {FOLDSEEK}")
    except Exception as e:
        log(f"ERROR: foldseek not found. Install it or set FOLDSEEK env var.")
        log(f"  Error: {e}")
        sys.exit(1)

    # Create directories
    DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    STRUCTURES_DIR.mkdir(parents=True, exist_ok=True)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    # Determine pairs to process
    if args.pair_ids:
        # Specific pairs from command line
        pairs = []
        for pid in args.pair_ids:
            if '_' in pid:
                g1, g2 = pid.split('_', 1)
                pairs.append((g1, g2))
            else:
                log(f"WARNING: Invalid pair format '{pid}', expected GENE1_GENE2")
        log(f"Processing {len(pairs)} pairs from command line")
    elif args.all:
        # All pairs from features CSV
        if not FEATURES_CSV.exists():
            log(f"ERROR: Features CSV not found at {FEATURES_CSV}")
            sys.exit(1)
        log(f"Loading all pairs from {FEATURES_CSV}...")
        df = pd.read_csv(FEATURES_CSV, usecols=['A1', 'A2'], low_memory=False)
        pairs = [(row['A1'], row['A2']) for _, row in df.iterrows() if pd.notna(row['A1']) and pd.notna(row['A2'])]
        log(f"  Found {len(pairs)} pairs")
    else:
        # Default: pairs from input file
        if not args.pairs.exists():
            log(f"ERROR: Pairs file not found at {args.pairs}")
            log(f"Create it or use --all to process all pairs")
            sys.exit(1)
        pairs = load_pairs_from_csv(args.pairs)
        log(f"Loaded {len(pairs)} pairs from {args.pairs}")

    if not pairs:
        log("No pairs to process!")
        sys.exit(1)

    # Connect to database
    conn = sqlite3.connect(DB_PATH)

    log("Creating tables...")
    create_tables(conn)

    log(f"Processing {len(pairs)} pairs...")
    proteins_cache = {}
    success = 0
    errors = 0

    for i, (gene_a, gene_b) in enumerate(pairs):
        log(f"\n[{i+1}/{len(pairs)}] {gene_a}_{gene_b}")
        try:
            if populate_pair(conn, gene_a, gene_b, proteins_cache):
                success += 1
            else:
                errors += 1
        except Exception as e:
            log(f"  ERROR: {e}")
            errors += 1

        # Commit every pair
        conn.commit()

    # Summary
    pairs_count = conn.execute("SELECT COUNT(*) FROM pairs").fetchone()[0]
    proteins_count = conn.execute("SELECT COUNT(*) FROM proteins").fetchone()[0]
    alignments_count = conn.execute("SELECT COUNT(*) FROM alignments").fetchone()[0]

    log(f"\n{'='*50}")
    log(f"Database summary:")
    log(f"  Pairs: {pairs_count}")
    log(f"  Proteins: {proteins_count}")
    log(f"  Alignments: {alignments_count}")
    log(f"  Success: {success}, Errors: {errors}")
    log(f"  Database: {DB_PATH}")
    if DB_PATH.exists():
        log(f"  Size: {DB_PATH.stat().st_size / 1024 / 1024:.1f} MB")

    conn.close()
    log("Done!")


if __name__ == "__main__":
    main()