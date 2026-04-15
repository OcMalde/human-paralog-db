#!/usr/bin/env python3
"""
export_static.py

Exports database contents to static JSON files for GitHub Pages / CDN hosting.

Usage:
    python scripts/export_static.py                    # Export all pairs + genes
    python scripts/export_static.py SMARCA2_SMARCA4    # Export specific pair(s)
    python scripts/export_static.py --clean            # Clean and re-export

Output structure:
    data/
    ├── index.json              # List of all pairs with metadata
    ├── family_index.json       # Gene -> pair IDs mapping
    ├── meta/
    │   └── boxplot_stats.json  # Global statistics for boxplots
    ├── genes/
    │   └── {GENE}.json         # Per-gene data (PDB, AM, domains, drugs, ...)
    └── pairs/
        └── GENE1_GENE2/
            ├── report.json     # Pair-specific report data (alignment, rects, domPairs)
            ├── summary.json    # Pair-level summary (conservation, SL, sim search, ...)
            └── pdb.json        # Combined aligned PDB only (pdb64_full)

Per-gene data (in data/genes/{GENE}.json, NOT duplicated per pair):
    - af_pdb64          Canonical AlphaFold PDB (unaligned)
    - pdbeComplexes     PDBe experimental structures with coordB64
    - bfactors          AlphaMissense per-residue scores
    - plddt             pLDDT per-residue scores
    - domains           All domain annotations (EBI, TED, disorder, cavities, DrugCLIP)
    - amMatrixRects     AM saturation matrix rects by normalisation mode
    - seq               Canonical protein sequence
    - info              Gene metadata (description, known drugs, essentiality, ...)
"""

import argparse
import base64
import gzip
import json
import shutil
import sqlite3
import sys
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional

# ============= Configuration =============
BASE_DIR = Path(__file__).parent.parent
DB_PATH = BASE_DIR / "data" / "pairs.db"
STRUCTURES_DIR = BASE_DIR / "data" / "structures"
OUTPUT_DIR = BASE_DIR / "data" / "pairs"
GENES_DIR = BASE_DIR / "data" / "genes"
META_DIR = BASE_DIR / "data" / "meta"

# Fields moved from report.json to per-gene files
PER_GENE_REPORT_FIELDS = {
    'bfactorsA', 'bfactorsB',
    'plddtA', 'plddtB',
    'domainsA', 'domainsB',
    'amMatrixA_rectsByMode', 'amMatrixB_rectsByMode',
    'seq1', 'seq2',
    'pdbeComplexes',
    'amModes',
}

# Fields moved from summary.json to per-gene files
PER_GENE_SUMMARY_FIELDS = {'gene1', 'gene2'}


def log(msg: str):
    print(f"[export] {msg}")


def get_db() -> sqlite3.Connection:
    """Get database connection."""
    if not DB_PATH.exists():
        raise FileNotFoundError(f"Database not found at {DB_PATH}")
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    return conn


# ============= Per-gene export =============

def export_gene(conn: sqlite3.Connection, gene: str, acc: str, genes_dir: Path) -> bool:
    """Export per-gene data to data/genes/{GENE}.json.

    Reads the canonical AF structure from proteins table and per-residue /
    domain / AM data from the first pair in report_data that contains this gene.
    """
    # Canonical AlphaFold structure (unaligned, chain A)
    prot = conn.execute("SELECT pdb_b64 FROM proteins WHERE acc = ?", (acc,)).fetchone()
    af_pdb64 = prot['pdb_b64'] if prot else None

    # All pairs for this gene, ordered for stability
    pairs = conn.execute("""
        SELECT pair_id, gene_a, gene_b, acc_a, acc_b FROM pairs
        WHERE gene_a = ? OR gene_b = ?
        ORDER BY pair_id
    """, (gene, gene)).fetchall()

    gene_out = {
        'gene': gene,
        'acc': acc,
        'af_pdb64': af_pdb64,
        'bfactors': [],
        'plddt': [],
        'domains': [],
        'amMatrixRects': {},
        'seq': '',
        'pdbeComplexes': [],
        'info': {},
    }

    for pair in pairs:
        report = conn.execute(
            "SELECT data_json, summary_json FROM report_data WHERE pair_id = ?",
            (pair['pair_id'],)
        ).fetchone()

        if not report or not report['data_json']:
            continue

        data = json.loads(report['data_json'])
        side = 'A' if pair['gene_a'] == gene else 'B'
        seq_key = 'seq1' if side == 'A' else 'seq2'

        gene_out['bfactors'] = data.get(f'bfactors{side}', [])
        gene_out['plddt'] = data.get(f'plddt{side}', [])
        gene_out['domains'] = data.get(f'domains{side}', [])
        gene_out['amMatrixRects'] = data.get(f'amMatrix{side}_rectsByMode', {})
        gene_out['seq'] = data.get(seq_key, '')
        gene_out['pdbeComplexes'] = [
            c for c in data.get('pdbeComplexes', [])
            if c.get('source_acc') == acc
        ]

        if report['summary_json']:
            summary = json.loads(report['summary_json'])
            gene_key = 'gene1' if side == 'A' else 'gene2'
            gene_out['info'] = summary.get(gene_key, {})

        break  # First pair with data is sufficient

    genes_dir.mkdir(parents=True, exist_ok=True)
    with open(genes_dir / f"{gene}.json", 'w') as f:
        json.dump(gene_out, f, separators=(',', ':'))

    return True


def export_all_genes(conn: sqlite3.Connection, genes_dir: Path):
    """Export a gene file for every unique gene in the database."""
    genes = conn.execute(
        "SELECT gene_a AS gene, acc_a AS acc FROM pairs "
        "UNION SELECT gene_b, acc_b FROM pairs ORDER BY gene"
    ).fetchall()

    log(f"Exporting {len(genes)} gene files...")
    ok = 0
    for g in genes:
        try:
            if export_gene(conn, g['gene'], g['acc'], genes_dir):
                ok += 1
        except Exception as e:
            log(f"  ERROR exporting gene {g['gene']}: {e}")

    log(f"Gene export complete: {ok}/{len(genes)} genes")


# ============= Per-pair export =============

def export_pair(conn: sqlite3.Connection, pair_id: str, output_dir: Path) -> bool:
    """Export a single pair to JSON files (per-gene fields stripped)."""

    pair = conn.execute(
        "SELECT * FROM pairs WHERE pair_id = ?", (pair_id,)
    ).fetchone()

    if not pair:
        log(f"  Pair {pair_id} not found in database")
        return False

    report = conn.execute(
        "SELECT data_json, summary_json FROM report_data WHERE pair_id = ?", (pair_id,)
    ).fetchone()

    if not report or not report['data_json']:
        log(f"  No report data for {pair_id}")
        return False

    pair_dir = output_dir / pair_id
    pair_dir.mkdir(parents=True, exist_ok=True)

    # report.json — strip per-gene fields (now in data/genes/{GENE}.json)
    report_data = json.loads(report['data_json'])
    for field in PER_GENE_REPORT_FIELDS:
        report_data.pop(field, None)
    with open(pair_dir / "report.json", 'w') as f:
        json.dump(report_data, f, separators=(',', ':'))

    # summary.json — strip per-gene fields
    if report['summary_json']:
        summary_data = json.loads(report['summary_json'])
        for field in PER_GENE_SUMMARY_FIELDS:
            summary_data.pop(field, None)
        with open(pair_dir / "summary.json", 'w') as f:
            json.dump(summary_data, f, separators=(',', ':'))

    # pdb.json — only the combined aligned structure (pdb64_a/b moved to gene files)
    pdb_variants = generate_pdb_variants(conn, pair_id, pair['acc_a'], pair['acc_b'])
    if pdb_variants:
        with open(pair_dir / "pdb.json", 'w') as f:
            json.dump(pdb_variants, f, separators=(',', ':'))

    return True


def generate_pdb_variants(conn: sqlite3.Connection, pair_id: str, acc_a: str, acc_b: str) -> Dict[str, str]:
    """Generate pdb.json — only pdb64_full (combined aligned structure).

    pdb64_a and pdb64_b (individual gene structures) are now stored in
    data/genes/{GENE}.json as af_pdb64 (canonical, unaligned from proteins table).
    """
    prot_a = conn.execute("SELECT pdb_b64 FROM proteins WHERE acc = ?", (acc_a,)).fetchone()
    prot_b = conn.execute("SELECT pdb_b64 FROM proteins WHERE acc = ?", (acc_b,)).fetchone()

    if not prot_a or not prot_b:
        return {}

    pdb_a = base64.b64decode(prot_a['pdb_b64']).decode('utf-8', errors='ignore')
    pdb_b = base64.b64decode(prot_b['pdb_b64']).decode('utf-8', errors='ignore')

    aln = conn.execute(
        "SELECT u_matrix, t_vector FROM alignments WHERE pair_id = ?", (pair_id,)
    ).fetchone()

    U = json.loads(aln['u_matrix']) if aln and aln['u_matrix'] else None
    T = json.loads(aln['t_vector']) if aln and aln['t_vector'] else None

    combined = create_aligned_pdb(pdb_a, pdb_b, U, T, True, True)
    return {'pdb64_full': base64.b64encode(combined.encode()).decode()}


def create_aligned_pdb(pdb_a: str, pdb_b: str, U: List[float], T: List[float],
                       include_chain_a: bool = True, include_chain_b: bool = True) -> str:
    """Create combined PDB with chain B transformed."""
    lines = []

    if include_chain_a:
        for line in pdb_a.splitlines():
            if line.startswith(("ATOM  ", "HETATM")):
                lines.append(line[:21] + "A" + line[22:])

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


# ============= Index files =============

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


def create_family_graph(conn: sqlite3.Connection, data_dir: Path):
    """Create family_graph.json for the homepage bubble chart.

    Groups genes into families (connected components of pairs), computes a
    human-readable label for each family, and serialises genes + pairs so the
    frontend can render both the bubble chart and the per-family force network.
    Reads gene metadata from data/genes/{GENE}.json (written by export_all_genes).
    """

    rows = conn.execute("""
        SELECT p.pair_id, p.gene_a, p.gene_b, a.fident
        FROM pairs p
        LEFT JOIN alignments a ON p.pair_id = a.pair_id
    """).fetchall()

    adjacency = defaultdict(set)
    pair_map = {}
    for r in rows:
        adjacency[r['gene_a']].add(r['gene_b'])
        adjacency[r['gene_b']].add(r['gene_a'])
        pair_map[r['pair_id']] = {
            'id': r['pair_id'],
            'gA': r['gene_a'],
            'gB': r['gene_b'],
            'fi': round(float(r['fident']), 4) if r['fident'] else None,
        }

    # Connected components → families
    visited = set()
    family_groups = []
    for gene in sorted(adjacency.keys()):
        if gene in visited:
            continue
        group = set()
        queue = [gene]
        while queue:
            g = queue.pop()
            if g in visited:
                continue
            visited.add(g)
            group.add(g)
            queue.extend(adjacency.get(g, set()) - visited)
        family_groups.append(sorted(group))

    # Load per-gene info — prefer gene files, fall back to pair summary files
    genes_dir = data_dir / 'genes'
    gene_info: dict = {}

    for pair_id in pair_map:
        for gkey in ('gA', 'gB'):
            sym = pair_map[pair_id][gkey]
            if sym in gene_info:
                continue
            # Try gene file first
            gene_file = genes_dir / f"{sym}.json"
            if gene_file.exists():
                try:
                    with open(gene_file) as f:
                        gd = json.load(f)
                    info = gd.get('info', {}) or {}
                    desc = info.get('description', {}) or {}
                    dl = len((desc.get('function') or '') + (desc.get('name') or ''))
                    nd = len(info.get('known_drugs') or [])
                    gene_info[sym] = {'dl': dl, 'nd': nd}
                    continue
                except Exception:
                    pass
            # Fall back to summary.json
            summary_path = data_dir / 'pairs' / pair_id / 'summary.json'
            if summary_path.exists():
                try:
                    with open(summary_path) as f:
                        s = json.load(f)
                    for sk in ('gene1', 'gene2'):
                        g = s.get(sk, {})
                        gsym = g.get('symbol')
                        if gsym and gsym not in gene_info:
                            desc = g.get('description', {}) or {}
                            dl = len((desc.get('function') or '') + (desc.get('name') or ''))
                            nd = len(g.get('known_drugs') or [])
                            gene_info[gsym] = {'dl': dl, 'nd': nd}
                except Exception:
                    pass

    def make_label(genes: list) -> str:
        prefixes = [re.match(r'^([A-Za-z]+)', g).group(1) for g in genes
                    if re.match(r'^([A-Za-z]+)', g)]
        counts = Counter(prefixes)
        common = [p for p, c in counts.most_common() if c > 1]
        covered = {g for g in genes if any(g.startswith(p) for p in common)}
        uncovered = sorted(
            [g for g in genes if g not in covered],
            key=lambda g: (-(gene_info.get(g, {}).get('dl', 0)),
                           -(gene_info.get(g, {}).get('nd', 0)))
        )
        parts = common[:5] + uncovered[:max(0, 5 - len(common))]
        return '/'.join(parts[:5]) if parts else genes[0]

    families = []
    for i, genes in enumerate(sorted(family_groups, key=lambda g: -len(g))):
        gene_set = set(genes)
        fam_pairs = [pd for pd in pair_map.values()
                     if pd['gA'] in gene_set and pd['gB'] in gene_set]
        fam_pairs.sort(key=lambda p: -(p['fi'] or 0))

        families.append({
            'id': i,
            'label': make_label(genes),
            'size': len(genes),
            'nPairs': len(fam_pairs),
            'genes': [{'s': g,
                       'dl': gene_info.get(g, {}).get('dl', 0),
                       'nd': gene_info.get(g, {}).get('nd', 0)}
                      for g in genes],
            'pairs': fam_pairs,
        })

    out_path = data_dir / 'family_graph.json'
    with open(out_path, 'w') as f:
        json.dump({'families': families}, f, separators=(',', ':'))
    log(f"Created family_graph.json with {len(families)} families")


# ============= Main =============

def main():
    parser = argparse.ArgumentParser(description='Export database to static JSON files')
    parser.add_argument('pair_ids', nargs='*', help='Specific pair IDs to export')
    parser.add_argument('--clean', action='store_true', help='Remove existing exports first')
    parser.add_argument('--skip-genes', action='store_true', help='Skip gene file export')
    args = parser.parse_args()

    if not DB_PATH.exists():
        log(f"ERROR: Database not found at {DB_PATH}")
        log("Run populate_db.py and generate_report_data.py first")
        sys.exit(1)

    conn = get_db()

    if args.clean:
        for d in [OUTPUT_DIR, GENES_DIR]:
            if d.exists():
                log(f"Cleaning {d}...")
                shutil.rmtree(d)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    GENES_DIR.mkdir(parents=True, exist_ok=True)
    META_DIR.mkdir(parents=True, exist_ok=True)

    # Export per-gene files first (pair exports depend on gene data being available
    # for family_graph; and gene export reads from DB, not from pair files)
    if not args.skip_genes:
        export_all_genes(conn, GENES_DIR)

    # Determine pairs to export
    if args.pair_ids:
        pair_ids = args.pair_ids
    else:
        pair_ids = [r[0] for r in conn.execute("SELECT pair_id FROM pairs").fetchall()]

    log(f"\nExporting {len(pair_ids)} pairs...")

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

    # Index files
    log("\nCreating index files...")
    create_index_file(conn, BASE_DIR / "data" / "index.json")
    create_family_index(conn, BASE_DIR / "data" / "family_index.json")
    create_search_index(conn, BASE_DIR / "data" / "search_index.json")
    create_family_graph(conn, BASE_DIR / "data")

    conn.close()

    log(f"\n{'='*50}")
    log(f"Export complete!")
    log(f"  Pairs: {success} ok, {errors} errors")
    log(f"  Output: {OUTPUT_DIR}")

    pairs_size = sum(f.stat().st_size for f in OUTPUT_DIR.rglob('*') if f.is_file())
    genes_size = sum(f.stat().st_size for f in GENES_DIR.rglob('*') if f.is_file())
    log(f"  data/pairs/ size: {pairs_size / 1024 / 1024:.1f} MB")
    log(f"  data/genes/ size: {genes_size / 1024 / 1024:.1f} MB")
    log(f"  Total data size:  {(pairs_size + genes_size) / 1024 / 1024:.1f} MB")


if __name__ == "__main__":
    main()
