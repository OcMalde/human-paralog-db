#!/usr/bin/env python3
"""
generate_local_report.py

Generate a self-contained HTML report for a specific paralog pair.
Reuses the DB codebase (report.html, app.js) so any future UI changes
are automatically picked up on next run.

Usage:
    python scripts/generate_local_report.py CDK4_CDK6
    python scripts/generate_local_report.py CDK4_CDK6 -o ~/reports/CDK4_CDK6.html
    python scripts/generate_local_report.py CDK4_CDK6 --use-cached

    # Full pipeline: populate DB + generate report data + PLMA + HTML
    # (no need to edit pairs.csv — pair stays local only)
    python scripts/generate_local_report.py ADAMTS12_ADAMTS7 --full
"""

import argparse
import json
import os
import re
import sqlite3
import subprocess
import sys
import urllib.request
from pathlib import Path

# Project paths
PROJECT_DIR = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = PROJECT_DIR / "scripts"
DATA_DIR = PROJECT_DIR / "data"
DB_PATH = DATA_DIR / "pairs.db"
LIB_CACHE_DIR = SCRIPTS_DIR / ".lib_cache"
OUTPUT_DIR = PROJECT_DIR / "output"

# CDN libraries to download and inline
CDN_LIBS = {
    "molstar.css": "https://unpkg.com/molstar@3.43.0/build/viewer/molstar.css",
    "chart.umd.min.js": "https://cdn.jsdelivr.net/npm/chart.js@4.4.1/dist/chart.umd.min.js",
    "openchemlib-full.js": "https://unpkg.com/openchemlib@8.13.0/dist/openchemlib-full.js",
    "molstar.js": "https://unpkg.com/molstar@3.43.0/build/viewer/molstar.js",
}


def log(msg):
    print(f"[local-report] {msg}", flush=True)


def download_and_cache(name, url):
    """Download a CDN library to local cache if not already cached."""
    cache_path = LIB_CACHE_DIR / name
    if cache_path.exists():
        return cache_path.read_text(encoding="utf-8", errors="replace")

    LIB_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    log(f"  Downloading {name}...")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        content = resp.read().decode("utf-8", errors="replace")
    cache_path.write_text(content, encoding="utf-8")
    return content


def get_db():
    DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(DB_PATH))
    conn.row_factory = sqlite3.Row
    return conn


def run_full_pipeline(pair_id, conn):
    """Run the full pipeline: populate DB + generate report data + extract PLMA.

    This allows generating local reports for pairs not in pairs.csv,
    keeping them out of the online DB.
    """
    gene_a, gene_b = pair_id.split("_", 1)

    # Step 1: Populate DB (structure fetch + foldseek alignment)
    pair = conn.execute("SELECT * FROM pairs WHERE pair_id = ?", (pair_id,)).fetchone()
    if not pair:
        log(f"Step 1/3: Populating DB for {gene_a} vs {gene_b}...")
        script = SCRIPTS_DIR / "populate_db.py"
        result = subprocess.run(
            [sys.executable, str(script), pair_id],
            cwd=str(PROJECT_DIR),
        )
        if result.returncode != 0:
            raise RuntimeError(f"populate_db.py failed for {pair_id}")
        # Reconnect to pick up new data (populate_db uses its own connection)
        conn.close()
        conn = get_db()
    else:
        log("Step 1/3: Pair already in DB (skipping populate)")

    # Step 2: Generate report data
    report = conn.execute(
        "SELECT data_json FROM report_data WHERE pair_id = ?", (pair_id,)
    ).fetchone()
    if not report or not report["data_json"]:
        log("Step 2/3: Generating report data...")
        sys.path.insert(0, str(SCRIPTS_DIR))
        from generate_report_data import generate_report_data as gen_report

        data_obj, summary_obj = gen_report(pair_id, conn)
        conn.execute(
            """INSERT OR REPLACE INTO report_data (pair_id, data_json, summary_json)
               VALUES (?, ?, ?)""",
            (pair_id, json.dumps(data_obj), json.dumps(summary_obj)),
        )
        conn.commit()
    else:
        log("Step 2/3: Report data already in DB (skipping)")

    # Step 3: Extract PLMA data
    plma_path = DATA_DIR / "pairs" / pair_id / "plma.json"
    if not plma_path.exists():
        log("Step 3/3: Extracting PLMA data...")
        script = SCRIPTS_DIR / "extract_plma.py"
        result = subprocess.run(
            [sys.executable, str(script), pair_id],
            cwd=str(PROJECT_DIR),
        )
        if result.returncode != 0:
            log("  WARNING: PLMA extraction failed (family data may not be available)")
    else:
        log("Step 3/3: PLMA data already exists (skipping)")

    return conn


def load_pair_data(conn, pair_id):
    """Load report data, summary, and PDB variants from the database."""
    # Check if pair exists
    pair = conn.execute("SELECT * FROM pairs WHERE pair_id = ?", (pair_id,)).fetchone()
    if not pair:
        raise ValueError(f"Pair {pair_id} not found in database")

    # Check for generated report data
    report = conn.execute(
        "SELECT data_json, summary_json FROM report_data WHERE pair_id = ?",
        (pair_id,),
    ).fetchone()

    if not report or not report["data_json"]:
        raise ValueError(
            f"No report data for {pair_id}. Run: python scripts/generate_report_data.py {pair_id}"
        )

    report_data = json.loads(report["data_json"])
    summary_data = json.loads(report["summary_json"]) if report["summary_json"] else {}

    # Generate PDB variants (reuse export_static logic)
    sys.path.insert(0, str(SCRIPTS_DIR))
    from export_static import generate_pdb_variants

    pdb_variants = generate_pdb_variants(conn, pair_id, pair["acc_a"], pair["acc_b"])

    return report_data, summary_data, pdb_variants


def load_plma_data(pair_id):
    """Load PLMA data from the exported JSON or generate it."""
    plma_path = DATA_DIR / "pairs" / pair_id / "plma.json"
    if plma_path.exists():
        return json.loads(plma_path.read_text())
    return None


def load_family_data():
    """Load full_families.json, family_index.json, and index.json."""
    full_families = None
    family_index = None
    index_data = []

    ff_path = DATA_DIR / "full_families.json"
    if ff_path.exists():
        full_families = json.loads(ff_path.read_text())

    fi_path = DATA_DIR / "family_index.json"
    if fi_path.exists():
        family_index = json.loads(fi_path.read_text())

    idx_path = DATA_DIR / "index.json"
    if idx_path.exists():
        index_data = json.loads(idx_path.read_text())

    return full_families, family_index, index_data


def assemble_html(
    report_data, summary_data, pdb_variants, plma_data,
    full_families, family_index, index_data,
):
    """Assemble a self-contained HTML file from the DB codebase."""

    # Read template files from the codebase
    html = (PROJECT_DIR / "report.html").read_text(encoding="utf-8")
    app_js = (PROJECT_DIR / "static" / "js" / "app.js").read_text(encoding="utf-8")

    # --- Step 1: Download and cache CDN libraries ---
    log("Loading libraries...")
    libs = {}
    for name, url in CDN_LIBS.items():
        libs[name] = download_and_cache(name, url)

    # --- Step 2: Build inline data script ---
    inline_data = {
        "DATA": report_data,
        "SUMMARY": summary_data,
        "PDB": pdb_variants or {},
        "PLMA": plma_data,
        "FULL_FAMILIES": full_families,
        "FAMILY_INDEX": family_index,
        "INDEX": index_data,
    }
    # Use separators to minimize size
    data_json = json.dumps(inline_data, separators=(",", ":"))
    data_script = f"<script>window.__INLINE__={data_json};</script>"

    # --- Step 3: Transform HTML ---

    # 3a: Replace Molstar CSS link with inline style
    html = html.replace(
        '<link rel="stylesheet" href="https://unpkg.com/molstar@3.43.0/build/viewer/molstar.css"/>',
        f'<style>{libs["molstar.css"]}</style>',
    )

    # 3b: Replace CDN script tags with inline scripts
    cdn_replacements = {
        '<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.1/dist/chart.umd.min.js"></script>':
            f'<script>{libs["chart.umd.min.js"]}</script>',
        '<script src="https://unpkg.com/openchemlib@8.13.0/dist/openchemlib-full.js"></script>':
            f'<script>{libs["openchemlib-full.js"]}</script>',
        '<script src="https://unpkg.com/molstar@3.43.0/build/viewer/molstar.js"></script>':
            f'<script>{libs["molstar.js"]}</script>',
    }
    for old, new in cdn_replacements.items():
        html = html.replace(old, new)

    # 3c: Replace app.js reference with data injection + inline app.js
    html = html.replace(
        '<script src="static/js/app.js"></script>',
        f'{data_script}\n<script>{app_js}</script>',
    )

    return html


def main():
    parser = argparse.ArgumentParser(
        description="Generate a self-contained HTML report for a paralog pair"
    )
    parser.add_argument("pair_id", help="Pair ID (e.g., CDK4_CDK6)")
    parser.add_argument(
        "-o", "--output", help="Output file path (default: output/<pair_id>.html)"
    )
    parser.add_argument(
        "--use-cached",
        action="store_true",
        help="Use cached report data from DB (skip regeneration)",
    )
    parser.add_argument(
        "--regenerate",
        action="store_true",
        help="Force regeneration of report data",
    )
    parser.add_argument(
        "--full",
        action="store_true",
        help="Run full pipeline (populate DB + report data + PLMA) for pairs not yet in DB",
    )
    args = parser.parse_args()

    pair_id = args.pair_id
    log(f"Generating self-contained report for {pair_id}...")

    # Connect to DB
    conn = get_db()

    # Full pipeline mode: populate + report data + PLMA in one go
    if args.full:
        conn = run_full_pipeline(pair_id, conn)
    else:
        # Check if pair exists in DB
        pair = conn.execute("SELECT * FROM pairs WHERE pair_id = ?", (pair_id,)).fetchone()
        if not pair:
            log(f"ERROR: Pair {pair_id} not found in database")
            log("  Hint: use --full to run the full pipeline (populate + report data + PLMA)")
            log("Available pairs:")
            for r in conn.execute("SELECT pair_id FROM pairs").fetchall():
                log(f"  {r[0]}")
            sys.exit(1)

        # Check if report data needs generation
        report = conn.execute(
            "SELECT data_json FROM report_data WHERE pair_id = ?", (pair_id,)
        ).fetchone()

        if not report or not report["data_json"] or args.regenerate:
            log("Generating report data (this may take a minute)...")
            sys.path.insert(0, str(SCRIPTS_DIR))
            from generate_report_data import generate_report_data as gen_report

            data_obj, summary_obj = gen_report(pair_id, conn)
            conn.execute(
                """INSERT OR REPLACE INTO report_data (pair_id, data_json, summary_json)
                   VALUES (?, ?, ?)""",
                (pair_id, json.dumps(data_obj), json.dumps(summary_obj)),
            )
            conn.commit()
            log("  Report data generated and cached in DB")

    # Load all data
    log("Loading pair data...")
    report_data, summary_data, pdb_variants = load_pair_data(conn, pair_id)

    # Load PLMA data
    plma_data = load_plma_data(pair_id)
    if not plma_data:
        log("  No PLMA data found (run extract_plma.py first if needed)")

    # Load family data
    full_families, family_index, index_data = load_family_data()

    conn.close()

    # Assemble HTML
    log("Assembling self-contained HTML...")
    html = assemble_html(
        report_data, summary_data, pdb_variants, plma_data,
        full_families, family_index, index_data,
    )

    # Write output
    if args.output:
        out_path = Path(args.output)
    else:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        out_path = OUTPUT_DIR / f"{pair_id}.html"

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html, encoding="utf-8")

    size_mb = out_path.stat().st_size / 1024 / 1024
    log(f"Written: {out_path} ({size_mb:.1f} MB)")
    log("Open in browser to view the report.")


if __name__ == "__main__":
    main()
