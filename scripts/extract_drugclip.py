#!/usr/bin/env python3
"""
extract_drugclip.py

Extract DrugCLIP/GenomeScreen pocket data from the (partial) zip download
for target UniProt IDs used in the paralog database.

Reads local file headers from the zip (works even with incomplete downloads),
extracts leader.csv, grid.in, and one representative PDB for each pocket.

Usage:
    python scripts/extract_drugclip.py
    python scripts/extract_drugclip.py --zip /path/to/screen_results.zip
"""

import argparse
import csv
import gzip
import io
import json
import math
import os
import struct
import sys
import zlib
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent
INPUT_DIR = PROJECT_DIR / "input"
DRUGCLIP_DIR = INPUT_DIR / "drugclip"

# Default zip path (HuggingFace partial download)
DEFAULT_ZIP = Path("/Users/olivierdennler/Documents/data/DrugCLIP/GenomeScreenDB/.huggingface/download")

# Target UniProt IDs from pairs.csv mapped to gene names
TARGET_UNIPROTS = {
    "P51531": "SMARCA2", "P51532": "SMARCA4",
    "O14497": "ARID1A", "Q8NFD5": "ARID1B",
    "Q68DW7": "STAG1", "Q8N3U4": "STAG2",
    "P06733": "ENO1", "P09104": "ENO2", "P13929": "ENO3", "A6NNW6": "ENO4",
    "Q92922": "SMARCC1", "Q8TAQ2": "SMARCC2",
    "P0CG47": "UBB", "P0CG48": "UBC",
    "P61326": "MAGOH", "Q96A72": "MAGOHB",
    "P23368": "ME2", "Q16798": "ME3",
}


def log(msg):
    print(f"[drugclip] {msg}", flush=True)


def find_zip_file(zip_dir):
    """Find the screen_results zip (complete or incomplete) in a directory."""
    zip_dir = Path(zip_dir)

    # Direct file
    if zip_dir.is_file():
        return zip_dir

    # Look for screen_results.zip or incomplete download
    candidates = []
    for f in zip_dir.rglob("*"):
        if f.is_file() and "screen_results" in f.name and f.stat().st_size > 1_000_000:
            candidates.append(f)

    if not candidates:
        return None

    # Prefer complete over incomplete
    for c in candidates:
        if c.suffix == ".zip":
            return c
    return candidates[0]


def iter_zip_entries(zip_path, target_ids):
    """Iterate over local file headers in a zip, yielding entries matching target IDs.

    Works with incomplete zip files by reading local headers sequentially.
    """
    target_set = set(target_ids)

    with open(zip_path, "rb") as f:
        offset = 0
        count = 0
        while True:
            f.seek(offset)
            sig = f.read(4)
            if sig != b"PK\x03\x04":
                break

            header = f.read(26)
            if len(header) < 26:
                break

            (ver, flags, method, mtime, mdate, crc,
             comp_size, uncomp_size, fname_len, extra_len) = struct.unpack(
                "<HHHHH III HH", header
            )

            fname_bytes = f.read(fname_len)
            if len(fname_bytes) < fname_len:
                break
            fname = fname_bytes.decode("utf-8", errors="replace")

            extra = f.read(extra_len)
            data_offset = f.tell()
            next_offset = data_offset + comp_size

            # Check if this entry matches any target
            matched_uid = None
            for uid in target_set:
                if uid in fname:
                    matched_uid = uid
                    break

            if matched_uid is not None:
                # Read and decompress data
                raw = f.read(comp_size)
                data = None
                if method == 0:  # stored
                    data = raw
                elif method == 8:  # deflate
                    try:
                        data = zlib.decompress(raw, -15)
                    except Exception:
                        data = None

                yield fname, matched_uid, data

            offset = next_offset
            count += 1
            if count % 20000 == 0:
                log(f"  Scanned {count} zip entries...")

        log(f"  Total zip entries scanned: {count}")


def parse_grid_in(text):
    """Parse a grid.in file and return grid center coordinates."""
    for line in text.strip().splitlines():
        if line.startswith("GRID_CENTER"):
            parts = line.split(None, 1)[1]
            coords = [float(x.strip().rstrip(",")) for x in parts.split(",")]
            return coords
    return None


def parse_leader_csv(data_bytes):
    """Parse leader.csv and return list of drug hits."""
    text = data_bytes.decode("utf-8", errors="replace")
    reader = csv.DictReader(io.StringIO(text))
    drugs = []
    for row in reader:
        try:
            drugs.append({
                "smiles": row.get("smiles", ""),
                "name": row.get("Name", ""),
                "oid": row.get("oid", ""),
                "score": float(row.get("score", 0)),
            })
        except (ValueError, KeyError):
            continue
    return drugs


def pocket_residues_from_pdb_and_grid(pdb_bytes, grid_center, radius=8.0):
    """Find residues within radius of grid center in a PDB file.

    Returns (min_resnum, max_resnum, list_of_resnums).
    """
    if pdb_bytes is None or grid_center is None:
        return None, None, []

    # Decompress if gzipped
    try:
        pdb_text = gzip.decompress(pdb_bytes).decode("utf-8", errors="replace")
    except Exception:
        pdb_text = pdb_bytes.decode("utf-8", errors="replace")

    gx, gy, gz = grid_center
    r2 = radius * radius
    pocket_residues = set()

    for line in pdb_text.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            resnum = int(line[22:26].strip())
            chain = line[21]
        except (ValueError, IndexError):
            continue

        if chain != "A":
            continue

        dx = x - gx
        dy = y - gy
        dz = z - gz
        if dx * dx + dy * dy + dz * dz <= r2:
            pocket_residues.add(resnum)

    if not pocket_residues:
        return None, None, []

    sorted_res = sorted(pocket_residues)
    return min(sorted_res), max(sorted_res), sorted_res


def main():
    parser = argparse.ArgumentParser(description="Extract DrugCLIP pocket data")
    parser.add_argument("--zip", default=str(DEFAULT_ZIP),
                        help="Path to screen_results.zip or directory containing it")
    args = parser.parse_args()

    zip_path = find_zip_file(args.zip)
    if not zip_path:
        log(f"ERROR: Could not find GenomeScreen zip in {args.zip}")
        sys.exit(1)

    log(f"Using zip: {zip_path} ({zip_path.stat().st_size / 1e9:.2f} GB)")

    # Collect data per pocket
    # Structure: pockets[uniprot_id][pocket_name] = {leader, grids, pdb}
    pockets = {}

    log("Scanning zip for target UniProt IDs...")
    for fname, uid, data in iter_zip_entries(zip_path, TARGET_UNIPROTS.keys()):
        if data is None:
            continue

        # Parse pocket folder name from path
        # e.g. screen_results/AF-P09104-F1-model_v4_0_pocket4/leader.csv
        parts = fname.split("/")
        if len(parts) < 3:
            continue

        pocket_folder = parts[1]  # AF-P09104-F1-model_v4_0_pocket4

        if uid not in pockets:
            pockets[uid] = {}
        if pocket_folder not in pockets[uid]:
            pockets[uid][pocket_folder] = {
                "leader": None,
                "grids": [],
                "pdb0": None,
            }

        pocket = pockets[uid][pocket_folder]
        basename = parts[-1] if len(parts) > 2 else ""

        if basename == "leader.csv":
            pocket["leader"] = data
        elif basename.endswith("_grid.in"):
            grid_text = data.decode("utf-8", errors="replace")
            coords = parse_grid_in(grid_text)
            if coords:
                pocket["grids"].append(coords)
        elif basename.endswith(".pdbgz") and "_0_complex" in basename:
            pocket["pdb0"] = data

    # Process and save
    DRUGCLIP_DIR.mkdir(parents=True, exist_ok=True)
    total_pockets = 0

    for uid, pocket_folders in sorted(pockets.items()):
        gene = TARGET_UNIPROTS.get(uid, uid)
        log(f"\n{gene} ({uid}): {len(pocket_folders)} pocket(s)")

        acc_dir = DRUGCLIP_DIR / uid
        acc_dir.mkdir(parents=True, exist_ok=True)

        pockets_data = []

        for pocket_name, pocket_data in sorted(pocket_folders.items()):
            # Parse leader.csv
            drugs = []
            if pocket_data["leader"]:
                drugs = parse_leader_csv(pocket_data["leader"])
                log(f"  {pocket_name}: {len(drugs)} drug hits")

            # Average grid center from all conformers
            grids = pocket_data["grids"]
            avg_grid = None
            if grids:
                avg_grid = [
                    sum(g[0] for g in grids) / len(grids),
                    sum(g[1] for g in grids) / len(grids),
                    sum(g[2] for g in grids) / len(grids),
                ]

            # Find pocket residues from PDB + grid center
            min_res, max_res, res_list = pocket_residues_from_pdb_and_grid(
                pocket_data["pdb0"], avg_grid
            )

            if min_res is None and avg_grid is not None:
                # Fallback: try with larger radius
                min_res, max_res, res_list = pocket_residues_from_pdb_and_grid(
                    pocket_data["pdb0"], avg_grid, radius=12.0
                )

            # Extract pocket number from folder name
            pocket_num = ""
            if "pocket" in pocket_name:
                pocket_num = pocket_name.split("pocket")[-1]

            pocket_info = {
                "name": pocket_name,
                "pocket_num": pocket_num,
                "grid_center": avg_grid,
                "start": min_res,
                "end": max_res,
                "residues": res_list,
                "n_residues": len(res_list),
                "top_drugs": drugs[:10],
                "n_total_hits": len(drugs),
                "top_score": drugs[0]["score"] if drugs else None,
            }
            pockets_data.append(pocket_info)
            total_pockets += 1

            log(f"    Pocket {pocket_num}: residues {min_res}-{max_res} ({len(res_list)} res), "
                f"grid={avg_grid and [round(c,1) for c in avg_grid]}, "
                f"top_score={pocket_info['top_score']}")

            # Save representative PDB
            if pocket_data["pdb0"]:
                pdb_path = acc_dir / f"pocket{pocket_num}_complex.pdb.gz"
                with open(pdb_path, "wb") as f:
                    # Already gzipped
                    f.write(pocket_data["pdb0"])

        # Save pockets.json for this accession
        with open(acc_dir / "pockets.json", "w") as f:
            json.dump(pockets_data, f, indent=2)

    log(f"\nDone! Extracted {total_pockets} pockets for {len(pockets)} genes")
    log(f"Output: {DRUGCLIP_DIR}")

    # Summary of missing genes
    found = set(pockets.keys())
    missing = set(TARGET_UNIPROTS.keys()) - found
    if missing:
        log(f"\nGenes NOT found in zip (may be in undownloaded portion):")
        for uid in sorted(missing):
            log(f"  {TARGET_UNIPROTS[uid]} ({uid})")


if __name__ == "__main__":
    main()
