#!/usr/bin/env python3
"""
download_drugclip.py

Download DrugCLIP/GenomeScreen pocket data from HuggingFace for specified
UniProt accessions. Uses HTTP range requests to extract only the needed files
from the remote ZIP archive without downloading the full ~27 GB file.

Dataset: https://huggingface.co/datasets/THU-ATOM/GenomeScreen

Usage:
    python scripts/download_drugclip.py P11802 Q00534
    python scripts/download_drugclip.py --all          # all accessions from pairs.csv
"""

import argparse
import csv
import io
import json
import math
import os
import struct
import sys
import time
import zlib
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
PROJECT_DIR = Path(__file__).resolve().parent.parent
INPUT_DIR = PROJECT_DIR / "input"
DRUGCLIP_DIR = INPUT_DIR / "drugclip"
STRUCTURES_DIR = PROJECT_DIR / "data" / "structures"

DATASET_REPO = "THU-ATOM/GenomeScreen"
ZIP_FILENAME = "screen_results.zip"

# Radius (Angstroms) around grid center to find pocket residues in AlphaFold PDB
POCKET_RADIUS = 10.0


def log(msg):
    print(f"[drugclip-dl] {msg}", flush=True)


# ---------------------------------------------------------------------------
# HTTP range-request helpers
# ---------------------------------------------------------------------------
def get_http_client():
    """Return an httpx client with redirect support. Falls back to urllib."""
    try:
        import httpx
        return httpx.Client(follow_redirects=True, timeout=60)
    except ImportError:
        return None


def range_get(client, url, start, end):
    """Fetch bytes [start, end] inclusive via HTTP Range request."""
    headers = {"Range": f"bytes={start}-{end}"}
    if client is not None:  # httpx
        resp = client.get(url, headers=headers)
        resp.raise_for_status()
        return resp.content
    else:
        import urllib.request
        req = urllib.request.Request(url, headers=headers)
        with urllib.request.urlopen(req, timeout=60) as resp:
            return resp.read()


# ---------------------------------------------------------------------------
# HuggingFace helpers
# ---------------------------------------------------------------------------
def resolve_hf_zip_url():
    """Get the direct URL and total size for screen_results.zip on HuggingFace."""
    try:
        from huggingface_hub import HfApi
        api = HfApi()
        info = api.dataset_info(DATASET_REPO, files_metadata=True)
        for sib in info.siblings:
            if sib.rfilename == ZIP_FILENAME:
                size = sib.size or (sib.lfs.size if sib.lfs else None)
                url = f"https://huggingface.co/datasets/{info.id}/resolve/main/{ZIP_FILENAME}"
                return url, size
    except Exception as e:
        log(f"huggingface_hub failed ({e}), falling back to API")

    # Fallback: use HuggingFace API directly
    import urllib.request
    api_url = f"https://huggingface.co/api/datasets/{DATASET_REPO}"
    req = urllib.request.Request(api_url)
    with urllib.request.urlopen(req, timeout=30) as resp:
        data = json.loads(resp.read())

    siblings = data.get("siblings", [])
    for sib in siblings:
        if sib.get("rfilename") == ZIP_FILENAME:
            size = sib.get("size") or sib.get("lfs", {}).get("size")
            repo_id = data.get("id", DATASET_REPO)
            url = f"https://huggingface.co/datasets/{repo_id}/resolve/main/{ZIP_FILENAME}"
            return url, size

    raise RuntimeError(f"Could not find {ZIP_FILENAME} in {DATASET_REPO}")


# ---------------------------------------------------------------------------
# Remote ZIP parsing
# ---------------------------------------------------------------------------
def read_zip64_eocd(client, url, total_size):
    """Read the ZIP / ZIP64 end-of-central-directory and return
    (cd_offset, cd_size, total_entries).
    """
    # Read last 200 KB to find EOCD records
    tail_size = min(200_000, total_size)
    tail_start = total_size - tail_size
    tail = range_get(client, url, tail_start, total_size - 1)

    # Find EOCD signature
    eocd_pos = tail.rfind(b"PK\x05\x06")
    if eocd_pos < 0:
        raise RuntimeError("Cannot find ZIP end-of-central-directory")

    # Parse standard EOCD
    (_, disk_num, disk_cd, num_disk, num_total,
     cd_size, cd_offset, comment_len) = struct.unpack_from(
        "<IHHHHIIH", tail, eocd_pos
    )

    # Check for ZIP64
    needs_zip64 = (num_total == 0xFFFF or cd_offset == 0xFFFFFFFF or
                   cd_size == 0xFFFFFFFF)

    if needs_zip64:
        # Find ZIP64 EOCD Locator
        loc_pos = tail.rfind(b"PK\x06\x07", 0, eocd_pos)
        if loc_pos < 0:
            raise RuntimeError("ZIP64 EOCD Locator not found")

        _, loc_disk, z64_eocd_abs, loc_disks = struct.unpack_from(
            "<IIQI", tail, loc_pos
        )

        # Read ZIP64 EOCD (may be in our tail already)
        z64_tail_off = z64_eocd_abs - tail_start
        if 0 <= z64_tail_off < len(tail) - 56:
            z64 = tail[z64_tail_off:]
        else:
            z64 = range_get(client, url, z64_eocd_abs, z64_eocd_abs + 100)

        if z64[:4] != b"PK\x06\x06":
            raise RuntimeError("Invalid ZIP64 EOCD signature")

        (_, rec_size, ver_made, ver_need, disk, disk_cd,
         entries_disk, entries_total, cd_size, cd_offset) = struct.unpack_from(
            "<I Q H H I I Q Q Q Q", z64
        )
        num_total = entries_total

    return cd_offset, cd_size, num_total


def parse_central_directory(cd_data, target_prefixes):
    """Parse the central directory and return entries matching any target prefix.

    Returns dict: {prefix: [entry_dict, ...]}
    """
    found = {t: [] for t in target_prefixes}
    pos = 0
    count = 0

    while pos + 46 <= len(cd_data):
        sig = struct.unpack_from("<I", cd_data, pos)[0]
        if sig != 0x02014b50:
            break

        (ver_made, ver_need, flags, method, mtime, mdate, crc,
         comp_size, uncomp_size, fname_len, extra_len, comment_len,
         disk_start, int_attr, ext_attr, local_offset) = struct.unpack_from(
            "<HHHHHH III HHH HH I I", cd_data, pos + 4
        )

        fname = cd_data[pos + 46: pos + 46 + fname_len].decode("utf-8", errors="replace")

        # Resolve ZIP64 extra field
        actual_offset = local_offset
        actual_comp = comp_size
        actual_uncomp = uncomp_size

        if uncomp_size == 0xFFFFFFFF or comp_size == 0xFFFFFFFF or local_offset == 0xFFFFFFFF:
            extra = cd_data[pos + 46 + fname_len: pos + 46 + fname_len + extra_len]
            epos = 0
            while epos + 4 <= len(extra):
                eid, esz = struct.unpack_from("<HH", extra, epos)
                if eid == 0x0001:
                    z64 = extra[epos + 4: epos + 4 + esz]
                    zp = 0
                    if uncomp_size == 0xFFFFFFFF and zp + 8 <= len(z64):
                        actual_uncomp = struct.unpack_from("<Q", z64, zp)[0]
                        zp += 8
                    if comp_size == 0xFFFFFFFF and zp + 8 <= len(z64):
                        actual_comp = struct.unpack_from("<Q", z64, zp)[0]
                        zp += 8
                    if local_offset == 0xFFFFFFFF and zp + 8 <= len(z64):
                        actual_offset = struct.unpack_from("<Q", z64, zp)[0]
                        zp += 8
                    break
                epos += 4 + esz

        for prefix in target_prefixes:
            if prefix in fname:
                found[prefix].append({
                    "name": fname,
                    "offset": actual_offset,
                    "comp_size": actual_comp,
                    "uncomp_size": actual_uncomp,
                    "method": method,
                })
                break

        pos += 46 + fname_len + extra_len + comment_len
        count += 1
        if count % 100000 == 0:
            log(f"  Parsed {count} central directory entries...")

    log(f"  Central directory: {count} total entries parsed")
    return found


def download_file_from_zip(client, url, entry):
    """Download and decompress a single file from the remote zip via its
    local file header offset."""
    offset = entry["offset"]

    # Read local file header (30 bytes + variable)
    local_header = range_get(client, url, offset, offset + 29)
    if local_header[:4] != b"PK\x03\x04":
        log(f"  WARNING: Bad local header for {entry['name']}")
        return None

    fname_len = struct.unpack_from("<H", local_header, 26)[0]
    extra_len = struct.unpack_from("<H", local_header, 28)[0]
    data_offset = offset + 30 + fname_len + extra_len

    comp_size = entry["comp_size"]
    if comp_size == 0:
        return b""

    raw = range_get(client, url, data_offset, data_offset + comp_size - 1)

    method = entry["method"]
    if method == 0:  # stored
        return raw
    elif method == 8:  # deflate
        try:
            return zlib.decompress(raw, -15)
        except zlib.error as e:
            log(f"  WARNING: Decompression failed for {entry['name']}: {e}")
            return None
    else:
        log(f"  WARNING: Unknown compression method {method} for {entry['name']}")
        return None


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------
def parse_grid_in(text):
    """Parse a grid.in file and return [x, y, z] grid center coordinates."""
    for line in text.strip().splitlines():
        if line.startswith("GRID_CENTER"):
            parts = line.split(None, 1)[1]
            coords = [float(x.strip().rstrip(",")) for x in parts.split(",")]
            return coords
    return None


def parse_leader_csv(data_bytes):
    """Parse leader.csv and return list of drug hit dicts."""
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


def find_pocket_residues(pdb_path, grid_center, radius=POCKET_RADIUS):
    """Find CA atoms within `radius` Angstroms of grid_center in an AlphaFold PDB.

    Returns (start_residue, end_residue, sorted_residue_list).
    """
    if grid_center is None or not pdb_path.exists():
        return None, None, []

    gx, gy, gz = grid_center
    r2 = radius * radius
    pocket_residues = set()

    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                resnum = int(line[22:26].strip())
            except (ValueError, IndexError):
                continue

            dx, dy, dz = x - gx, y - gy, z - gz
            if dx * dx + dy * dy + dz * dz <= r2:
                pocket_residues.add(resnum)

    if not pocket_residues:
        return None, None, []

    sorted_res = sorted(pocket_residues)
    return sorted_res[0], sorted_res[-1], sorted_res


# ---------------------------------------------------------------------------
# Main logic
# ---------------------------------------------------------------------------
def process_accession(client, url, cd_entries, accession, pdb_path):
    """Download and process all pockets for one UniProt accession.

    Returns list of pocket dicts ready for JSON serialization.
    """
    prefix = f"AF-{accession}"
    entries = cd_entries.get(prefix, [])
    if not entries:
        log(f"  No entries found for {accession}")
        return []

    # Group entries by pocket folder
    folders = {}
    for e in entries:
        parts = e["name"].split("/")
        if len(parts) < 2:
            continue
        folder = parts[1]
        if folder not in folders:
            folders[folder] = []
        folders[folder].append(e)

    log(f"  {accession}: {len(folders)} pocket folder(s), {len(entries)} total files")

    pockets_data = []
    for folder_name in sorted(folders.keys()):
        folder_entries = folders[folder_name]

        # Identify file types
        leader_entries = [e for e in folder_entries if e["name"].endswith("leader.csv")]
        grid_entries = [e for e in folder_entries if e["name"].endswith("_grid.in")]

        # Download leader.csv
        drugs = []
        if leader_entries:
            data = download_file_from_zip(client, url, leader_entries[0])
            if data:
                drugs = parse_leader_csv(data)

        # Download and parse grid.in files -> collect grid centers
        grid_centers = []
        for ge in grid_entries:
            data = download_file_from_zip(client, url, ge)
            if data:
                text = data.decode("utf-8", errors="replace")
                coords = parse_grid_in(text)
                if coords:
                    grid_centers.append(coords)

        # Average grid center
        avg_grid = None
        if grid_centers:
            n = len(grid_centers)
            avg_grid = [
                round(sum(g[0] for g in grid_centers) / n, 3),
                round(sum(g[1] for g in grid_centers) / n, 3),
                round(sum(g[2] for g in grid_centers) / n, 3),
            ]

        # Extract pocket number from folder name
        # Formats: AF-P11802-F1-model_v4_0_pocket4  or  AF-P11802-F1-model_v4_0_4
        pocket_label = folder_name.split("_")[-1]  # "pocket4" or "4" or "0"
        if pocket_label.startswith("pocket"):
            pocket_num_str = pocket_label.replace("pocket", "")
        else:
            pocket_num_str = pocket_label

        try:
            pocket_num = int(pocket_num_str)
        except ValueError:
            pocket_num = 0

        # Find pocket residues from AlphaFold PDB + grid center
        start_res, end_res, res_list = find_pocket_residues(pdb_path, avg_grid)

        pocket_info = {
            "pocket_num": pocket_num,
            "name": folder_name,
            "start": start_res,
            "end": end_res,
            "grid_center": avg_grid,
            "top_drugs": drugs[:10],
            "n_total_hits": len(drugs),
            "top_score": round(drugs[0]["score"], 3) if drugs else None,
            "n_pocket_residues": len(res_list),
            "residues": res_list,
        }
        pockets_data.append(pocket_info)

        log(f"    {folder_name}: pocket_num={pocket_num}, "
            f"residues {start_res}-{end_res} ({len(res_list)} CA), "
            f"grid={avg_grid}, "
            f"hits={len(drugs)}, top_score={pocket_info['top_score']}")

    # Sort by pocket_num
    pockets_data.sort(key=lambda p: p["pocket_num"])
    return pockets_data


def main():
    parser = argparse.ArgumentParser(
        description="Download DrugCLIP/GenomeScreen pocket data from HuggingFace"
    )
    parser.add_argument(
        "accessions", nargs="*",
        help="UniProt accession(s) to download (e.g. P11802 Q00534)"
    )
    parser.add_argument(
        "--all", action="store_true",
        help="Download for all accessions found in pairs.csv"
    )
    args = parser.parse_args()

    if args.all:
        pairs_csv = PROJECT_DIR / "input" / "pairs.csv"
        if not pairs_csv.exists():
            log(f"ERROR: {pairs_csv} not found")
            sys.exit(1)
        accessions = set()
        with open(pairs_csv) as f:
            reader = csv.DictReader(f)
            for row in reader:
                for col in ["accession_a", "accession_b"]:
                    if col in row and row[col]:
                        accessions.add(row[col])
        accessions = sorted(accessions)
        log(f"Found {len(accessions)} accessions from pairs.csv")
    elif args.accessions:
        accessions = args.accessions
    else:
        parser.print_help()
        sys.exit(1)

    # Step 1: Resolve HuggingFace ZIP URL
    log("Resolving HuggingFace dataset URL...")
    zip_url, total_size = resolve_hf_zip_url()
    log(f"ZIP URL: {zip_url}")
    log(f"ZIP size: {total_size / 1e9:.2f} GB")

    # Step 2: Open HTTP client
    client = get_http_client()
    if client is None:
        log("WARNING: httpx not available, using urllib (slower)")

    # Step 3: Read ZIP central directory
    log("Reading ZIP end-of-central-directory...")
    cd_offset, cd_size, n_entries = read_zip64_eocd(client, zip_url, total_size)
    log(f"Central directory: {n_entries} entries, {cd_size / 1e6:.1f} MB at offset {cd_offset}")

    log("Downloading central directory (this may take a moment)...")
    t0 = time.time()
    chunk_size = 10_000_000
    cd_data = bytearray()
    for chunk_start in range(0, cd_size, chunk_size):
        chunk_end = min(chunk_start + chunk_size, cd_size) - 1
        abs_start = cd_offset + chunk_start
        abs_end = cd_offset + chunk_end
        chunk = range_get(client, zip_url, abs_start, abs_end)
        cd_data.extend(chunk)
    log(f"Central directory downloaded: {len(cd_data) / 1e6:.1f} MB in {time.time() - t0:.0f}s")

    # Step 4: Parse central directory for target accessions
    target_prefixes = [f"AF-{acc}" for acc in accessions]
    log(f"Searching for {len(target_prefixes)} accession(s) in central directory...")
    cd_entries = parse_central_directory(cd_data, target_prefixes)

    # Free central directory memory
    del cd_data

    # Step 5: Process each accession
    DRUGCLIP_DIR.mkdir(parents=True, exist_ok=True)
    results_summary = []

    for accession in accessions:
        prefix = f"AF-{accession}"
        n_files = len(cd_entries.get(prefix, []))
        log(f"\nProcessing {accession} ({n_files} files in ZIP)...")

        # Determine PDB path
        pdb_path = STRUCTURES_DIR / f"AF-{accession}-F1-AM.pdb"
        if not pdb_path.exists():
            log(f"  WARNING: AlphaFold PDB not found at {pdb_path}")
            log(f"           Pocket residue ranges will be unavailable")

        pockets = process_accession(client, zip_url, cd_entries, accession, pdb_path)

        if not pockets:
            log(f"  No pockets found for {accession}")
            continue

        # Save pockets.json
        acc_dir = DRUGCLIP_DIR / accession
        acc_dir.mkdir(parents=True, exist_ok=True)
        out_path = acc_dir / "pockets.json"
        with open(out_path, "w") as f:
            json.dump(pockets, f, indent=2)

        log(f"  Saved {len(pockets)} pocket(s) to {out_path}")
        results_summary.append((accession, len(pockets)))

    # Close HTTP client
    if hasattr(client, "close"):
        client.close()

    # Summary
    log(f"\n{'='*60}")
    log(f"DONE - Results summary:")
    for acc, n_pockets in results_summary:
        log(f"  {acc}: {n_pockets} pocket(s) -> input/drugclip/{acc}/pockets.json")

    missing = [acc for acc in accessions
               if not any(acc == r[0] for r in results_summary)]
    if missing:
        log(f"\nAccessions with NO pockets in GenomeScreen:")
        for acc in missing:
            log(f"  {acc}")

    log(f"Output directory: {DRUGCLIP_DIR}")


if __name__ == "__main__":
    main()
