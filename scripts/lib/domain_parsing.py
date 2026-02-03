"""Domain and annotation parsing functions."""

import csv
import json
import re
import tarfile
from typing import Any, Dict, List, Optional, Tuple

from .config import CAVITY_INDEX_PATH, CAVITY_DATA_DIR, log


def parse_ted_domains(ted_json: Optional[Dict]) -> List[Dict[str, Any]]:
    """Parse TED domain annotations from /summary endpoint."""
    if not ted_json:
        return []
    
    # Handle different response formats
    rows = []
    if isinstance(ted_json, dict):
        rows = ted_json.get("results") or ted_json.get("data") or []
    elif isinstance(ted_json, list):
        rows = ted_json
    
    out = []
    for i, r in enumerate(rows):
        # Modern format: explicit start/end
        if "start" in r and "end" in r:
            try:
                s = int(r["start"])
                e = int(r["end"])
                if e > s:
                    name = (r.get("name") or r.get("cath_label") or "").strip()
                    if not name or name == "-":
                        name = f"TED{i+1}"
                    out.append({
                        "start": s,
                        "end": e,
                        "label": name,
                        "name": name,
                        "type": "TED",
                        "uid": f"TED:{s}-{e}:{name}"
                    })
            except:
                continue
        # Older "chopping" format
        elif r.get("chopping"):
            chop = str(r.get("chopping", "")).strip()
            if "-" not in chop:
                continue
            try:
                for j, part in enumerate(chop.split(",")):
                    part = part.strip()
                    if "-" not in part:
                        continue
                    s, e = map(int, part.split("-", 1))
                    if e <= s:
                        continue
                    lab = (r.get("cath_label") or r.get("ted_id") or "").strip()
                    if not lab or lab == "-":
                        lab = f"TED{i+1}_{j+1}"
                    out.append({
                        "start": s,
                        "end": e,
                        "label": lab,
                        "name": lab,
                        "type": "TED",
                        "uid": f"TED:{s}-{e}:{lab}"
                    })
            except:
                continue
    return out


def parse_ebi_domains(ebi_json: Optional[Dict], dtype: str) -> List[Dict[str, Any]]:
    """Parse EBI/UniProt domain annotations."""
    if not ebi_json:
        return []
    features = ebi_json.get("features") or (ebi_json if isinstance(ebi_json, list) else [])
    out = []
    for f in features:
        cat = (f.get("category") or "").upper()
        if cat != "DOMAINS_AND_SITES":
            continue
        desc = f.get("description") or f.get("type") or ""
        ftype = (f.get("type") or "").lower()
        is_disorder = "disorder" in desc.lower() or ftype == "disorder"
        if (dtype == "Disordered") != is_disorder:
            continue
        s = f.get("begin") or (f.get("location", {}).get("begin"))
        e = f.get("end") or (f.get("location", {}).get("end"))
        if s and e:
            try:
                out.append({
                    "start": int(s), "end": int(e),
                    "label": desc or f.get("type") or "Domain",
                    "type": "Disordered" if is_disorder else "Domain",
                    "uid": f"EBI:{s}-{e}",
                    "raw_type": f.get("type") or ""
                })
            except:
                pass
    return out


def parse_cavities(cav_json: Optional[List]) -> List[Dict[str, Any]]:
    """Parse cavity annotations from database JSON."""
    if not cav_json or not isinstance(cav_json, list):
        return []
    return [
        {
            "start": c.get("start"), 
            "end": c.get("end"), 
            "label": c.get("label", "Cavity"), 
            "type": "Cavity", 
            "uid": f"CAV:{c.get('start')}-{c.get('end')}",
            "drug_score": c.get("drug_score"),
            "druggability": c.get("druggability"),
        } 
        for c in cav_json if c.get("start") and c.get("end")
    ]


# ============= Local Cavity Loading from AF-all_cavities =============

def _load_cavity_rows_for_acc(acc: str) -> List[Dict[str, str]]:
    """
    Return all index rows for a given UniProt acc from AF-all_cavities.idx.
    Primary match: uniprot_id == acc
    Fallback: cavity_name startswith 'AF-{acc}-'
    """
    rows = []
    if not CAVITY_INDEX_PATH.exists():
        return rows
    
    try:
        with CAVITY_INDEX_PATH.open(newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                if row.get("uniprot_id") == acc:
                    rows.append(row)
    except Exception as e:
        pass
    
    if not rows:
        try:
            with CAVITY_INDEX_PATH.open(newline="", encoding="utf-8") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                key = f"AF-{acc}-"
                for row in reader:
                    c = row.get("cavity_name") or ""
                    if c.startswith(key):
                        rows.append(row)
        except:
            pass
    
    return rows


def _cavity_range_from_tar(acc: str, cavity_name: str) -> Optional[Tuple[int, int]]:
    """
    Find residue range for a cavity by scanning the PDB inside the corresponding
    AF-{acc}-F1-model_v1_cavity_result_*.tar.gz.
    """
    if not CAVITY_DATA_DIR.exists():
        return None
    
    cav_id = None
    m = re.search(r"_C(\d+)$", cavity_name)
    if m:
        cav_id = m.group(1)
    
    candidates = []
    if cav_id:
        pattern = f"AF-{acc}-F1-model_v1_cavity_result_{cav_id}.tar.gz"
        candidates = list(CAVITY_DATA_DIR.glob(pattern))
    
    if not candidates:
        pattern = f"AF-{acc}-F1-model_v1_cavity_result_*.tar.gz"
        candidates = list(CAVITY_DATA_DIR.glob(pattern))
    
    for tar_path in sorted(candidates):
        try:
            with tarfile.open(tar_path, "r:gz") as tar:
                pdb_member = None
                for mem in tar.getmembers():
                    if mem.name.lower().endswith(".pdb"):
                        pdb_member = mem
                        break
                if pdb_member is None:
                    continue
                
                fobj = tar.extractfile(pdb_member)
                if fobj is None:
                    continue
                
                seen = set()
                for bline in fobj:
                    try:
                        ln = bline.decode("ascii", "ignore")
                    except:
                        continue
                    if not ln.startswith("ATOM"):
                        continue
                    try:
                        res = int(ln[22:26])
                        seen.add(res)
                    except:
                        continue
                
                if not seen:
                    continue
                
                return min(seen), max(seen)
        except Exception as e:
            pass
    
    return None


def load_cavities_for_acc(acc: str) -> List[Dict[str, Any]]:
    """
    Load cavity annotations from local AF-all_cavities data.
    Returns domain-like dicts with start/end, druggability, etc.
    """
    rows = _load_cavity_rows_for_acc(acc)
    out = []
    
    for row in rows:
        cav_name = row.get("cavity_name") or ""
        if not cav_name:
            continue
        
        rng = _cavity_range_from_tar(acc, cav_name)
        if not rng:
            continue
        s, e = rng
        if e <= s:
            continue
        
        drug_score = row.get("drug_score")
        druggability = row.get("druggability")
        
        # Create label with druggability info
        label = cav_name.split("_")[-1] if "_" in cav_name else cav_name  # e.g., "C1"
        if druggability:
            label = f"{label} ({druggability})"
        
        out.append({
            "start": s,
            "end": e,
            "name": cav_name,
            "label": label,
            "type": "Cavity",
            "uid": f"CAV:{s}-{e}:{cav_name}",
            "drug_score": drug_score,
            "druggability": druggability,
            "vacant_volume": row.get("vacant_volume"),
            "surface_area": row.get("surface_area"),
        })
    
    return out


def rects_on_alignment(
    doms: List[Dict], 
    which: str, 
    color: str, 
    aligned_cols: List, 
    alnLen: int
) -> List[Dict[str, Any]]:
    """Map domain coordinates to alignment column coordinates."""
    if not doms:
        return []
    
    pos_to_col = {}
    for col, qpos, tpos in aligned_cols:
        pos = qpos if which == 'A' else tpos
        if pos is not None:
            pos_to_col[pos] = col
    
    out = []
    for d in doms:
        s, e = d.get("start"), d.get("end")
        if s is None or e is None:
            continue
        try:
            s, e = int(s), int(e)
        except:
            continue
        if e < s:
            s, e = e, s
        
        cols = [pos_to_col[p] for p in range(s, e + 1) if p in pos_to_col]
        if not cols:
            continue
        cols.sort()
        
        # Group into contiguous segments
        seg_start = cols[0]
        prev = cols[0]
        
        for c in cols[1:] + [None]:
            if c is None or c != prev + 1:
                a = max(1, min(alnLen, seg_start))
                b = max(1, min(alnLen, prev))
                if b >= a:
                    out.append({
                        "start": a, "end": b,
                        "color": d.get("color") or color,
                        "label": d.get("label") or "",
                        "id": d.get("uid") or f"{s}-{e}",
                        "druggability": d.get("druggability"),  # Pass through for cavities
                    })
                if c is not None:
                    seg_start = c
            if c is not None:
                prev = c
    return out


def get_annotation_from_db(conn, acc: str, source: str) -> Optional[Dict]:
    """Get annotation JSON from database."""
    row = conn.execute(
        "SELECT data_json FROM annotations WHERE acc=? AND source=?",
        (acc, source)
    ).fetchone()
    if row and row["data_json"]:
        try:
            return json.loads(row["data_json"])
        except:
            pass
    return None
