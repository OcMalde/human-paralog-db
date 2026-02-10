"""API fetching functions for PDBe and AlphaMissense."""

import base64
import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any, Dict, List, Optional, Tuple

import requests

from .config import PDBE_BASE, AM_HOTSPOT_API, CACHE_DIR, AA_ORDER, log


# ============= HTTP Helpers =============

def http_get_json(url: str, timeout: int = 30) -> Optional[Dict]:
    """GET JSON from URL with error handling."""
    try:
        resp = requests.get(url, timeout=timeout, headers={"Accept": "application/json"})
        if resp.ok:
            return resp.json()
    except:
        pass
    return None


def http_get_binary(url: str, timeout: int = 60) -> Optional[bytes]:
    """GET binary content from URL."""
    try:
        resp = requests.get(url, timeout=timeout)
        if resp.ok:
            return resp.content
    except:
        pass
    return None


# ============= AlphaMissense Hotspot Fetching =============

def _fetch_single_hotspot(acc: str, resi: int) -> Tuple[int, Optional[Dict]]:
    """Fetch a single hotspot - for threading."""
    url = f"{AM_HOTSPOT_API}?uid={acc}&resi={resi}"
    try:
        resp = requests.get(url, timeout=20, headers={"Accept": "application/json"})
        if resp.ok:
            return resi, resp.json()
    except:
        pass
    return resi, None


def fetch_am_hotspots(acc: str, prot_len: int, threads: int = 8) -> List[Optional[Dict[str, Any]]]:
    """Fetch AlphaMissense hotspots for all residues, with caching and threading."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_file = CACHE_DIR / f"am_hotspots_{acc}.json"
    
    # Load existing cache
    cache = {}
    if cache_file.exists():
        try:
            cache = json.loads(cache_file.read_text(encoding='utf-8'))
            cache = {str(k): v for k, v in cache.items()}
        except:
            cache = {}
    
    # Find missing residues
    missing = [i for i in range(1, prot_len + 1) if str(i) not in cache]
    
    if missing:
        log(f"  Fetching {len(missing)}/{prot_len} AM hotspots for {acc} ({threads} threads)...")
        
        # Use threading for faster fetching
        fetched = 0
        with ThreadPoolExecutor(max_workers=threads) as executor:
            futures = {executor.submit(_fetch_single_hotspot, acc, resi): resi for resi in missing}
            
            for future in as_completed(futures):
                resi, data = future.result()
                cache[str(resi)] = data
                fetched += 1
                
                if fetched % 100 == 0:
                    log(f"    {fetched}/{len(missing)} fetched...")
        
        # Save cache
        try:
            cache_file.write_text(json.dumps(cache), encoding='utf-8')
            log(f"  Cached AM hotspots for {acc}")
        except Exception as e:
            log(f"  Failed to save cache: {e}")
    else:
        log(f"  AM hotspots for {acc}: {len(cache)} cached (skipped)")
    
    return [cache.get(str(i)) for i in range(1, prot_len + 1)]


def build_am_matrix_scores(
    hotspots: List[Optional[Dict[str, Any]]],
    pos_by_col: Dict[int, Optional[int]],
    aln_len: int,
) -> Dict[str, List[Optional[float]]]:
    """Build per-AA substitution score arrays on the alignment axis."""
    scores_by_aa = {aa: [None] * aln_len for aa in AA_ORDER}
    if not hotspots or aln_len <= 0:
        return scores_by_aa
    
    def parse_aa_list(val):
        if not val:
            return set()
        if isinstance(val, (list, tuple)):
            return set(str(v).strip() for v in val)
        # Handle format like "15:A,D,E,G,H,..." - strip count prefix
        val_str = str(val).replace(' ', '')
        if ':' in val_str:
            val_str = val_str.split(':', 1)[1]  # Take part after colon
        return set(val_str.split(','))
    
    def score_for_substitution(hp, aa):
        if not hp:
            return None
        benign = parse_aa_list(hp.get('benign_all') or hp.get('benign'))
        ambiguous = parse_aa_list(hp.get('ambiguous_all') or hp.get('ambiguous'))
        pathogenic = parse_aa_list(hp.get('pathogenic_all') or hp.get('pathogenic'))
        
        if aa in pathogenic:
            return 0.90
        if aa in ambiguous:
            return 0.55
        if aa in benign:
            return 0.15
        return None
    
    for col in range(1, aln_len + 1):
        resi = pos_by_col.get(col)
        hp = None
        if isinstance(resi, int) and 1 <= resi <= len(hotspots):
            hp = hotspots[resi - 1]
        for aa in AA_ORDER:
            scores_by_aa[aa][col - 1] = score_for_substitution(hp, aa)
    
    return scores_by_aa


# ============= PDBe Complex Fetching (with caching) =============

# Global cache for PDBe complexes (persists across pairs in same run)
_pdbe_complex_cache: Dict[str, Dict] = {}


def fetch_pdbe_best_structures(acc: str, max_results: int = 3) -> List[Dict[str, Any]]:
    """Get best PDB structures for a UniProt accession."""
    url = f"{PDBE_BASE}/mappings/best_structures/{acc.upper()}"
    data = http_get_json(url)
    if not data:
        return []
    
    rows = data.get(acc.upper()) or data.get(acc) or []
    if not isinstance(rows, list):
        return []
    
    def score(r):
        cov = float(r.get("coverage") or 0)
        res = float(r.get("resolution") or 9999)
        return (cov, -res)
    
    rows = sorted(rows, key=score, reverse=True)[:max_results]
    
    results = []
    for r in rows:
        results.append({
            "pdb_id": (r.get("pdb_id") or "").upper(),
            "chain_id": r.get("chain_id") or "",
            "coverage": float(r.get("coverage") or 0),
            "resolution": float(r.get("resolution")) if r.get("resolution") else None,
            "start": r.get("start"),
            "end": r.get("end"),
        })
    return results


def download_pdb_structure(pdb_id: str) -> Optional[Tuple[str, str]]:
    """Download PDB structure, return (format, base64_content). Uses file cache."""
    pid = pdb_id.lower()
    
    # Check file cache
    cache_file = CACHE_DIR / f"pdb_{pid}.json"
    if cache_file.exists():
        try:
            cached = json.loads(cache_file.read_text())
            return (cached["format"], cached["b64"])
        except:
            pass
    
    sources = [
        (f"https://files.rcsb.org/download/{pid}.pdb", "pdb"),
        (f"https://www.ebi.ac.uk/pdbe/entry-files/download/pdb{pid}.ent", "pdb"),
        (f"https://files.rcsb.org/download/{pid}.cif", "mmcif"),
        (f"https://www.ebi.ac.uk/pdbe/entry-files/download/{pid}.cif", "mmcif"),
    ]
    
    for url, fmt in sources:
        content = http_get_binary(url)
        if content:
            text = content.decode('utf-8', errors='ignore')[:100]
            if fmt == "pdb" and ("ATOM" in text or "HEADER" in text):
                b64 = base64.b64encode(content).decode('ascii')
                # Cache to file
                try:
                    CACHE_DIR.mkdir(parents=True, exist_ok=True)
                    cache_file.write_text(json.dumps({"format": fmt, "b64": b64}))
                except:
                    pass
                return (fmt, b64)
            if fmt == "mmcif" and ("_atom_site" in text or "data_" in text):
                b64 = base64.b64encode(content).decode('ascii')
                try:
                    CACHE_DIR.mkdir(parents=True, exist_ok=True)
                    cache_file.write_text(json.dumps({"format": fmt, "b64": b64}))
                except:
                    pass
                return (fmt, b64)
    
    return None


def fetch_pdbe_ligands(pdb_id: str) -> Tuple[Dict[str, int], str]:
    """Get ligand info for a PDB entry."""
    url = f"{PDBE_BASE}/pdb/entry/ligand_monomers/{pdb_id.lower()}"
    data = http_get_json(url)
    if not data:
        return {}, ""
    
    entry = data.get(pdb_id.lower()) or data.get(pdb_id.upper()) or []
    counts = {}
    for lig in entry:
        chem_id = lig.get("chem_comp_id") or ""
        if chem_id and chem_id.upper() not in ("HOH", "WAT"):
            counts[chem_id] = counts.get(chem_id, 0) + 1
    
    if counts:
        summary = ", ".join(f"{k}x{v}" if v > 1 else k for k, v in sorted(counts.items()))
    else:
        summary = ""
    
    return counts, summary


def fetch_pdbe_title(pdb_id: str) -> str:
    """Fetch the title/description of a PDB entry from PDBe summary API."""
    url = f"{PDBE_BASE}/pdb/entry/summary/{pdb_id.lower()}"
    data = http_get_json(url)
    if not data:
        return ""
    entry = data.get(pdb_id.lower(), [])
    if entry and isinstance(entry, list) and len(entry) > 0:
        return entry[0].get("title", "")
    return ""


def fetch_pdbe_complexes(acc_a: str, acc_b: str, max_per_protein: int = 3) -> List[Dict[str, Any]]:
    """Fetch PDBe complexes for both proteins. Caches downloaded structures."""
    global _pdbe_complex_cache
    
    complexes = []
    seen_pdb_ids = set()
    
    for acc in [acc_a, acc_b]:
        if not acc:
            continue
        
        structures = fetch_pdbe_best_structures(acc, max_per_protein)
        for s in structures:
            pdb_id = s["pdb_id"]
            if not pdb_id or pdb_id in seen_pdb_ids:
                continue
            seen_pdb_ids.add(pdb_id)
            
            # Check memory cache first
            if pdb_id in _pdbe_complex_cache:
                cached = _pdbe_complex_cache[pdb_id].copy()
                cached["source_acc"] = acc  # Update source
                complexes.append(cached)
                log(f"    PDBe {pdb_id}: cached")
                continue
            
            dl = download_pdb_structure(pdb_id)
            if not dl:
                log(f"    Failed to download {pdb_id}")
                continue
            
            coord_fmt, coord_b64 = dl
            ligand_counts, ligand_summary = fetch_pdbe_ligands(pdb_id)
            title = fetch_pdbe_title(pdb_id)

            complex_entry = {
                "source_acc": acc,
                "pdb_id": pdb_id,
                "chain_id": s["chain_id"],
                "chains": [s["chain_id"]] if s["chain_id"] else [],
                "coverage": s["coverage"],
                "resolution": s["resolution"],
                "coordFormat": coord_fmt,
                "coordB64": coord_b64,
                "ligandCounts": ligand_counts,
                "ligandSummary": ligand_summary,
                "title": title,
                "uniprot_start": s["start"],
                "uniprot_end": s["end"],
            }
            
            # Cache in memory
            _pdbe_complex_cache[pdb_id] = complex_entry.copy()
            
            complexes.append(complex_entry)
            log(f"    Added PDBe: {pdb_id} for {acc} (ligands: {ligand_summary or 'none'})")
    
    return complexes
