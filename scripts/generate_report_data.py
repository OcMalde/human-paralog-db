#!/usr/bin/env python3
"""
generate_report_data.py

Pre-generates all report data for pairs in the database.
Uses modular library for clean separation of concerns.

Usage:
    python scripts/generate_report_data.py [pair_id]
    python scripts/generate_report_data.py --force  # Re-generate all
    python scripts/generate_report_data.py --threads 4  # Parallel processing
"""

import argparse
import base64
import json
import os
import sqlite3
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

# Add scripts directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from lib.config import DB_PATH, CACHE_DIR, AM_MODES, log, STRUCTURES_DIR
from lib.data_loading import (
    get_pair_row, get_chromosome_info, load_essential_genes,
    get_conservation_percentiles, get_boxplot_data, build_ppi_network_info
)
from lib.api_fetching import (
    fetch_am_hotspots, build_am_matrix_scores, fetch_pdbe_complexes
)
from lib.am_processing import compute_am_rects_all_modes
from lib.domain_parsing import (
    parse_ted_domains, parse_ebi_domains, parse_cavities,
    rects_on_alignment, get_annotation_from_db, load_cavities_for_acc
)
from lib.seq_align import needleman_wunsch, compute_identity, alignment_to_columns

# Foldseek binary
FOLDSEEK = os.environ.get("FOLDSEEK", "foldseek")


def slice_pdb_range(pdb_text: str, start: int, end: int, chain_id: str = "A") -> str:
    """Extract a residue range from PDB text."""
    lines = []
    for ln in pdb_text.splitlines():
        if ln.startswith("ATOM") or ln.startswith("HETATM"):
            try:
                res_num = int(ln[22:26].strip())
                if start <= res_num <= end:
                    # Update chain ID
                    ln = ln[:21] + chain_id + ln[22:]
                    lines.append(ln)
            except:
                pass
    return "\n".join(lines) + "\n"


def extend_range_for_foldseek(start: int, end: int, prot_len: int, target_min: int = 50) -> Tuple[int, int]:
    """Extend range to meet minimum size for foldseek."""
    current = end - start + 1
    if current >= target_min:
        return start, end
    
    needed = target_min - current
    extend_each = needed // 2 + 1
    new_start = max(1, start - extend_each)
    new_end = min(prot_len, end + extend_each)
    return new_start, new_end


def run_foldseek_domain(pdb1_text: str, pdb2_text: str) -> Optional[Tuple[str, str, float, float, List[float], List[float]]]:
    """Run foldseek on domain PDB texts, return (qaln, taln, tm, fident, U, T)."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpd = Path(tmpdir)
        f1 = tmpd / "domA.pdb"
        f2 = tmpd / "domB.pdb"
        f1.write_text(pdb1_text)
        f2.write_text(pdb2_text)
        
        out_tsv = tmpd / "result.tsv"
        
        args = [
            FOLDSEEK, "easy-search",
            str(f1), str(f2),
            str(out_tsv), str(tmpd / "tmp"),
            "-a",
            "--max-seqs", "1",
            "--format-output", "query,target,fident,alntmscore,qaln,taln,u,t",
        ]
        
        try:
            result = subprocess.run(args, capture_output=True, timeout=60)
            if result.returncode != 0:
                return None
            
            if not out_tsv.exists() or out_tsv.stat().st_size == 0:
                return None
            
            line = out_tsv.read_text().strip().split('\n')[0]
            cols = line.split('\t')
            
            if len(cols) < 8:
                return None
            
            fident = float(cols[2]) if cols[2] else None
            tm = float(cols[3]) if cols[3] else None
            qaln = cols[4]
            taln = cols[5]
            U = [float(x) for x in cols[6].split(",")] if cols[6] else None
            T = [float(x) for x in cols[7].split(",")] if cols[7] else None
            
            return qaln, taln, tm, fident, U, T
        except:
            return None


def apply_transform_to_pdb(pdb_text: str, U: List[float], T: List[float]) -> str:
    """Apply rotation (U) and translation (T) to PDB coordinates."""
    if not U or not T or len(U) != 9 or len(T) != 3:
        return pdb_text
    
    lines = []
    for ln in pdb_text.splitlines():
        if ln.startswith(("ATOM", "HETATM")):
            try:
                x = float(ln[30:38])
                y = float(ln[38:46])
                z = float(ln[46:54])
                
                # Apply rotation and translation
                nx = U[0]*x + U[1]*y + U[2]*z + T[0]
                ny = U[3]*x + U[4]*y + U[5]*z + T[1]
                nz = U[6]*x + U[7]*y + U[8]*z + T[2]
                
                ln = f"{ln[:30]}{nx:8.3f}{ny:8.3f}{nz:8.3f}{ln[54:]}"
            except:
                pass
        lines.append(ln)
    return "\n".join(lines)


def compute_dam_percent(qaln: str, taln: str, startA: int, startB: int, 
                         bfA: List[float], bfB: List[float], thresh: float = 0.30) -> Optional[float]:
    """Compute difference in AlphaMissense score for aligned positions."""
    diffs = []
    qi = startA - 1
    ti = startB - 1
    
    for qc, tc in zip(qaln, taln):
        if qc not in "-.":
            qi += 1
        if tc not in "-.":
            ti += 1
        
        if qc not in "-." and tc not in "-.":
            try:
                vA = bfA[qi - 1] if 0 <= qi - 1 < len(bfA) else None
                vB = bfB[ti - 1] if 0 <= ti - 1 < len(bfB) else None
                if vA is not None and vB is not None:
                    diffs.append(abs(vA - vB))
            except:
                pass
    
    if not diffs:
        return None
    
    above_thresh = sum(1 for d in diffs if d >= thresh)
    return 100.0 * above_thresh / len(diffs)


def generate_domain_alignments(
    domsA: List[Dict], domsB: List[Dict],
    pdb_a_text: str, pdb_b_text: str,
    lenA: int, lenB: int,
    bfA: List[float], bfB: List[float],
    max_doms: int = 6
) -> List[Dict[str, Any]]:
    """Generate domain×domain sub-alignments using foldseek."""
    
    def dom_candidates(doms: List[Dict], maxn: int) -> List[Dict]:
        """Select top domains by size."""
        cands = []
        for d in doms:
            if d.get("type") not in ("Domain", "TED"):
                continue
            s, e = d.get("start", 0), d.get("end", 0)
            if e > s:
                cands.append((e - s, id(d), d))  # Use id() as tiebreaker to avoid dict comparison
        cands.sort(key=lambda x: x[0], reverse=True)
        return [x[2] for x in cands[:maxn]]
    
    candA = dom_candidates(domsA, max_doms)
    candB = dom_candidates(domsB, max_doms)
    
    if not candA or not candB:
        return []
    
    log(f"  Domain×domain: {len(candA)} × {len(candB)} candidates")
    
    domPairs = []
    for da in candA:
        for db in candB:
            qa, qb = int(da["start"]), int(da["end"])
            ta, tb = int(db["start"]), int(db["end"])
            
            qa_ext, qb_ext = extend_range_for_foldseek(qa, qb, lenA, target_min=50)
            ta_ext, tb_ext = extend_range_for_foldseek(ta, tb, lenB, target_min=50)
            
            A_txt = slice_pdb_range(pdb_a_text, qa_ext, qb_ext, "A")
            B_txt = slice_pdb_range(pdb_b_text, ta_ext, tb_ext, "B")
            
            if len(A_txt.splitlines()) < 10 or len(B_txt.splitlines()) < 10:
                continue
            
            result = run_foldseek_domain(A_txt, B_txt)
            if not result:
                continue
            
            q_aln, t_aln, tm, fid, U, T = result
            
            # Apply transform to B
            B_transformed = apply_transform_to_pdb(B_txt, U, T) if U and T else B_txt
            sub_pdb = A_txt + B_transformed + "END\n"
            sub_b64 = base64.b64encode(sub_pdb.encode("utf-8")).decode("ascii")
            
            dam_pct = compute_dam_percent(q_aln, t_aln, qa_ext, ta_ext, bfA, bfB, thresh=0.30)
            
            dom_label_A = da.get("label") or da.get("name") or "A-dom"
            dom_label_B = db.get("label") or db.get("name") or "B-dom"
            
            extA = (qa_ext != qa or qb_ext != qb)
            extB = (ta_ext != ta or tb_ext != tb)
            
            if extA:
                dom_label_A = f"{dom_label_A} [ext:{qa_ext}-{qb_ext}]"
            if extB:
                dom_label_B = f"{dom_label_B} [ext:{ta_ext}-{tb_ext}]"
            
            domPairs.append({
                "Aname": dom_label_A,
                "Arng": f"{qa}-{qb}",
                "ArngExt": f"{qa_ext}-{qb_ext}" if extA else None,
                "Bname": dom_label_B,
                "Brng": f"{ta}-{tb}",
                "BrngExt": f"{ta_ext}-{tb_ext}" if extB else None,
                "tm": tm,
                "fident": (fid * 100.0 if fid is not None and fid <= 1.0 else fid),
                "damPct": dam_pct,
                "pdb64": sub_b64,
            })
    
    # Sort by fident descending
    domPairs.sort(key=lambda r: (
        r["fident"] if r["fident"] is not None else -1,
        r["tm"] if r["tm"] is not None else -1
    ), reverse=True)
    
    log(f"  Generated {len(domPairs)} domain×domain sub-alignments")
    return domPairs


def get_pdb_text_for_acc(conn: sqlite3.Connection, acc: str) -> Optional[str]:
    """Get PDB text for a protein from database or file."""
    # Try database first
    prot = conn.execute("SELECT pdb_b64 FROM proteins WHERE acc = ?", (acc,)).fetchone()
    if prot and prot["pdb_b64"]:
        try:
            return base64.b64decode(prot["pdb_b64"]).decode("utf-8", errors="ignore")
        except:
            pass
    
    # Try file system
    for path in [
        STRUCTURES_DIR / f"AF-{acc}-F1-AM.pdb",
        STRUCTURES_DIR / f"AF-{acc}-F1-AM_v4.pdb",
    ]:
        if path.exists():
            return path.read_text()
    
    return None


# Standard amino acid 3-letter to 1-letter mapping
AA3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}


def extract_sequence_from_pdb(pdb_text: str) -> str:
    """Extract amino acid sequence from PDB text."""
    if not pdb_text:
        return ""

    residues = {}
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            res_name = line[17:20].strip()
            res_num = int(line[22:26].strip())
            if res_name in AA3TO1 and res_num not in residues:
                residues[res_num] = AA3TO1[res_name]

    if not residues:
        return ""

    # Build sequence from sorted residue numbers
    seq = []
    for i in range(min(residues.keys()), max(residues.keys()) + 1):
        seq.append(residues.get(i, 'X'))

    return ''.join(seq)


def compute_sequence_alignment(seq1: str, seq2: str) -> Optional[Dict[str, Any]]:
    """Compute sequence alignment and return data for client."""
    if not seq1 or not seq2:
        return None

    try:
        aln1, aln2, score = needleman_wunsch(seq1, seq2)
        identity = compute_identity(aln1, aln2)
        aligned_cols = alignment_to_columns(aln1, aln2)

        return {
            'qaln': aln1,
            'taln': aln2,
            'aligned_cols': aligned_cols,
            'identity': identity,
            'score': score,
        }
    except Exception as e:
        log(f"  WARNING: Sequence alignment failed: {e}")
        return None


def generate_report_data(pair_id: str, conn: sqlite3.Connection) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Generate complete DATA_OBJ and SUMMARY for a pair."""
    
    conn.row_factory = sqlite3.Row
    pair = conn.execute("SELECT * FROM pairs WHERE pair_id = ?", (pair_id,)).fetchone()
    if not pair:
        raise ValueError(f"Pair {pair_id} not found in database")
    
    gene_a = pair["gene_a"]
    gene_b = pair["gene_b"]
    acc_a = pair["acc_a"]
    acc_b = pair["acc_b"]
    
    log(f"  {gene_a} ({acc_a}) vs {gene_b} ({acc_b})")
    
    # Get alignment
    aln = conn.execute("SELECT * FROM alignments WHERE pair_id = ?", (pair_id,)).fetchone()
    qaln = aln["qaln"] if aln else ""
    taln = aln["taln"] if aln else ""
    tm_score = float(aln["tm_score"]) if aln and aln["tm_score"] else None
    fident = float(aln["fident"]) if aln and aln["fident"] else None
    
    alnLen = len(qaln)
    log(f"  Alignment: {alnLen} cols, TM={tm_score}")
    
    # Get proteins
    prot_a = conn.execute("SELECT * FROM proteins WHERE acc = ?", (acc_a,)).fetchone()
    prot_b = conn.execute("SELECT * FROM proteins WHERE acc = ?", (acc_b,)).fetchone()
    
    lenA = prot_a["length"] if prot_a else 0
    lenB = prot_b["length"] if prot_b else 0
    log(f"  Lengths: {gene_a}={lenA}, {gene_b}={lenB}")
    
    # Get bfactors (AM values)
    bfA = json.loads(prot_a["bfactors_json"]) if prot_a and prot_a["bfactors_json"] else []
    bfB = json.loads(prot_b["bfactors_json"]) if prot_b and prot_b["bfactors_json"] else []
    
    # Build column position maps
    qpos_by_col = {}
    tpos_by_col = {}
    aligned_cols = []
    qi = ti = 0
    for col in range(1, alnLen + 1):
        qc = qaln[col-1] if col-1 < len(qaln) else '-'
        tc = taln[col-1] if col-1 < len(taln) else '-'
        qpos = tpos = None
        if qc not in "-.":
            qi += 1
            qpos = qi
        if tc not in "-.":
            ti += 1
            tpos = ti
        qpos_by_col[col] = qpos
        tpos_by_col[col] = tpos
        aligned_cols.append([col, qpos, tpos])
    
    # Get annotations from database
    ted_a = get_annotation_from_db(conn, acc_a, "ted")
    ted_b = get_annotation_from_db(conn, acc_b, "ted")
    ebi_a = get_annotation_from_db(conn, acc_a, "ebi")
    ebi_b = get_annotation_from_db(conn, acc_b, "ebi")
    
    # Parse domains
    doms_ted_a = parse_ted_domains(ted_a)
    doms_ted_b = parse_ted_domains(ted_b)
    doms_ebi_a = parse_ebi_domains(ebi_a, "Domain")
    doms_ebi_b = parse_ebi_domains(ebi_b, "Domain")
    disorder_a = parse_ebi_domains(ebi_a, "Disordered")
    disorder_b = parse_ebi_domains(ebi_b, "Disordered")
    
    # Load cavities from local AF-all_cavities data (if available)
    cavities_a = load_cavities_for_acc(acc_a)
    cavities_b = load_cavities_for_acc(acc_b)
    
    # Combine all domains for domain×domain alignments
    domsA_all = doms_ebi_a + doms_ted_a + disorder_a + cavities_a
    domsB_all = doms_ebi_b + doms_ted_b + disorder_b + cavities_b
    
    log(f"  Domains: {gene_a}={len(domsA_all)} (TED:{len(doms_ted_a)}, EBI:{len(doms_ebi_a)}, Cav:{len(cavities_a)})")
    log(f"  Domains: {gene_b}={len(domsB_all)} (TED:{len(doms_ted_b)}, EBI:{len(doms_ebi_b)}, Cav:{len(cavities_b)})")
    
    # Compute alignment rects for domains
    # IMPORTANT: domA/B_rects should ONLY include EBI domains, NOT TED (TED has its own track)
    domA_rects = rects_on_alignment([d for d in doms_ebi_a if d.get("type") == "Domain"], 'A', "#2ca02c", aligned_cols, alnLen)
    domB_rects = rects_on_alignment([d for d in doms_ebi_b if d.get("type") == "Domain"], 'B', "#2ca02c", aligned_cols, alnLen)
    disorderA_rects = rects_on_alignment(disorder_a, 'A', "#9467bd", aligned_cols, alnLen)
    disorderB_rects = rects_on_alignment(disorder_b, 'B', "#9467bd", aligned_cols, alnLen)
    # TED domains have their own separate track
    tedA_rects = rects_on_alignment(doms_ted_a, 'A', "#d62728", aligned_cols, alnLen)
    tedB_rects = rects_on_alignment(doms_ted_b, 'B', "#d62728", aligned_cols, alnLen)
    # Cavities have their own track
    cavA_rects = rects_on_alignment(cavities_a, 'A', "#ff7d45", aligned_cols, alnLen)
    cavB_rects = rects_on_alignment(cavities_b, 'B', "#ff7d45", aligned_cols, alnLen)
    
    # Fetch AM hotspots for substitution matrix
    log(f"  Fetching AM hotspots...")
    hotspotsA = fetch_am_hotspots(acc_a, lenA) if lenA > 0 else []
    hotspotsB = fetch_am_hotspots(acc_b, lenB) if lenB > 0 else []
    
    matA_scores_raw = build_am_matrix_scores(hotspotsA, qpos_by_col, alnLen)
    matB_scores_raw = build_am_matrix_scores(hotspotsB, tpos_by_col, alnLen)
    
    # Compute AM rects for all modes
    am_rects = compute_am_rects_all_modes(
        bfA, bfB, aligned_cols, alnLen, matA_scores_raw, matB_scores_raw
    )
    
    # Fetch PDBe complexes
    log(f"  Fetching PDBe complexes...")
    pdbe_complexes = fetch_pdbe_complexes(acc_a, acc_b)
    log(f"  Found {len(pdbe_complexes)} PDBe complexes")
    
    # Generate domain×domain alignments
    pdb_a_text = get_pdb_text_for_acc(conn, acc_a)
    pdb_b_text = get_pdb_text_for_acc(conn, acc_b)
    
    if pdb_a_text and pdb_b_text and (domsA_all or domsB_all):
        domPairs = generate_domain_alignments(
            domsA_all, domsB_all,
            pdb_a_text, pdb_b_text,
            lenA, lenB, bfA, bfB
        )
    else:
        domPairs = []
        log(f"  Skipping domain×domain (no PDB or no domains)")

    # Extract full sequences from PDB and compute sequence alignment
    seq1 = extract_sequence_from_pdb(pdb_a_text) if pdb_a_text else ""
    seq2 = extract_sequence_from_pdb(pdb_b_text) if pdb_b_text else ""
    log(f"  Sequences: {gene_a}={len(seq1)} aa, {gene_b}={len(seq2)} aa")

    seq_alignment = None
    if seq1 and seq2:
        log(f"  Computing sequence alignment...")
        seq_alignment = compute_sequence_alignment(seq1, seq2)
        if seq_alignment:
            log(f"  Sequence alignment: {len(seq_alignment['qaln'])} cols, identity={seq_alignment['identity']*100:.1f}%")

    # Build DATA_OBJ
    DATA_OBJ = {
        "PAIR": pair_id,
        "g1": gene_a, "g2": gene_b,
        "a1": acc_a, "a2": acc_b,
        "qaln": qaln, "taln": taln,
        "tm": tm_score,
        "fident": fident,
        "lenA": lenA,
        "lenB": lenB,
        "domainsA": domsA_all,
        "domainsB": domsB_all,
        "domA_alnRects": domA_rects,
        "domB_alnRects": domB_rects,
        "disorderA_alnRects": disorderA_rects,
        "disorderB_alnRects": disorderB_rects,
        "tedA_alnRects": tedA_rects,
        "tedB_alnRects": tedB_rects,
        "cavA_alnRects": cavA_rects,
        "cavB_alnRects": cavB_rects,
        "qposByCol": qpos_by_col,
        "tposByCol": tpos_by_col,
        "amModes": AM_MODES,
        "bfactorsA": bfA,
        "bfactorsB": bfB,
        "domPairs": domPairs,
        "pdbeComplexes": pdbe_complexes,
        "seq1": seq1,
        "seq2": seq2,
        "seqAlignment": seq_alignment,
        **am_rects,  # Unpack all AM rect data
    }
    
    # Build SUMMARY
    pair_row = get_pair_row(pair_id)
    essential_genes = load_essential_genes()
    
    if pair_row is not None:
        log(f"  Found pair in features CSV")
    else:
        log(f"  WARNING: Pair not in CSV - summary incomplete")
    
    gene1_info = {
        'symbol': gene_a,
        'uniprot': acc_a,
        'is_essential': gene_a in essential_genes,
        'chromosome': get_chromosome_info(gene_a),
    }
    
    gene2_info = {
        'symbol': gene_b,
        'uniprot': acc_b,
        'is_essential': gene_b in essential_genes,
        'chromosome': get_chromosome_info(gene_b),
    }
    
    log(f"  {gene_a}: essential={gene1_info['is_essential']}, chr={gene1_info['chromosome'].get('chromosome', 'NA')}")
    log(f"  {gene_b}: essential={gene2_info['is_essential']}, chr={gene2_info['chromosome'].get('chromosome', 'NA')}")
    
    pair_info = {}
    if pair_row is not None:
        def safe_bool(val):
            if pd.isna(val): return None
            if isinstance(val, bool): return val
            if isinstance(val, str): return val.lower() in ('true', '1', 'yes')
            return bool(val)
        
        def safe_int(val):
            if pd.isna(val): return None
            try: return int(val)
            except: return None
        
        pair_info = {
            'wgd': safe_bool(pair_row.get('WGD')),
            'family_size': safe_int(pair_row.get('family_size')),
            'closest': safe_bool(pair_row.get('closest')),
            'same_chr': safe_bool(pair_row.get('same_chr')),
            'interact_bioplex': safe_bool(pair_row.get('interact')),
            'n_shared_ppi': safe_int(pair_row.get('n_shared_ppi')) or 0,
            'ppi_network': build_ppi_network_info(pair_row),
        }
    
    conservation = get_conservation_percentiles(pair_row)
    log(f"  Conservation metrics: {list(conservation.keys())}")
    
    boxplots = {}
    for metric_col, metric_info in conservation.items():
        boxplots[metric_col] = get_boxplot_data(metric_col, metric_info.get('value'))
    
    SUMMARY = {
        'gene1': gene1_info,
        'gene2': gene2_info,
        'pair': pair_info,
        'conservation': conservation,
        'boxplots': boxplots,
    }
    
    return DATA_OBJ, SUMMARY


def store_report_data(conn: sqlite3.Connection, pair_id: str, data_obj: Dict, summary: Dict):
    """Store generated report data in database."""
    conn.execute("""
        CREATE TABLE IF NOT EXISTS report_data (
            pair_id TEXT PRIMARY KEY,
            data_json TEXT,
            summary_json TEXT,
            generated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    """)
    
    conn.execute("""
        INSERT OR REPLACE INTO report_data (pair_id, data_json, summary_json)
        VALUES (?, ?, ?)
    """, (pair_id, json.dumps(data_obj), json.dumps(summary)))
    conn.commit()


def main():
    parser = argparse.ArgumentParser(description='Generate report data for pairs')
    parser.add_argument('pair_id', nargs='?', help='Specific pair to generate')
    parser.add_argument('--force', action='store_true', help='Re-generate all (ignore existing)')
    parser.add_argument('--am-threads', type=int, default=8, help='Threads for AM hotspot fetching')
    args = parser.parse_args()
    
    if not DB_PATH.exists():
        print(f"ERROR: Database not found at {DB_PATH}")
        sys.exit(1)
    
    # Create cache directory
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    
    # Ensure report_data table exists
    conn.execute("""
        CREATE TABLE IF NOT EXISTS report_data (
            pair_id TEXT PRIMARY KEY,
            data_json TEXT,
            summary_json TEXT,
            generated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    """)
    
    if args.pair_id:
        pairs = [args.pair_id]
    else:
        pairs = [r[0] for r in conn.execute("SELECT pair_id FROM pairs").fetchall()]
    
    # Check which are already done (unless --force)
    if not args.force and not args.pair_id:
        existing = set(r[0] for r in conn.execute("SELECT pair_id FROM report_data").fetchall())
        to_process = [p for p in pairs if p not in existing]
        log(f"Found {len(pairs)} pairs, {len(existing)} already processed, {len(to_process)} remaining")
        pairs = to_process
    
    if not pairs:
        log("Nothing to process!")
        conn.close()
        return
    
    log(f"Processing {len(pairs)} pair(s) (AM threads: {args.am_threads})...")
    
    success = 0
    errors = 0
    
    for i, pair_id in enumerate(pairs):
        try:
            log(f"[{i+1}/{len(pairs)}] {pair_id}")
            data_obj, summary = generate_report_data(pair_id, conn)
            store_report_data(conn, pair_id, data_obj, summary)
            success += 1
            log(f"  Done!")
        except Exception as e:
            errors += 1
            log(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()
    
    conn.close()
    log(f"All done! Success: {success}, Errors: {errors}")


if __name__ == "__main__":
    main()
