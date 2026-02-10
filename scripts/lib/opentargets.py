"""OpenTargets known drug loading for report generation."""

import csv
import glob
import re
from pathlib import Path
from typing import Any, Dict, List, Optional

from .config import INPUT_DIR, GENE_ID_MAP_CSV, log

KNOWN_DRUG_DIR = INPUT_DIR / "opentargets" / "known_drug"

# Cache: ensembl_id -> list of drug dicts
_drug_cache: Optional[Dict[str, List[Dict[str, Any]]]] = None
# Cache: uniprot_acc -> ensembl_id
_acc_to_ens: Optional[Dict[str, str]] = None


def _build_acc_to_ens_map() -> Dict[str, str]:
    """Build UniProt accession -> Ensembl ID mapping from gene_ids CSV."""
    global _acc_to_ens
    if _acc_to_ens is not None:
        return _acc_to_ens

    _acc_to_ens = {}
    if not GENE_ID_MAP_CSV.exists():
        return _acc_to_ens

    with open(GENE_ID_MAP_CSV) as f:
        reader = csv.DictReader(f)
        for row in reader:
            ens_id = row.get('ensembl_id', '').strip()
            acc_str = row.get('UniProt_Acc', '')
            if not ens_id or not acc_str:
                continue
            # Format: "[P11802, Q96BE9]" or "[Q00534]" - parse with regex
            accs = re.findall(r'[A-Z][A-Z0-9]{4,9}', acc_str)
            for acc in accs:
                _acc_to_ens[acc.strip()] = ens_id

    return _acc_to_ens


def _load_all_known_drugs() -> Dict[str, List[Dict[str, Any]]]:
    """Load all known drugs from Parquet files, indexed by target Ensembl ID."""
    global _drug_cache
    if _drug_cache is not None:
        return _drug_cache

    _drug_cache = {}
    parquet_files = sorted(glob.glob(str(KNOWN_DRUG_DIR / "*.parquet")))
    if not parquet_files:
        log("  OpenTargets: no known_drug parquet files found")
        return _drug_cache

    try:
        import pyarrow.parquet as pq
    except ImportError:
        log("  OpenTargets: pyarrow not available, skipping known drugs")
        return _drug_cache

    for f in parquet_files:
        table = pq.read_table(f, columns=[
            'drugId', 'targetId', 'prefName', 'drugType',
            'mechanismOfAction', 'phase', 'tradeNames',
        ])
        df = table.to_pandas()
        for _, row in df.iterrows():
            tid = row.get('targetId', '')
            if not tid:
                continue
            if tid not in _drug_cache:
                _drug_cache[tid] = []
            _drug_cache[tid].append({
                'drugId': row.get('drugId', ''),
                'name': row.get('prefName', ''),
                'type': row.get('drugType', ''),
                'moa': row.get('mechanismOfAction', ''),
                'phase': int(row['phase']) if row.get('phase') and row['phase'] == row['phase'] else None,
                'tradeNames': list(row['tradeNames']) if hasattr(row.get('tradeNames', None), '__iter__') and not isinstance(row.get('tradeNames'), str) else [],
            })

    log(f"  OpenTargets: loaded {sum(len(v) for v in _drug_cache.values())} drug entries for {len(_drug_cache)} targets")
    return _drug_cache


def get_known_drugs(ensembl_id: str) -> List[Dict[str, Any]]:
    """Get known drugs for a target by Ensembl gene ID.

    Returns deduplicated list of drugs sorted by phase (highest first).
    """
    cache = _load_all_known_drugs()
    raw = cache.get(ensembl_id, [])
    if not raw:
        return []

    # Deduplicate by drugId, keep highest phase
    seen = {}
    for d in raw:
        did = d['drugId']
        if did not in seen or (d.get('phase') or 0) > (seen[did].get('phase') or 0):
            seen[did] = d

    drugs = sorted(seen.values(), key=lambda x: -(x.get('phase') or 0))
    return drugs


def get_known_drugs_for_acc(acc: str) -> List[Dict[str, Any]]:
    """Get known drugs for a target by UniProt accession.

    Maps to Ensembl ID first, then looks up known drugs.
    """
    acc_map = _build_acc_to_ens_map()
    ens_id = acc_map.get(acc)
    if not ens_id:
        return []
    return get_known_drugs(ens_id)
