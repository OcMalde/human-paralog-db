"""Data loading functions for CSV files."""

import ast
from typing import Any, Dict, List, Optional, Set

import numpy as np
import pandas as pd

from .config import (
    FEATURES_CSV, GENE_LOC_CSV, ESSENTIAL_CSV, GENE_ID_MAP_CSV, log
)

# Caches
_features_df_cache = None
_gene_loc_cache = None
_essential_genes_cache = None
_gene_id_map_cache = None


def load_features_df() -> pd.DataFrame:
    global _features_df_cache
    if _features_df_cache is None:
        if FEATURES_CSV.exists():
            _features_df_cache = pd.read_csv(FEATURES_CSV, low_memory=False)
            log(f"Loaded {len(_features_df_cache)} pairs from features CSV")
        else:
            log(f"WARNING: Features CSV not found at {FEATURES_CSV}")
            _features_df_cache = pd.DataFrame()
    return _features_df_cache


def load_gene_locations() -> pd.DataFrame:
    global _gene_loc_cache
    if _gene_loc_cache is None:
        if GENE_LOC_CSV.exists():
            # Try comma first (notebook default), then tab if that fails
            try:
                _gene_loc_cache = pd.read_csv(GENE_LOC_CSV, dtype=str)
                if _gene_loc_cache.empty or len(_gene_loc_cache.columns) < 2:
                    _gene_loc_cache = pd.read_csv(GENE_LOC_CSV, sep='\t', dtype=str)
            except:
                _gene_loc_cache = pd.read_csv(GENE_LOC_CSV, sep='\t', dtype=str)
            log(f"Loaded {len(_gene_loc_cache)} gene locations, columns: {list(_gene_loc_cache.columns)}")
        else:
            _gene_loc_cache = pd.DataFrame()
    return _gene_loc_cache


def load_essential_genes() -> Set[str]:
    global _essential_genes_cache
    if _essential_genes_cache is None:
        if ESSENTIAL_CSV.exists():
            df = pd.read_csv(ESSENTIAL_CSV)
            # Column is 'Essentials' in Depmap CSV, format: "GENE1 (12345)"
            col = None
            for c in ['Essentials', 'gene', 'Gene', df.columns[0]]:
                if c in df.columns:
                    col = c
                    break
            if col:
                # Extract just the gene symbol (before space)
                _essential_genes_cache = set(
                    str(g).split(' ')[0] for g in df[col].dropna()
                )
                log(f"Loaded {len(_essential_genes_cache)} essential genes")
            else:
                _essential_genes_cache = set()
        else:
            _essential_genes_cache = set()
    return _essential_genes_cache


def load_gene_id_map() -> Dict[str, str]:
    global _gene_id_map_cache
    if _gene_id_map_cache is None:
        _gene_id_map_cache = {}
        if GENE_ID_MAP_CSV.exists():
            df = pd.read_csv(GENE_ID_MAP_CSV, dtype=str)
            entrez_col = symbol_col = None
            for col in df.columns:
                cl = col.lower()
                if 'entrez' in cl or 'ncbi' in cl:
                    entrez_col = col
                elif 'symbol' in cl or 'gene_name' in cl.replace(' ', '_'):
                    symbol_col = col
            if entrez_col and symbol_col:
                for _, row in df.iterrows():
                    eid_raw = str(row.get(entrez_col, '')).strip()
                    sym = str(row.get(symbol_col, '')).strip()
                    if not sym or sym == 'nan':
                        continue
                    # Parse Entrez ID - handle list format like [7105]
                    eid_raw = eid_raw.strip('[]').strip()
                    if not eid_raw or eid_raw == 'nan':
                        continue
                    # Handle comma-separated lists (take first)
                    for eid in eid_raw.split(','):
                        eid = eid.strip().strip("'\"")
                        if eid and eid != 'nan':
                            try:
                                eid = str(int(float(eid)))
                                _gene_id_map_cache[eid] = sym
                            except ValueError:
                                pass
            log(f"Loaded {len(_gene_id_map_cache)} gene ID mappings")
    return _gene_id_map_cache


def get_pair_row(pair_id: str) -> Optional[pd.Series]:
    df = load_features_df()
    if df.empty:
        return None
    for col in ['sorted_gene_pair', 'pair', 'pair_id', 'Pair', 'PAIR']:
        if col in df.columns:
            match = df[df[col] == pair_id]
            if not match.empty:
                return match.iloc[0]
    return None


def get_chromosome_info(gene_symbol: str) -> Dict[str, Any]:
    gene_locs = load_gene_locations()
    if gene_locs.empty:
        return {}
    
    match = None
    
    # Try matching by Gene name first
    if 'Gene name' in gene_locs.columns:
        match = gene_locs[gene_locs['Gene name'] == gene_symbol]
    
    # Try alternative column names
    if (match is None or match.empty):
        for col in ['gene_name', 'symbol', 'Symbol', 'SYMBOL', 'gene']:
            if col in gene_locs.columns:
                match = gene_locs[gene_locs[col] == gene_symbol]
                if match is not None and not match.empty:
                    break
    
    if match is None or match.empty:
        return {}
    
    row = match.iloc[0]
    
    # Try multiple possible column names for chromosome info
    def get_col(row, *names):
        for name in names:
            if name in row.index:
                val = row.get(name)
                if pd.notna(val) and str(val).strip():
                    return str(val)
        return 'NA'
    
    return {
        'chromosome': get_col(row, 'Chromosome/scaffold name', 'chromosome', 'chr', 'Chromosome'),
        'start': get_col(row, 'Gene start (bp)', 'start', 'Start', 'gene_start'),
        'end': get_col(row, 'Gene end (bp)', 'end', 'End', 'gene_end'),
        'strand': get_col(row, 'Strand', 'strand'),
    }


def compute_percentile(value: float, series: pd.Series) -> float:
    if pd.isna(value) or series.empty:
        return 50.0
    valid = series.dropna()
    if valid.empty:
        return 50.0
    return float((np.sum(valid < value) / len(valid)) * 100)


def get_conservation_percentiles(pair_row: Optional[pd.Series]) -> Dict[str, Dict[str, Any]]:
    df = load_features_df()
    if df.empty or pair_row is None:
        return {}
    
    metrics = {
        'min_sequence_identity': {'label': 'Sequence Identity', 'higher_is_more_conserved': True},
        'alntmscore_struct': {'label': 'TM-score (Structure)', 'higher_is_more_conserved': True},
        'fident_struct': {'label': 'Foldseek Identity', 'higher_is_more_conserved': True},
        'lddt_struct': {'label': 'lDDT Score', 'higher_is_more_conserved': True},
        'esm2_mean_of_residue_tokens_cosine': {'label': 'ESM2 Distance', 'higher_is_more_conserved': False},
        'ProtT5_per-protein_cosine': {'label': 'ProtT5 Distance', 'higher_is_more_conserved': False},
    }
    
    results = {}
    for col, info in metrics.items():
        if col in pair_row.index and col in df.columns:
            val = pair_row[col]
            pct = compute_percentile(val, df[col])
            display_pct = float(pct) if info['higher_is_more_conserved'] else float(100 - pct)
            display_pct = max(0.0, min(100.0, display_pct))
            results[col] = {
                'value': float(val) if pd.notna(val) else None,
                'percentile': display_pct,
                'raw_percentile': float(pct),
                'label': info['label'],
                'higher_is_more_conserved': info['higher_is_more_conserved'],
                'radar_value': display_pct,
                'direction_hint': 'Higher = more conserved' if info['higher_is_more_conserved'] else 'Lower = more conserved'
            }
    return results


def get_boxplot_data(metric_col: str, pair_value: Optional[float]) -> Dict[str, Any]:
    df = load_features_df()
    if df.empty or metric_col not in df.columns:
        return {}
    
    series = df[metric_col].dropna()
    if series.empty:
        return {}
    
    q1 = float(series.quantile(0.25))
    q2 = float(series.quantile(0.50))
    q3 = float(series.quantile(0.75))
    iqr = q3 - q1
    
    return {
        'min': float(series.min()),
        'max': float(series.max()),
        'q1': q1,
        'median': q2,
        'q3': q3,
        'whisker_low': max(float(series.min()), q1 - 1.5 * iqr),
        'whisker_high': min(float(series.max()), q3 + 1.5 * iqr),
        'pair_value': float(pair_value) if pair_value is not None and pd.notna(pair_value) else None,
        'n': int(len(series)),
    }


def parse_ppi_field(raw_value) -> Set[str]:
    if raw_value is None or (isinstance(raw_value, float) and pd.isna(raw_value)):
        return set()
    
    if isinstance(raw_value, (list, set, tuple)):
        values = list(raw_value)
    else:
        text = str(raw_value).strip()
        if not text:
            return set()
        try:
            parsed = ast.literal_eval(text)
            if isinstance(parsed, (list, set, tuple)):
                values = list(parsed)
            elif parsed is None:
                return set()
            else:
                values = [parsed]
        except:
            stripped = text.strip('{}[] ')
            if not stripped:
                return set()
            values = [frag.strip() for frag in stripped.split(',') if frag.strip()]
    
    partners = set()
    for val in values or []:
        if val is None:
            continue
        if isinstance(val, (int, np.integer)):
            partners.add(str(int(val)))
            continue
        sval = str(val).strip().strip("'\"")
        if not sval:
            continue
        try:
            partners.add(str(int(float(sval))))
        except:
            partners.add(sval)
    return partners


def enrich_partner_ids(partners: Set[str]) -> List[Dict[str, Optional[str]]]:
    if not partners:
        return []
    mapping = load_gene_id_map()
    
    result = []
    for pid in sorted(partners):
        symbol = mapping.get(pid, '')
        display = f"{symbol} ({pid})" if symbol and symbol != pid else pid
        result.append({'id': pid, 'symbol': symbol or None, 'display': display})
    return result


def build_ppi_network_info(pair_row: Optional[pd.Series]) -> Dict[str, Any]:
    if pair_row is None:
        return {}
    
    a1 = parse_ppi_field(pair_row.get('A1_ppi'))
    a2 = parse_ppi_field(pair_row.get('A2_ppi'))
    if not a1 and not a2:
        return {}
    
    shared = a1 & a2
    unique_a = a1 - shared
    unique_b = a2 - shared
    
    return {
        'shared': enrich_partner_ids(shared),
        'unique_gene1': enrich_partner_ids(unique_a),
        'unique_gene2': enrich_partner_ids(unique_b),
        'counts': {
            'shared': len(shared),
            'unique_gene1': len(unique_a),
            'unique_gene2': len(unique_b),
            'total_gene1': len(a1),
            'total_gene2': len(a2),
        }
    }
