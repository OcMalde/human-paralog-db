"""Data loading functions for CSV files."""

import ast
import time
from typing import Any, Dict, List, Optional, Set

import numpy as np
import pandas as pd
import requests

from .config import (
    FEATURES_CSV, GENE_LOC_CSV, ESSENTIAL_CSV, GENE_ID_MAP_CSV, SL_DATA_CSV, log
)

# Caches
_features_df_cache = None
_gene_loc_cache = None
_essential_genes_cache = None
_gene_id_map_cache = None
_sl_data_cache = None


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


def load_sl_data() -> pd.DataFrame:
    global _sl_data_cache
    if _sl_data_cache is None:
        if SL_DATA_CSV.exists():
            _sl_data_cache = pd.read_csv(SL_DATA_CSV, low_memory=False)
            log(f"Loaded {len(_sl_data_cache)} pairs from SL data CSV")
        else:
            log(f"WARNING: SL data CSV not found at {SL_DATA_CSV}")
            _sl_data_cache = pd.DataFrame()
    return _sl_data_cache


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


def compute_rank_relative_percentile(
    value: float,
    target_rank: float,
    df: pd.DataFrame,
    metric_col: str,
    rank_col: str,
    bucket_size: int = 10
) -> float:
    """Compute percentile of a metric within pairs of similar rank.

    Uses a bucket approach: pairs where rank is within ±bucket_size of target_rank.
    """
    if pd.isna(value) or pd.isna(target_rank) or df.empty:
        return 50.0

    if metric_col not in df.columns or rank_col not in df.columns:
        return 50.0

    # Filter to pairs with similar rank (±bucket_size)
    rank_min = max(0, target_rank - bucket_size)
    rank_max = target_rank + bucket_size

    similar_rank_mask = (df[rank_col] >= rank_min) & (df[rank_col] <= rank_max)
    similar_df = df[similar_rank_mask]

    if similar_df.empty or len(similar_df) < 5:
        # Fall back to overall percentile if bucket too small
        return compute_percentile(value, df[metric_col])

    return compute_percentile(value, similar_df[metric_col])


def get_similarity_search_percentiles(pair_row: Optional[pd.Series]) -> Dict[str, Dict[str, Any]]:
    """Get similarity search metrics with percentiles for visualization.

    All 6 metrics (3 structure + 3 sequence) are pair-level comparisons.
    - rank: Position in search results (lower = more similar)
    - selfSP: Number of human proteins between A and B (lower = closer)
    - taxid: Number of species between A and B (lower = closer)

    For selfSP and taxid, also computes rank-relative percentiles
    (comparing against pairs with similar rank).

    Percentile interpretation: low value = empty bar (close paralogs)
    """
    df = load_features_df()
    if df.empty or pair_row is None:
        return {}

    results = {}

    # Process each suffix (struct and seq)
    for suffix in ['_struct', '_seq']:
        rank_col = 'rank' + suffix
        selfsp_col = 'selfSP' + suffix
        taxid_col = 'taxid' + suffix

        # Get rank value for this pair (needed for rank-relative percentiles)
        rank_val = pair_row.get(rank_col) if rank_col in pair_row.index else None

        # Process rank metric
        if rank_col in pair_row.index and rank_col in df.columns:
            val = pair_row[rank_col]
            pct = compute_percentile(val, df[rank_col])
            # For rank: low value = good (close paralogs), so percentile shows position
            # Low rank = low percentile (empty bar = close)
            results[rank_col] = {
                'value': float(val) if pd.notna(val) else None,
                'percentile': float(pct),  # Low value = low percentile = empty bar
                'raw_percentile': float(pct),
                'label': rank_col,
                'higher_is_better': False,
                'radar_value': float(100 - pct),  # Inverted for radar (higher = better)
            }

        # Process selfSP metric (with rank-relative percentile)
        if selfsp_col in pair_row.index and selfsp_col in df.columns:
            val = pair_row[selfsp_col]
            pct_overall = compute_percentile(val, df[selfsp_col])
            pct_rank_rel = compute_rank_relative_percentile(
                val, rank_val, df, selfsp_col, rank_col
            ) if pd.notna(rank_val) else pct_overall

            # For selfSP: low value = good (fewer humans between), low percentile = empty bar
            results[selfsp_col] = {
                'value': float(val) if pd.notna(val) else None,
                'percentile': float(pct_overall),  # Overall percentile
                'percentile_rank_relative': float(pct_rank_rel),  # Vs similar rank pairs
                'raw_percentile': float(pct_overall),
                'label': selfsp_col,
                'higher_is_better': False,  # Lower = closer paralogs
                'radar_value': float(100 - pct_overall),  # Inverted for radar
            }

        # Process taxid metric (with rank-relative percentile)
        if taxid_col in pair_row.index and taxid_col in df.columns:
            val = pair_row[taxid_col]
            pct_overall = compute_percentile(val, df[taxid_col])
            pct_rank_rel = compute_rank_relative_percentile(
                val, rank_val, df, taxid_col, rank_col
            ) if pd.notna(rank_val) else pct_overall

            # For taxid: low value = good (fewer species between), low percentile = empty bar
            results[taxid_col] = {
                'value': float(val) if pd.notna(val) else None,
                'percentile': float(pct_overall),  # Overall percentile
                'percentile_rank_relative': float(pct_rank_rel),  # Vs similar rank pairs
                'raw_percentile': float(pct_overall),
                'label': taxid_col,
                'higher_is_better': False,  # Lower = closer paralogs
                'radar_value': float(100 - pct_overall),  # Inverted for radar
            }

    return results


def get_family_feature_percentiles(pair_row: Optional[pd.Series]) -> Dict[str, Dict[str, Any]]:
    """Get family feature metrics with percentiles for radar chart.

    Uses ratio (r_*) columns for proper comparison across pairs.
    """
    df = load_features_df()
    if df.empty or pair_row is None:
        return {}

    results = {}

    # Use ratio columns (r_*) for proper comparison - exact column names as labels
    metrics = {
        'rmean_shared_aa_withFamily': {'higher_is_better': True},
        'rmean_shared_aa_pairExclusive': {'higher_is_better': True},
        'rmean_shared_aa_onlyWithFamily': {'higher_is_better': True},
        'clustalo_r_shared_aa_withFamily': {'higher_is_better': True},
        'clustalo_r_shared_aa_pairExclusive': {'higher_is_better': True},
        'clustalo_r_shared_aa_onlyWithFamily': {'higher_is_better': True},
        'clustalo_r_sum_specific': {'higher_is_better': True},
    }

    for col, info in metrics.items():
        if col in pair_row.index and col in df.columns:
            val = pair_row[col]
            pct = compute_percentile(val, df[col])
            display_pct = float(pct) if info['higher_is_better'] else float(100 - pct)
            display_pct = max(0.0, min(100.0, display_pct))
            results[col] = {
                'value': float(val) if pd.notna(val) else None,
                'percentile': display_pct,
                'raw_percentile': float(pct),
                'label': col,  # Use exact column name
                'higher_is_better': info['higher_is_better'],
                'radar_value': display_pct,
            }

    return results


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


def get_sl_row(pair_id: str) -> Optional[pd.Series]:
    """Get a row from the SL data by pair ID."""
    sl_df = load_sl_data()
    if sl_df.empty:
        return None
    for col in ['sorted_gene_pair', 'pair', 'pair_id', 'Pair', 'PAIR']:
        if col in sl_df.columns:
            match = sl_df[sl_df[col] == pair_id]
            if not match.empty:
                return match.iloc[0]
    return None


def get_sl_functional_overlap(pair_id: str, pair_row: Optional[pd.Series]) -> Dict[str, Any]:
    """Get synthetic lethality and functional overlap data for a pair.

    Combines:
    - SL flag (SL_consensus only - lenient is similar)
    - Only screens with hits (value > 0)
    - GO similarity scores from features CSV (BPO, CCO, MFO)
    """
    result = {
        'is_sl': None,
        'sl_screens': [],  # Only screens with hits
        'go_similarity': {},
    }

    # Get SL data
    sl_row = get_sl_row(pair_id)
    if sl_row is not None:
        # Only use SL_consensus flag (lenient is redundant)
        if 'SL_consensus' in sl_row.index:
            val = sl_row['SL_consensus']
            result['is_sl'] = bool(val) if pd.notna(val) else None

        # Only include screens where pair was found SL (value > 0)
        screen_cols = {
            'n_SL_thompson': 'Thompson',
            'n_SL_dede': 'Dede',
            'n_SL_parrish': 'Parrish',
            'n_SL_chymera': 'Chymera',
            'n_SL_ito': 'Ito',
            'n_SL_TCGA_DepMap': 'TCGA/DepMap',
        }
        for col, label in screen_cols.items():
            if col in sl_row.index:
                val = sl_row[col]
                if pd.notna(val) and int(val) > 0:
                    result['sl_screens'].append(label)

    # Get GO similarity from features CSV
    if pair_row is not None:
        df = load_features_df()
        go_metrics = {
            'BPO': {'label': 'Biological Process'},
            'CCO': {'label': 'Cellular Component'},
            'MFO': {'label': 'Molecular Function'},
        }
        for col, info in go_metrics.items():
            if col in pair_row.index and col in df.columns:
                val = pair_row[col]
                pct = compute_percentile(val, df[col])
                result['go_similarity'][col] = {
                    'label': info['label'],
                    'value': float(val) if pd.notna(val) else None,
                    'percentile': float(pct),
                }

    return result


# UniProt description cache (persists across calls)
_uniprot_cache: Dict[str, Dict[str, str]] = {}


def fetch_uniprot_description(accession: str) -> Dict[str, str]:
    """Fetch protein description from UniProt API.

    Returns dict with 'name' (recommended name) and 'function' (function description).
    """
    global _uniprot_cache

    if accession in _uniprot_cache:
        return _uniprot_cache[accession]

    result = {'name': '', 'function': ''}

    try:
        url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
        resp = requests.get(url, timeout=10)
        if resp.status_code == 200:
            data = resp.json()

            # Get recommended protein name
            prot_desc = data.get('proteinDescription', {})
            rec_name = prot_desc.get('recommendedName', {})
            if rec_name:
                full_name = rec_name.get('fullName', {})
                if isinstance(full_name, dict):
                    result['name'] = full_name.get('value', '')
                elif isinstance(full_name, str):
                    result['name'] = full_name

            # Get function from comments
            comments = data.get('comments', [])
            for comment in comments:
                if comment.get('commentType') == 'FUNCTION':
                    texts = comment.get('texts', [])
                    if texts:
                        result['function'] = texts[0].get('value', '')
                    break

        # Rate limit - UniProt allows ~10 req/sec
        time.sleep(0.1)

    except Exception as e:
        log(f"Error fetching UniProt data for {accession}: {e}")

    _uniprot_cache[accession] = result
    return result


def get_gene_descriptions(acc_a: str, acc_b: str) -> Dict[str, Dict[str, str]]:
    """Get protein descriptions for both genes in a pair."""
    return {
        'gene_a': fetch_uniprot_description(acc_a) if acc_a else {'name': '', 'function': ''},
        'gene_b': fetch_uniprot_description(acc_b) if acc_b else {'name': '', 'function': ''},
    }
