"""AlphaMissense processing - normalization, coloring, and rect generation.

This module provides:
1. Normalization functions (raw, percentile, minmax, zscore)
2. Color banding for AM scores (grey→red for pathogenicity)
3. Rectangle generation for SVG tracks
4. Full AM matrix computation for all 20 amino acids

FIXES:
- Empty data (None) shows as white (#ffffff) to distinguish from actual data
- Score of 0 or low benign shows as light grey (#dddddd)
- V row included in AA_ORDER
- Bin summary properly computed as average of available data
"""

import math
import bisect
from typing import Any, Dict, List, Optional, Tuple

# All 20 standard amino acids - V is at position 13
AA_ORDER = ['K','R','H','E','D','N','Q','T','S','C','G','A','V','L','I','M','P','Y','F','W']
AM_MODES = ["raw", "percentile", "minmax", "zscore"]

# Color for missing/no-data regions (white to distinguish from grey benign scores)
NO_DATA_COLOR = "#ffffff"


def make_normalizer(values: List[Optional[float]], mode: str):
    """Create a normalization function based on mode."""
    xs = [float(v) for v in values if isinstance(v, (int, float))]

    if not xs or mode == "raw":
        return lambda v: v if isinstance(v, (int, float)) else None

    if mode == "minmax":
        lo, hi = min(xs), max(xs)
        if hi <= lo:
            return lambda v: 0.5 if isinstance(v, (int, float)) else None
        return lambda v: (v - lo) / (hi - lo) if isinstance(v, (int, float)) else None

    if mode == "percentile":
        xs_sorted = sorted(xs)
        n = len(xs_sorted)
        def f(v):
            if not isinstance(v, (int, float)):
                return None
            idx = bisect.bisect_right(xs_sorted, v)
            return idx / n if n > 1 else 0.5
        return f

    if mode == "zscore":
        mu = sum(xs) / len(xs)
        var = sum((x - mu) ** 2 for x in xs) / len(xs)
        sigma = math.sqrt(var) if var > 0 else 1e-8
        def f(v):
            if not isinstance(v, (int, float)):
                return None
            z = (v - mu) / sigma
            return 1.0 / (1.0 + math.exp(-z))  # Sigmoid to 0-1
        return f

    return lambda v: v if isinstance(v, (int, float)) else None


def band_color_am(v: Optional[float]) -> Optional[str]:
    """
    Color for AM score: grey (benign) → red (pathogenic).

    Thresholds:
    - None → None (will be handled separately as no-data)
    - < 0.2 → #dddddd (grey/benign)
    - < 0.4 → #bbbbbb (light grey)
    - < 0.7 → #ff7d45 (orange)
    - >= 0.7 → #d62728 (red/pathogenic)
    """
    if v is None:
        return None
    if v < 0.2:
        return "#dddddd"  # Benign - grey
    if v < 0.4:
        return "#bbbbbb"  # Low
    if v < 0.7:
        return "#ff7d45"  # Medium - orange
    return "#d62728"      # High - red


def dam_color(v: Optional[float]) -> Optional[str]:
    """
    Color for delta-AM (difference between paralogs).
    Low difference = grey, high difference = red.
    """
    if v is None:
        return None
    if v <= 0.10:
        return "#dddddd"
    if v <= 0.30:
        return "#bbbbbb"
    if v <= 0.50:
        return "#ffdb13"  # Yellow
    if v <= 0.70:
        return "#ff7d45"  # Orange
    return "#d62728"      # Red


def rects_from_values(
    vals: List[Optional[float]],
    aln_len: int,
    color_fn=None,
    show_no_data: bool = True
) -> List[Dict[str, Any]]:
    """
    Convert value array to colored rectangles for tracks.

    Args:
        vals: List of values (None means no data)
        aln_len: Alignment length
        color_fn: Function to convert value to color
        show_no_data: If True, show dark grey for None values
    """
    if color_fn is None:
        color_fn = band_color_am

    out = []
    cur_color = None
    seg_start = None

    for i, v in enumerate(vals, start=1):
        c = color_fn(v)

        # Handle no-data case
        if c is None:
            if show_no_data:
                c = NO_DATA_COLOR
            else:
                # Close current segment if any
                if cur_color is not None and seg_start is not None:
                    out.append({"start": seg_start, "end": i - 1, "color": cur_color})
                    cur_color = None
                    seg_start = None
                continue

        if cur_color is None:
            cur_color = c
            seg_start = i
        elif c != cur_color:
            out.append({"start": seg_start, "end": i - 1, "color": cur_color})
            cur_color = c
            seg_start = i

    if cur_color is not None and seg_start is not None:
        out.append({"start": seg_start, "end": aln_len, "color": cur_color})

    return out


def compute_am_rects_all_modes(
    bfA: List[float],
    bfB: List[float],
    aligned_cols: List,  # [[col, qpos, tpos], ...]
    aln_len: int,
    mat_a_scores_raw: Dict[str, List[Optional[float]]],
    mat_b_scores_raw: Dict[str, List[Optional[float]]],
) -> Dict[str, Any]:
    """
    Compute all AM-related rectangles for all normalization modes.

    This function matches the call signature in generate_report_data.py.

    Returns dict with all the rect data needed for visualization.

    IMPORTANT COLOR DISTINCTION:
    - No data available → white (#ffffff)
    - Data exists, benign (0-0.2) → light grey (#dddddd)
    """
    result = {
        "amModes": AM_MODES,
        "amAlignedARectsByMode": {},
        "amAlignedBRectsByMode": {},
        "damAlignedRectsByMode": {},
        "amMatrixA_rectsByMode": {},
        "amMatrixB_rectsByMode": {},
    }

    # Build qpos/tpos maps from aligned_cols
    qpos_by_col = {}
    tpos_by_col = {}
    for item in aligned_cols:
        col, qpos, tpos = item[0], item[1], item[2]
        qpos_by_col[col] = qpos
        tpos_by_col[col] = tpos

    for mode in AM_MODES:
        normA = make_normalizer(bfA, mode)
        normB = make_normalizer(bfB, mode)

        # AM on alignment axis (summary tracks from bfactors)
        # These show the per-residue AM score
        amA_vals = [None] * aln_len
        amB_vals = [None] * aln_len

        for col in range(1, aln_len + 1):
            qpos = qpos_by_col.get(col)
            tpos = tpos_by_col.get(col)

            if qpos is not None and 1 <= qpos <= len(bfA):
                amA_vals[col - 1] = normA(bfA[qpos - 1])
            if tpos is not None and 1 <= tpos <= len(bfB):
                amB_vals[col - 1] = normB(bfB[tpos - 1])

        # Show no-data as dark grey for summary tracks
        result["amAlignedARectsByMode"][mode] = rects_from_values(amA_vals, aln_len, band_color_am, show_no_data=True)
        result["amAlignedBRectsByMode"][mode] = rects_from_values(amB_vals, aln_len, band_color_am, show_no_data=True)

        # ΔAM per alignment column (only where both are aligned)
        dam_vals = [None] * aln_len
        for col, qpos, tpos in aligned_cols:
            if qpos and tpos and 1 <= qpos <= len(bfA) and 1 <= tpos <= len(bfB):
                va = normA(bfA[qpos - 1])
                vb = normB(bfB[tpos - 1])
                if isinstance(va, (int, float)) and isinstance(vb, (int, float)):
                    dam_vals[col - 1] = abs(va - vb)

        result["damAlignedRectsByMode"][mode] = rects_from_values(dam_vals, aln_len, dam_color, show_no_data=True)

        # Substitution matrix rects - all 20 AAs including V
        matA_rects_mode = {}
        matB_rects_mode = {}

        # Collect all matrix scores for normalization
        all_mat_a_values = []
        all_mat_b_values = []
        for aa in AA_ORDER:
            all_mat_a_values.extend([v for v in (mat_a_scores_raw.get(aa) or []) if v is not None])
            all_mat_b_values.extend([v for v in (mat_b_scores_raw.get(aa) or []) if v is not None])

        normMatA = make_normalizer(all_mat_a_values, mode)
        normMatB = make_normalizer(all_mat_b_values, mode)

        for aa in AA_ORDER:
            arrA = mat_a_scores_raw.get(aa) or []
            arrB = mat_b_scores_raw.get(aa) or []

            # Ensure arrays are correct length
            while len(arrA) < aln_len:
                arrA.append(None)
            while len(arrB) < aln_len:
                arrB.append(None)

            valsA = [normMatA(v) if isinstance(v, (int, float)) else None for v in arrA[:aln_len]]
            valsB = [normMatB(v) if isinstance(v, (int, float)) else None for v in arrB[:aln_len]]

            # For matrix rows, do NOT show no-data (leave as gaps/transparent)
            matA_rects_mode[aa] = rects_from_values(valsA, aln_len, band_color_am, show_no_data=False)
            matB_rects_mode[aa] = rects_from_values(valsB, aln_len, band_color_am, show_no_data=False)

        result["amMatrixA_rectsByMode"][mode] = matA_rects_mode
        result["amMatrixB_rectsByMode"][mode] = matB_rects_mode

    return result
