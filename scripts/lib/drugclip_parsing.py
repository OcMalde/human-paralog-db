"""
drugclip_parsing.py

Load DrugCLIP/GenomeScreen pocket data from extracted JSON files.
Returns domain-compatible dicts for integration with the report pipeline.
"""

import json
from pathlib import Path
from typing import List, Dict, Any

from .config import DRUGCLIP_DIR, log


def load_drugclip_pockets(acc: str) -> List[Dict[str, Any]]:
    """Load DrugCLIP pockets for a UniProt accession.

    Returns list of domain-compatible dicts with DrugCLIP-specific fields.
    """
    pockets_file = DRUGCLIP_DIR / acc / "pockets.json"
    if not pockets_file.exists():
        return []

    with open(pockets_file) as f:
        pockets_data = json.load(f)

    domains = []
    for p in pockets_data:
        if p.get("start") is None or p.get("end") is None:
            continue

        pocket_num = p.get("pocket_num", "?")
        top_score = p.get("top_score")
        score_str = f" score={top_score:.1f}" if top_score else ""

        domain = {
            "type": "DrugCLIP",
            "label": f"Pocket {pocket_num} (DrugCLIP)",
            "start": p["start"],
            "end": p["end"],
            "uid": f"DC:{p['start']}-{p['end']}:{p.get('name', '')}",
            "raw_type": "DrugCLIP",
            # DrugCLIP-specific fields
            "grid_center": p.get("grid_center"),
            "residues": p.get("residues", []),
            "top_drugs": p.get("top_drugs", []),
            "n_total_hits": p.get("n_total_hits", 0),
            "top_score": top_score,
            "pocket_num": pocket_num,
            "pocket_name": p.get("name", ""),
        }
        domains.append(domain)

    return domains
