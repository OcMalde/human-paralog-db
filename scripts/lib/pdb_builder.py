"""PDB file building with custom B-factors for coloring.

This module generates PDB variants with B-factors set to different values
for visualization in Molstar (AM scores, pLDDT, delta AM, etc.)
"""

import base64
from typing import List, Optional, Dict, Any


def parse_pdb_atoms(pdb_text: str) -> List[Dict[str, Any]]:
    """Parse ATOM/HETATM lines from PDB text into structured records."""
    atoms = []
    for line in pdb_text.splitlines():
        if line.startswith(('ATOM', 'HETATM')):
            try:
                atom = {
                    'type': line[0:6].strip(),
                    'serial': int(line[6:11].strip()),
                    'name': line[12:16].strip(),
                    'altLoc': line[16:17],
                    'resName': line[17:20].strip(),
                    'chainID': line[21:22],
                    'resSeq': int(line[22:26].strip()),
                    'iCode': line[26:27],
                    'x': float(line[30:38].strip()),
                    'y': float(line[38:46].strip()),
                    'z': float(line[46:54].strip()),
                    'occupancy': float(line[54:60].strip()) if len(line) > 54 and line[54:60].strip() else 1.0,
                    'tempFactor': float(line[60:66].strip()) if len(line) > 60 and line[60:66].strip() else 0.0,
                    'element': line[76:78].strip() if len(line) > 76 else '',
                    'charge': line[78:80].strip() if len(line) > 78 else '',
                }
                atoms.append(atom)
            except (ValueError, IndexError):
                pass
    return atoms


def format_pdb_atom(atom: Dict[str, Any]) -> str:
    """Format an atom record back to PDB format."""
    return (
        f"{atom['type']:<6s}"
        f"{atom['serial']:>5d} "
        f"{atom['name']:>4s}"
        f"{atom['altLoc']:1s}"
        f"{atom['resName']:>3s} "
        f"{atom['chainID']:1s}"
        f"{atom['resSeq']:>4d}"
        f"{atom['iCode']:1s}   "
        f"{atom['x']:>8.3f}"
        f"{atom['y']:>8.3f}"
        f"{atom['z']:>8.3f}"
        f"{atom['occupancy']:>6.2f}"
        f"{atom['tempFactor']:>6.2f}"
        f"          "
        f"{atom['element']:>2s}"
        f"{atom['charge']:>2s}"
    )


def set_bfactors_by_residue(
    pdb_text: str,
    bfactor_map: Dict[tuple, float],  # (chainID, resSeq) -> bfactor
    default_bfactor: float = 50.0
) -> str:
    """Set B-factors for each residue based on a mapping.

    Args:
        pdb_text: Original PDB text
        bfactor_map: Dict mapping (chainID, resSeq) to B-factor value
        default_bfactor: Default B-factor for unmapped residues

    Returns:
        Modified PDB text with updated B-factors
    """
    lines = []
    for line in pdb_text.splitlines():
        if line.startswith(('ATOM', 'HETATM')):
            try:
                chain = line[21:22]
                resSeq = int(line[22:26].strip())
                key = (chain, resSeq)
                bfactor = bfactor_map.get(key, default_bfactor)

                # Update B-factor in line
                line = line[:60] + f"{bfactor:>6.2f}" + line[66:]
            except (ValueError, IndexError):
                pass
        lines.append(line)

    return '\n'.join(lines)


def filter_chains(pdb_text: str, chains_to_keep: List[str]) -> str:
    """Keep only specified chains in PDB.

    Args:
        pdb_text: Original PDB text
        chains_to_keep: List of chain IDs to keep (e.g., ['A'] or ['B'])

    Returns:
        Filtered PDB text
    """
    lines = []
    chains_set = set(chains_to_keep)

    for line in pdb_text.splitlines():
        if line.startswith(('ATOM', 'HETATM')):
            chain = line[21:22]
            if chain in chains_set:
                lines.append(line)
        elif line.startswith(('MODEL', 'ENDMDL', 'END', 'TER')):
            lines.append(line)
        elif not line.startswith(('ATOM', 'HETATM', 'CONECT', 'MASTER')):
            # Keep header lines, SEQRES, etc.
            lines.append(line)

    return '\n'.join(lines)


def build_am_bfactor_map(
    aligned_cols: List[tuple],  # [(col, qpos, tpos), ...]
    bfactors_a: List[float],
    bfactors_b: List[float],
    chain_a_id: str = 'A',
    chain_b_id: str = 'B'
) -> Dict[tuple, float]:
    """Build B-factor map from AM scores for aligned residues.

    Args:
        aligned_cols: List of (col, qpos, tpos) tuples
        bfactors_a: AM scores for protein A (0-1 scale, or raw)
        bfactors_b: AM scores for protein B
        chain_a_id: Chain ID for protein A
        chain_b_id: Chain ID for protein B

    Returns:
        Dict mapping (chainID, resSeq) to B-factor (scaled 0-100 for visualization)
    """
    bfactor_map = {}

    for col, qpos, tpos in aligned_cols:
        if qpos and 1 <= qpos <= len(bfactors_a):
            val = bfactors_a[qpos - 1]
            if isinstance(val, (int, float)):
                # Scale to 0-100 range for Molstar uncertainty coloring
                scaled = val * 100 if val <= 1 else val
                bfactor_map[(chain_a_id, qpos)] = scaled

        if tpos and 1 <= tpos <= len(bfactors_b):
            val = bfactors_b[tpos - 1]
            if isinstance(val, (int, float)):
                scaled = val * 100 if val <= 1 else val
                bfactor_map[(chain_b_id, tpos)] = scaled

    return bfactor_map


def build_plddt_bfactor_map(
    aligned_cols: List[tuple],
    plddt_a: List[float],
    plddt_b: List[float],
    chain_a_id: str = 'A',
    chain_b_id: str = 'B'
) -> Dict[tuple, float]:
    """Build B-factor map from pLDDT scores (already 0-100 scale)."""
    bfactor_map = {}

    for col, qpos, tpos in aligned_cols:
        if qpos and 1 <= qpos <= len(plddt_a):
            val = plddt_a[qpos - 1]
            if isinstance(val, (int, float)):
                bfactor_map[(chain_a_id, qpos)] = val

        if tpos and 1 <= tpos <= len(plddt_b):
            val = plddt_b[tpos - 1]
            if isinstance(val, (int, float)):
                bfactor_map[(chain_b_id, tpos)] = val

    return bfactor_map


def pdb_to_base64(pdb_text: str) -> str:
    """Convert PDB text to base64 string."""
    return base64.b64encode(pdb_text.encode('utf-8')).decode('ascii')
