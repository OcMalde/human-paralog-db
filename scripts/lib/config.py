"""Configuration and constants for report generation."""

from pathlib import Path

# Paths - relative to working directory (repo root)
BASE_DIR = Path.cwd()
DB_PATH = BASE_DIR / "data" / "pairs.db"
STRUCTURES_DIR = BASE_DIR / "data" / "structures"
CACHE_DIR = BASE_DIR / "data" / "cache"

# Input files (reference data)
INPUT_DIR = BASE_DIR / "input"
FEATURES_CSV = INPUT_DIR / "ens111_human_allFeatures.csv"
GENE_LOC_CSV = INPUT_DIR / "mart_ens105_gene_loc.txt"
ESSENTIAL_CSV = INPUT_DIR / "CRISPRInferredCommonEssentials.csv"
GENE_ID_MAP_CSV = INPUT_DIR / "all_genes_ids.csv"

# Cavity data paths (optional - HegeLab AF-all_cavities dump)
CAVITY_INDEX_PATH = INPUT_DIR / "AF-all_cavities" / "AF-all_cavities.idx"
CAVITY_DATA_DIR = INPUT_DIR / "AF-all_cavities" / "all_cavities_data"

# API endpoints
PDBE_BASE = "https://www.ebi.ac.uk/pdbe/api"
AM_HOTSPOT_API = "https://alphamissense.hegelab.org/hotspotapi"

# AM normalization modes
AM_MODES = ["raw", "percentile", "minmax", "zscore"]
AA_ORDER = ['K','R','H','E','D','N','Q','T','S','C','G','A','V','L','I','M','P','Y','F','W']


def log(msg: str):
    print(f"[generate] {msg}")