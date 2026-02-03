# Human Paralog Database

Interactive explorer for human paralog gene pairs.

**Live:** https://ocmalde.github.io/human-paralog-db/

## Data Structure

```
data/
├── index.json              # Pair index with metadata
├── family_index.json       # Gene -> pair_id mapping
├── search_index.json       # Search terms
└── pairs/
    └── GENE1_GENE2/
        ├── report.json     # Full report data
        ├── summary.json    # Summary statistics
        └── pdb.json        # PDB structure variants
```