# Human Paralog Database

An interactive explorer for human paralog gene pairs with structural alignments, AlphaMissense pathogenicity data, and functional annotations.

**Live Demo:** https://ocmalde.github.io/human-paralog-db/

## Features

- **Structural alignment visualization** using Foldseek and Mol* viewer
- **AlphaMissense pathogenicity** scores with multiple normalization modes
- **Protein-protein interaction** networks (shared vs unique partners)
- **Domain annotations** from Pfam, CATH, and TED databases
- **Family navigation** to explore related paralog pairs
- **Summary statistics** with radar charts and boxplot comparisons

## Quick Start

### View the Explorer

Simply visit the GitHub Pages site or open `index.html` in a browser after cloning.

### Add New Pairs

1. Add gene pairs to `input/pairs.csv`:
   ```csv
   A1,A2
   GENE1,GENE2
   GENE3,GENE4
   ```

2. Run the pipeline:
   ```bash
   # Populate database with structures and alignments
   python scripts/populate_db.py

   # Generate report data (features, domains, PPIs)
   python scripts/generate_report_data.py

   # Export to static JSON files
   python scripts/export_static.py
   ```

3. Commit and push to update the live site.

## Project Structure

```
human-paralog-db/
├── index.html              # Main entry point (pair selector)
├── report.html             # Detailed pair report viewer
├── static/
│   ├── css/style.css       # Styles
│   └── js/app.js           # Report visualization code
├── data/
│   ├── index.json          # Pair index with metadata
│   ├── family_index.json   # Gene -> pair_id mapping
│   ├── search_index.json   # Search terms index
│   └── pairs/
│       └── GENE1_GENE2/
│           ├── report.json  # Full report data
│           ├── summary.json # Summary statistics
│           └── pdb.json     # PDB structure variants
├── input/
│   ├── pairs.csv                      # Input gene pairs
│   ├── all_genes_ids.csv              # Gene -> UniProt mapping
│   ├── ens111_human_allFeatures.csv   # Ensembl features
│   ├── mart_ens105_gene_loc.txt       # Gene locations
│   └── CRISPRInferredCommonEssentials.csv  # Essential genes
└── scripts/
    ├── populate_db.py         # Build SQLite DB with structures
    ├── generate_report_data.py # Generate report JSON
    ├── export_static.py       # Export to static files
    └── lib/
        ├── config.py          # Configuration
        ├── fetch_features.py  # Feature fetching
        ├── structure_utils.py # PDB utilities
        └── ...
```

## Data Pipeline

### 1. populate_db.py

Builds the SQLite database with:
- Gene symbol to UniProt accession mapping
- AlphaFold structure downloads
- Foldseek structural alignments
- Domain annotations (Pfam, CATH, TED)

### 2. generate_report_data.py

Generates report data including:
- AlphaMissense pathogenicity tracks
- Protein-protein interactions (BioGRID)
- Conservation scores
- Summary statistics

### 3. export_static.py

Exports everything to static JSON files for GitHub Pages hosting.

## Configuration

### Data Source (app.js)

```javascript
// For GitHub Pages (default)
const DATA_BASE = './data';

// For Bunny CDN (large datasets)
const DATA_BASE = 'https://your-zone.b-cdn.net';
```

## Input Data Requirements

The following files are needed in `input/`:

| File | Description | Source |
|------|-------------|--------|
| `pairs.csv` | Gene pairs to analyze | User-defined |
| `all_genes_ids.csv` | Gene -> UniProt mapping | Ensembl BioMart |
| `ens111_human_allFeatures.csv` | Domain/feature annotations | Ensembl |
| `mart_ens105_gene_loc.txt` | Chromosome locations | Ensembl BioMart |
| `CRISPRInferredCommonEssentials.csv` | Essential genes | DepMap |

## Dependencies

### Python

```bash
pip install pandas requests biopython
```

### External Tools

- **Foldseek** - Structural alignment (must be in PATH)
  ```bash
  # Install via conda
  conda install -c conda-forge -c bioconda foldseek
  ```

## Scaling to 100K+ Pairs

For large datasets:

1. **Use Bunny CDN** for data hosting (~$10/month for 1TB)
2. **Update `DATA_BASE`** in app.js to your CDN URL
3. **Batch processing**: Run scripts with specific pair IDs
   ```bash
   python scripts/export_static.py PAIR1 PAIR2 PAIR3
   ```

## API Endpoints (for reference)

The static site replaces these former API endpoints:

| Original API | Static File |
|--------------|-------------|
| `/api/pairs` | `data/index.json` |
| `/api/report/{pair}` | `data/pairs/{pair}/report.json` |
| `/api/summary/{pair}` | `data/pairs/{pair}/summary.json` |
| `/api/pdb/{pair}` | `data/pairs/{pair}/pdb.json` |
| `/api/family/{gene}` | `data/family_index.json` |

## License

MIT License

## Credits

- **AlphaFold** structures from DeepMind/EBI
- **AlphaMissense** pathogenicity from Google DeepMind
- **Foldseek** alignment from Steinegger Lab
- **Mol\*** viewer from PDBe
- **BioGRID** for protein interactions