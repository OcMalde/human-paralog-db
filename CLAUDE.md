# CLAUDE.md — Instructions for Claude Code in this repo

## Adding new pairs to the ONLINE database

> **Online DB** = the GitHub Pages site (data/pairs/, data/index.json, static/data/).
> **Offline reports** = generate_local_report.py → standalone HTML. These are separate.

### Full pipeline (MUST follow in order)

```
1. Edit input/pairs.csv          — add new pairs (A1,A2 format)

2. python scripts/populate_db.py [GENE1_GENE2 ...]
   — fetches AF structures, runs foldseek, populates pairs.db

3. python scripts/run_plma.py [GENE1_GENE2 ...]
   — REQUIRED for family alignment section in the report
   — Needs Docker + container 4fa044ff5162 (fuse-phylotree image, has paloma-D)
   — If Docker is not running: ASK THE USER to start Docker before proceeding
   — Do NOT skip this step; missing plma.json → "PLMA not available" in the UI

4. /Users/olivierdennler/miniconda/bin/python3 scripts/generate_report_data.py [GENE1_GENE2]
   — MUST use conda Python (not system Python) — system Python lacks pyarrow → 0 drugs
   — MUST run from repo root: cd /Users/olivierdennler/Documents/data/human-paralog-db
   — Runs for all new pairs if no pair_id given (skips existing unless --force)
   — Known drug count = 0 for BOTH genes is suspicious → check:
       a) used conda Python?  b) parquet files in input/opentargets/known_drug/?
       c) gene in input/all_genes_ids.csv with Ensembl ID?

5. /Users/olivierdennler/miniconda/bin/python3 scripts/export_static.py
   — Exports per-gene files: data/genes/{GENE}.json (PDB, AM, domains, drugs, ...)
   — Exports per-pair files: data/pairs/GENE1_GENE2/{report,summary,pdb,plma}.json
     (pair files no longer contain per-gene data — browser loads both in parallel)
   — Updates data/index.json, data/family_index.json, data/search_index.json

6. git add data/pairs/ data/genes/ data/index.json data/family_index.json data/search_index.json input/pairs.csv
   git commit && git push
```

### Key facts to remember

- **Docker container**: `4fa044ff5162` (fuse-phylotree image) — has paloma-D at `/paloma-D`
- **conda Python**: `/Users/olivierdennler/miniconda/bin/python3` — required for pyarrow/drug data
- **Repo root**: `/Users/olivierdennler/Documents/data/human-paralog-db` — all scripts must run from here
- **pairs.db**: SQLite at `data/pairs.db` — source of truth; static JSONs are derived from it
- **data/genes/**: Per-gene JSON files — shared across all pairs, avoid duplication at scale
- **data/pairs/**: Per-pair JSON files — alignment, rects, domPairs, conservation, SL; NO per-gene data
- **PLMA agraph caching**: `data/cache/plma/family_NNN/` — if paloma-D segfaults, check if agraph is already in the Docker container (`docker exec CONTAINER ls /family_*`) and copy it out manually with `docker cp`

### Common mistakes to avoid

| Mistake | Symptom | Fix |
|---|---|---|
| System Python instead of conda | All known_drugs = 0 | Use `/miniconda/bin/python3` |
| Wrong working directory | "No such file" errors | `cd human-paralog-db` first |
| Skip run_plma.py | "PLMA data not available" in UI | Run step 3, check Docker is running |
| export_static.py overwrites file-level fixes | Drug data reverts to 0 | Always fix in DB via generate_report_data.py, not files |
| Adding pairs without re-exporting gene files | Missing gene data in browser | export_static.py always re-exports gene files; no --skip-genes needed for normal runs |
| paloma-D segfault | "ERROR: paloma-D did not produce .agraph" | Check if agraph is in container; copy out with docker cp |
