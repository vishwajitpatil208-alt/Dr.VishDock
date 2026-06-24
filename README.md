# Dr. VishDOCK

**Advanced Automated Protein-Ligand Docking Pipeline**

A single-command, end-to-end pipeline that implements all five stages of
structure-based docking (Protein Modeling → Fixing → Ligand Prep → Docking →
Visualization & Report) with one robust Python script — featuring a
**pro-level hacker/gamer terminal UI** powered by Rich.

---

## 🚀 Quick Start

```bash
# 1. Create + activate the conda environment (one-time, ≈ 3.5 GB, 5–10 min)
conda env create -f environment.yml
conda activate docking-pipeline

# 2. Run your first docking job — PDB ID + SMILES, one line, no extra setup
python docking_pipeline.py --protein 1HSG --ligand "CC(=O)Oc1ccccc1C(=O)O" -o results/
```

That's it. `1HSG` (HIV-1 protease) is fetched from RCSB, fixed and minimized,
the SMILES (aspirin) is embedded in 3D and converted to PDBQT, AutoDock Vina
docks it, and PyMOL renders `best_pose.png` / `interactions.png` /
`overview.png` plus a one-page `report.pdf` — all inside `results/`.

```
results/
└── run_1HSG__ligand/
    ├── best_pose.png       ← hero shot with H-bonds (yellow) + hydrophobic contacts (pink)
    ├── interactions.png
    ├── overview.png
    ├── report.pdf
    └── results.json
```

See [Run Commands](#-run-commands) below for CSV-batch and folder-batch modes.

---

## What's new in v2.1.0

| Upgrade | Details |
|---------|---------|
| **Rebranded: Dr. VishDOCK** | New ASCII banner, all UI text and CLI help updated. |
| **Hardened ligand format support** | `.pdb`, `.sdf`, `.mol2`, `.mol`, `.smi`, and inline SMILES all go through a **two-strategy fallback chain**: RDKit's native reader first, then Open Babel bond-order perception if RDKit fails or returns zero bonds (the common failure mode for `.pdb` ligand snippets lacking CONECT records). |
| **Salt/fragment stripping** | Multi-fragment ligand inputs (e.g. a SMILES salt form `"CC(=O)O.[Na+]"`, or a PDB/MOL2 snippet that picked up a stray ion) automatically keep only the largest fragment - the real ligand. |
| **Fixed a path/SMILES ambiguity bug** | A mistyped or missing path like `ligand.pdb` previously fell through silently and was misclassified as the literal SMILES string `"ligand.pdb"`, producing a baffling "Invalid SMILES" error. Now raises a clear `FileNotFoundError` immediately. Verified safe against legitimate stereo-SMILES containing `/` or `\` (e.g. `C/C=C/C`). |
| **CSV mode hardened** | Empty/malformed CSVs, missing required columns, blank cells, and missing files now all produce clear `SystemExit` messages instead of raw tracebacks. Relative paths in a CSV resolve against the CSV's own folder if not found relative to the current working directory - so a CSV can be run from anywhere. A `smiles` column is accepted as an alias for `ligand`. |
| **Folder mode hardened** | Validates `--protein-dir`/`--ligand-dir` actually exist before scanning. Automatically filters out unsupported/junk files (`.DS_Store`, `README.txt`, etc.) with a warning, rather than crashing on them mid-run. |
| **PDBQT validation** | Every ligand-PDBQT-writing strategy (Meeko / Open Babel pybel / obabel CLI) is now followed by a sanity check that the output file actually contains atom records before being trusted - a "successful" but empty file no longer silently reaches Vina. |

---

## What's new in v2.0.2

| Patch | Details |
|-------|---------|
| **AlphaFold DB now uses the official API** | Calls `https://alphafold.ebi.ac.uk/api/prediction/{accession}` which returns the *current* `pdbUrl` regardless of model version. Static URL fallback walks `model_v6 → v5 → v4 → v3 → v2`. Fixes 404s introduced by the AF DB v6 release. |
| **Long-sequence handling** | If AF DB has no model for a UniProt ID *and* the sequence exceeds ESMFold's ~400-aa public limit, the job fails with a clear message ("supply a PDB ID or local .pdb instead") rather than crashing with a cryptic 413. |
| **PyMOL `one_letter` bug fixed** | Residue labels now use `pymol.stored.aa3to1`, the canonical namespace visible inside `cmd.label()` expressions. Includes all 20 AA + common non-standard (MSE, HID/HIE/HIP, CYX…). Wrapped in its own try/except so labels never kill the render. |
| **Clean console output** | Full tracebacks now go to `pipeline.log` and `ERROR.txt` only; the console gets a single one-line failure summary so the rich progress bar stays intact. |
| **Bundled `test_data/` folder** | New `test_data/setup_test_data.py` downloads three small PDBs and generates four ligands across all supported formats, so CSV + folder modes have a one-command reproducible test bed. |
| **Updated `example_jobs.csv`** | Replaced the EGFR/P00533 entry (which was breaking on AF DB v4 URLs) with safe PDB-ID-only entries. |

---

## ⚡ Quick Test of CSV + Folder Modes (recommended first run)

```bash
conda activate docking-pipeline

# Step 1 — populate test_data/ (one-time, ~30 sec)
cd test_data
python setup_test_data.py

# Step 2 — CSV mode
python ../docking_pipeline.py --csv example_csv_jobs.csv --parallel 3 -o ../csv_test/

# Step 3 — folder mode, 1-to-1 pairing
python ../docking_pipeline.py \
    --protein-dir proteins/ --ligand-dir ligands/ \
    --pair-mode zip --parallel 3 -o ../folder_zip_test/

# Step 4 — folder mode, cross-docking
python ../docking_pipeline.py \
    --protein-dir proteins/ --ligand-dir ligands/ \
    --pair-mode cross --parallel 4 -o ../folder_cross_test/
```

See [`test_data/README.md`](test_data/README.md) for details.

---

## What's new in v2.0 (original UI overhaul)

| Upgrade | Details |
|---------|---------|
| **Dr. VishDOCK cyberpunk banner** | Bold double-edged ASCII art splash on launch, gradient green→cyan. |
| **Live progress UI** | Animated spinner, color-coded stage markers (⟦01⟧ … ⟦05⟧), per-job briefing panels, neon-green completion ticks. |
| **Final results table** | Beautiful Rich table with ΔG, pose counts, H-bond counts, hydrophobic counts, status — ending in a `MISSION COMPLETE` panel. |
| **`best_pose.png`** | Hero shot of the top-scoring pose with element-coloured ligand sticks, ray-traced outlines, 300 dpi. |
| **H-bonds visualised** | Drawn as **yellow dashed lines with distance labels** between the exact donor/acceptor atom pairs. |
| **Hydrophobic contacts visualised** | Drawn as **pink/salmon dashed lines** between ligand-C and pocket-C atoms. |
| **Labeled pocket residues** | Each interacting residue is labeled (one-letter code + residue number) at its Cα. |
| **`interactions.png`** | Tight zoomed close-up of the interaction network only. |
| **`overview.png`** | Whole-protein context view. |
| **`--no-banner`** flag | For headless / CI environments where ANSI escapes would clutter logs. |

---

## Features

| Stage | Capability |
|-------|-----------|
| **Stage 1** — Protein Modeling | Accepts **PDB ID** (RCSB), **UniProt ID** (AlphaFold DB / ESMFold fallback), **FASTA file or raw sequence**, or local `.pdb` file. Auto-detects the input type. |
| **Stage 2** — Protein Fixing | PDBFixer cleanup (waters removed, missing residues/atoms rebuilt, hydrogens at pH 7.4) + OpenMM Amber14/GBN2 energy minimization + Gasteiger-charged PDBQT receptor (Meeko → Open Babel fallback). |
| **Stage 3** — Ligand Prep | Accepts **SMILES string**, `.smi`, `.sdf`, `.mol2`, `.mol`, or `.pdb`. 3D embed (ETKDGv3) + MMFF94 optimization + Meeko PDBQT with proper torsions. |
| **Stage 4** — Docking | **AutoDock Vina** via Python API (CLI fallback). Auto binding-site detection: user grid → co-crystal HETATM center → blind protein COM. Multi-pose ranking. |
| **Stage 5** — Visualization | **PyMOL** (open-source) publication-quality PNG of best pose + interaction residues highlighted, plus distance-based H-bond / hydrophobic interaction analysis, plus a single-page **PDF report** per run via ReportLab. |
| **Batch / Parallel** | CSV-driven batch (`--csv`), folder × folder cross-docking (`--protein-dir`/`--ligand-dir`), N parallel workers (`--parallel N`), aggregated `summary.csv` + `summary.json`. |

---

## Installation

```bash
# 1. Create the conda environment (≈ 3.5 GB, takes 5–10 min)
conda env create -f environment.yml

# 2. Activate it
conda activate docking-pipeline

# 3. Verify
python docking_pipeline.py --version
python docking_pipeline.py --help
```

> On Apple Silicon, replace `conda` with `mamba` for a faster solve, or use
> `CONDA_SUBDIR=osx-64 conda env create -f environment.yml` if any package is
> not yet built for arm64.

---

## ▶ Run Commands

### 1. Single protein + single ligand

```bash
# PDB ID + SMILES (short -o flag)
python docking_pipeline.py --protein 1HSG --ligand "CC(=O)Oc1ccccc1C(=O)O" -o results/

# Same thing, long-form flags
python docking_pipeline.py \
    --protein 1HSG \
    --ligand "CC(=O)Oc1ccccc1C(=O)O" \
    --output results/

# UniProt ID (auto-fetches AlphaFold) + SDF file
python docking_pipeline.py \
    --protein P00533 \
    --ligand erlotinib.sdf \
    --output results/

# Local PDB file + MOL2 file, user-defined grid box
python docking_pipeline.py \
    --protein receptor.pdb \
    --ligand ligand.mol2 \
    --center 12.3 -4.5 27.8 \
    --size 22 22 22 \
    --output results/

# FASTA sequence (ESMFold) + SMILES, skip minimization for speed
python docking_pipeline.py \
    --protein "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLS..." \
    --ligand "c1ccccc1O" \
    --no-minimize \
    --output results/
```

### 2. CSV-driven batch (recommended for many jobs)

Create a CSV with at least `protein,ligand` columns (optionally `name,cx,cy,cz,sx,sy,sz`):

```csv
protein,ligand,name
1HSG,CC(=O)Oc1ccccc1C(=O)O,1HSG_aspirin
P00533,erlotinib.sdf,EGFR_erlotinib
6LU7,nirmatrelvir.mol2,SARS_Mpro
```

```bash
python docking_pipeline.py \
    --csv example_jobs.csv \
    --parallel 4 \
    --output csv_results/
```

### 3. Folder × Folder batch (cross-docking)

```bash
# All-vs-all (default)
python docking_pipeline.py \
    --protein-dir proteins/ \
    --ligand-dir ligands/ \
    --pair-mode cross \
    --parallel 8 \
    --output cross_results/

# 1-to-1 zipped pairing (sorted file order)
python docking_pipeline.py \
    --protein-dir proteins/ \
    --ligand-dir ligands/ \
    --pair-mode zip \
    --parallel 4 \
    --output paired_results/
```

---

## Output Layout

```
results/
├── pipeline.log              # top-level log
├── summary.csv               # aggregated results across ALL jobs
├── summary.json
└── run_<protein>__<ligand>/
    ├── pipeline.log              # per-job log
    ├── 01_protein_raw.pdb        # original download
    ├── 02_protein_fixed.pdb      # PDBFixer output
    ├── 03_protein_minimized.pdb  # OpenMM-minimized
    ├── 04_receptor.pdbqt         # ready for Vina
    ├── 05_ligand_3d.sdf          # 3D-embedded ligand
    ├── 06_ligand.pdbqt
    ├── 07_docked.pdbqt           # all Vina poses
    ├── 08_docked.sdf             # poses in SDF
    ├── 09_complex.pdb            # protein + best pose
    ├── best_pose.png             # ★ hero shot (H-bonds, hydrophobic, labels)
    ├── interactions.png          # close-up of binding interactions
    ├── overview.png              # whole-protein context
    ├── results.json              # all metrics
    └── report.pdf                # one-page summary (all 3 images embedded)
```

---

## CLI Reference (key flags)

| Flag | Description |
|------|-------------|
| `--protein <X>` | PDB ID / UniProt ID / FASTA / file |
| `--ligand <Y>`  | SMILES / .sdf / .mol2 / .pdb / .mol / .smi |
| `--csv FILE`    | CSV batch mode |
| `--protein-dir DIR --ligand-dir DIR` | Folder batch mode |
| `--pair-mode {cross,zip}` | Folder pairing strategy |
| `--center X Y Z`  | Manual grid center (Å) |
| `--size SX SY SZ` | Manual grid size (Å) |
| `--exhaustiveness N` | Vina exhaustiveness (default 16) |
| `--num-poses N`   | Poses per ligand (default 9) |
| `--no-minimize`   | Skip OpenMM minimization (faster) |
| `--parallel N`    | Run N jobs in parallel |
| `--cpu-per-job N` | Threads per Vina job (0 = all) |
| `-o, --output DIR`| Output directory |
| `-v, --verbose`   | DEBUG-level logging |

---

## Quick Smoke Test

```bash
conda activate docking-pipeline
python docking_pipeline.py \
    --protein 1HSG \
    --ligand "CC(=O)Oc1ccccc1C(=O)O" \
    --output smoke_test/ -v
```

When it finishes you should have `smoke_test/run_1HSG__1HSG/report.pdf` and a
non-empty `smoke_test/summary.csv`.

---

## Robustness Notes

* Every stage has **multiple fallback paths** (e.g. Meeko → Open Babel CLI; Vina
  Python API → Vina binary; AlphaFold DB → ESMFold).
* Per-job failures **do not crash the batch** — they are logged and reported
  in `summary.csv` with the error string.
* Heavy modules (RDKit, OpenMM, Meeko, Vina, PyMOL) are imported **lazily**
  inside their classes, so `--help` works on a bare interpreter and the script
  remains importable for testing.
* All outputs are deterministic where possible (RDKit ETKDG seeded; Vina
  exhaustiveness configurable).
