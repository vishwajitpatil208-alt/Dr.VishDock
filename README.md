# MedVolt.Ai
Fully automated Protein Modeling, Structure Preparation, AutoDock Vina Docking, Visualization, and PDF Reporting pipeline.


# Advanced Automated Protein-Ligand Docking Pipeline

**Medvolt Tech Private Limited — Computational Chemistry Assignment**

A single-command, end-to-end pipeline that implements all five stages of the
assignment (Protein Modeling → Fixing → Ligand Prep → Docking → Visualization &
Report) with one robust Python script — now with powered by Rich.

---

## What's new in v2.0

| Upgrade | Details |
|---------|---------|
| **MedVolt.Ai cyberpunk banner** | Bold double-edged ASCII art splash on launch, gradient green→cyan. |
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

###  protein + ligand

```bash
# PDB ID + SMILES
python docking_pipeline.py \
    --protein 1HSG \
    --ligand "CC(=O)Oc1ccccc1C(=O)O" \
    --output results/

```



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
| `--protein <X>` | PDB ID 
| `--ligand <Y>`  | SMILES 
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
