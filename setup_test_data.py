#!/usr/bin/env python3
"""
============================================================================
  Dr. VishDOCK — Test Data Setup
  --------------------------------------------------------------------------
  Run once to populate the test_data/ folder with:
    * proteins/   ← three small PDB structures from RCSB (bound-ligand)
    * ligands/    ← four ligands in different formats (.sdf, .mol2, .smi, .pdb)
    * example_csv_jobs.csv   ← safe CSV using only PDB IDs (no UniProt issues)

  After this, you can immediately run:
    python ../docking_pipeline.py --csv example_csv_jobs.csv -o ../csv_test/
    python ../docking_pipeline.py --protein-dir proteins/ \
                                  --ligand-dir  ligands/  \
                                  --pair-mode cross --parallel 4 \
                                  -o ../folder_test/
============================================================================
"""
from __future__ import annotations
import sys
import urllib.request
from pathlib import Path

HERE  = Path(__file__).parent.resolve()
PROT  = HERE / "proteins"
LIG   = HERE / "ligands"
PROT.mkdir(exist_ok=True)
LIG.mkdir(exist_ok=True)

GREEN = "\033[1;92m"; CYAN = "\033[1;96m"; RED = "\033[1;91m"
DIM   = "\033[2m";    OFF  = "\033[0m"


def say(prefix: str, msg: str, color: str = GREEN) -> None:
    print(f"  {color}{prefix}{OFF}  {msg}")


print()
print(f"{CYAN}╔══════════════════════════════════════════════════════════════╗{OFF}")
print(f"{CYAN}║       Dr. VishDOCK  ─  Test Data Setup                       ║{OFF}")
print(f"{CYAN}╚══════════════════════════════════════════════════════════════╝{OFF}")
print()


# ─── Step 1: download proteins ─────────────────────────────────────────────
PDBS = [
    ("1HSG", "HIV-1 protease + indinavir         (small,  ~99 res)"),
    ("3HTB", "β-glucosidase + nojirimycin        (medium, ~470 res)"),
    ("4LDE", "β2-adrenergic GPCR + ligand        (membrane receptor)"),
]
print(f"{CYAN}[1/3]{OFF}  Downloading test proteins from RCSB ...")
for pid, desc in PDBS:
    dest = PROT / f"{pid.lower()}.pdb"
    if dest.exists() and dest.stat().st_size > 1000:
        say("✓", f"{pid:6}  already present  {DIM}({desc}){OFF}")
        continue
    url = f"https://files.rcsb.org/download/{pid}.pdb"
    try:
        print(f"  ↓  fetching {pid}  ({desc})...", end="", flush=True)
        urllib.request.urlretrieve(url, dest)
        print(f"  {GREEN}done{OFF}  ({dest.stat().st_size // 1024} KB)")
    except Exception as e:
        print(f"  {RED}FAIL{OFF}: {e}")
        sys.exit(1)
print()


# ─── Step 2: generate ligand files via RDKit ───────────────────────────────
print(f"{CYAN}[2/3]{OFF}  Generating ligand files ...")
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    print(f"  {RED}RDKit not found.{OFF}  Activate the conda env first:")
    print(f"        conda activate docking-pipeline")
    print(f"        python {Path(__file__).name}")
    sys.exit(1)


LIGANDS = {
    "aspirin":   ("CC(=O)Oc1ccccc1C(=O)O",                "anti-inflammatory"),
    "ibuprofen": ("CC(C)Cc1ccc(C(C)C(O)=O)cc1",           "NSAID"),
    "caffeine":  ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C",         "stimulant"),
    "phenol":    ("Oc1ccccc1",                            "small aromatic"),
}


def build_mol(smi: str):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    AllChem.EmbedMolecule(mol, params)
    try:
        AllChem.MMFFOptimizeMolecule(mol, mmffVariant="MMFF94", maxIters=500)
    except Exception:
        AllChem.UFFOptimizeMolecule(mol, maxIters=500)
    return mol


# write SDF for all four
for name, (smi, desc) in LIGANDS.items():
    mol = build_mol(smi)
    if mol is None:
        say("✗", f"could not build {name}", RED)
        continue
    mol.SetProp("_Name", name)
    sdf = LIG / f"{name}.sdf"
    with Chem.SDWriter(str(sdf)) as w:
        w.write(mol)
    say("✓", f"{sdf.name:18}  {DIM}{desc}{OFF}")

# Convert caffeine.sdf → caffeine.mol2 via Open Babel (if available)
try:
    from openbabel import pybel
    src = next(pybel.readfile("sdf", str(LIG / "caffeine.sdf")))
    src.write("mol2", str(LIG / "caffeine.mol2"), overwrite=True)
    say("✓", f"caffeine.mol2       {DIM}(MOL2 via Open Babel){OFF}")
    # also produce a PDB-format ligand
    src.write("pdb", str(LIG / "caffeine.pdb"), overwrite=True)
    say("✓", f"caffeine.pdb        {DIM}(PDB ligand file){OFF}")
except Exception as e:
    print(f"  {RED}!{OFF}  Open Babel conversion skipped: {e}")

# Plain SMILES file with one molecule
(LIG / "ethanol.smi").write_text("CCO ethanol\n")
say("✓", "ethanol.smi         (plain SMILES file)")
print()


# ─── Step 3: example CSV ───────────────────────────────────────────────────
print(f"{CYAN}[3/3]{OFF}  Writing example CSV (PDB IDs only - guaranteed to work) ...")
csv_path = HERE / "example_csv_jobs.csv"
csv_path.write_text(
    "protein,ligand,name\n"
    "1HSG,CC(=O)Oc1ccccc1C(=O)O,1HSG_aspirin\n"
    "3HTB,CC(C)Cc1ccc(C(C)C(O)=O)cc1,3HTB_ibuprofen\n"
    "4LDE,CN1C=NC2=C1C(=O)N(C(=O)N2C)C,4LDE_caffeine\n"
)
say("✓", csv_path.name)
print()


# ─── Summary ──────────────────────────────────────────────────────────────
print(f"{GREEN}╔══════════════════════════════════════════════════════════════╗{OFF}")
print(f"{GREEN}║                    SETUP COMPLETE                            ║{OFF}")
print(f"{GREEN}╚══════════════════════════════════════════════════════════════╝{OFF}")
print()
print(f"  Proteins ({len(list(PROT.glob('*.pdb')))} files):")
for p in sorted(PROT.glob("*.pdb")):
    print(f"    • {p.relative_to(HERE)}")
print(f"\n  Ligands ({len(list(LIG.iterdir()))} files):")
for l in sorted(LIG.iterdir()):
    print(f"    • {l.relative_to(HERE)}")
print()
print(f"{CYAN}Try one of these commands next:{OFF}\n")
print(f"  {GREEN}#{OFF} CSV mode (3 jobs)")
print(f"  python ../docking_pipeline.py --csv {csv_path.name} \\")
print(f"      --parallel 3 -o ../csv_test/\n")
print(f"  {GREEN}#{OFF} Folder mode - 1-to-1 pairing (3 jobs)")
print(f"  python ../docking_pipeline.py \\")
print(f"      --protein-dir proteins/ --ligand-dir ligands/ \\")
print(f"      --pair-mode zip --parallel 3 -o ../folder_zip_test/\n")
print(f"  {GREEN}#{OFF} Folder mode - cross-docking (3 proteins x N ligands)")
print(f"  python ../docking_pipeline.py \\")
print(f"      --protein-dir proteins/ --ligand-dir ligands/ \\")
print(f"      --pair-mode cross --parallel 4 -o ../folder_cross_test/")
print()
