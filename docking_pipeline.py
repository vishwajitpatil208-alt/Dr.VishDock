#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
  ADVANCED AUTOMATED PROTEIN-LIGAND DOCKING PIPELINE
  ------------------------------------------------------------------------------
  Medvolt Tech Private Limited - Computational Chemistry Assignment
  ------------------------------------------------------------------------------
  Implements the full PDF workflow (5 stages) end-to-end with one command:

    Stage 1 - Protein retrieval / modeling
              accepts: RCSB PDB ID | UniProt ID | FASTA sequence/file | .pdb file
    Stage 2 - Protein fixing & preparation
              PDBFixer  -> missing atoms/residues, remove waters, add H's
              OpenMM    -> energy minimization (steric clash relief)
              Meeko / Open Babel -> Gasteiger charges + .pdbqt receptor
    Stage 3 - Ligand preparation
              accepts: SMILES | .smi | .sdf | .mol2 | .pdb | .mol
              RDKit / Open Babel -> 3D embedding, MMFF94 optimization
              Meeko / Open Babel -> .pdbqt ligand (with torsions)
    Stage 4 - Docking
              AutoDock Vina (python API or binary fallback)
              Auto binding-site detection (co-crystal HETATM or center-of-mass)
              Multiple poses, scored & ranked
    Stage 5 - Visualization & Report
              PyMOL (open-source) -> publication-quality complex images
              H-bond / hydrophobic interaction analysis
              ReportLab -> single PDF report per run + global summary.csv

  Modes:
    * Single        :   --protein <X> --ligand <Y>
    * CSV-driven    :   --csv pairs.csv             (parallel across rows)
    * Folder x Folder:  --protein-dir P/ --ligand-dir L/   (parallel pairs)

  Author : Generated for Medvolt assignment
  Python : 3.10+
================================================================================
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
import traceback
import urllib.error
import urllib.request
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

warnings.filterwarnings("ignore")

# ----------------------------------------------------------------------------
#  Third-party imports - imported lazily where possible so the script can at
#  least show --help even if a dependency is missing.
# ----------------------------------------------------------------------------
try:
    import numpy as np
except ImportError:                                                  # pragma: no cover
    np = None

try:
    import pandas as pd
except ImportError:                                                  # pragma: no cover
    pd = None

try:
    from tqdm import tqdm
except ImportError:                                                  # pragma: no cover
    def tqdm(x, **_kw):  # minimal stand-in
        return x

# Rich - for the hacker/gamer terminal UI. All UI is optional; pipeline keeps
# running if rich is missing, falling back to plain text.
try:
    from rich.console import Console, Group
    from rich.panel  import Panel
    from rich.text   import Text
    from rich.table  import Table
    from rich.align  import Align
    from rich.live   import Live
    from rich.progress import (
        Progress, SpinnerColumn, BarColumn, TextColumn,
        TimeElapsedColumn, MofNCompleteColumn,
    )
    from rich.rule    import Rule
    from rich.columns import Columns
    from rich import box as _rbox
    from rich.theme  import Theme
    RICH = True
    _theme = Theme({
        "ok":     "bold bright_green",
        "warn":   "bold yellow",
        "err":    "bold red",
        "info":   "bold cyan",
        "dim":    "grey50",
        "neon":   "bold spring_green1",
        "accent": "bold magenta",
        "label":  "bold cyan",
        "value":  "bright_white",
    })
    console = Console(theme=_theme, highlight=False)
except ImportError:                                                  # pragma: no cover
    RICH = False
    console = None


# ============================================================================
#  GLOBAL CONFIG
# ============================================================================

PIPELINE_VERSION = "2.0.1"
DEFAULT_EXHAUSTIVENESS = 16
DEFAULT_NUM_POSES = 9
DEFAULT_BOX_PADDING = 5.0          # Å around binding site
DEFAULT_BLIND_BOX_SIZE = 25.0      # Å when no pocket info available
DEFAULT_GRID_SPACING = 0.375       # Å, Vina default

RCSB_URL   = "https://files.rcsb.org/download/{pdb_id}.pdb"
AF_URL     = "https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v4.pdb"
UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{uid}.fasta"
ESMFOLD_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"

PDB_ID_RE       = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
UNIPROT_ID_RE   = re.compile(
    r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$"
)
SMILES_CHARS_RE = re.compile(r"^[A-Za-z0-9@+\-\[\]\(\)=#\\/\.%\*:]+$")


# ============================================================================
#  HACKER / GAMER TERMINAL UI
# ============================================================================

# Embedded ASCII art (no pyfiglet dependency required)
_MEDVOLT_BANNER = r"""
███╗   ███╗███████╗██████╗ ██╗   ██╗ ██████╗ ██╗  ████████╗     █████╗ ██╗
████╗ ████║██╔════╝██╔══██╗██║   ██║██╔═══██╗██║  ╚══██╔══╝    ██╔══██╗██║
██╔████╔██║█████╗  ██║  ██║██║   ██║██║   ██║██║     ██║   ██╗ ███████║██║
██║╚██╔╝██║██╔══╝  ██║  ██║╚██╗ ██╔╝██║   ██║██║     ██║   ╚═╝ ██╔══██║██║
██║ ╚═╝ ██║███████╗██████╔╝ ╚████╔╝ ╚██████╔╝███████╗ ██║       ██║  ██║██║
╚═╝     ╚═╝╚══════╝╚═════╝   ╚═══╝   ╚═════╝ ╚══════╝ ╚═╝       ╚═╝  ╚═╝╚═╝
"""

_TAGLINE = "  P R O T E I N • L I G A N D • D O C K I N G • S U I T E  "

# Stage glyph + colour map for status output
_STAGE_META = {
    "stage1": ("⟦ 01 ⟧", "PROTEIN  ACQUISITION  ", "bright_cyan"),
    "stage2": ("⟦ 02 ⟧", "STRUCTURE  REFINEMENT ", "bright_green"),
    "stage3": ("⟦ 03 ⟧", "LIGAND  PREPARATION   ", "bright_magenta"),
    "stage4": ("⟦ 04 ⟧", "DOCKING  EXECUTION    ", "bright_yellow"),
    "stage5": ("⟦ 05 ⟧", "VISUAL  &  REPORT     ", "bright_red"),
}


class HackerUI:
    """Pro-level cyberpunk terminal UI built on rich. Safe no-op if rich missing."""

    def __init__(self, enabled: bool = True):
        self.enabled = enabled and RICH and console is not None
        self.c = console if self.enabled else None

    # -- helpers -------------------------------------------------------------
    def _say(self, *parts):
        if not self.enabled:
            text = " ".join(str(p) for p in parts if not isinstance(p, dict))
            print(text)
            return
        self.c.print(*parts)

    # -- banner --------------------------------------------------------------
    def banner(self):
        if not self.enabled:
            print("=" * 72)
            print("  M E D V O L T . A I   --  Protein-Ligand Docking Suite")
            print(f"  v{PIPELINE_VERSION}")
            print("=" * 72)
            return

        gradient_lines = []
        colors = ["spring_green1", "spring_green3", "green1",
                  "bright_green", "green3", "cyan1", "bright_cyan"]
        for i, line in enumerate(_MEDVOLT_BANNER.strip("\n").splitlines()):
            gradient_lines.append(Text(line, style=f"bold {colors[i % len(colors)]}"))
        body = Group(
            *gradient_lines,
            Text(""),
            Align.center(Text(_TAGLINE, style="bold black on bright_green")),
            Text(""),
            Align.center(Text(
                f"⚡ v{PIPELINE_VERSION}  •  ⟨ Computational Chemistry Engine ⟩  •  © MedVolt Tech",
                style="dim bright_cyan",
            )),
        )
        self.c.print(Panel(
            Align.center(body),
            border_style="bright_green",
            box=_rbox.DOUBLE_EDGE,
            padding=(1, 2),
            title="[bold black on bright_green] SYSTEM ONLINE [/]",
            subtitle="[bold black on cyan] AUTH: ROOT [/]",
        ))

    # -- boot sequence -------------------------------------------------------
    def boot_sequence(self, output_dir: str, n_jobs: int, parallel: int):
        if not self.enabled:
            print(f"Output dir: {output_dir}")
            print(f"Jobs: {n_jobs}  |  Workers: {parallel}")
            return
        lines = [
            ("[bright_green]>[/] initializing docking subsystems ........ [ok]OK[/]"),
            ("[bright_green]>[/] loading toolchain (rdkit / vina / openmm) [ok]OK[/]"),
            ("[bright_green]>[/] establishing remote endpoints (RCSB / AF / ESM) [ok]OK[/]"),
            ("[bright_green]>[/] forging output directory ............... [ok]OK[/]"),
        ]
        for ln in lines:
            self.c.print(ln)
        cfg = Table.grid(padding=(0, 2))
        cfg.add_column(style="label", justify="right")
        cfg.add_column(style="value")
        cfg.add_row("◆ OUTPUT_DIR",     f"[neon]{output_dir}[/]")
        cfg.add_row("◆ JOBS_QUEUED",    f"[neon]{n_jobs}[/]")
        cfg.add_row("◆ WORKERS",        f"[neon]{parallel}[/]")
        cfg.add_row("◆ SESSION_TIME",   f"[neon]{datetime.now():%Y-%m-%d %H:%M:%S}[/]")
        self.c.print(Panel(cfg, title="[bold]◤ MISSION BRIEFING ◥[/]",
                           border_style="cyan", box=_rbox.HEAVY))

    # -- per-job header ------------------------------------------------------
    def job_header(self, idx: int, total: int, name: str, protein: str, ligand: str):
        if not self.enabled:
            print(f"\n>>> JOB [{idx}/{total}]  {name}")
            return
        t = Table.grid(padding=(0, 2))
        t.add_column(style="label", justify="right")
        t.add_column(style="value")
        t.add_row("⟨ JOB ⟩",     f"[neon]{idx}/{total}[/]   [accent]{name}[/]")
        t.add_row("⟨ PROTEIN ⟩", f"[bright_white]{protein[:80]}[/]")
        t.add_row("⟨ LIGAND ⟩",  f"[bright_white]{ligand[:80]}[/]")
        self.c.print(Panel(t, border_style="bright_magenta", box=_rbox.HEAVY,
                           title="[bold black on bright_magenta] ▶ EXECUTING ◀ [/]"))

    # -- stage indicator -----------------------------------------------------
    def stage(self, key: str, msg: str = ""):
        if not self.enabled:
            print(f"  [{key}] {msg}")
            return
        glyph, label, colr = _STAGE_META.get(key, ("⟦??⟧", key, "white"))
        self.c.print(
            f"  [bold {colr}]{glyph}[/]  [bold {colr}]{label}[/]  [dim]{msg}[/]"
        )

    # -- status lines --------------------------------------------------------
    def ok(self, msg: str):
        if self.enabled:
            self.c.print(f"    [ok]✔[/]  [bright_white]{msg}[/]")
        else:
            print(f"    [OK] {msg}")

    def warn(self, msg: str):
        if self.enabled:
            self.c.print(f"    [warn]▲[/]  [yellow]{msg}[/]")
        else:
            print(f"    [WARN] {msg}")

    def fail(self, msg: str):
        if self.enabled:
            self.c.print(f"    [err]✖[/]  [red]{msg}[/]")
        else:
            print(f"    [FAIL] {msg}")

    def divider(self):
        if self.enabled:
            self.c.print(Rule(style="bright_green"))
        else:
            print("-" * 72)

    # -- final summary table -------------------------------------------------
    def summary(self, results: List[Dict[str, Any]], out_dir: str):
        if not self.enabled:
            n_ok = sum(1 for r in results if r.get("success"))
            print(f"\nDONE. {n_ok}/{len(results)} jobs succeeded. Results: {out_dir}")
            return
        n_ok   = sum(1 for r in results if r.get("success"))
        n_fail = len(results) - n_ok

        t = Table(
            title="[bold bright_green]◤ DOCKING RESULTS ◥[/]",
            box=_rbox.DOUBLE_EDGE, border_style="bright_green",
            header_style="bold black on bright_green",
            title_style="bold bright_green",
            show_lines=False, expand=True,
        )
        t.add_column("#",        style="dim",           justify="right", no_wrap=True)
        t.add_column("JOB",      style="bold magenta",  no_wrap=False)
        t.add_column("STATUS",   justify="center")
        t.add_column("ΔG (kcal/mol)", style="bold bright_green", justify="right")
        t.add_column("POSES",    justify="right", style="cyan")
        t.add_column("H-BONDS",  justify="right", style="yellow")
        t.add_column("HYDRO",    justify="right", style="magenta")
        t.add_column("SITE",     style="dim")

        for i, r in enumerate(results, 1):
            status = "[ok]✔ OK[/]" if r.get("success") else "[err]✖ FAIL[/]"
            best = r.get("best_affinity")
            best_str = f"{best:.2f}" if isinstance(best, (int, float)) else "—"
            t.add_row(
                str(i),
                str(r.get("name", "?"))[:30],
                status,
                best_str,
                str(len(r.get("affinities", []) or [])),
                str(r.get("n_hbonds", 0)),
                str(r.get("n_hydrophobic", 0)),
                str(r.get("site_source", "")),
            )
        self.c.print(t)

        verdict_style = "ok" if n_fail == 0 else ("warn" if n_ok > 0 else "err")
        verdict_text  = "MISSION COMPLETE" if n_fail == 0 else (
            "MISSION PARTIAL" if n_ok > 0 else "MISSION FAILED"
        )
        self.c.print(Panel(
            Align.center(Text(
                f"  ►  {verdict_text}  ◄        {n_ok} success   /   {n_fail} failed   /   {len(results)} total",
                style=f"bold {verdict_style}",
            )),
            border_style="bright_green", box=_rbox.DOUBLE,
        ))
        self.c.print(f"  [dim]all artifacts saved to[/]  [neon]{out_dir}[/]")


# ============================================================================
#  LOGGING
# ============================================================================

def setup_logger(out_dir: Path, name: str = "docking", verbose: bool = True) -> logging.Logger:
    """Console + file logger; safe to call multiple times."""
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    logger.handlers.clear()

    fmt = logging.Formatter(
        "[%(asctime)s] [%(levelname)-7s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    # When rich UI is active we route status through the UI itself; keep stream
    # handler at WARNING so the progress bar/panels are not interrupted.
    sh.setLevel(logging.WARNING if RICH else logging.INFO)
    logger.addHandler(sh)

    out_dir.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(out_dir / "pipeline.log", mode="a")
    fh.setFormatter(fmt)
    fh.setLevel(logging.DEBUG)
    logger.addHandler(fh)

    logger.propagate = False
    return logger


# ============================================================================
#  UTILITIES
# ============================================================================

def safe_name(s: str) -> str:
    """Make a string safe for use as a filename."""
    s = str(s).strip()
    s = re.sub(r"[^A-Za-z0-9_.\-]+", "_", s)
    return s[:64] if len(s) > 64 else s


def http_download(url: str, dest: Path, retries: int = 3, timeout: int = 60) -> bool:
    """Robust HTTP download with retry."""
    last_err: Optional[Exception] = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "DockingPipeline/1.0"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = resp.read()
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(data)
            return True
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError) as e:
            last_err = e
    raise RuntimeError(f"Download failed after {retries} attempts: {url}  ({last_err})")


def http_post(url: str, payload: Union[str, bytes], retries: int = 3, timeout: int = 300) -> bytes:
    """Robust HTTP POST (used for ESMFold)."""
    data = payload.encode() if isinstance(payload, str) else payload
    last_err: Optional[Exception] = None
    for _ in range(retries):
        try:
            req = urllib.request.Request(
                url, data=data,
                headers={"User-Agent": "DockingPipeline/1.0", "Content-Type": "text/plain"},
            )
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except Exception as e:
            last_err = e
    raise RuntimeError(f"POST failed: {url}  ({last_err})")


def which(cmd: str) -> Optional[str]:
    return shutil.which(cmd)


def run_cmd(cmd: List[str], cwd: Optional[Path] = None, timeout: int = 600) -> Tuple[int, str, str]:
    """Run a subprocess, capture stdout/stderr."""
    try:
        p = subprocess.run(
            cmd, cwd=str(cwd) if cwd else None,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            timeout=timeout, check=False, text=True,
        )
        return p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired as e:
        return 124, "", f"Timeout: {e}"


# ============================================================================
#  INPUT TYPE DETECTION
# ============================================================================

def detect_protein_input(value: str) -> str:
    """
    Return one of: 'pdb_file', 'pdb_id', 'uniprot_id',
                   'fasta_file', 'fasta_seq'
    """
    p = Path(value)
    if p.exists() and p.is_file():
        suf = p.suffix.lower()
        if suf in {".pdb", ".ent"}:
            return "pdb_file"
        if suf in {".fasta", ".fa", ".faa", ".seq", ".txt"}:
            return "fasta_file"
        # Try sniff content
        try:
            head = p.read_text(errors="ignore")[:2048]
            if head.startswith(">") or all(c.isalpha() or c.isspace() for c in head[:200]):
                return "fasta_file"
            if "ATOM" in head or "HEADER" in head:
                return "pdb_file"
        except Exception:
            pass
        raise ValueError(f"Cannot determine file type: {value}")

    s = value.strip().upper()
    if PDB_ID_RE.match(s):
        return "pdb_id"
    if UNIPROT_ID_RE.match(s):
        return "uniprot_id"
    # Looks like an amino-acid sequence?
    if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWYBXZUO\s]+", s) and len(s.replace(" ", "")) >= 10:
        return "fasta_seq"
    raise ValueError(f"Could not classify protein input: {value!r}")


def detect_ligand_input(value: str) -> str:
    """Return one of: 'sdf', 'mol2', 'pdb', 'mol', 'smi_file', 'smiles'."""
    p = Path(value)
    if p.exists() and p.is_file():
        suf = p.suffix.lower().lstrip(".")
        if suf in {"sdf", "mol2", "pdb", "mol"}:
            return suf
        if suf in {"smi", "smiles"}:
            return "smi_file"
        raise ValueError(f"Unsupported ligand file type: {value}")
    # treat as SMILES
    if SMILES_CHARS_RE.match(value.strip()):
        return "smiles"
    raise ValueError(f"Could not classify ligand input: {value!r}")


# ============================================================================
#  STAGE 1+2 : PROTEIN PREPARATION
# ============================================================================

class ProteinPreparer:
    """Resolves any protein input to a fully-prepared receptor.pdbqt."""

    def __init__(self, log: logging.Logger, work_dir: Path, minimize: bool = True):
        self.log = log
        self.work = work_dir
        self.work.mkdir(parents=True, exist_ok=True)
        self.minimize = minimize
        self.cocrystal_center: Optional[Tuple[float, float, float]] = None
        self.cocrystal_resname: Optional[str] = None

    # ---- Stage 1 -----------------------------------------------------------

    def resolve(self, value: str) -> Path:
        """Return a raw .pdb file for any supported input."""
        kind = detect_protein_input(value)
        self.log.info(f"Protein input classified as: {kind}")
        raw = self.work / "01_protein_raw.pdb"

        if kind == "pdb_file":
            shutil.copy(value, raw)

        elif kind == "pdb_id":
            url = RCSB_URL.format(pdb_id=value.upper())
            self.log.info(f"Downloading from RCSB: {url}")
            http_download(url, raw)

        elif kind == "uniprot_id":
            url = AF_URL.format(uid=value.upper())
            self.log.info(f"Trying AlphaFold DB: {url}")
            try:
                http_download(url, raw)
            except Exception as e:
                self.log.warning(f"AlphaFold DB miss ({e}); falling back to ESMFold")
                fasta_url = UNIPROT_FASTA.format(uid=value.upper())
                fasta_path = self.work / "uniprot.fasta"
                http_download(fasta_url, fasta_path)
                seq = self._read_fasta(fasta_path)
                self._fold_esm(seq, raw)

        elif kind == "fasta_file":
            seq = self._read_fasta(Path(value))
            self._fold_esm(seq, raw)

        elif kind == "fasta_seq":
            self._fold_esm(value.strip().replace(" ", ""), raw)

        else:
            raise ValueError(f"Unhandled protein kind: {kind}")

        # snag co-crystal ligand info for binding-site fallback
        try:
            self._detect_cocrystal(raw)
        except Exception as e:
            self.log.debug(f"Co-crystal detection skipped: {e}")

        return raw

    # ---- helpers -----------------------------------------------------------

    @staticmethod
    def _read_fasta(path: Path) -> str:
        text = Path(path).read_text(errors="ignore")
        if text.startswith(">"):
            lines = [ln.strip() for ln in text.splitlines() if not ln.startswith(">")]
            return "".join(lines)
        return "".join(text.split())

    def _fold_esm(self, seq: str, dest: Path) -> None:
        if len(seq) > 400:
            self.log.warning(f"Sequence length {len(seq)} > 400; ESMFold public API may reject it.")
        self.log.info(f"Folding sequence (len={len(seq)}) with ESMFold public API...")
        data = http_post(ESMFOLD_URL, seq, timeout=600)
        if b"ATOM" not in data[:5000]:
            raise RuntimeError("ESMFold returned no PDB. Try a shorter sequence or local model.")
        dest.write_bytes(data)

    def _detect_cocrystal(self, pdb_path: Path) -> None:
        """Scan HETATM (excluding waters/common ions) for the largest ligand."""
        skip = {"HOH", "WAT", "DOD", "NA", "CL", "K", "MG", "CA", "ZN",
                "MN", "FE", "CU", "SO4", "PO4", "GOL", "EDO", "PEG", "ACT", "DMS"}
        groups: Dict[str, List[Tuple[float, float, float]]] = {}
        for line in pdb_path.read_text(errors="ignore").splitlines():
            if line.startswith("HETATM"):
                resname = line[17:20].strip()
                if resname in skip:
                    continue
                try:
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                except ValueError:
                    continue
                groups.setdefault(resname, []).append((x, y, z))
        if not groups:
            return
        # Pick the group with the most heavy atoms
        resname, coords = max(groups.items(), key=lambda kv: len(kv[1]))
        if len(coords) < 5:
            return
        cx = sum(c[0] for c in coords) / len(coords)
        cy = sum(c[1] for c in coords) / len(coords)
        cz = sum(c[2] for c in coords) / len(coords)
        self.cocrystal_center = (cx, cy, cz)
        self.cocrystal_resname = resname
        self.log.info(f"Co-crystal ligand detected: {resname}  center=({cx:.2f},{cy:.2f},{cz:.2f})")

    # ---- Stage 2 -----------------------------------------------------------

    def fix(self, raw_pdb: Path) -> Path:
        """Remove waters, add missing atoms/residues/hydrogens via PDBFixer."""
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile

        self.log.info("Running PDBFixer (missing residues/atoms, hydrogens, strip waters)...")
        fixer = PDBFixer(filename=str(raw_pdb))
        fixer.findMissingResidues()
        # drop terminal missing residues (PDBFixer convention)
        chains = list(fixer.topology.chains())
        keys_to_del = []
        for key in list(fixer.missingResidues.keys()):
            chain_idx, ins_idx = key
            chain = chains[chain_idx]
            nres = len(list(chain.residues()))
            if ins_idx == 0 or ins_idx == nres:
                keys_to_del.append(key)
        for k in keys_to_del:
            del fixer.missingResidues[k]

        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(keepWater=False)
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(pH=7.4)

        fixed = self.work / "02_protein_fixed.pdb"
        with open(fixed, "w") as fh:
            PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
        self.log.info(f"Fixed protein written: {fixed.name}")
        return fixed

    def minimize_structure(self, fixed_pdb: Path) -> Path:
        """Light OpenMM energy minimization to relieve steric clashes."""
        if not self.minimize:
            return fixed_pdb
        try:
            from openmm import app, unit, LangevinIntegrator, Platform
            from openmm.app import ForceField, Modeller, PDBFile, Simulation
        except ImportError:
            self.log.warning("OpenMM not installed; skipping minimization.")
            return fixed_pdb

        self.log.info("Energy minimization (OpenMM, amber14, implicit solvent)...")
        try:
            pdb = PDBFile(str(fixed_pdb))
            ff = ForceField("amber14-all.xml", "implicit/gbn2.xml")
            modeller = Modeller(pdb.topology, pdb.positions)
            modeller.addHydrogens(ff, pH=7.4)
            system = ff.createSystem(
                modeller.topology, nonbondedMethod=app.NoCutoff,
                constraints=app.HBonds,
            )
            integrator = LangevinIntegrator(
                300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds
            )
            try:
                platform = Platform.getPlatformByName("CPU")
                sim = Simulation(modeller.topology, system, integrator, platform)
            except Exception:
                sim = Simulation(modeller.topology, system, integrator)
            sim.context.setPositions(modeller.positions)
            sim.minimizeEnergy(maxIterations=500)
            state = sim.context.getState(getPositions=True)
            out = self.work / "03_protein_minimized.pdb"
            with open(out, "w") as fh:
                PDBFile.writeFile(modeller.topology, state.getPositions(), fh, keepIds=True)
            self.log.info(f"Minimized protein written: {out.name}")
            return out
        except Exception as e:
            self.log.warning(f"Minimization failed ({e}); using PDBFixer output.")
            return fixed_pdb

    def to_pdbqt(self, prepared_pdb: Path) -> Path:
        """Convert receptor to .pdbqt (Gasteiger charges)."""
        out = self.work / "04_receptor.pdbqt"

        # Strategy 1: Meeko's receptor preparation if available
        if which("mk_prepare_receptor.py"):
            rc, _so, se = run_cmd([
                "mk_prepare_receptor.py", "--read_pdb", str(prepared_pdb),
                "-o", str(out.with_suffix("")), "-p",
            ])
            if rc == 0 and out.exists():
                self.log.info("Receptor PDBQT via mk_prepare_receptor.py")
                return out
            self.log.debug(f"mk_prepare_receptor.py failed (rc={rc}): {se[:200]}")

        # Strategy 2: Open Babel python API
        try:
            from openbabel import pybel
            mol = next(pybel.readfile("pdb", str(prepared_pdb)))
            mol.OBMol.AddPolarHydrogens()
            mol.calccharges("gasteiger")
            mol.write("pdbqt", str(out), overwrite=True, opt={"r": True, "x": True})
            # ensure file looks like a receptor (no ROOT/TORSDOF lines)
            self._sanitize_receptor_pdbqt(out)
            self.log.info("Receptor PDBQT via Open Babel (pybel)")
            return out
        except Exception as e:
            self.log.debug(f"pybel route failed: {e}")

        # Strategy 3: obabel CLI
        if which("obabel"):
            rc, _so, se = run_cmd([
                "obabel", str(prepared_pdb), "-O", str(out),
                "-xr", "--partialcharge", "gasteiger",
            ])
            if rc == 0 and out.exists():
                self._sanitize_receptor_pdbqt(out)
                self.log.info("Receptor PDBQT via obabel CLI")
                return out

        raise RuntimeError("Could not produce receptor PDBQT (install meeko or open-babel).")

    @staticmethod
    def _sanitize_receptor_pdbqt(p: Path) -> None:
        """Strip ROOT/ENDROOT/TORSDOF/BRANCH lines that obabel sometimes writes for receptors."""
        bad = ("ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF")
        lines = p.read_text().splitlines()
        keep = [ln for ln in lines if not ln.strip().startswith(bad)]
        p.write_text("\n".join(keep) + "\n")


# ============================================================================
#  BINDING SITE DETECTION
# ============================================================================

class BindingSiteDetector:
    """Decides the docking grid box."""

    def __init__(self, log: logging.Logger):
        self.log = log

    def determine(
        self,
        receptor_pdb: Path,
        user_center: Optional[Tuple[float, float, float]] = None,
        user_size: Optional[Tuple[float, float, float]] = None,
        cocrystal_center: Optional[Tuple[float, float, float]] = None,
        padding: float = DEFAULT_BOX_PADDING,
    ) -> Tuple[Tuple[float, float, float], Tuple[float, float, float], str]:
        """
        Return (center, size, source_description).
        Priority: user -> co-crystal -> protein-COM (blind).
        """
        if user_center and user_size:
            self.log.info(f"Using user-supplied grid: center={user_center} size={user_size}")
            return user_center, user_size, "user-specified"

        if cocrystal_center:
            size = (20.0, 20.0, 20.0)
            self.log.info(f"Using co-crystal ligand center: {cocrystal_center}, size={size}")
            return cocrystal_center, size, "co-crystal ligand"

        # Blind docking: protein center-of-mass, large box
        coords: List[Tuple[float, float, float]] = []
        for line in receptor_pdb.read_text(errors="ignore").splitlines():
            if line.startswith(("ATOM", "HETATM")):
                try:
                    coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                except ValueError:
                    continue
        if not coords:
            raise RuntimeError("No atoms parsed from receptor.")
        cx = sum(c[0] for c in coords) / len(coords)
        cy = sum(c[1] for c in coords) / len(coords)
        cz = sum(c[2] for c in coords) / len(coords)
        xs = [c[0] for c in coords]
        ys = [c[1] for c in coords]
        zs = [c[2] for c in coords]
        # Use bounding box but cap at a reasonable size (Vina max 126 Å)
        sx = min(max(max(xs) - min(xs) + padding, 25.0), 60.0)
        sy = min(max(max(ys) - min(ys) + padding, 25.0), 60.0)
        sz = min(max(max(zs) - min(zs) + padding, 25.0), 60.0)
        self.log.info(
            f"Blind docking: center=({cx:.2f},{cy:.2f},{cz:.2f}) size=({sx:.1f},{sy:.1f},{sz:.1f})"
        )
        return (cx, cy, cz), (sx, sy, sz), "blind (protein COM)"


# ============================================================================
#  STAGE 3 : LIGAND PREPARATION
# ============================================================================

class LigandPreparer:
    """SMILES / SDF / MOL2 / PDB / MOL -> 3D, optimized, .pdbqt."""

    def __init__(self, log: logging.Logger, work_dir: Path):
        self.log = log
        self.work = work_dir
        self.work.mkdir(parents=True, exist_ok=True)

    def prepare(self, value: str, name: str = "ligand") -> Tuple[Path, Path]:
        """Return (ligand_3d_sdf, ligand_pdbqt)."""
        kind = detect_ligand_input(value)
        self.log.info(f"Ligand input classified as: {kind}")
        sdf_out = self.work / "05_ligand_3d.sdf"
        pdbqt_out = self.work / "06_ligand.pdbqt"

        # Build / read RDKit Mol
        if kind == "smiles":
            mol = self._from_smiles(value)
        elif kind == "smi_file":
            text = Path(value).read_text().strip().splitlines()[0].split()[0]
            mol = self._from_smiles(text)
        else:
            mol = self._from_file(Path(value), kind)

        # Generate 3D conformer & MMFF optimize
        mol = self._embed_and_optimize(mol)

        # Write SDF
        self._write_sdf(mol, sdf_out, name)

        # Convert SDF -> PDBQT (Meeko preferred, OpenBabel fallback)
        self._sdf_to_pdbqt(sdf_out, pdbqt_out)
        return sdf_out, pdbqt_out

    # ---- builders ----------------------------------------------------------

    def _from_smiles(self, smi: str):
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smi!r}")
        mol = Chem.AddHs(mol)
        return mol

    def _from_file(self, path: Path, kind: str):
        from rdkit import Chem
        if kind == "sdf":
            suppl = Chem.SDMolSupplier(str(path), removeHs=False)
            mols = [m for m in suppl if m is not None]
            if not mols:
                raise ValueError(f"No molecules in SDF: {path}")
            mol = mols[0]
        elif kind == "mol":
            mol = Chem.MolFromMolFile(str(path), removeHs=False)
        elif kind == "pdb":
            mol = Chem.MolFromPDBFile(str(path), removeHs=False, sanitize=True)
        elif kind == "mol2":
            try:
                mol = Chem.MolFromMol2File(str(path), removeHs=False, sanitize=True)
            except Exception:
                mol = None
            if mol is None:
                # convert via OpenBabel to SDF then read
                from openbabel import pybel
                obmol = next(pybel.readfile("mol2", str(path)))
                tmp_sdf = self.work / "_tmp_lig.sdf"
                obmol.write("sdf", str(tmp_sdf), overwrite=True)
                mol = Chem.MolFromMolFile(str(tmp_sdf), removeHs=False)
        else:
            raise ValueError(f"Unknown ligand kind: {kind}")

        if mol is None:
            raise ValueError(f"Could not parse ligand: {path}")
        mol = Chem.AddHs(mol, addCoords=True)
        return mol

    def _embed_and_optimize(self, mol):
        from rdkit.Chem import AllChem
        # if conformer already present (from file), just optimize
        if mol.GetNumConformers() == 0:
            params = AllChem.ETKDGv3()
            params.randomSeed = 0xC0FFEE
            if AllChem.EmbedMolecule(mol, params) != 0:
                # fallback random coords
                AllChem.EmbedMolecule(mol, useRandomCoords=True)
        try:
            AllChem.MMFFOptimizeMolecule(mol, mmffVariant="MMFF94", maxIters=500)
        except Exception:
            try:
                AllChem.UFFOptimizeMolecule(mol, maxIters=500)
            except Exception:
                self.log.warning("Force-field optimization failed; using raw coordinates.")
        return mol

    def _write_sdf(self, mol, path: Path, name: str) -> None:
        from rdkit import Chem
        mol.SetProp("_Name", safe_name(name))
        with Chem.SDWriter(str(path)) as w:
            w.write(mol)

    def _sdf_to_pdbqt(self, sdf_in: Path, pdbqt_out: Path) -> None:
        # Try Meeko first (preferred for Vina; handles torsions correctly)
        try:
            from meeko import MoleculePreparation, PDBQTWriterLegacy  # type: ignore
            from rdkit import Chem
            mol = next(Chem.SDMolSupplier(str(sdf_in), removeHs=False))
            prep = MoleculePreparation()
            try:
                # newer meeko API (>=0.5)
                mol_setups = prep.prepare(mol)
                pdbqt_text, _, _ = PDBQTWriterLegacy.write_string(mol_setups[0])
            except Exception:
                # older meeko API
                prep.prepare(mol)
                pdbqt_text = prep.write_pdbqt_string()
            pdbqt_out.write_text(pdbqt_text)
            self.log.info("Ligand PDBQT via Meeko")
            return
        except Exception as e:
            self.log.debug(f"Meeko ligand prep failed: {e}; falling back to OpenBabel.")

        # OpenBabel fallback
        try:
            from openbabel import pybel
            mol = next(pybel.readfile("sdf", str(sdf_in)))
            mol.OBMol.AddPolarHydrogens()
            mol.calccharges("gasteiger")
            mol.write("pdbqt", str(pdbqt_out), overwrite=True)
            self.log.info("Ligand PDBQT via Open Babel")
            return
        except Exception as e:
            self.log.debug(f"pybel SDF->PDBQT failed: {e}")

        if which("obabel"):
            rc, _, se = run_cmd([
                "obabel", str(sdf_in), "-O", str(pdbqt_out),
                "--partialcharge", "gasteiger", "-h",
            ])
            if rc == 0 and pdbqt_out.exists():
                self.log.info("Ligand PDBQT via obabel CLI")
                return

        raise RuntimeError("Could not produce ligand PDBQT (install meeko or open-babel).")


# ============================================================================
#  STAGE 4 : DOCKING (AutoDock Vina)
# ============================================================================

@dataclass
class DockingResult:
    poses_pdbqt: Path
    log_text: str
    affinities: List[float] = field(default_factory=list)  # kcal/mol per pose
    best_affinity: Optional[float] = None
    center: Tuple[float, float, float] = (0, 0, 0)
    size: Tuple[float, float, float] = (0, 0, 0)
    site_source: str = ""
    engine: str = ""


class DockingEngine:
    """AutoDock Vina via Python API or CLI."""

    def __init__(self, log: logging.Logger, work_dir: Path,
                 exhaustiveness: int = DEFAULT_EXHAUSTIVENESS,
                 num_poses: int = DEFAULT_NUM_POSES,
                 cpu: int = 0):
        self.log = log
        self.work = work_dir
        self.exhaustiveness = exhaustiveness
        self.num_poses = num_poses
        self.cpu = cpu  # 0 = all

    def run(self, receptor: Path, ligand: Path,
            center: Tuple[float, float, float],
            size: Tuple[float, float, float]) -> DockingResult:
        out = self.work / "07_docked.pdbqt"
        result = DockingResult(
            poses_pdbqt=out, log_text="",
            center=center, size=size,
        )

        # Try Python API
        try:
            from vina import Vina  # type: ignore
            self.log.info("Running AutoDock Vina (Python API)...")
            v = Vina(sf_name="vina", cpu=self.cpu, verbosity=1)
            v.set_receptor(str(receptor))
            v.set_ligand_from_file(str(ligand))
            v.compute_vina_maps(center=list(center), box_size=list(size))
            v.dock(exhaustiveness=self.exhaustiveness, n_poses=self.num_poses)
            v.write_poses(str(out), n_poses=self.num_poses, overwrite=True)
            try:
                energies = v.energies(n_poses=self.num_poses)
                # energies: list of [affinity, ...]
                result.affinities = [float(row[0]) for row in energies]
            except Exception:
                result.affinities = self._parse_affinities(out)
            result.engine = "vina-python"
            result.log_text = f"Vina Python API; exhaustiveness={self.exhaustiveness}"
        except Exception as e:
            self.log.warning(f"Vina Python API unavailable ({e}); trying CLI.")
            result = self._run_cli(receptor, ligand, center, size, out, result)

        if result.affinities:
            result.best_affinity = min(result.affinities)
            self.log.info(
                f"Docking complete. Best affinity: {result.best_affinity:.2f} kcal/mol "
                f"({len(result.affinities)} poses)"
            )
        return result

    # ---- CLI fallback ------------------------------------------------------

    def _run_cli(self, receptor, ligand, center, size, out, result) -> DockingResult:
        if not which("vina"):
            raise RuntimeError("AutoDock Vina is not installed (neither python nor CLI).")
        cmd = [
            "vina",
            "--receptor", str(receptor),
            "--ligand", str(ligand),
            "--out", str(out),
            "--center_x", f"{center[0]:.3f}",
            "--center_y", f"{center[1]:.3f}",
            "--center_z", f"{center[2]:.3f}",
            "--size_x", f"{size[0]:.3f}",
            "--size_y", f"{size[1]:.3f}",
            "--size_z", f"{size[2]:.3f}",
            "--exhaustiveness", str(self.exhaustiveness),
            "--num_modes", str(self.num_poses),
        ]
        if self.cpu > 0:
            cmd += ["--cpu", str(self.cpu)]
        rc, so, se = run_cmd(cmd, timeout=3600)
        if rc != 0:
            raise RuntimeError(f"Vina CLI failed: {se[:500]}")
        result.log_text = so
        result.affinities = self._parse_affinities(out)
        result.engine = "vina-cli"
        return result

    @staticmethod
    def _parse_affinities(pdbqt: Path) -> List[float]:
        """Parse 'REMARK VINA RESULT' lines."""
        out: List[float] = []
        if not pdbqt.exists():
            return out
        for line in pdbqt.read_text().splitlines():
            if line.startswith("REMARK VINA RESULT"):
                parts = line.split()
                # REMARK VINA RESULT:  -8.234  0.000  0.000
                try:
                    out.append(float(parts[3]))
                except (IndexError, ValueError):
                    pass
        return out


# ============================================================================
#  POSE POST-PROCESSING + INTERACTION ANALYSIS
# ============================================================================

class PoseProcessor:
    """Convert PDBQT poses to SDF + build protein-ligand complex PDB."""

    def __init__(self, log: logging.Logger, work_dir: Path):
        self.log = log
        self.work = work_dir

    def split_and_convert(self, poses_pdbqt: Path) -> Path:
        """Convert multi-model PDBQT into a single SDF (best pose first)."""
        out_sdf = self.work / "08_docked.sdf"
        # Try Meeko
        try:
            from meeko import PDBQTMolecule, RDKitMolCreate  # type: ignore
            from rdkit import Chem
            pmol = PDBQTMolecule.from_file(str(poses_pdbqt), skip_typing=True)
            rdmols = RDKitMolCreate.from_pdbqt_mol(pmol)
            with Chem.SDWriter(str(out_sdf)) as w:
                for m in rdmols:
                    if m is not None:
                        w.write(m)
            return out_sdf
        except Exception as e:
            self.log.debug(f"Meeko pose->SDF failed: {e}")

        # Fallback: openbabel
        try:
            from openbabel import pybel
            with open(out_sdf, "w") as fh:
                for mol in pybel.readfile("pdbqt", str(poses_pdbqt)):
                    fh.write(mol.write("sdf"))
            return out_sdf
        except Exception as e:
            self.log.warning(f"Could not convert poses to SDF: {e}")
            return poses_pdbqt

    def make_complex(self, receptor_pdb: Path, ligand_pdbqt: Path) -> Path:
        """Concatenate receptor + best ligand pose into a single PDB."""
        out = self.work / "09_complex.pdb"
        # Extract first MODEL from ligand PDBQT
        lig_atoms: List[str] = []
        in_model = False
        first_model_done = False
        for line in ligand_pdbqt.read_text().splitlines():
            if line.startswith("MODEL"):
                if first_model_done:
                    break
                in_model = True
                continue
            if line.startswith("ENDMDL"):
                first_model_done = True
                in_model = False
                continue
            if in_model and line.startswith(("ATOM", "HETATM")):
                # convert PDBQT ATOM -> PDB HETATM, strip Vina columns >66
                pdb_line = "HETATM" + line[6:66]
                # Force residue name LIG
                pdb_line = pdb_line[:17] + "LIG" + pdb_line[20:]
                if len(pdb_line) < 80:
                    pdb_line = pdb_line + " " * (80 - len(pdb_line))
                lig_atoms.append(pdb_line)

        rec_lines = []
        for line in receptor_pdb.read_text().splitlines():
            if line.startswith(("ATOM", "HETATM", "TER")):
                rec_lines.append(line)

        with open(out, "w") as fh:
            for ln in rec_lines:
                fh.write(ln + "\n")
            fh.write("TER\n")
            for ln in lig_atoms:
                fh.write(ln + "\n")
            fh.write("END\n")
        return out


class InteractionAnalyzer:
    """Lightweight distance-based H-bond / hydrophobic detection."""

    HBOND_DONORS    = {"N", "O", "S"}    # very approximate
    HBOND_ACCEPTORS = {"N", "O", "F", "S"}
    HBOND_MAX = 3.5   # Å
    HYDROPHOBIC_MAX = 4.5
    HYDROPHOBIC_ATOMS = {"C"}

    def __init__(self, log: logging.Logger):
        self.log = log

    def analyze(self, complex_pdb: Path) -> Dict[str, Any]:
        prot_atoms: List[Tuple[str, str, int, np.ndarray]] = []
        lig_atoms:  List[Tuple[str, np.ndarray]] = []
        for line in complex_pdb.read_text().splitlines():
            if not line.startswith(("ATOM", "HETATM")):
                continue
            try:
                elem = line[76:78].strip() or line[12:14].strip()[0]
                resname = line[17:20].strip()
                resnum = int(line[22:26])
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                coord = np.array([x, y, z])
            except (ValueError, IndexError):
                continue
            if resname == "LIG":
                lig_atoms.append((elem.upper(), coord))
            else:
                if elem.upper() == "H":
                    continue
                prot_atoms.append((elem.upper(), resname, resnum, coord))

        hbonds: List[Dict[str, Any]] = []
        hydro:  List[Dict[str, Any]] = []

        for (lel, lc) in lig_atoms:
            for (pel, pres, pnum, pc) in prot_atoms:
                d = float(np.linalg.norm(lc - pc))
                if (lel in self.HBOND_DONORS and pel in self.HBOND_ACCEPTORS) or \
                   (lel in self.HBOND_ACCEPTORS and pel in self.HBOND_DONORS):
                    if d <= self.HBOND_MAX:
                        hbonds.append({
                            "residue": f"{pres}{pnum}",
                            "protein_atom": pel,
                            "ligand_atom": lel,
                            "distance": round(d, 2),
                        })
                if lel in self.HYDROPHOBIC_ATOMS and pel in self.HYDROPHOBIC_ATOMS and d <= self.HYDROPHOBIC_MAX:
                    hydro.append({
                        "residue": f"{pres}{pnum}",
                        "distance": round(d, 2),
                    })

        # de-duplicate hydrophobic by residue (keep closest)
        seen: Dict[str, Dict[str, Any]] = {}
        for h in hydro:
            r = h["residue"]
            if r not in seen or h["distance"] < seen[r]["distance"]:
                seen[r] = h
        hydro = list(seen.values())

        # interaction residues (unique)
        residues = sorted({h["residue"] for h in hbonds} | {h["residue"] for h in hydro})

        self.log.info(f"Interactions: {len(hbonds)} H-bond candidates, {len(hydro)} hydrophobic contacts")
        return {
            "hbonds": hbonds[:30],
            "hydrophobic": hydro[:30],
            "interaction_residues": residues,
        }


# ============================================================================
#  STAGE 5 : VISUALIZATION
# ============================================================================

class Visualizer:
    """
    Render publication-quality images using PyMOL (open source).

    Produces THREE images per job:
      * best_pose.png    - hero shot of the best-scoring pose with H-bonds &
                           hydrophobic contacts highlighted, labels, sticks.
      * interactions.png - tight close-up of the interaction network.
      * overview.png     - whole-protein cartoon + ligand context.
    """

    def __init__(self, log: logging.Logger, work_dir: Path):
        self.log = log
        self.work = work_dir

    def render(
        self,
        complex_pdb: Path,
        interaction_residues: List[str],
        hbonds: List[Dict[str, Any]],
        hydrophobic: List[Dict[str, Any]],
    ) -> Dict[str, Optional[Path]]:
        out_best     = self.work / "best_pose.png"
        out_inter    = self.work / "interactions.png"
        out_over     = self.work / "overview.png"
        results: Dict[str, Optional[Path]] = {
            "best_pose": None, "interactions": None, "overview": None,
        }

        try:
            import pymol  # type: ignore
            from pymol import cmd
        except Exception:
            self.log.warning("PyMOL not available; skipping image rendering.")
            return results

        try:
            pymol.finish_launching(["pymol", "-qc"])     # quiet, no GUI
            cmd.reinitialize()
            cmd.load(str(complex_pdb), "cplx")

            # ---- shared scene setup ------------------------------------------------
            cmd.bg_color("white")
            cmd.hide("everything")
            cmd.show("cartoon", "polymer")
            cmd.color("gray80", "polymer")
            cmd.set("cartoon_transparency", 0.55)

            # Ligand: prominent magenta sticks with element-coloured heteroatoms
            cmd.show("sticks", "resn LIG")
            cmd.color("magenta", "resn LIG and elem C")
            cmd.util.cnc("resn LIG")
            cmd.set("stick_radius", 0.20, "resn LIG")

            # Pocket residues: cyan sticks, side chains only
            pocket_sel_parts: List[str] = []
            for r in interaction_residues or []:
                m = re.match(r"([A-Z]{3})(\d+)", r)
                if m:
                    pocket_sel_parts.append(f"(resn {m.group(1)} and resi {m.group(2)})")
            if pocket_sel_parts:
                cmd.select("pocket", " or ".join(pocket_sel_parts))
                cmd.show("sticks", "pocket and not (name C+N+O+CA+H*)")
                cmd.color("palecyan", "pocket and elem C")
                cmd.util.cnc("pocket")
                cmd.set("stick_radius", 0.14, "pocket")

                # Inject a 3-letter→1-letter map into PyMOL's `stored` namespace,
                # which IS visible inside cmd.label() expressions (the default
                # namespace there does NOT include arbitrary Python globals).
                try:
                    from pymol import stored                       # type: ignore
                    stored.aa3to1 = {
                        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
                        "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
                        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
                        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
                        # common non-standard
                        "MSE": "M", "SEC": "U", "PYL": "O", "HID": "H", "HIE": "H",
                        "HIP": "H", "CYX": "C", "CYM": "C", "ASH": "D", "GLH": "E",
                        "LYN": "K",
                    }
                    cmd.label(
                        "pocket and name CA",
                        "'%s%s' % (stored.aa3to1.get(resn, resn), resi)",
                    )
                    cmd.set("label_color", "black")
                    cmd.set("label_size", 14)
                    cmd.set("label_font_id", 7)
                    cmd.set("label_outline_color", "white")
                except Exception as _le:
                    # Labels are nice-to-have; never let them kill the render
                    self.log.debug(f"Residue labeling skipped: {_le}")

            # ---- H-bonds: yellow dashes with distance labels -----------------------
            n_hb_drawn = 0
            for i, h in enumerate(hbonds[:15]):
                m = re.match(r"([A-Z]{3})(\d+)", h["residue"])
                if not m:
                    continue
                resn, resi = m.group(1), m.group(2)
                # Pair the named atom-pair; if elements are ambiguous, use any matching elem
                lel = h.get("ligand_atom", "")
                pel = h.get("protein_atom", "")
                sel_lig  = f"resn LIG and elem {lel}" if lel else "resn LIG"
                sel_prot = f"resn {resn} and resi {resi} and elem {pel}" if pel else \
                           f"resn {resn} and resi {resi}"
                obj = f"hb_{i:02d}"
                try:
                    # mode=0 = simple distance; we filter atoms ourselves
                    cmd.distance(obj, sel_lig, sel_prot, cutoff=h.get("distance", 3.5) + 0.05, mode=0)
                    cmd.color("yellow", obj)
                    cmd.set("dash_width", 3.5, obj)
                    cmd.set("dash_radius", 0.06, obj)
                    cmd.set("dash_gap", 0.30, obj)
                    cmd.set("dash_length", 0.30, obj)
                    n_hb_drawn += 1
                except Exception:
                    pass
            cmd.set("dash_color", "yellow", "hb_*")

            # ---- Hydrophobic contacts: purple/orange dashes ------------------------
            n_hy_drawn = 0
            for i, h in enumerate(hydrophobic[:15]):
                m = re.match(r"([A-Z]{3})(\d+)", h["residue"])
                if not m:
                    continue
                resn, resi = m.group(1), m.group(2)
                sel_lig  = "resn LIG and elem C"
                sel_prot = f"resn {resn} and resi {resi} and elem C and sidechain"
                obj = f"hy_{i:02d}"
                try:
                    cmd.distance(obj, sel_lig, sel_prot, cutoff=h.get("distance", 4.5) + 0.05, mode=0)
                    cmd.color("deepsalmon", obj)   # orange-pink
                    cmd.set("dash_width", 2.2, obj)
                    cmd.set("dash_radius", 0.04, obj)
                    cmd.set("dash_gap", 0.50, obj)
                    cmd.set("dash_length", 0.20, obj)
                    cmd.hide("labels", obj)        # too noisy when many drawn
                    n_hy_drawn += 1
                except Exception:
                    pass

            # ---- Global ray-tracing quality settings ------------------------------
            cmd.set("ray_opaque_background", 1)
            cmd.set("ray_shadows", 0)
            cmd.set("antialias", 2)
            cmd.set("ambient", 0.30)
            cmd.set("specular", 0.25)
            cmd.set("ray_trace_mode", 1)           # outlines
            cmd.set("ray_trace_color", "black")

            # =====================================================================
            # IMAGE 1 -- best_pose.png : hero shot of the pose
            # =====================================================================
            cmd.orient("resn LIG")
            cmd.zoom("resn LIG", buffer=6)
            cmd.ray(1600, 1200)
            cmd.png(str(out_best), dpi=300)
            results["best_pose"] = out_best
            self.log.info(f"Saved: {out_best.name}  (H-bonds drawn: {n_hb_drawn}, hydrophobic: {n_hy_drawn})")

            # =====================================================================
            # IMAGE 2 -- interactions.png : tight close-up
            # =====================================================================
            cmd.zoom("resn LIG", buffer=3)
            cmd.set("cartoon_transparency", 0.75)
            cmd.ray(1600, 1200)
            cmd.png(str(out_inter), dpi=300)
            results["interactions"] = out_inter
            self.log.info(f"Saved: {out_inter.name}")

            # =====================================================================
            # IMAGE 3 -- overview.png : whole protein context
            # =====================================================================
            cmd.hide("labels")
            cmd.set("cartoon_transparency", 0.20)
            cmd.color("lightblue", "polymer")
            cmd.zoom("polymer", buffer=2)
            cmd.ray(1600, 1200)
            cmd.png(str(out_over), dpi=200)
            results["overview"] = out_over
            self.log.info(f"Saved: {out_over.name}")

            cmd.delete("all")
            return results

        except Exception as e:
            self.log.warning(f"PyMOL render failed: {e}")
            return results


# ============================================================================
#  REPORT GENERATION
# ============================================================================

class ReportGenerator:
    """Single PDF summary per run."""

    def __init__(self, log: logging.Logger, work_dir: Path):
        self.log = log
        self.work = work_dir

    def write(self, info: Dict[str, Any]) -> Path:
        out = self.work / "report.pdf"
        try:
            from reportlab.lib.pagesizes import A4
            from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
            from reportlab.lib.units import cm
            from reportlab.lib import colors
            from reportlab.platypus import (
                SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, PageBreak
            )
        except ImportError:
            self.log.warning("reportlab not installed; writing text report instead.")
            return self._text_report(info)

        styles = getSampleStyleSheet()
        title_style = ParagraphStyle(
            "title", parent=styles["Heading1"], fontSize=18, textColor=colors.HexColor("#1f4e79"),
            spaceAfter=12, alignment=1,
        )
        h2 = ParagraphStyle(
            "h2", parent=styles["Heading2"], fontSize=13, textColor=colors.HexColor("#1f4e79"),
            spaceBefore=10, spaceAfter=6,
        )
        body = styles["BodyText"]

        doc = SimpleDocTemplate(str(out), pagesize=A4,
                                leftMargin=2*cm, rightMargin=2*cm,
                                topMargin=1.5*cm, bottomMargin=1.5*cm)
        story = []
        story.append(Paragraph("Protein-Ligand Docking Report", title_style))
        story.append(Paragraph(f"<i>Generated: {datetime.now():%Y-%m-%d %H:%M:%S}</i>", body))
        story.append(Spacer(1, 0.4*cm))

        # ---- Protein info ----
        story.append(Paragraph("Protein Information", h2))
        prot_tbl = [
            ["Input", str(info.get("protein_input", ""))],
            ["Resolved type", info.get("protein_kind", "")],
            ["Raw file", info.get("protein_raw", "")],
            ["Prepared receptor (PDBQT)", info.get("receptor_pdbqt", "")],
        ]
        if info.get("cocrystal_resname"):
            prot_tbl.append(["Co-crystal ligand", info["cocrystal_resname"]])
        story.append(self._table(prot_tbl))
        story.append(Spacer(1, 0.3*cm))

        # ---- Ligand info ----
        story.append(Paragraph("Ligand Information", h2))
        story.append(self._table([
            ["Input", str(info.get("ligand_input", ""))[:80]],
            ["Resolved type", info.get("ligand_kind", "")],
            ["3D SDF",      info.get("ligand_sdf", "")],
            ["PDBQT",       info.get("ligand_pdbqt", "")],
        ]))
        story.append(Spacer(1, 0.3*cm))

        # ---- Docking parameters ----
        story.append(Paragraph("Docking Setup", h2))
        c = info.get("box_center", (0, 0, 0))
        s = info.get("box_size",   (0, 0, 0))
        story.append(self._table([
            ["Engine",          info.get("engine", "")],
            ["Binding site",    info.get("site_source", "")],
            ["Grid center (Å)", f"({c[0]:.2f}, {c[1]:.2f}, {c[2]:.2f})"],
            ["Grid size (Å)",   f"({s[0]:.1f}, {s[1]:.1f}, {s[2]:.1f})"],
            ["Exhaustiveness",  info.get("exhaustiveness", "")],
            ["Poses requested", info.get("num_poses", "")],
        ]))
        story.append(Spacer(1, 0.3*cm))

        # ---- Docking results ----
        story.append(Paragraph("Docking Scores (kcal/mol)", h2))
        aff = info.get("affinities", []) or []
        if aff:
            rows = [["Pose", "Affinity (kcal/mol)"]]
            for i, a in enumerate(aff, 1):
                rows.append([str(i), f"{a:.3f}"])
            story.append(self._table(rows, header=True))
            best = info.get("best_affinity")
            if best is not None:
                story.append(Spacer(1, 0.2*cm))
                story.append(Paragraph(
                    f"<b>Best binding affinity: {best:.3f} kcal/mol</b>", body
                ))
        else:
            story.append(Paragraph("No affinities were recovered.", body))
        story.append(Spacer(1, 0.3*cm))

        # ---- Interactions ----
        story.append(Paragraph("Key Interactions", h2))
        hbonds = info.get("hbonds", []) or []
        if hbonds:
            rows = [["Residue", "Protein atom", "Ligand atom", "Distance (Å)"]]
            for h in hbonds[:15]:
                rows.append([h["residue"], h["protein_atom"], h["ligand_atom"], f"{h['distance']:.2f}"])
            story.append(Paragraph("<b>H-bond candidates:</b>", body))
            story.append(self._table(rows, header=True))
            story.append(Spacer(1, 0.2*cm))
        hydro = info.get("hydrophobic", []) or []
        if hydro:
            rows = [["Residue", "Distance (Å)"]]
            for h in hydro[:15]:
                rows.append([h["residue"], f"{h['distance']:.2f}"])
            story.append(Paragraph("<b>Hydrophobic contacts:</b>", body))
            story.append(self._table(rows, header=True))

        # ---- Images ----
        for img_key, caption in (
            ("image_best_pose",    "Best Docked Pose  (H-bonds = yellow dashes, hydrophobic = pink dashes)"),
            ("image_interactions", "Close-up of the Interaction Network"),
            ("image_overview",     "Whole-protein Context"),
        ):
            img = info.get(img_key) or (info.get("image") if img_key == "image_best_pose" else "")
            if img and Path(img).exists():
                story.append(PageBreak())
                story.append(Paragraph(caption, h2))
                try:
                    im = Image(str(img), width=16*cm, height=12*cm, kind="proportional")
                    story.append(im)
                except Exception:
                    pass

        doc.build(story)
        self.log.info(f"PDF report written: {out}")
        return out

    @staticmethod
    def _table(data, header=False):
        from reportlab.platypus import Table, TableStyle
        from reportlab.lib import colors
        from reportlab.lib.units import cm
        t = Table(data, hAlign="LEFT")
        style = [
            ("FONT", (0, 0), (-1, -1), "Helvetica", 9),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("GRID", (0, 0), (-1, -1), 0.3, colors.HexColor("#aaaaaa")),
            ("LEFTPADDING", (0, 0), (-1, -1), 4),
            ("RIGHTPADDING", (0, 0), (-1, -1), 4),
        ]
        if header:
            style += [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1f4e79")),
                ("TEXTCOLOR",  (0, 0), (-1, 0), colors.white),
                ("FONT", (0, 0), (-1, 0), "Helvetica-Bold", 9),
            ]
        else:
            style.append(("BACKGROUND", (0, 0), (0, -1), colors.HexColor("#eef2f6")))
            style.append(("FONT", (0, 0), (0, -1), "Helvetica-Bold", 9))
        t.setStyle(TableStyle(style))
        return t

    def _text_report(self, info):
        out = self.work / "report.txt"
        with open(out, "w") as fh:
            fh.write("DOCKING REPORT\n==============\n\n")
            for k, v in info.items():
                fh.write(f"{k}: {v}\n")
        return out


# ============================================================================
#  PIPELINE ORCHESTRATOR
# ============================================================================

@dataclass
class JobSpec:
    protein: str
    ligand: str
    name: str = ""
    center: Optional[Tuple[float, float, float]] = None
    size: Optional[Tuple[float, float, float]] = None

    def auto_name(self) -> str:
        if self.name:
            return safe_name(self.name)
        p = Path(self.protein).stem if Path(self.protein).exists() else safe_name(self.protein)
        l = Path(self.ligand).stem  if Path(self.ligand).exists()  else "ligand"
        return safe_name(f"{p}__{l}")


@dataclass
class JobResult:
    name: str
    success: bool
    best_affinity: Optional[float] = None
    affinities: List[float] = field(default_factory=list)
    site_source: str = ""
    n_hbonds: int = 0
    n_hydrophobic: int = 0
    interaction_residues: List[str] = field(default_factory=list)
    output_dir: str = ""
    report: str = ""
    error: str = ""


def _run_one_job(args_blob: Dict[str, Any]) -> Dict[str, Any]:
    """Top-level worker function (must be picklable for multiprocessing)."""
    spec       = JobSpec(**args_blob["spec"])
    out_root   = Path(args_blob["out_root"])
    config     = args_blob["config"]

    job_dir = out_root / f"run_{spec.auto_name()}"
    job_dir.mkdir(parents=True, exist_ok=True)
    log = setup_logger(job_dir, name=spec.auto_name(), verbose=config["verbose"])
    log.info("=" * 72)
    log.info(f"JOB: {spec.auto_name()}")
    log.info(f"  protein = {spec.protein}")
    log.info(f"  ligand  = {spec.ligand}")
    log.info("=" * 72)

    info: Dict[str, Any] = {
        "protein_input": spec.protein,
        "ligand_input":  spec.ligand,
        "name":          spec.auto_name(),
    }
    result = JobResult(name=spec.auto_name(), success=False, output_dir=str(job_dir))

    try:
        # ---- STAGE 1+2 ----
        prep = ProteinPreparer(log, job_dir, minimize=config["minimize"])
        raw  = prep.resolve(spec.protein)
        info["protein_kind"]      = detect_protein_input(spec.protein)
        info["protein_raw"]       = str(raw)
        info["cocrystal_resname"] = prep.cocrystal_resname or ""
        fixed     = prep.fix(raw)
        minimized = prep.minimize_structure(fixed)
        receptor  = prep.to_pdbqt(minimized)
        info["receptor_pdbqt"] = str(receptor)

        # ---- STAGE 3 ----
        lig = LigandPreparer(log, job_dir)
        sdf, lpdbqt = lig.prepare(spec.ligand, name=spec.auto_name())
        info["ligand_kind"]  = detect_ligand_input(spec.ligand)
        info["ligand_sdf"]   = str(sdf)
        info["ligand_pdbqt"] = str(lpdbqt)

        # ---- STAGE 4 ----
        site = BindingSiteDetector(log)
        center, size, site_src = site.determine(
            receptor_pdb=minimized,
            user_center=spec.center,
            user_size=spec.size,
            cocrystal_center=prep.cocrystal_center,
        )
        engine = DockingEngine(
            log, job_dir,
            exhaustiveness=config["exhaustiveness"],
            num_poses=config["num_poses"],
            cpu=config["cpu_per_job"],
        )
        dock = engine.run(receptor, lpdbqt, center, size)
        info.update({
            "engine":         dock.engine,
            "exhaustiveness": config["exhaustiveness"],
            "num_poses":      config["num_poses"],
            "box_center":     center,
            "box_size":       size,
            "site_source":    site_src,
            "affinities":     dock.affinities,
            "best_affinity":  dock.best_affinity,
        })

        # ---- Pose post-process + interactions ----
        post = PoseProcessor(log, job_dir)
        _    = post.split_and_convert(dock.poses_pdbqt)
        cplx = post.make_complex(minimized, dock.poses_pdbqt)
        ia   = InteractionAnalyzer(log).analyze(cplx)
        info.update(ia)

        # ---- STAGE 5 ----
        viz = Visualizer(log, job_dir)
        imgs = viz.render(
            cplx,
            ia["interaction_residues"],
            ia["hbonds"],
            ia["hydrophobic"],
        )
        info["image"]              = str(imgs.get("best_pose") or "")
        info["image_best_pose"]    = str(imgs.get("best_pose") or "")
        info["image_interactions"] = str(imgs.get("interactions") or "")
        info["image_overview"]     = str(imgs.get("overview") or "")

        # ---- Report ----
        report_path = ReportGenerator(log, job_dir).write(info)

        # JSON dump
        (job_dir / "results.json").write_text(json.dumps(info, indent=2, default=str))

        result.success              = True
        result.best_affinity        = dock.best_affinity
        result.affinities           = dock.affinities
        result.site_source          = site_src
        result.n_hbonds             = len(ia["hbonds"])
        result.n_hydrophobic        = len(ia["hydrophobic"])
        result.interaction_residues = ia["interaction_residues"]
        result.report               = str(report_path)
        log.info(f"JOB SUCCESS - best affinity {dock.best_affinity}")

    except Exception as e:
        tb = traceback.format_exc()
        log.error(f"JOB FAILED: {e}\n{tb}")
        result.error = f"{e}"
        (job_dir / "ERROR.txt").write_text(f"{e}\n\n{tb}")

    return asdict(result)


class Pipeline:
    """High-level orchestrator that handles single / CSV / folder modes."""

    def __init__(self, args: argparse.Namespace, log: logging.Logger, ui: Optional["HackerUI"] = None):
        self.args = args
        self.log  = log
        self.ui   = ui or HackerUI(enabled=False)
        self.out_root = Path(args.output).resolve()
        self.out_root.mkdir(parents=True, exist_ok=True)

    # ---- input parsers -----------------------------------------------------

    def build_jobs(self) -> List[JobSpec]:
        a = self.args
        jobs: List[JobSpec] = []

        if a.csv:
            return self._jobs_from_csv(Path(a.csv))

        if a.protein_dir or a.ligand_dir:
            return self._jobs_from_dirs(a.protein_dir, a.ligand_dir, a.pair_mode)

        if not a.protein or not a.ligand:
            raise SystemExit(
                "ERROR: provide either --protein & --ligand, --csv, or --protein-dir & --ligand-dir"
            )

        center = tuple(a.center) if a.center else None
        size   = tuple(a.size)   if a.size   else None
        jobs.append(JobSpec(
            protein=a.protein, ligand=a.ligand, name=a.name or "",
            center=center, size=size,
        ))
        return jobs

    def _jobs_from_csv(self, csv_path: Path) -> List[JobSpec]:
        if pd is None:
            raise RuntimeError("pandas is required for --csv mode.")
        df = pd.read_csv(csv_path)
        cols = {c.lower().strip(): c for c in df.columns}
        if "protein" not in cols or "ligand" not in cols:
            raise SystemExit("CSV must contain at least 'protein' and 'ligand' columns.")
        jobs: List[JobSpec] = []
        for i, row in df.iterrows():
            name = ""
            if "name" in cols:
                v = row[cols["name"]]
                name = "" if pd.isna(v) else str(v)
            center = None
            size = None
            if all(k in cols for k in ("cx", "cy", "cz")):
                try:
                    center = (float(row[cols["cx"]]), float(row[cols["cy"]]), float(row[cols["cz"]]))
                except (TypeError, ValueError):
                    center = None
            if all(k in cols for k in ("sx", "sy", "sz")):
                try:
                    size = (float(row[cols["sx"]]), float(row[cols["sy"]]), float(row[cols["sz"]]))
                except (TypeError, ValueError):
                    size = None
            jobs.append(JobSpec(
                protein=str(row[cols["protein"]]).strip(),
                ligand=str(row[cols["ligand"]]).strip(),
                name=name, center=center, size=size,
            ))
        self.log.info(f"Loaded {len(jobs)} jobs from CSV {csv_path}")
        return jobs

    def _jobs_from_dirs(self, protein_dir: Optional[str], ligand_dir: Optional[str], pair_mode: str) -> List[JobSpec]:
        if not protein_dir or not ligand_dir:
            raise SystemExit("Folder mode requires BOTH --protein-dir and --ligand-dir")
        prots = sorted([p for p in Path(protein_dir).iterdir() if p.is_file()])
        ligs  = sorted([l for l in Path(ligand_dir).iterdir()  if l.is_file()])
        if not prots or not ligs:
            raise SystemExit("Empty protein or ligand directory.")
        jobs: List[JobSpec] = []

        if pair_mode == "zip":
            n = min(len(prots), len(ligs))
            if len(prots) != len(ligs):
                self.log.warning(f"--pair-mode zip: protein/ligand count mismatch ({len(prots)}/{len(ligs)}); using first {n}")
            for p, l in zip(prots[:n], ligs[:n]):
                jobs.append(JobSpec(protein=str(p), ligand=str(l), name=f"{p.stem}__{l.stem}"))
        else:  # cross
            for p in prots:
                for l in ligs:
                    jobs.append(JobSpec(protein=str(p), ligand=str(l), name=f"{p.stem}__{l.stem}"))
        self.log.info(
            f"Folder mode ({pair_mode}): {len(prots)} proteins x {len(ligs)} ligands -> {len(jobs)} jobs"
        )
        return jobs

    # ---- execution ---------------------------------------------------------

    def run(self) -> int:
        jobs = self.build_jobs()
        n_parallel = max(1, int(self.args.parallel))
        self.log.info(f"Total jobs to run: {len(jobs)}")

        # ── Mission briefing panel (hacker UI) ─────────────────────────────
        self.ui.boot_sequence(str(self.out_root), len(jobs), n_parallel)

        config = {
            "exhaustiveness": self.args.exhaustiveness,
            "num_poses":      self.args.num_poses,
            "cpu_per_job":    self.args.cpu_per_job,
            "minimize":       not self.args.no_minimize,
            "verbose":        self.args.verbose,
        }
        results: List[Dict[str, Any]] = []

        # Build the rich progress bar (or fall back to tqdm)
        use_rich = self.ui.enabled and RICH
        if use_rich:
            progress = Progress(
                SpinnerColumn(style="bright_green", spinner_name="dots12"),
                TextColumn("[bold bright_green]▶[/]"),
                TextColumn("[bold magenta]{task.description}[/]"),
                BarColumn(bar_width=None,
                          complete_style="bright_green",
                          finished_style="bold bright_green",
                          pulse_style="bright_cyan"),
                MofNCompleteColumn(),
                TextColumn("[cyan]•[/]"),
                TimeElapsedColumn(),
                console=console, transient=False, expand=True,
            )
        else:
            progress = None

        # ── Sequential mode ───────────────────────────────────────────────
        if n_parallel == 1 or len(jobs) == 1:
            if progress:
                with progress:
                    task = progress.add_task("DOCKING IN PROGRESS", total=len(jobs))
                    for idx, spec in enumerate(jobs, 1):
                        self.ui.divider()
                        self.ui.job_header(idx, len(jobs), spec.auto_name(),
                                           spec.protein, spec.ligand)
                        self.ui.stage("stage1", f"resolving protein → {spec.protein[:60]}")
                        self.ui.stage("stage2", "PDBFixer + OpenMM minimization")
                        self.ui.stage("stage3", f"preparing ligand → {spec.ligand[:60]}")
                        self.ui.stage("stage4", f"AutoDock Vina (exhaustiveness={config['exhaustiveness']})")
                        self.ui.stage("stage5", "PyMOL render + interaction analysis + PDF report")
                        res = _run_one_job({
                            "spec": asdict(spec),
                            "out_root": str(self.out_root),
                            "config": config,
                        })
                        results.append(res)
                        if res.get("success"):
                            ba = res.get("best_affinity")
                            self.ui.ok(
                                f"{res['name']}  →  ΔG = {ba:.2f} kcal/mol  "
                                f"(H-bonds: {res.get('n_hbonds',0)}, "
                                f"hydrophobic: {res.get('n_hydrophobic',0)})"
                            )
                        else:
                            self.ui.fail(f"{res.get('name','?')}  →  {res.get('error','unknown error')}")
                        progress.update(task, advance=1)
            else:
                for spec in tqdm(jobs, desc="Docking"):
                    res = _run_one_job({
                        "spec": asdict(spec),
                        "out_root": str(self.out_root),
                        "config": config,
                    })
                    results.append(res)

        # ── Parallel mode ─────────────────────────────────────────────────
        else:
            self.log.info(f"Running with {n_parallel} parallel workers")
            if progress:
                with progress:
                    task = progress.add_task(
                        f"DOCKING ({n_parallel} parallel workers)", total=len(jobs))
                    with ProcessPoolExecutor(max_workers=n_parallel) as ex:
                        futs = [ex.submit(_run_one_job, {
                            "spec": asdict(spec),
                            "out_root": str(self.out_root),
                            "config": config,
                        }) for spec in jobs]
                        for f in as_completed(futs):
                            try:
                                res = f.result()
                            except Exception as e:
                                self.log.error(f"Worker died: {e}")
                                res = {"name": "?", "success": False, "error": str(e)}
                            results.append(res)
                            if res.get("success"):
                                ba = res.get("best_affinity")
                                ba_str = f"{ba:.2f}" if isinstance(ba, (int, float)) else "—"
                                self.ui.ok(f"{res.get('name','?')}  →  ΔG = {ba_str}")
                            else:
                                self.ui.fail(f"{res.get('name','?')}  →  {res.get('error','err')[:60]}")
                            progress.update(task, advance=1)
            else:
                with ProcessPoolExecutor(max_workers=n_parallel) as ex:
                    futs = [ex.submit(_run_one_job, {
                        "spec": asdict(spec),
                        "out_root": str(self.out_root),
                        "config": config,
                    }) for spec in jobs]
                    for f in tqdm(as_completed(futs), total=len(futs), desc="Docking"):
                        try:
                            results.append(f.result())
                        except Exception as e:
                            self.log.error(f"Worker died: {e}")
                            results.append({"name": "?", "success": False, "error": str(e)})

        # ── Write artifacts + summary ──────────────────────────────────────
        self._write_summary(results)
        n_ok = sum(1 for r in results if r.get("success"))
        self.log.info(f"DONE. {n_ok}/{len(results)} jobs succeeded.")
        self.ui.summary(results, str(self.out_root))
        return 0 if n_ok == len(results) else 1

    def _write_summary(self, results: List[Dict[str, Any]]) -> None:
        # CSV
        summary = self.out_root / "summary.csv"
        try:
            import csv
            with open(summary, "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow([
                    "name", "success", "best_affinity_kcal_per_mol", "n_poses",
                    "site_source", "n_hbonds", "n_hydrophobic",
                    "interaction_residues", "output_dir", "report", "error",
                ])
                for r in results:
                    w.writerow([
                        r.get("name", ""),
                        r.get("success", False),
                        r.get("best_affinity", ""),
                        len(r.get("affinities", []) or []),
                        r.get("site_source", ""),
                        r.get("n_hbonds", 0),
                        r.get("n_hydrophobic", 0),
                        "; ".join(r.get("interaction_residues", []) or []),
                        r.get("output_dir", ""),
                        r.get("report", ""),
                        r.get("error", ""),
                    ])
            self.log.info(f"Global summary: {summary}")
        except Exception as e:
            self.log.warning(f"Could not write summary.csv: {e}")

        (self.out_root / "summary.json").write_text(json.dumps(results, indent=2, default=str))


# ============================================================================
#  CLI
# ============================================================================

def build_cli() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="docking_pipeline.py",
        description=(
            "Advanced Automated Protein-Ligand Docking Pipeline. "
            "Implements Stages 1-5 from the Medvolt assignment."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples
--------
  # Single run with PDB ID + SMILES
  python docking_pipeline.py --protein 1HSG --ligand "CC(=O)Oc1ccccc1C(=O)O" -o out/

  # UniProt ID (auto AlphaFold) + SDF
  python docking_pipeline.py --protein P00533 --ligand erlotinib.sdf -o out/

  # CSV-driven batch (4 parallel workers)
  python docking_pipeline.py --csv jobs.csv -o out/ --parallel 4

  # Folder x Folder cross-docking
  python docking_pipeline.py --protein-dir proteins/ --ligand-dir ligands/ \\
                             --pair-mode cross --parallel 8 -o out/
""",
    )
    g_in = p.add_argument_group("INPUT (choose one mode)")
    g_in.add_argument("--protein", help="PDB ID | UniProt ID | FASTA file | FASTA sequence | .pdb file")
    g_in.add_argument("--ligand",  help="SMILES | .sdf | .mol2 | .pdb | .mol | .smi file")
    g_in.add_argument("--csv",     help="CSV file with columns: protein,ligand[,name,cx,cy,cz,sx,sy,sz]")
    g_in.add_argument("--protein-dir", help="Directory of protein files for batch mode")
    g_in.add_argument("--ligand-dir",  help="Directory of ligand files for batch mode")
    g_in.add_argument("--pair-mode", choices=["cross", "zip"], default="cross",
                      help="cross = all-vs-all (default); zip = pair files by sorted order")
    g_in.add_argument("--name", default="", help="Custom run name (single mode)")

    g_site = p.add_argument_group("BINDING SITE (optional, otherwise auto)")
    g_site.add_argument("--center", nargs=3, type=float, metavar=("X", "Y", "Z"),
                        help="Grid box center (Å)")
    g_site.add_argument("--size",   nargs=3, type=float, metavar=("SX", "SY", "SZ"),
                        help="Grid box size (Å)")

    g_dock = p.add_argument_group("DOCKING")
    g_dock.add_argument("--exhaustiveness", type=int, default=DEFAULT_EXHAUSTIVENESS,
                        help=f"Vina exhaustiveness (default {DEFAULT_EXHAUSTIVENESS})")
    g_dock.add_argument("--num-poses", type=int, default=DEFAULT_NUM_POSES,
                        help=f"Number of poses (default {DEFAULT_NUM_POSES})")
    g_dock.add_argument("--no-minimize", action="store_true",
                        help="Skip OpenMM energy minimization step")

    g_exec = p.add_argument_group("EXECUTION")
    g_exec.add_argument("-o", "--output", default="docking_results",
                        help="Output directory (default: docking_results)")
    g_exec.add_argument("--parallel", type=int, default=1,
                        help="Number of parallel jobs (default 1)")
    g_exec.add_argument("--cpu-per-job", type=int, default=0,
                        help="CPU threads per Vina job (0 = all available)")
    g_exec.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    g_exec.add_argument("--no-banner", action="store_true",
                        help="Disable the hacker-style banner & rich UI (plain text only)")
    g_exec.add_argument("--version", action="version", version=f"%(prog)s {PIPELINE_VERSION}")

    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_cli()
    args = parser.parse_args(argv)

    out = Path(args.output).resolve()
    out.mkdir(parents=True, exist_ok=True)

    # Build UI before logging so we control std-out
    ui = HackerUI(enabled=not args.no_banner)
    ui.banner()

    log = setup_logger(out, name="pipeline", verbose=args.verbose)
    log.info(f"=== Docking Pipeline v{PIPELINE_VERSION} ===")
    log.info(f"Output dir: {out}")

    try:
        pipe = Pipeline(args, log, ui)
        return pipe.run()
    except KeyboardInterrupt:
        log.error("Interrupted by user.")
        ui.fail("Interrupted by user.")
        return 130
    except SystemExit:
        raise
    except Exception as e:
        log.error(f"FATAL: {e}\n{traceback.format_exc()}")
        ui.fail(f"FATAL: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
