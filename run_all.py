"""
Orchestrator: runs both the NMDAR and neurodegeneration pipelines sequentially.

Usage:
    python run_all.py

Global settings are declared in ``CFG`` below.  Adjust them before running to
control both pipelines at once.  Per-pipeline defaults in ``run_NMDAR.py`` and
``run_neurodegen.py`` are overridden by the values here.

Artefacts are written to:
    ./output/NMDAR/images/          — static PNG (and optionally PDF) plots
    ./output/NMDAR/animations/      — animated GIFs
    ./output/NMDAR/snoopy/          — Snoopy native files
    ./output/Neurodegeneration/images/     — static plots
    ./output/Neurodegeneration/animations/ — animated GIFs
"""

from __future__ import annotations

import os
import sys

_ROOT = os.path.dirname(os.path.abspath(__file__))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)


# ---------------------------------------------------------------------------
# Global settings — override before running, or import this dict from
# another script to pass the same settings into both pipelines.
# ---------------------------------------------------------------------------
CFG = {
    # ---- Animation quality ----
    "fps": 2,
    "n_frames": 60,

    # ---- Output format (mutually exclusive) ----
    # False → save static figures as PNG only; True → save as PDF only
    "save_pdf_plots": False,

    # ---- Simulation duration (seconds) ----
    "t_max_nmdar": 8.0,
    "t_max_neuro": 10.0,

    # ---- RNG seeds ---- 
    "seed_baseline": 42,    # for the baseline (healthy) condition
    "seed_perturbation": 7, # for the perturbation (disease) condition; keep different from baseline to avoid identical trajectories

    # ---- PSD-95 subbranch (REF51: PDZ anchor bind transition) ----
    # False omits the PSD-95 places and T_PSD95 transition from NMDAR figures.
    "include_PSD95_bind": False,

    # ---- Realistic synapse side-panel on NMDAR Petri net graphs ----
    # True  → attaches a narrow biological synapse cross-section image to the
    #         right of every NMDAR Petri net diagram (static + animation).
    # False → omit the side panel.
    "inc_cell_layer_sim": True,
}


def _make_dirs(root: str) -> None:
    for pathway in ("NMDAR", "Neurodegeneration"):
        for sub in ("images", "animations", "charlie"):
            os.makedirs(os.path.join(root, "output", pathway, sub), exist_ok=True)
    # NMDAR additionally emits Snoopy native files
    os.makedirs(os.path.join(root, "output", "NMDAR", "snoopy"), exist_ok=True)


def main() -> None:
    _make_dirs(_ROOT)

    print("=" * 72)
    print("  Pipeline 1/2: Anti-NMDAR encephalitis Petri Net study")
    print("=" * 72)
    import run_NMDAR as nmdar
    nmdar.main(cfg=CFG)

    print("\n" + "=" * 72)
    print("  Pipeline 2/2: Downstream neurodegeneration pathway")
    print("=" * 72)
    import run_neurodegen as neuro
    neuro.main(cfg=CFG)

    print("\n" + "=" * 72)
    print("  Both pipelines complete.")
    print(f"  NMDAR images              -> {os.path.join(_ROOT, 'output', 'NMDAR', 'images')}")
    print(f"  NMDAR animations          -> {os.path.join(_ROOT, 'output', 'NMDAR', 'animations')}")
    print(f"  NMDAR Snoopy files        -> {os.path.join(_ROOT, 'output', 'NMDAR', 'snoopy')}")
    print(f"  NMDAR Charlie APNN        -> {os.path.join(_ROOT, 'output', 'NMDAR', 'charlie')}")
    print(f"  Neurodegeneration images  -> {os.path.join(_ROOT, 'output', 'Neurodegeneration', 'images')}")
    print(f"  Neurodegeneration anims   -> {os.path.join(_ROOT, 'output', 'Neurodegeneration', 'animations')}")
    print(f"  Neurodegeneration Charlie -> {os.path.join(_ROOT, 'output', 'Neurodegeneration', 'charlie')}")
    print("=" * 72)

if __name__ == "__main__":
    main()
