"""
Generates the downstream Ca2+ -> CASP3 / neurodegeneration figures.

Usage (from project root):
    python run_neurodegen.py

Output files written to ./output/Neurodegeneration/:
    neurodegeneration_pn.png:        static side-by-side Petri net diagram
    neurodegeneration_comparison.png:        ensemble bar chart + absolute difference
    neurodegeneration_trajectory.png:        static trajectory (CASP3, Apoptosis, cumulative ND)
    neurodegeneration_trajectory.gif:        animated version of the above (3x2 grid, cumulative spans bottom row)
    neurodegeneration_trajectory_5_1.gif:    animated version, 5x1 single-column layout
    neurodegeneration_downstream_pn.gif:     animated side-by-side Petri nets + cumulative ND
"""

from __future__ import annotations

import os
import sys
from typing import Any, Dict, Optional

_ROOT = os.path.dirname(os.path.abspath(__file__))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np
from src.compare_neurodegen import (
    animate_downstream_petri_nets,
    render_comparison,
    render_petri_nets,
    render_trajectory_comparison,
    render_trajectory_gif,
)
from src.neurodegeneration import build_downstream_baseline, build_downstream_perturbation
from src.simulation import gillespie
from src.charlie_export import export_apnn, export_andl, export_legend

_OUT_BASE = os.path.join(_ROOT, "output", "Neurodegeneration")

# Failed render steps are collected and reported at the end instead of aborting
# the whole run (so one rendering error no longer leaves only the Charlie files).
_FAILED: list = []


def _step(label: str, fn, *args, **kwargs) -> None:
    """Run one render step; on error print a full traceback and keep going."""
    print(label, flush=True)
    try:
        fn(*args, **kwargs)
    except Exception:
        import traceback
        print(f"  !! FAILED: {label}", flush=True)
        traceback.print_exc()
        _FAILED.append(label)


# ---------------------------------------------------------------------------
# Default settings (used when run standalone)
# ---------------------------------------------------------------------------
_DEFAULT_CFG = {
    "fps": 8,
    "n_frames": 60,
    "save_pdf_plots": False,   # False -> save PNG only; True -> save PDF only
    "t_max_neuro": 10.0,
    "seed_baseline": 42,
    "seed_perturbation": 7,
}


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main(cfg: Optional[Dict[str, Any]] = None) -> None:
    if cfg is None:
        cfg = _DEFAULT_CFG

    fps = cfg.get("fps", 8)
    n_frames = cfg.get("n_frames", 60)
    pdf_only = cfg.get("save_pdf_plots", False)
    t_max = cfg.get("t_max_neuro", 10.0)
    seed_b = cfg.get("seed_baseline", 42)
    seed_p = cfg.get("seed_perturbation", 7)

    img_dir = os.path.join(_OUT_BASE, "images")
    anim_dir = os.path.join(_OUT_BASE, "animations")
    charlie_dir = os.path.join(_OUT_BASE, "charlie")
    for d in (img_dir, anim_dir, charlie_dir):
        os.makedirs(d, exist_ok=True)

    # ---- Charlie APNN export (before any simulation mutates tokens) ----
    print("[0/5] Exporting Charlie APNN files (downstream baseline + perturbation) ...")
    for _b, _tag in ((build_downstream_baseline, "baseline"),
                     (build_downstream_perturbation, "perturbation")):
        _m = _b()
        export_apnn(_m, os.path.join(charlie_dir, f"neurodegeneration_{_tag}.apnn"))
        export_andl(_m, os.path.join(charlie_dir, f"neurodegeneration_{_tag}.andl"))
        export_legend(_m, os.path.join(charlie_dir, f"neurodegeneration_{_tag}_legend.tsv"))

    # ---- Pipeline steps ----
    _step("[1/5] Rendering static Petri net diagrams ...",
          render_petri_nets,
          os.path.join(img_dir, "neurodegeneration_pn.png"),
          pdf_only=pdf_only)

    _step("[2/5] Running ensemble simulation and rendering comparison ...",
          render_comparison,
          os.path.join(img_dir, "neurodegeneration_comparison.png"),
          pdf_only=pdf_only)

    print("[3/5] Running single-trace simulations for trajectory plots ...")
    base = build_downstream_baseline()
    pert = build_downstream_perturbation()
    trace_b = gillespie(base, t_max=t_max, rng=np.random.default_rng(seed_b))
    trace_p = gillespie(pert, t_max=t_max, rng=np.random.default_rng(seed_p))
    print(f"  Baseline: {len(trace_b)} events  |  Perturbation: {len(trace_p)} events")

    _step("[3/5] Rendering static trajectory comparison ...",
          render_trajectory_comparison, trace_b, trace_p,
          out_path=os.path.join(img_dir, "neurodegeneration_trajectory.png"),
          t_max=t_max, pdf_only=pdf_only)

    _step("[4/5] Rendering animated trajectory GIF (3x2) ...",
          render_trajectory_gif, trace_b, trace_p,
          out_gif=os.path.join(anim_dir, "neurodegeneration_trajectory.gif"),
          t_max=t_max, n_frames=80, fps=12, layout="3x2")

    _step("[4/5] Rendering animated trajectory GIF (5x1) ...",
          render_trajectory_gif, trace_b, trace_p,
          out_gif=os.path.join(anim_dir, "neurodegeneration_trajectory_5_1.gif"),
          t_max=t_max, n_frames=80, fps=12, layout="5x1")

    print("[5/5] Rendering animated downstream Petri net GIF ...")
    # Re-build fresh models (the trajectory animations mutated token counts)
    base2 = build_downstream_baseline()
    pert2 = build_downstream_perturbation()
    trace_b2 = gillespie(base2, t_max=t_max, rng=np.random.default_rng(seed_b))
    trace_p2 = gillespie(pert2, t_max=t_max, rng=np.random.default_rng(seed_p))
    _step("[5/5] Animating downstream Petri nets ...",
          animate_downstream_petri_nets,
          base2, pert2, trace_b2, trace_p2,
          out_gif=os.path.join(anim_dir, "neurodegeneration_downstream_pn.gif"),
          t_max=t_max, fps=fps, n_frames=n_frames,
          title="Ca\u00b2\u207a \u2192 CASP3 / Neurodegeneration: Downstream Petri Net Token Flow")

    print("\nDone. Output written to ./output/Neurodegeneration/")
    if _FAILED:
        print("\n  !! The following steps failed (see tracebacks above):",
              flush=True)
        for label in _FAILED:
            print(f"     - {label}")
        print(flush=True)
    else:
        print("  All steps completed successfully.\n", flush=True)


if __name__ == "__main__":
    main()
