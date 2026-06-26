"""
End-to-end driver for the anti-NMDAR encephalitis Petri-Net study.

Produces these artefacts in ./output/NMDAR/:

    Static diagrams (images/ subdir):
        1. NMDAR_baseline_pn.png   / .pdf
        2. NMDAR_perturbation_pn.png / .pdf
        3. ca_influx_comparison.png / .pdf
        4. NMDAR_PN.png / .pdf
    Token-flow animations (animations/ subdir):
        5. NMDAR_baseline_animation.gif
        6. NMDAR_perturbation_animation.gif
        7. NMDAR_baseline_combined.gif
        8. NMDAR_perturbation_combined.gif
        9. ca_influx_comparison.gif
       10. NMDAR_PN.gif
    Snoopy native files (snoopy/ subdir):
       11. baseline.spn     / baseline.cpn
       12. perturbation.spn / perturbation.cpn
    Charlie APNN files (charlie/ subdir):
       13. NMDAR_baseline.apnn     + .andl + _legend.tsv
       14. NMDAR_perturbation.apnn + .andl + _legend.tsv
"""

from __future__ import annotations

import os
import sys
import time
from typing import Any, Dict, Optional

import numpy as np

from src.compare_ca_influx import animated_comparison, static_comparison
from src.models import build_baseline, build_perturbation
from src.simulation import gillespie
from src.snoopy_export import export_pair
from src.charlie_export import export_apnn, export_andl, export_legend
from src.visualization import (
    animate,
    animate_combined,
    animate_nmdar_pn_comparison,
    draw_nmdar_pn_comparison,
    draw_static,
    set_vertical_scale,
)

_ROOT = os.path.dirname(os.path.abspath(__file__))
_OUT_BASE = os.path.join(_ROOT, "output", "NMDAR")

# ---------------------------------------------------------------------------
# Default settings (used when run standalone, i.e. ``python run_NMDAR.py``)
# ---------------------------------------------------------------------------
_DEFAULT_CFG = {
    "fps": 2,
    "n_frames": 60,
    "save_pdf_plots": False,   # False -> save PNG only; True -> save PDF only
    "t_max_nmdar": 8.0,
    "seed_baseline": 42,
    "seed_perturbation": 7,
    "inc_cell_layer_sim": True,  # attach realistic synapse side-panel to NMDAR PN graphs
}


# ---------------------------------------------------------------------------
# Visual layout tuning (affects the rendered images only)
# ---------------------------------------------------------------------------
# Two knobs for stretching the whole diagram apart.  This is purely a visual
# experiment: edit the numbers, re-run, and compare the images to find the
# spacing that reads best.  They are deliberately NOT part of the config dict.
#
#   * plasm_increase_mul_hl -- HORIZONTAL factor.  Every node (all layers) is
#     fanned apart about the rightmost column, so columns stay vertically
#     aligned, the right edge is fixed, and the leftmost nodes move the most
#     (toward the left).
#   * plasm_increase_mul_vl -- VERTICAL factor.  Every node AND every anatomical
#     band is stretched about the diagram mid-line, making the layers thicker.
#
# Set BOTH to 1.0 to return to the original, un-stretched diagram (no other
# edits needed).  Try e.g. 1.3, 1.5, 1.8 ...
plasm_increase_mul_hl = 1.5
plasm_increase_mul_vl = 1.5

# Vertical-spread anchor (diagram mid-height).  Must match the value handed to
# ``set_vertical_scale`` so the bands and the nodes stretch about the same line.
_V_ANCHOR = 3.8


def _apply_plasma_spacing(pn, hl: float, vl: float) -> None:
    """Stretch the whole diagram apart for the static/animated images.

    Horizontal: every node is fanned out by ``hl`` about the rightmost column,
    so columns stay vertically aligned and the right edge is fixed while the
    leftmost nodes move the most (toward the left).
    Vertical: every node is spread by ``vl`` about ``_V_ANCHOR`` -- the same
    transform ``set_vertical_scale`` applies to the compartment bands, so the
    bands grow thicker in lockstep.  The attached subnet block keeps its fixed
    offset from its anchor place, so it just shifts along unscaled.

    Mutates each node's ``position`` in place; a no-op when both factors are 1.
    Call only after the Snoopy/Charlie exports so those keep the grid coords.
    """
    if hl == 1.0 and vl == 1.0:
        return

    nodes = list(pn.places.values()) + list(pn.transitions.values())
    x_anchor = max(n.position[0] for n in nodes)   # rightmost column stays fixed

    for n in nodes:
        x, y = n.position
        x = x_anchor - (x_anchor - x) * hl             # horizontal spread (leftward)
        y = _V_ANCHOR + (y - _V_ANCHOR) * vl           # vertical spread (thicker layers)
        n.position = (x, y)


def _log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# Steps that fail are recorded here and reported in a summary at the end, so a
# single rendering error on one machine no longer silently aborts the whole run
# (previously a crash right after the Charlie/Snoopy export left only those
# files behind, with no clear reason why the images never appeared).
_FAILED: list = []


def _step(label: str, fn, *args, **kwargs) -> None:
    """Run one render step; on error print a full traceback and keep going."""
    _log(label)
    try:
        fn(*args, **kwargs)
    except Exception:
        import traceback
        print(f"  !! FAILED: {label}", flush=True)
        traceback.print_exc()
        _FAILED.append(label)


def main(cfg: Optional[Dict[str, Any]] = None) -> None:
    if cfg is None:
        cfg = _DEFAULT_CFG

    fps = cfg.get("fps", 8)
    n_frames = cfg.get("n_frames", 60)
    pdf_only = cfg.get("save_pdf_plots", False)
    t_max = cfg.get("t_max_nmdar", 8.0)
    seed_b = cfg.get("seed_baseline", 42)
    seed_p = cfg.get("seed_perturbation", 7)
    include_psd95 = cfg.get("include_PSD95_bind", False)
    inc_cell = cfg.get("inc_cell_layer_sim", True)

    img_dir = os.path.join(_OUT_BASE, "images")
    anim_dir = os.path.join(_OUT_BASE, "animations")
    snoopy_dir = os.path.join(_OUT_BASE, "snoopy")
    charlie_dir = os.path.join(_OUT_BASE, "charlie")
    for d in (img_dir, anim_dir, snoopy_dir, charlie_dir):
        os.makedirs(d, exist_ok=True)

    rng = np.random.default_rng(seed=seed_b)

    # ---- Build models ----
    _log("Building Petri-Net models"
         + ("" if include_psd95 else " (PSD-95 bind subbranch hidden)"))
    base = build_baseline(include_psd95_bind=include_psd95)
    pert = build_perturbation(include_psd95_bind=include_psd95)

    # ---- Snoopy native files (.spn + .cpn per model) ----
    _log("Exporting Snoopy files: Baseline (.spn + .cpn)")
    _step("Snoopy baseline", export_pair, base,
          os.path.join(snoopy_dir, "baseline"))
    _log("Exporting Snoopy files: Perturbation (.spn + .cpn)")
    _step("Snoopy perturbation", export_pair, pert,
          os.path.join(snoopy_dir, "perturbation"))

    # ---- Charlie-readable APNN files (structure for Charlie's analyses) ----
    _log("Exporting Charlie APNN files: Baseline + Perturbation (.apnn)")
    for _m, _tag in ((base, "baseline"), (pert, "perturbation")):
        _c = export_apnn(_m, os.path.join(charlie_dir, f"NMDAR_{_tag}.apnn"))
        export_andl(_m, os.path.join(charlie_dir, f"NMDAR_{_tag}.andl"))
        export_legend(_m, os.path.join(charlie_dir, f"NMDAR_{_tag}_legend.tsv"))
        _log(f"  NMDAR_{_tag} (.apnn + .andl): {_c['places']} places, "
             f"{_c['transitions']} transitions, {_c['arcs']} arcs")

    # ---- Visual spacing (images only; exports above keep the grid coords) ----
    set_vertical_scale(plasm_increase_mul_vl, _V_ANCHOR)   # stretch the bands too
    if plasm_increase_mul_hl != 1.0 or plasm_increase_mul_vl != 1.0:
        _log(f"Applying spacing factors H={plasm_increase_mul_hl} "
             f"V={plasm_increase_mul_vl} to images")
        for _m in (base, pert):
            _apply_plasma_spacing(_m, plasm_increase_mul_hl, plasm_increase_mul_vl)

    # ---- Static figures ----
    _step("Rendering Baseline static figure", draw_static, base,
          os.path.join(img_dir, "NMDAR_baseline_pn.png"),
          pdf_only=pdf_only, inc_cell_layer_sim=inc_cell)
    _step("Rendering Perturbation static figure", draw_static, pert,
          os.path.join(img_dir, "NMDAR_perturbation_pn.png"),
          pdf_only=pdf_only, inc_cell_layer_sim=inc_cell)

    # ---- Run Gillespie simulations ----
    _log("Running Gillespie SSA: Baseline")
    trace_b = gillespie(base, t_max=t_max, rng=rng)
    _log(f"  -> {len(trace_b)} events")

    _log("Running Gillespie SSA: Perturbation")
    trace_p = gillespie(pert, t_max=t_max,
                        rng=np.random.default_rng(seed=seed_p))
    _log(f"  -> {len(trace_p)} events")

    # ---- Token-flow GIFs ----
    _step("Animating Baseline token flow", animate, base, trace_b,
          out_gif=os.path.join(anim_dir, "NMDAR_baseline_animation.gif"),
          fps=fps, n_frames=n_frames,
          title="Baseline NMDAR signalling - Petri Net token flow",
          inc_cell_layer_sim=inc_cell)

    _step("Animating Perturbation token flow", animate, pert, trace_p,
          out_gif=os.path.join(anim_dir, "NMDAR_perturbation_animation.gif"),
          fps=fps, n_frames=n_frames,
          title="Anti-NMDAR encephalitis - Petri Net token flow",
          inc_cell_layer_sim=inc_cell)

    # ---- Combined GIFs: Petri net + dynamic curves side-by-side ----
    _step("Animating Baseline combined (Petri net + curves)",
          animate_combined, base, trace_b,
          out_gif=os.path.join(anim_dir, "NMDAR_baseline_combined.gif"),
          t_max=t_max, fps=fps, n_frames=n_frames,
          title="Baseline: Petri Net token flow + Ca\u00b2\u207a curves",
          inc_cell_layer_sim=inc_cell)

    _step("Animating Perturbation combined (Petri net + curves)",
          animate_combined, pert, trace_p,
          out_gif=os.path.join(anim_dir, "NMDAR_perturbation_combined.gif"),
          t_max=t_max, fps=fps, n_frames=n_frames,
          title="Perturbation: Petri Net token flow + Ca\u00b2\u207a curves",
          inc_cell_layer_sim=inc_cell)

    # ---- Ca2+ comparison ----
    _step("Rendering Ca2+ influx comparison (static)",
          static_comparison, trace_b, trace_p,
          out_path=os.path.join(img_dir, "ca_influx_comparison.png"),
          t_max=t_max, pdf_only=pdf_only)

    _step("Rendering Ca2+ influx comparison (animated)",
          animated_comparison, trace_b, trace_p,
          out_gif=os.path.join(anim_dir, "ca_influx_comparison.gif"),
          t_max=t_max, n_frames=80, fps=14)

    # ---- NMDAR side-by-side Petri net + cumulative Ca2+ ----
    _step("Rendering NMDAR PN comparison (static)",
          draw_nmdar_pn_comparison, base, pert, trace_b, trace_p,
          out_path=os.path.join(img_dir, "NMDAR_PN.png"),
          t_max=t_max, pdf_only=pdf_only,
          inc_cell_layer_sim=inc_cell)

    _step("Rendering NMDAR PN comparison (animated)",
          animate_nmdar_pn_comparison, base, pert, trace_b, trace_p,
          out_gif=os.path.join(anim_dir, "NMDAR_PN.gif"),
          t_max=t_max, fps=fps, n_frames=n_frames,
          title="Anti-NMDAR Encephalitis: Petri Net Models and Cumulative Ca\u00b2\u207a Influx",
          inc_cell_layer_sim=inc_cell)

    # ---- Summary ----
    _log("Done. Artefacts written to ./output/NMDAR/")
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
