"""
Ca2+ influx comparison between Baseline and Perturbation.

The Gillespie traces from `simulation.gillespie` are sampled onto a common
time grid and turned into:
  * a static side-by-side PNG comparing instantaneous and cumulative influx
  * an animated GIF in which the two curves grow in real time
"""

from __future__ import annotations

from typing import Dict, List, Tuple

import matplotlib
matplotlib.use("Agg")          # file-only backend: pipeline saves, never shows
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

from .simulation import cumulative_ca_influx, sample_trace


def _curves(trace_b, trace_p, t_max: float, n: int = 600):
    t_grid = np.linspace(0, t_max, n)
    cb = sample_trace(trace_b, t_grid, "Ca_in")
    cp = sample_trace(trace_p, t_grid, "Ca_in")
    tb, cumb = cumulative_ca_influx(trace_b)
    tp, cump = cumulative_ca_influx(trace_p)
    # interp cumulative onto t_grid
    cumb_g = np.interp(t_grid, tb, cumb)
    cump_g = np.interp(t_grid, tp, cump)
    return t_grid, cb, cp, cumb_g, cump_g


def static_comparison(trace_b, trace_p, out_path: str, t_max: float,
                      pdf_only: bool = False) -> None:
    t, cb, cp, cumb, cump = _curves(trace_b, trace_p, t_max)
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), dpi=140, sharex=True)

    ax = axes[0]
    ax.plot(t, cb, color="#2e7d32", lw=2.0, label="Baseline (healthy)")
    ax.plot(t, cp, color="#c62828", lw=2.0, label="Perturbation (anti-NMDAR)")
    ax.fill_between(t, cb, alpha=0.18, color="#2e7d32")
    ax.fill_between(t, cp, alpha=0.18, color="#c62828")
    ax.set_ylabel("Cytosolic Ca²⁺ tokens", fontsize=11)
    ax.set_title("Instantaneous cytosolic Ca²⁺ around the NMDAR",
                 fontsize=12, weight="bold")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper left", framealpha=0.9, fontsize=10)

    ax = axes[1]
    ax.plot(t, cumb, color="#2e7d32", lw=2.2, label="Baseline cumulative influx")
    ax.plot(t, cump, color="#c62828", lw=2.2, label="Perturbation cumulative influx")
    ax.set_xlabel("Simulation time (s)", fontsize=11)
    ax.set_ylabel("Σ Ca²⁺ through open NMDAR", fontsize=11)
    ax.set_title("Cumulative Ca²⁺ influx (T_Ca_in firings)",
                 fontsize=12, weight="bold")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper left", framealpha=0.9, fontsize=10)

    # annotation: ratio at t_end
    if cumb[-1] > 0:
        ratio = cump[-1] / cumb[-1]
        axes[1].annotate(
            f"Perturbation / Baseline at t = {t_max:.1f}s : "
            f"{cump[-1]:.0f} / {cumb[-1]:.0f}  =  {ratio*100:.1f} %",
            xy=(t[-1], cump[-1]), xytext=(t[-1]*0.55, cump[-1] + (cumb[-1]-cump[-1])*0.5),
            fontsize=10, color="#222", weight="bold",
            arrowprops=dict(arrowstyle="->", color="#666"),
            bbox=dict(facecolor="#fff9c4", edgecolor="#888", pad=4),
        )

    fig.suptitle("Ca²⁺ influx: Baseline vs. anti-NMDAR encephalitis",
                 fontsize=14, weight="bold", y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    import os
    if pdf_only:
        base, _ = os.path.splitext(out_path)
        fig.savefig(base + ".pdf", facecolor="white", bbox_inches="tight")
    else:
        fig.savefig(out_path, dpi=150, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def animated_comparison(trace_b, trace_p, out_gif: str, t_max: float,
                        n_frames: int = 80, fps: int = 12) -> None:
    t, cb, cp, cumb, cump = _curves(trace_b, trace_p, t_max, n=n_frames)

    fig, axes = plt.subplots(2, 1, figsize=(11, 7), dpi=120, sharex=True)
    fig.suptitle("Ca²⁺ influx: Baseline vs. anti-NMDAR encephalitis (animated)",
                 fontsize=13, weight="bold", y=0.98)

    ax_top, ax_bot = axes
    for ax in axes:
        ax.grid(True, alpha=0.25)

    ax_top.set_xlim(0, t_max)
    ax_top.set_ylim(0, max(cb.max(), cp.max()) * 1.15 + 0.5)
    ax_top.set_ylabel("Cytosolic Ca²⁺ tokens")
    ax_top.set_title("Instantaneous cytosolic Ca²⁺",
                     fontsize=11, weight="bold")

    ax_bot.set_xlim(0, t_max)
    ax_bot.set_ylim(0, max(cumb.max(), cump.max()) * 1.15 + 1.0)
    ax_bot.set_xlabel("Simulation time (s)")
    ax_bot.set_ylabel("Σ Ca²⁺ influx events")
    ax_bot.set_title("Cumulative Ca²⁺ influx",
                     fontsize=11, weight="bold")

    line_cb,   = ax_top.plot([], [], color="#2e7d32", lw=2.0, label="Baseline")
    line_cp,   = ax_top.plot([], [], color="#c62828", lw=2.0, label="Perturbation")
    line_cumb, = ax_bot.plot([], [], color="#2e7d32", lw=2.2, label="Baseline")
    line_cump, = ax_bot.plot([], [], color="#c62828", lw=2.2, label="Perturbation")
    head_b = ax_top.scatter([], [], color="#2e7d32", s=40, zorder=5)
    head_p = ax_top.scatter([], [], color="#c62828", s=40, zorder=5)
    ax_top.legend(loc="upper left", framealpha=0.9)
    ax_bot.legend(loc="upper left", framealpha=0.9)

    hud = fig.text(0.985, 0.015, "", fontsize=10, ha="right",
                   family="monospace",
                   bbox=dict(facecolor="#fff9c4", edgecolor="#888", pad=4))

    def update(k):
        line_cb.set_data(t[:k+1],   cb[:k+1])
        line_cp.set_data(t[:k+1],   cp[:k+1])
        line_cumb.set_data(t[:k+1], cumb[:k+1])
        line_cump.set_data(t[:k+1], cump[:k+1])
        head_b.set_offsets([[t[k], cb[k]]])
        head_p.set_offsets([[t[k], cp[k]]])
        ratio = (cump[k] / cumb[k]) * 100 if cumb[k] > 0 else 0.0
        hud.set_text(
            f"t={t[k]:5.2f}s  Σ baseline={cumb[k]:5.0f}  "
            f"Σ perturb={cump[k]:5.0f}  ratio={ratio:5.1f}%"
        )
        return line_cb, line_cp, line_cumb, line_cump, head_b, head_p

    anim = FuncAnimation(fig, update, frames=len(t),
                         interval=1000 / fps, blit=False)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    anim.save(out_gif, writer=PillowWriter(fps=fps), dpi=110)
    plt.close(fig)
