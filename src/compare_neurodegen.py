"""
Comparison figures for the downstream Ca2+ -> CASP3 / neurodegeneration model.

Entry points (imported by ``run_neurodegen.py``):

    render_petri_nets            static side-by-side Petri net diagram
    animate_downstream_petri_nets animated side-by-side Petri nets + cumulative ND curve
    render_comparison            ensemble bar chart + absolute difference
    render_trajectory_comparison static trajectory (CASP3, Apoptosis, cumulative ND)
    render_trajectory_gif        animated version of the trajectory
"""

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.lines import Line2D
from matplotlib.animation import FuncAnimation, PillowWriter
import matplotlib.gridspec as gridspec

from src.simulation import gillespie, sample_trace, cumulative_ca_influx
from src.neurodegeneration import build_downstream_baseline, build_downstream_perturbation
from src.visualization import _place_edge, _box_edge, _BIDIR_RAD

# ---------------------------------------------------------------------------
# Compartment colours / bands for the downstream layout
# ---------------------------------------------------------------------------
COMP_COLORS = {
    "cytosol":   "#caffbf",
    "endosome":  "#ffadad",
}

# Y range has been extended upward to y=7.0 so the subnet block above
# Ca_in (at y ~ 5.3) is fully visible and not clipped.
COMP_BANDS = [
    ("NMDAR input interface",             4.0,  7.0, "#f5f0e8"),
    ("Cytosol / Signalling",             -0.5,  4.0, "#e9fbe0"),
    ("Caspase cascade",                  -4.0, -0.5, "#ffd6a5"),
    ("Apoptosis / Neurodegeneration",    -6.2, -4.0, "#ffadad"),
]

PLACE_R = 0.42

# Subnet block style (hierarchical link marker) — matches src/visualization.py
_SUBNET_FACE    = "#dfe6f0"
_SUBNET_EDGE    = "#34495e"
_SUBNET_SHADOW  = "#c4cfdd"
_SUBNET_W       = 2.25
_SUBNET_H       = 1.45

# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------


def _net_variant_tag(pn) -> str:
    """Return "(baseline)"/"(perturbation)" for a known downstream net, else ""."""
    name = (getattr(pn, "name", "") or "").lower()
    if "perturb" in name:
        return "(perturbation)"
    if "baseline" in name:
        return "(baseline)"
    return ""


def _draw_subnet_block(ax, place, pn=None):
    """Draw a stacked coarse-node block beside *place* signalling a linked subnet.

    When *pn* is a baseline/perturbation net, a "(baseline)"/"(perturbation)"
    line is appended so the block states which variant of the subnet it links to.
    """
    px, py = place.position
    dx, dy = getattr(place, "subnet_offset", (0.0, -2.3))
    bx, by = px + dx, py + dy

    label = place.subnet_link
    tag = _net_variant_tag(pn)
    if tag:
        label = f"{label}\n{tag}"

    n_lines = label.count("\n") + 1
    w = _SUBNET_W
    h = max(_SUBNET_H, 0.34 * n_lines + 0.40)   # grow to fit the line count

    # Dashed connector: place -> block
    ax.add_patch(FancyArrowPatch(
        (px, py + PLACE_R), (bx, by - h / 2),
        arrowstyle="-|>", color=_SUBNET_EDGE, linewidth=1.5,
        mutation_scale=13, linestyle=(0, (4, 2)),
        shrinkA=2, shrinkB=2, zorder=3,
        connectionstyle="arc3,rad=0.0",
    ))

    # Shadow sheet (stacked behind)
    ax.add_patch(FancyBboxPatch(
        (bx - w / 2 + 0.16, by - h / 2 - 0.16), w, h,
        boxstyle="round,pad=0.02,rounding_size=0.08",
        facecolor=_SUBNET_SHADOW, edgecolor=_SUBNET_EDGE,
        linewidth=1.0, alpha=0.95, zorder=4))

    # Body — double outline (Snoopy coarse-node convention)
    ax.add_patch(FancyBboxPatch(
        (bx - w / 2, by - h / 2), w, h,
        boxstyle="round,pad=0.02,rounding_size=0.08",
        facecolor=_SUBNET_FACE, edgecolor=_SUBNET_EDGE,
        linewidth=1.9, zorder=5))
    ax.add_patch(FancyBboxPatch(
        (bx - w / 2 + 0.11, by - h / 2 + 0.11), w - 0.22, h - 0.22,
        boxstyle="round,pad=0.02,rounding_size=0.06",
        facecolor="none", edgecolor=_SUBNET_EDGE,
        linewidth=0.9, zorder=6))

    ax.text(bx, by, label,
            ha="center", va="center", fontsize=8.0, weight="bold",
            color="#1a1a1a", zorder=7)


def _draw_pn(pn, ax, title="", firing=None):
    """Render a downstream Petri net onto *ax*.

    Parameters
    ----------
    pn : PetriNet
    ax : matplotlib.Axes
    title : str
        Optional subplot title.
    firing : str or None
        Transition name to highlight as currently firing (yellow).
    """
    ax.set_aspect("equal")
    ax.set_axis_off()

    all_x = ([p.position[0] for p in pn.places.values()] +
             [t.position[0] for t in pn.transitions.values()])
    x_min, x_max = min(all_x) - 1.4, max(all_x) + 1.8
    y_min, y_max = -6.5, 7.0

    # --- Compartment bands ---
    for label, y0, y1, col in COMP_BANDS:
        ax.add_patch(mpatches.Rectangle(
            (x_min, y0), x_max - x_min, y1 - y0,
            facecolor=col, edgecolor="none", alpha=0.5, zorder=0))
        ax.text(x_max - 0.1, y1 - 0.15, label,
                fontsize=7.5, style="italic", color="#444",
                ha="right", va="top", zorder=1, weight="bold")
    for y in (-0.5, -4.0):
        ax.add_line(Line2D([x_min, x_max], [y, y],
                           color="#7a5a98", lw=0.8, alpha=0.6,
                           zorder=1, ls="--"))
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    p_pos = {n: p.position for n, p in pn.places.items()}

    # --- Arcs ---
    # Box size must match the transition boxes drawn below (0.52 x 0.85).
    _bw, _bh = 0.52, 0.85
    for t_name, t in pn.transitions.items():
        tp = t.position
        is_firing = (firing == t_name)
        in_places = {a.place for a in t.inputs}
        out_places = {a.place for a in t.outputs}
        bidir = in_places & out_places
        for arc in t.inputs:
            pp = p_pos[arc.place]
            colour = "#c0392b" if is_firing else "#444"
            lw     = 2.0 if is_firing else 1.0
            rad    = _BIDIR_RAD if arc.place in bidir else 0.0
            ax.add_patch(FancyArrowPatch(
                _place_edge(pp, tp, PLACE_R), _box_edge(tp, pp, _bw, _bh),
                arrowstyle="-|>", color=colour, lw=lw,
                mutation_scale=10, connectionstyle=f"arc3,rad={rad}",
                shrinkA=0.0, shrinkB=0.0, zorder=2))
            if arc.weight > 1:
                mx, my = (pp[0] + tp[0]) / 2, (pp[1] + tp[1]) / 2
                ax.text(mx, my + 0.1, str(arc.weight),
                        fontsize=7, color="#444", ha="center", va="center",
                        weight="bold")
        for arc in t.outputs:
            pp = p_pos[arc.place]
            colour = "#e67e22" if is_firing else "#4a148c"
            lw     = 2.0 if is_firing else 1.0
            rad    = _BIDIR_RAD if arc.place in bidir else 0.0
            ax.add_patch(FancyArrowPatch(
                _box_edge(tp, pp, _bw, _bh), _place_edge(pp, tp, PLACE_R),
                arrowstyle="-|>", color=colour, lw=lw,
                mutation_scale=10, connectionstyle=f"arc3,rad={rad}",
                shrinkA=0.0, shrinkB=0.0, zorder=2))

    # --- Places (including subnet blocks) ---
    for n, p in pn.places.items():
        col = COMP_COLORS.get(p.compartment, "#e0e0e0")
        ax.add_patch(mpatches.Circle(
            p.position, radius=PLACE_R,
            facecolor=col, edgecolor="#1a1a1a", lw=1.3, zorder=4))
        ax.text(p.position[0], p.position[1], str(p.tokens),
                ha="center", va="center", fontsize=10,
                weight="bold", zorder=6)

        # Label placement respects per-place `label_side`
        side = getattr(p, "label_side", "below")
        px, py = p.position
        gap = PLACE_R + 0.15
        if side == "above":
            lx, ly, ha, va = px, py + gap, "center", "bottom"
        elif side == "right":
            lx, ly, ha, va = px + gap, py, "left", "center"
        elif side == "left":
            lx, ly, ha, va = px - gap, py, "right", "center"
        else:  # below
            lx, ly, ha, va = px, py - gap, "center", "top"
        ax.text(lx, ly, p.label, ha=ha, va=va, fontsize=6.5,
                color="#1a1a1a", zorder=5,
                bbox=dict(facecolor="white", edgecolor="#888", alpha=0.9,
                          pad=1.2, boxstyle="round,pad=0.15,rounding_size=0.1"))

        if getattr(p, "subnet_link", ""):
            _draw_subnet_block(ax, p, pn)

    # --- Transitions ---
    for n, t in pn.transitions.items():
        x, y = t.position
        w, h = 0.52, 0.85
        is_firing = (firing == n)
        face = "#fdd835" if is_firing else "#2c3e50"
        edge = "#c0392b" if is_firing else "#1a1a1a"
        text_colour = "#1a1a1a" if is_firing else "white"
        ax.add_patch(FancyBboxPatch(
            (x - w / 2, y - h / 2), w, h,
            boxstyle="round,pad=0.04",
            facecolor=face, edgecolor=edge, lw=1.2, zorder=5))
        ax.text(x, y + 0.10, t.label, ha="center", va="center",
                fontsize=5.8, color=text_colour, zorder=6, weight="bold")
        ax.text(x, y - 0.22, f"k={t.rate:.2f}", ha="center", va="center",
                fontsize=4.8, color="#aed6f1", zorder=6)

    if title:
        ax.set_title(title, fontsize=11, weight="bold", pad=8)


def _make_legend_handles():
    """Return a list of legend handles for the downstream Petri net diagram."""
    return [
        mpatches.Patch(color="#caffbf", label="Cytosol"),
        mpatches.Patch(color="#ffadad", label="Mitochondria / Endosome"),
        Line2D([0], [0], color="#444",    lw=2, label="Input arc"),
        Line2D([0], [0], color="#4a148c", lw=2, label="Output arc"),
        mpatches.Patch(facecolor="#2c3e50", edgecolor="black",
                       label="Transition"),
        mpatches.Patch(facecolor="#fdd835", edgecolor="#c0392b",
                       label="Firing transition"),
        mpatches.Patch(facecolor=_SUBNET_FACE, edgecolor=_SUBNET_EDGE,
                       linewidth=1.5, label="Attached block (coarse node)"),
    ]


# ---------------------------------------------------------------------------
# Figure 1: static Petri net diagram
# ---------------------------------------------------------------------------
def render_petri_nets(out_path: str, pdf_only: bool = False) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(18, 11), dpi=140)
    fig.suptitle(
        "Ca2+ -> CASP3 / Neurodegeneration: Downstream Petri Net",
        fontsize=14, weight="bold", y=0.99)

    _draw_pn(build_downstream_baseline(),     axes[0], "Baseline (healthy)")
    _draw_pn(build_downstream_perturbation(), axes[1],
             "Perturbation (anti-NMDAR encephalitis)")

    fig.legend(handles=_make_legend_handles(),
               loc="lower center", ncol=4,
               fontsize=8.5, framealpha=0.92, bbox_to_anchor=(0.5, 0.0))

    fig.tight_layout(rect=[0, 0.04, 1, 0.97])
    import os as _os
    if pdf_only:
        base, _ = _os.path.splitext(out_path)
        out_path = base + ".pdf"
        fig.savefig(out_path, facecolor="white", bbox_inches="tight")
    else:
        fig.savefig(out_path, dpi=150, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Figure 2: animated side-by-side Petri nets + cumulative ND curve
# ---------------------------------------------------------------------------
def animate_downstream_petri_nets(
    base, pert, trace_b, trace_p,
    out_gif: str, t_max: float,
    fps: int = 8, n_frames: int = 60,
    title: str = "",
) -> None:
    """Animated GIF: baseline (left) + perturbation (right) downstream Petri
    nets with a synchronised cumulative neurodegeneration curve underneath.

    Token counts, firing highlights, and the cumulative curve all advance
    in lockstep across ``n_frames``.
    """
    # Sample trace indices for keyframes
    idxs = np.linspace(0, len(trace_b) - 1, n_frames, dtype=int)
    idxs = sorted(set(idxs.tolist()))

    # Cumulative neurodegeneration data
    def _nd_cum(trace):
        ts, cum = [0.0], [0.0]
        total = 0
        for t, tr, _ in trace[1:]:
            if tr == "T_Neurodegeneration":
                total += 1
            ts.append(t)
            cum.append(total)
        return np.array(ts), np.array(cum)

    tc_b, nd_b = _nd_cum(trace_b)
    tc_p, nd_p = _nd_cum(trace_p)
    t_grid = np.linspace(0, t_max, n_frames)
    nd_b_g = np.interp(t_grid, tc_b, nd_b)
    nd_p_g = np.interp(t_grid, tc_p, nd_p)

    # Build figure
    fig = plt.figure(figsize=(18, 14), dpi=120)
    gs = fig.add_gridspec(2, 2, height_ratios=[3, 1.2], hspace=0.25,
                          left=0.02, right=0.98, top=0.93, bottom=0.07)

    ax_base = fig.add_subplot(gs[0, 0])
    ax_pert = fig.add_subplot(gs[0, 1])
    ax_cum  = fig.add_subplot(gs[1, :])

    for ax in (ax_base, ax_pert):
        ax.set_axis_off()

    # Extents for each PN (pre-computed)
    def _extent(pn):
        xs = [p.position[0] for p in pn.places.values()] + \
             [t.position[0] for t in pn.transitions.values()]
        return min(xs) - 1.4, max(xs) + 1.8, -6.5, 7.0

    xmb, xMb, ymb, yMb = _extent(base)
    xmp, xMp, ymp, yMp = _extent(pert)

    for ax, (xm, xM, ym, yM) in zip(
        (ax_base, ax_pert), ((xmb, xMb, ymb, yMb), (xmp, xMp, ymp, yMp))
    ):
        ax.set_xlim(xm, xM)
        ax.set_ylim(ym, yM)
        ax.set_aspect("equal")

    # Cumulative curve setup
    y_max_cum = max(nd_b.max(), nd_p.max()) * 1.15 + 1.0
    ax_cum.set_xlim(0, t_max)
    ax_cum.set_ylim(0, y_max_cum)
    ax_cum.set_xlabel("Simulation time (s)", fontsize=11)
    ax_cum.set_ylabel(r"$\Sigma$ Neurodegeneration firings", fontsize=11)
    ax_cum.set_title("Cumulative neurodegeneration", fontsize=12, weight="bold")
    ax_cum.grid(True, alpha=0.25)

    line_b,   = ax_cum.plot([], [], color="#2e7d32", lw=2.2, label="Baseline")
    line_p,   = ax_cum.plot([], [], color="#c62828", lw=2.2, label="Perturbation")
    head_b    = ax_cum.scatter([], [], color="#2e7d32", s=60, zorder=5)
    head_p    = ax_cum.scatter([], [], color="#c62828", s=60, zorder=5)
    ax_cum.legend(loc="upper left", framealpha=0.9, fontsize=10)

    if title:
        fig.suptitle(title, fontsize=14, weight="bold", y=0.97)

    hud = fig.text(0.50, 0.01, "", fontsize=8, family="monospace", ha="center",
                   bbox=dict(facecolor="#fff9c4", edgecolor="#888", pad=4))

    def update(frame_k):
        k = frame_k

        # ---- Baseline Petri net ----
        ax_base.clear()
        i = idxs[min(k, len(idxs) - 1)]
        t_i, tr_i, marking = trace_b[i]
        for n, v in marking.items():
            base.places[n].tokens = v
        _draw_pn(base, ax_base, "Baseline (healthy)",
                 firing=(tr_i if tr_i not in ("<init>", "<end>") else None))

        # ---- Perturbation Petri net ----
        ax_pert.clear()
        tp_idx = min(range(len(trace_p)), key=lambda j: abs(trace_p[j][0] - t_i))
        tp_i, trp_i, marking_p = trace_p[tp_idx]
        for n, v in marking_p.items():
            pert.places[n].tokens = v
        _draw_pn(pert, ax_pert, "Perturbation (anti-NMDAR)",
                 firing=(trp_i if trp_i not in ("<init>", "<end>") else None))

        # ---- Cumulative curve ----
        line_b.set_data(t_grid[:k+1], nd_b_g[:k+1])
        line_p.set_data(t_grid[:k+1], nd_p_g[:k+1])
        head_b.set_offsets([[t_grid[k], nd_b_g[k]]])
        head_p.set_offsets([[t_grid[k], nd_p_g[k]]])

        ratio = (nd_p_g[k] / nd_b_g[k]) * 100 if nd_b_g[k] > 0 else 0.0
        hud.set_text(
            f"t = {t_grid[k]:5.2f}s    "
            f"Baseline cumulative ND: {nd_b_g[k]:5.0f}    "
            f"Perturbation cumulative ND: {nd_p_g[k]:5.0f}    "
            f"Ratio: {ratio:5.1f}%"
        )
        return []

    anim = FuncAnimation(fig, update, frames=len(idxs),
                         interval=1000 / fps, blit=False)
    anim.save(out_gif, writer=PillowWriter(fps=fps), dpi=120)
    plt.close(fig)
    print(f"  Saved: {out_gif}")


# ---------------------------------------------------------------------------
# Figure 3: ensemble comparison (bar chart + fold-change)
# ---------------------------------------------------------------------------

# Ensemble comparison (30 SSA runs): mean final levels and absolute difference (perturbation - baseline) 
# of downstream markers (CHOP, ROS, CASP3, Apoptosis, Neurodegeneration) at t=10s. 
# Demonstrates elevated pro-apoptotic signalling in the anti-NMDAR perturbation model.

MARKERS = ["CHOP", "ROS", "CASP3", "Apoptosis",
           "T_Apoptosis", "T_Neurodegeneration"]
LABELS  = ["CHOP", "ROS", "CASP3\n(active)", "Apoptosis\n(tokens)",
           "Apoptosis\n(firings)", "Neuro-\ndegeneration\n(firings)"]


def _run_ensemble(builder, n_runs: int = 30, t_max: float = 10.0) -> dict:
    results = {m: [] for m in MARKERS}
    for seed in range(n_runs):
        pn    = builder()
        trace = gillespie(pn, t_max=t_max,
                          rng=np.random.default_rng(seed))
        final = trace[-1][2]
        for m in ["CHOP", "ROS", "CASP3", "Apoptosis"]:
            results[m].append(final.get(m, 0))
        results["T_Apoptosis"].append(
            sum(1 for _, t, _ in trace if t == "T_Apoptosis"))
        results["T_Neurodegeneration"].append(
            sum(1 for _, t, _ in trace if t == "T_Neurodegeneration"))
    return results


def render_comparison(out_path: str, n_runs: int = 100,
                      pdf_only: bool = False) -> None:
    print(f"  Running {n_runs}-run ensemble (baseline) ...")
    r_b = _run_ensemble(build_downstream_baseline,     n_runs)
    print(f"  Running {n_runs}-run ensemble (perturbation) ...")
    r_p = _run_ensemble(build_downstream_perturbation, n_runs)

    mb  = np.array([np.mean(r_b[k]) for k in MARKERS])
    mp  = np.array([np.mean(r_p[k]) for k in MARKERS])
    eb  = np.array([np.std(r_b[k])  for k in MARKERS])
    ep  = np.array([np.std(r_p[k])  for k in MARKERS])
    # Absolute difference and propagated SD
    diff    = mp - mb
    diff_sd = np.sqrt(eb**2 + ep**2)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), dpi=130)
    fig.suptitle(
        f"Downstream marker accumulation: Baseline vs. Perturbation "
        f"({n_runs} SSA runs)",
        fontsize=13, weight="bold")

    x = np.arange(len(MARKERS))
    w = 0.34
    ax = axes[0]
    ax.bar(x - w / 2, mb, w, yerr=eb, color="#2e7d32", alpha=0.82,
           label="Baseline",     capsize=4)
    ax.bar(x + w / 2, mp, w, yerr=ep, color="#c62828", alpha=0.82,
           label="Perturbation", capsize=4)
    ax.set_xticks(x)
    ax.set_xticklabels(LABELS, fontsize=8.5)
    ax.set_ylabel("Mean token count / firing count at t = 10 s")
    ax.set_title("Mean ± SD across SSA runs", fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(axis="y", alpha=0.3)

    ax2 = axes[1]
    colors = ["#c62828" if d >= 0 else "#2e7d32" for d in diff]
    ax2.barh(LABELS, diff, xerr=diff_sd, color=colors, alpha=0.80, capsize=4)
    ax2.axvline(0, color="black", lw=1.4, ls="--", label="No change (Δ = 0)")
    ax2.set_xlabel("Change in mean (Perturbation - Baseline)")
    ax2.set_title("Absolute difference: Perturbation vs. Baseline", fontsize=10)
    ax2.legend(fontsize=9)
    ax2.grid(axis="x", alpha=0.3)
    for i, d in enumerate(diff):
        sign = "+" if d >= 0 else ""
        ax2.text(d + (0.02 if d >= 0 else -0.02), i,
                 f"{sign}{d:.2f}",
                 va="center", ha="left" if d >= 0 else "right",
                 fontsize=8.5, color="#222")

    fig.tight_layout()
    if pdf_only:
        import os as _os
        base, _ = _os.path.splitext(out_path)
        out_path = base + ".pdf"
        fig.savefig(out_path, facecolor="white", bbox_inches="tight")
    else:
        fig.savefig(out_path, dpi=150, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


def curves(trace_b, trace_p, t_max, n=600):

    t_grid = np.linspace(0, t_max, n)

    baseline = {
        "ROS": sample_trace(trace_b, t_grid, "ROS"),
        "CHOP": sample_trace(trace_b, t_grid, "CHOP"),
        "CASP3": sample_trace(trace_b, t_grid, "CASP3"),
        "Apoptosis": sample_trace(trace_b, t_grid, "Apoptosis"),
    }

    perturbed = {
        "ROS": sample_trace(trace_p, t_grid, "ROS"),
        "CHOP": sample_trace(trace_p, t_grid, "CHOP"),
        "CASP3": sample_trace(trace_p, t_grid, "CASP3"),
        "Apoptosis": sample_trace(trace_p, t_grid, "Apoptosis"),
    }

    return t_grid, baseline, perturbed


def _cumulative(trace, transition):
    """Cumulative firing count of a transition over time."""
    ts, cum = [0.0], [0.0]
    total = 0
    for t, tr, _ in trace[1:]:
        if tr == transition:
            total += 1
        ts.append(t)
        cum.append(total)
    return np.array(ts), np.array(cum)


_COL_B = "#2e7d32"
_COL_P = "#c62828"


# ---------------------------------------------------------------------------
# Figure 4: static trajectory comparison  (mirrors ca_influx_comparison.png)
# ---------------------------------------------------------------------------
def render_trajectory_comparison(trace_b, trace_p, out_path: str,
                                  t_max: float, pdf_only: bool = False) -> None:
    tg, base_d, pert_d = curves(trace_b, trace_p, t_max, n=600)

    tb, nd_b = _cumulative(trace_b, "T_Neurodegeneration")
    tp, nd_p = _cumulative(trace_p, "T_Neurodegeneration")
    nd_b_g = np.interp(tg, tb, nd_b)
    nd_p_g = np.interp(tg, tp, nd_p)

    panels = [
        ("CASP3",     "CASP3 tokens",    "Instantaneous active CASP3"),
        ("CHOP",      "CHOP tokens",     "Instantaneous CHOP (ER stress)"),
        ("ROS",       "ROS tokens",      "Instantaneous ROS"),
        ("Apoptosis", "Apoptosis tokens","Instantaneous apoptosis commitment"),
    ]

    fig = plt.figure(figsize=(14, 12), dpi=140)
    gs = gridspec.GridSpec(3, 2, figure=fig, hspace=0.45, wspace=0.3)

    fig.suptitle(
        "Neurodegeneration signalling: Baseline vs. anti-NMDAR encephalitis",
        fontsize=14, weight="bold", y=0.99)

    axes_grid = [
        fig.add_subplot(gs[0, 0]),  # CASP3
        fig.add_subplot(gs[0, 1]),  # CHOP
        fig.add_subplot(gs[1, 0]),  # ROS
        fig.add_subplot(gs[1, 1]),  # Apoptosis
    ]
    ax_cumul = fig.add_subplot(gs[2, :])  # cumulative, full width

    for ax, (key, ylabel, title) in zip(axes_grid, panels):
        ax.plot(tg, base_d[key], color=_COL_B, lw=2.0, label="Baseline (healthy)")
        ax.plot(tg, pert_d[key], color=_COL_P, lw=2.0, label="Perturbation (anti-NMDAR)")
        ax.fill_between(tg, base_d[key], alpha=0.18, color=_COL_B)
        ax.fill_between(tg, pert_d[key], alpha=0.18, color=_COL_P)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_title(title, fontsize=11, weight="bold")
        ax.legend(loc="upper left", framealpha=0.9, fontsize=9)
        ax.set_xlabel("Simulation time (s)", fontsize=10)
        ax.grid(True, alpha=0.25)

    # Cumulative neurodegeneration panel
    ax_cumul.plot(tg, nd_b_g, color=_COL_B, lw=2.2, label="Baseline cumulative")
    ax_cumul.plot(tg, nd_p_g, color=_COL_P, lw=2.2, label="Perturbation cumulative")
    ax_cumul.fill_between(tg, nd_b_g, alpha=0.18, color=_COL_B)
    ax_cumul.fill_between(tg, nd_p_g, alpha=0.18, color=_COL_P)
    ax_cumul.set_xlabel("Simulation time (s)", fontsize=11)
    ax_cumul.set_ylabel("Σ Neurodegeneration\nfirings", fontsize=11)
    ax_cumul.set_title("Cumulative neurodegeneration (T_Neurodegeneration firings)",
                       fontsize=12, weight="bold")
    ax_cumul.legend(loc="upper left", framealpha=0.9, fontsize=10)
    ax_cumul.grid(True, alpha=0.25)

    if nd_b_g[-1] > 0:
        ratio = nd_p_g[-1] / nd_b_g[-1]
        ax_cumul.annotate(
            f"Perturbation / Baseline at t = {t_max:.1f}s : "
            f"{nd_p_g[-1]:.0f} / {nd_b_g[-1]:.0f}  =  {ratio*100:.1f} %",
            xy=(tg[-1], nd_p_g[-1]),
            xytext=(tg[-1] * 0.45,
                    nd_p_g[-1] + max(nd_b_g[-1] - nd_p_g[-1], 0.5) * 0.4),
            fontsize=10, color="#222", weight="bold",
            arrowprops=dict(arrowstyle="->", color="#666"),
            bbox=dict(facecolor="#fff9c4", edgecolor="#888", pad=4),
        )

    fig.subplots_adjust(top=0.95, bottom=0.07, left=0.08, right=0.97,  hspace=0.45, wspace=0.3)
    import os as _os
    if pdf_only:
        base, _ = _os.path.splitext(out_path)
        out_path = base + ".pdf"
        fig.savefig(out_path, facecolor="white", bbox_inches="tight")
    else:
        fig.savefig(out_path, dpi=150, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Figure 5: animated trajectory GIF  (mirrors ca_influx_comparison.gif)
# ---------------------------------------------------------------------------
def render_trajectory_gif(trace_b, trace_p, out_gif: str, t_max: float,
                           n_frames: int = 80, fps: int = 12,
                           layout: str = "3x2") -> None:
    """Animated trajectory GIF.

    ``layout="5x1"`` stacks the five panels in a single column.
    ``layout="3x2"`` mirrors the static PNG (``render_trajectory_comparison``):
    CASP3 / CHOP / ROS / Apoptosis in a 2x2 grid with the cumulative
    neurodegeneration panel spanning the full width of the bottom row.
    """
    from matplotlib.animation import FuncAnimation, PillowWriter

    if layout not in ("5x1", "3x2"):
        raise ValueError(f"layout must be '5x1' or '3x2', got {layout!r}")

    tg, base_d, pert_d = curves(trace_b, trace_p, t_max, n=n_frames)

    tb, nd_b = _cumulative(trace_b, "T_Neurodegeneration")
    tp, nd_p = _cumulative(trace_p, "T_Neurodegeneration")
    nd_b_g = np.interp(tg, tb, nd_b)
    nd_p_g = np.interp(tg, tp, nd_p)

    anim_data_b = [base_d["CASP3"], base_d["CHOP"], base_d["ROS"], base_d["Apoptosis"], nd_b_g]
    anim_data_p = [pert_d["CASP3"], pert_d["CHOP"], pert_d["ROS"], pert_d["Apoptosis"], nd_p_g]
    titles  = ["Instantaneous active CASP3",
            "Instantaneous CHOP (ER stress)",
            "Instantaneous ROS",
            "Instantaneous apoptosis commitment",
            "Cumulative neurodegeneration firings"
            ]
    ylabels = ["CASP3 tokens", "CHOP tokens", "ROS tokens","Apoptosis tokens", "Σ Neurodegeneration"]

    if layout == "5x1":
        fig, axes_arr = plt.subplots(5, 1, figsize=(11, 16), dpi=140, sharex=True)
        axes = list(axes_arr)
    else:  # 3x2: 2x2 grid + full-width cumulative row (matches the static PNG)
        fig = plt.figure(figsize=(14, 12), dpi=140)
        gs = gridspec.GridSpec(3, 2, figure=fig, hspace=0.45, wspace=0.3)
        axes = [
            fig.add_subplot(gs[0, 0]),  # CASP3
            fig.add_subplot(gs[0, 1]),  # CHOP
            fig.add_subplot(gs[1, 0]),  # ROS
            fig.add_subplot(gs[1, 1]),  # Apoptosis
            fig.add_subplot(gs[2, :]),  # cumulative, full width
        ]

    fig.suptitle(
        "Neurodegeneration signalling: Baseline vs. anti-NMDAR encephalitis (animated)",
        fontsize=13, weight="bold", y=0.98)

    for ax, title, ylabel, db, dp in zip(
            axes, titles, ylabels, anim_data_b, anim_data_p):
        ax.set_xlim(0, t_max)
        ax.set_ylim(0, max(db.max(), dp.max()) * 1.2 + 0.5)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_title(title, fontsize=11, weight="bold")
        ax.grid(True, alpha=0.25)
        if layout == "3x2":
            ax.set_xlabel("Simulation time (s)", fontsize=10)
    if layout == "5x1":
        axes[2].set_xlabel("Simulation time (s)", fontsize=11)

    lines_b, lines_p, heads_b, heads_p = [], [], [], []
    for i, ax in enumerate(axes):
        lb, = ax.plot([], [], color=_COL_B, lw=2.0,
                      label="Baseline"     if i == 0 else "_")
        lp, = ax.plot([], [], color=_COL_P, lw=2.0,
                      label="Perturbation" if i == 0 else "_")
        hb  = ax.scatter([], [], color=_COL_B, s=40, zorder=5)
        hp  = ax.scatter([], [], color=_COL_P, s=40, zorder=5)
        lines_b.append(lb); lines_p.append(lp)
        heads_b.append(hb); heads_p.append(hp)
    axes[0].legend(loc="upper left", framealpha=0.9, fontsize=9)

    hud = fig.text(0.985, 0.012, "", fontsize=8, ha="right",
                   family="monospace",
                   bbox=dict(facecolor="#fff9c4", edgecolor="#888", pad=4))

    def update(k):
        for i in range(5):
            lines_b[i].set_data(tg[:k+1], anim_data_b[i][:k+1])
            lines_p[i].set_data(tg[:k+1], anim_data_p[i][:k+1])
            heads_b[i].set_offsets([[tg[k], anim_data_b[i][k]]])
            heads_p[i].set_offsets([[tg[k], anim_data_p[i][k]]])
        hud.set_text(
            f"t={tg[k]:5.2f}s  "
            f"CASP3 B={anim_data_b[0][k]:.0f} P={anim_data_p[0][k]:.0f}  "
            f"CHOP B={anim_data_b[1][k]:.0f} P={anim_data_p[1][k]:.0f}  "
            f"ROS B={anim_data_b[2][k]:.0f} P={anim_data_p[2][k]:.0f}  "
            f"Apoptosis B={anim_data_b[3][k]:.0f} P={anim_data_p[3][k]:.0f}  "
            f"Σ Neurodegeneration B={anim_data_b[4][k]:.0f} P={anim_data_p[4][k]:.0f}"
        )
        return lines_b + lines_p

    anim = FuncAnimation(fig, update, frames=n_frames,
                         interval=1000 / fps, blit=False)
    if layout == "5x1":
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    else:  # 3x2: keep the gridspec spacing (tight_layout breaks the spanning row)
        fig.subplots_adjust(top=0.92, bottom=0.06, left=0.07, right=0.97,
                            hspace=0.45, wspace=0.3)
    anim.save(out_gif, writer=PillowWriter(fps=fps), dpi=100)
    plt.close(fig)
    print(f"  Saved: {out_gif}")
