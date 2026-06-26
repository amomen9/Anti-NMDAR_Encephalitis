"""
Matplotlib renderer for `PetriNet` objects.

Two entry points:
    draw_static(pn, ..., out_path)  -> single annotated PNG
    animate(pn, trace, out_gif)     -> animated GIF showing token flow

Visual conventions
------------------
* Places       : filled circles, colour by `compartment`, token count printed
                 in the centre, name printed below.
* Transitions  : filled rounded rectangles, label centred, mass-action rate
                 annotation in italics on the side.
* Input arcs   : dark grey arrow Place -> Transition
* Output arcs  : indigo arrow Transition -> Place
* Compartment  : pastel rectangles spanning the canvas as anatomical bands.
"""

from __future__ import annotations

import math
import os
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")          # file-only backend: pipeline saves, never shows
import matplotlib.image as mpimg
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.lines import Line2D

from .petri_net import COMPARTMENT_COLORS, PetriNet
from .simulation import cumulative_ca_influx, sample_trace


# ---------------------------------------------------------------------------
# Anatomical compartment bands (y-range, label, colour)
# ---------------------------------------------------------------------------
COMPARTMENT_BANDS = [
    ("Pre-synaptic terminal",                  8.75, 11.0, "#fff9c4"),
    ("Extracellular space (synaptic cleft)",   6.25,  8.75, "#dff1ff"),
    ("Plasma membrane",                        3.75,  6.25, "#ece4f3"),
    ("Post-synaptic density (PSD)",            1.25,  3.75, "#fff1dc"),
    ("Cytosol",                               -1.25,  1.25, "#e9fbe0"),
    ("Endosomal / lysosomal",                 -3.6,  -1.25, "#ffe1e1"),
]


# ---------------------------------------------------------------------------
# Synapse cell-layer side panel
# ---------------------------------------------------------------------------
_CELL_LAYER_IMAGE_PATH: str = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "presentation", "synapse_candidates", "synapse_sidepanel.png",
)

# Fraction-from-TOP of the artwork at which each anatomical layer begins/ends,
# in the SAME top->bottom order as COMPARTMENT_BANDS.  Calibrated by visual
# inspection of synapse_sidepanel.png (the slight discrepancy from the curved
# layer borders in the artwork is intentionally ignored).  Consecutive pairs
# give the [top, bottom] image slice for each band, so each painted layer is
# stretched/compressed onto its matching net band rather than the image being
# cut at the (wrong) proportional positions.  Must hold len(BANDS)+1 entries.
#
#   0.000  top of pre-synaptic vesicles
#   0.180  pre-synaptic terminal | synaptic cleft (blue)
#   0.340  synaptic cleft        | plasma membrane (orange band)
#   0.420  plasma membrane       | post-synaptic density (orange scaffold band)
#   0.505  PSD                   | cytosol  (orange cobblestone ends, green begins)
#   0.785  cytosol               | endosomal/lysosomal (green ends, pink begins)
#   1.000  bottom of endosomal
_CELL_LAYER_IMG_EDGES: List[float] = [0.0, 0.18, 0.34, 0.42, 0.505, 0.785, 1.0]

# Side-panel shape/placement (all figure-fraction unless noted).
_CELL_PANEL_GAP = 0.006     # horizontal gap between the net and the panel
_CELL_PANEL_ASPECT = 6.5    # rendered height:width of the panel.  The source
                            # art is ~4:1; a larger value renders it narrower
                            # AND taller, as requested.


def _add_cell_layer_sidepanel(fig, ax_main) -> None:
    """Place the realistic synapse cross-section image as a narrow side panel
    to the right of *ax_main*.

    The artwork is sliced at its own anatomical layer boundaries
    (``_CELL_LAYER_IMG_EDGES``) and each layer is drawn into the figure-space
    rectangle occupied by its matching compartment band, so every image layer
    ends up exactly side-by-side with -- and the same height as -- its band on
    the net.

    Placement uses the live data->figure transform rather than
    ``ax_main.get_position()``: with ``aspect="equal"`` the diagram is drawn in
    a centred inset of its allocated box, so the allocated box is the wrong
    reference (it made the panel mis-sized and let it overlap the diagram).
    Transforming real data points instead pins the panel to the net's *drawn*
    edge and band heights, so the panel never overlaps and is exactly as tall
    as the anatomical stack.

    If the image file is missing, or the layer table is mis-sized, the call is a
    no-op.
    """
    if not os.path.isfile(_CELL_LAYER_IMAGE_PATH):
        return
    if len(_CELL_LAYER_IMG_EDGES) != len(COMPARTMENT_BANDS) + 1:
        return                                          # mis-calibrated: skip

    img = mpimg.imread(_CELL_LAYER_IMAGE_PATH)          # shape (H, W, 4)
    img_h = img.shape[0]

    # Finalise the transforms so the equal-aspect inset is resolved, then map
    # data coordinates -> figure fraction.
    fig.canvas.draw()
    inv = fig.transFigure.inverted()

    def _data_to_fig_y(y: float) -> float:
        return inv.transform(ax_main.transData.transform((0.0, y)))[1]

    def _data_to_fig_x(x: float) -> float:
        return inv.transform(ax_main.transData.transform((x, 0.0)))[0]

    # The compartment bands span x_min..x_max in data (see _draw_compartments);
    # the panel sits just right of that drawn edge.
    _, x_max = ax_main.get_xlim()
    panel_left = _data_to_fig_x(x_max) + _CELL_PANEL_GAP

    # Panel height = drawn height of the full anatomical stack (top of first band
    # to bottom of last band) -- i.e. exactly the net's height.
    top_fy = _data_to_fig_y(_vy(COMPARTMENT_BANDS[0][2]))    # y1 of first band
    bot_fy = _data_to_fig_y(_vy(COMPARTMENT_BANDS[-1][1]))   # y0 of last band
    panel_h = abs(top_fy - bot_fy)

    # Width from the requested rendered aspect ratio (narrower/taller than the
    # source art).  Convert through the figure's inch dimensions so the ratio is
    # honoured in real space, not in (anisotropic) figure fractions.
    fig_w_in, fig_h_in = fig.get_size_inches()
    panel_w = panel_h * (fig_h_in / fig_w_in) / _CELL_PANEL_ASPECT

    for i, (label, y0, y1, _colour) in enumerate(COMPARTMENT_BANDS):
        # Figure-space rectangle actually occupied by this band on the net.
        fy0 = _data_to_fig_y(_vy(y0))                   # band bottom
        fy1 = _data_to_fig_y(_vy(y1))                   # band top
        strip_bot = min(fy0, fy1)
        strip_h = abs(fy1 - fy0)

        # Image slice at the artwork's own layer boundaries (top-down).
        iy0 = int(round(_CELL_LAYER_IMG_EDGES[i] * img_h))
        iy1 = int(round(_CELL_LAYER_IMG_EDGES[i + 1] * img_h))
        if iy1 <= iy0 or strip_h <= 0:
            continue
        strip = img[iy0:iy1, :, :]                      # full-width strip

        ax_strip = fig.add_axes(
            [panel_left, strip_bot, panel_w, strip_h],
            zorder=15, frameon=False,
        )
        ax_strip.set_axis_off()
        ax_strip.imshow(strip, aspect="auto", interpolation="bilinear")


# ---------------------------------------------------------------------------
# Vertical layout scale (set by the driver to make the anatomical bands thicker)
# ---------------------------------------------------------------------------
# The driver (run_NMDAR.py) can stretch the diagram vertically by setting this
# scale via ``set_vertical_scale``.  It is applied to the compartment bands, the
# membrane bilayer hint, and every view's y-limits, so the bands grow in
# lockstep with the node positions (which the driver scales with the SAME factor
# and anchor).  Identity by default, so all other diagrams are unaffected.
_V_SCALE = 1.0
_V_ANCHOR = 3.8        # mid-height of the default view (-3.7 .. 11.3)


def set_vertical_scale(scale: float = 1.0, anchor: float = 3.8) -> None:
    """Set the active vertical layout scale (see ``_V_SCALE``)."""
    global _V_SCALE, _V_ANCHOR
    _V_SCALE, _V_ANCHOR = scale, anchor


def _vy(y: float) -> float:
    """Map a base y-coordinate through the active vertical layout scale."""
    return _V_ANCHOR + (y - _V_ANCHOR) * _V_SCALE


PLACE_RADIUS = 0.50

# ---------------------------------------------------------------------------
# Arc geometry
# ---------------------------------------------------------------------------
# Transition box size (must match the defaults in `_draw_transition`, used here
# to clip arc endpoints to the box border so arrowheads never hide under it).
_TBOX_W, _TBOX_H = 0.55, 0.95
# Gap (data units) left between a node's border and the arc end, so the
# arrowhead sits clearly outside the place circle / transition box.
_ARC_MARGIN = 0.07
# Curvature for the two arcs of a read-arc / catalyst pair (a place that is both
# pre- and post-place of the same transition). Both arcs of the pair receive the
# SAME rad; with matplotlib's arc3 that bows them to opposite physical sides, so
# they separate instead of overlapping. 0.0 would leave them on top of each other.
_BIDIR_RAD = 0.22


def _unit(dx: float, dy: float) -> Tuple[float, float]:
    d = math.hypot(dx, dy)
    return (0.0, 0.0) if d < 1e-9 else (dx / d, dy / d)


def _place_edge(center, toward, radius: float = PLACE_RADIUS,
                margin: float = _ARC_MARGIN) -> Tuple[float, float]:
    """Point on a place circle (+margin) in the direction of `toward`."""
    ux, uy = _unit(toward[0] - center[0], toward[1] - center[1])
    return (center[0] + (radius + margin) * ux, center[1] + (radius + margin) * uy)


def _box_edge(center, toward, w: float = _TBOX_W, h: float = _TBOX_H,
              margin: float = _ARC_MARGIN) -> Tuple[float, float]:
    """Point on a transition box border (+margin) in the direction of `toward`."""
    dx, dy = toward[0] - center[0], toward[1] - center[1]
    if abs(dx) < 1e-9 and abs(dy) < 1e-9:
        return (center[0], center[1])
    sx = (w / 2.0) / abs(dx) if abs(dx) > 1e-9 else math.inf
    sy = (h / 2.0) / abs(dy) if abs(dy) > 1e-9 else math.inf
    s = min(sx, sy)
    ex, ey = center[0] + s * dx, center[1] + s * dy
    ux, uy = _unit(dx, dy)
    return (ex + margin * ux, ey + margin * uy)

# Style for the "attached subnet" coarse-node block (hierarchical link marker).
_SUBNET_FACE = "#dfe6f0"      # pale slate-blue body
_SUBNET_EDGE = "#34495e"      # dark slate outline
_SUBNET_SHADOW = "#c4cfdd"    # offset "stacked sheets" sheet behind the body
_SUBNET_W = 3.7
_SUBNET_H = 1.25

# Transition-label font. The label is centred inside the box, so a long line
# (e.g. "endocytosis", "degradation") spills past the border at a fixed size.
# We auto-shrink: lines wider than what fits the box are scaled down
# proportionally, with a floor so they stay legible. `_TLABEL_CHARS_PER_W` is
# how many characters fill one data-unit of box width at the base size, so the
# fit stays correct if the box width `w` is changed.
_TLABEL_FS = 6.4              # base size (unchanged for short labels)
_TLABEL_FS_MIN = 4.6          # never shrink below this
_TLABEL_CHARS_PER_W = 14.5    # chars per unit of box width at the base size


# ---------------------------------------------------------------------------
def _draw_compartments(ax, x_min: float, x_max: float) -> None:
    for label, y0, y1, colour in COMPARTMENT_BANDS:
        sy0, sy1 = _vy(y0), _vy(y1)
        ax.add_patch(
            mpatches.Rectangle(
                (x_min, sy0), x_max - x_min, sy1 - sy0,
                facecolor=colour, edgecolor="none", alpha=0.55, zorder=0,
            )
        )
        ax.text(x_max - 0.15, sy1 - 0.18, label,
                fontsize=8.5, style="italic", color="#444",
                ha="right", va="top", zorder=1, weight="bold")
    # plasma-membrane bilayer hint
    for y in (3.75, 6.25):
        sy = _vy(y)
        ax.add_line(Line2D([x_min, x_max], [sy, sy],
                           color="#7a5a98", linewidth=0.8,
                           alpha=0.6, zorder=1, linestyle="--"))


def _arc(ax, p0: Tuple[float, float], p1: Tuple[float, float],
         colour: str, label: str = "", rad: float = 0.0,
         lw: float = 1.2) -> None:
    # Endpoints arrive already clipped to the node borders (see _draw_arcs), so
    # shrink is ~0: the arrowhead tip lands exactly at the supplied p1.
    arrow = FancyArrowPatch(
        p0, p1,
        arrowstyle="-|>",
        color=colour, linewidth=lw, mutation_scale=12,
        connectionstyle=f"arc3,rad={rad}",
        shrinkA=0.0, shrinkB=0.0, zorder=2,
    )
    ax.add_patch(arrow)
    if label:
        mx, my = (p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2
        # Nudge the label onto the bowed side of a curved arc (perpendicular to
        # the chord) so it sits on the line rather than between overlapping pair.
        if abs(rad) > 1e-9:
            ux, uy = _unit(p1[0] - p0[0], p1[1] - p0[1])
            px, py = -uy, ux                      # left-hand normal
            off = rad * math.hypot(p1[0] - p0[0], p1[1] - p0[1]) * 0.5
            mx, my = mx + px * off, my + py * off
        ax.text(mx, my + 0.10, label, fontsize=7,
                color=colour, alpha=0.95, zorder=3,
                ha="center", va="center", weight="bold",
                bbox=dict(facecolor="white", edgecolor="none",
                          pad=0.5, alpha=0.85))


def _net_variant_tag(pn) -> str:
    """Return "(baseline)"/"(perturbation)" for a known NMDAR net, else ""."""
    name = (getattr(pn, "name", "") or "").lower()
    if "perturb" in name:
        return "(perturbation)"
    if "baseline" in name:
        return "(baseline)"
    return ""


def _draw_subnet_block(ax, place, pn=None) -> None:
    """Draw a Snoopy-style 'coarse node' beside `place`.

    Signals that the subnet named by ``place.subnet_link`` is logically attached
    at this place but is rendered in a separate diagram (not expanded inline).
    Rendered as a double-bordered, stacked block joined to the place by a dashed
    connector arrow (token flow continues into the hidden subnet).  When `pn` is
    a baseline/perturbation net, a "(baseline)"/"(perturbation)" line is appended
    so the block states which variant of the subnet it links to.
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

    # Dashed connector: place -> block (logical link, distinct from token arcs).
    ax.add_patch(FancyArrowPatch(
        (px, py - PLACE_RADIUS), (bx, by + h / 2),
        arrowstyle="-|>", color=_SUBNET_EDGE, linewidth=1.5,
        mutation_scale=13, linestyle=(0, (4, 2)),
        shrinkA=2, shrinkB=2, zorder=3,
        connectionstyle="arc3,rad=0.0",
    ))

    # "Stacked sheets" sheet offset behind the body -> reads as a collapsed net.
    ax.add_patch(FancyBboxPatch(
        (bx - w / 2 + 0.16, by - h / 2 - 0.16), w, h,
        boxstyle="round,pad=0.02,rounding_size=0.08",
        facecolor=_SUBNET_SHADOW, edgecolor=_SUBNET_EDGE,
        linewidth=1.0, alpha=0.95, zorder=4))

    # Body with a doubled outline = Snoopy coarse/macro-node convention.
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


def _draw_place(ax, place, pn=None) -> None:
    colour = COMPARTMENT_COLORS.get(place.compartment, "#e0e0e0")
    circ = mpatches.Circle(place.position, radius=PLACE_RADIUS,
                           facecolor=colour, edgecolor="#1a1a1a",
                           linewidth=1.4, zorder=4)
    ax.add_patch(circ)
    ax.text(place.position[0], place.position[1], str(place.tokens),
            ha="center", va="center", fontsize=12, weight="bold",
            color="#1a1a1a", zorder=6)
    # label placed beside the place; default below, but per-place `label_side`
    # can move it (e.g. to "right") when a neighbouring transition would
    # otherwise hide a below-label.
    side = getattr(place, "label_side", "below")
    px, py = place.position
    gap = PLACE_RADIUS + 0.18
    if side == "above":
        lx, ly, ha, va = px, py + gap, "center", "bottom"
    elif side == "right":
        lx, ly, ha, va = px + gap, py, "left", "center"
    elif side == "left":
        lx, ly, ha, va = px - gap, py, "right", "center"
    else:  # below
        lx, ly, ha, va = px, py - gap, "center", "top"
    ax.text(lx, ly, place.label, ha=ha, va=va, fontsize=7.2,
            color="#1a1a1a", zorder=5,
            bbox=dict(facecolor="white", edgecolor="#888",
                      alpha=0.92, pad=1.4,
                      boxstyle="round,pad=0.18,rounding_size=0.10"))
    if getattr(place, "subnet_link", ""):
        _draw_subnet_block(ax, place, pn)


def _draw_transition(ax, transition, *, firing: bool = False,
                     w: float = 0.55, h: float = 0.95) -> None:
    x, y = transition.position
    face = "#fdd835" if firing else "#2c3e50"
    edge = "#c0392b" if firing else "#1a1a1a"
    text_colour = "#1a1a1a" if firing else "white"
    box = FancyBboxPatch(
        (x - w/2, y - h/2), w, h,
        boxstyle="round,pad=0.04,rounding_size=0.10",
        facecolor=face, edgecolor=edge, linewidth=1.4, zorder=5,
    )
    ax.add_patch(box)
    # Auto-fit the label: shrink only when the longest line is wider than the
    # box, so short labels keep the base size and wide words stay inside.
    longest = max((len(s) for s in transition.label.split("\n")), default=1)
    fit_chars = _TLABEL_CHARS_PER_W * w
    label_fs = _TLABEL_FS
    if longest > fit_chars:
        label_fs = max(_TLABEL_FS * fit_chars / longest, _TLABEL_FS_MIN)
    ax.text(x, y, transition.label, ha="center", va="center",
            fontsize=label_fs, color=text_colour, zorder=6, weight="bold")
    # rate annotation - direction controlled per-transition via annotation_side
    if transition.annotation:
        side = getattr(transition, "annotation_side", "below")
        if side == "above":
            ax_pos = (x, y + h / 2 + 0.30)
            ha, va = "center", "bottom"
        elif side == "left":
            ax_pos = (x - w / 2 - 0.10, y)
            ha, va = "right", "center"
        elif side == "right":
            ax_pos = (x + w / 2 + 0.10, y)
            ha, va = "left", "center"
        else:  # below (default)
            ax_pos = (x, y - h / 2 - 0.30)
            ha, va = "center", "top"
        ax.text(ax_pos[0], ax_pos[1], transition.annotation,
                ha=ha, va=va, fontsize=6.2, style="italic",
                color="#222", zorder=5,
                bbox=dict(facecolor="white", edgecolor="#aaa",
                          alpha=0.93, pad=1.3,
                          boxstyle="round,pad=0.18,rounding_size=0.08"))


def _draw_legend(ax, x_min: float, y_min: float) -> None:
    handles = [
        mpatches.Patch(color=c, label=k.capitalize())
        for k, c in COMPARTMENT_COLORS.items()
    ]
    handles += [
        Line2D([0], [0], color="#444",   lw=2, label="Input arc"),
        Line2D([0], [0], color="#4a148c", lw=2, label="Output arc"),
        mpatches.Patch(facecolor="#2c3e50", edgecolor="black",
                       label="Transition"),
        mpatches.Patch(facecolor="#fdd835", edgecolor="#c0392b",
                       label="Firing transition"),
        mpatches.Patch(facecolor=_SUBNET_FACE, edgecolor=_SUBNET_EDGE,
                       linewidth=1.5, label="Attached subnet (coarse node)"),
    ]
    ax.legend(handles=handles, loc="lower left",
              bbox_to_anchor=(0.0, 0.0),
              fontsize=7.5, framealpha=0.92, ncol=2,
              title="Legend", title_fontsize=8.5)


# ---------------------------------------------------------------------------
def _figure_for(pn: PetriNet, figsize: Tuple[float, float] = (22, 13)):
    fig, ax = plt.subplots(figsize=figsize, dpi=130)
    xs = [p.position[0] for p in pn.places.values()] + \
         [t.position[0] for t in pn.transitions.values()]
    ys = [p.position[1] for p in pn.places.values()] + \
         [t.position[1] for t in pn.transitions.values()]
    # Right margin: compartment-band labels are right-aligned at x_max - 0.15
    # (see _draw_compartments), so this gap also sets how far they sit from the
    # rightmost node column. It is wide enough to keep them clear of that column
    # even in the perturbation model, where the X[8] disease lane (e.g. the
    # "Bivalent crosslink" transition at x=19) would otherwise sit under the
    # "Post-synaptic density (PSD)" band label. The figure is Y-limited, so a
    # wider x-range here only shifts the labels right -- it doesn't rescale the
    # diagram. Baseline has no node that far right, so it just gains right
    # whitespace.
    x_min, x_max = min(xs) - 1.7, max(xs) + 3.6
    y_min, y_max = _vy(-3.7), _vy(11.3)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_aspect("equal")
    ax.set_axis_off()
    _draw_compartments(ax, x_min, x_max)
    return fig, ax, (x_min, x_max, y_min, y_max)


def _draw_arcs(ax, pn: PetriNet, firing: Optional[str] = None) -> None:
    """Draw all arcs, clipped to node borders and with catalyst pairs bowed apart.

    For every transition, a place that is *both* a pre- and post-place (a
    read-arc / catalyst, e.g. NMDAR_open in T_Ca_in) gets its two arcs drawn with
    the same `rad`, which arc3 renders on opposite sides so they no longer
    coincide. Endpoints are clipped to the place circle and transition box so the
    arrowheads sit just outside the nodes instead of under them.
    """
    for tname, t in pn.transitions.items():
        in_places = {a.place for a in t.inputs}
        out_places = {a.place for a in t.outputs}
        bidir = in_places & out_places
        firing_now = (firing == tname)

        for arc in t.inputs:                       # place -> transition
            pc = pn.places[arc.place].position
            tc = t.position
            p0 = _place_edge(pc, tc)
            p1 = _box_edge(tc, pc)
            colour = "#c0392b" if firing_now else "#3b3b3b"
            lw = 2.0 if firing_now else 1.1
            _arc(ax, p0, p1, colour=colour, lw=lw,
                 label=(f"\u00d7{arc.weight}" if arc.weight > 1 else ""),
                 rad=(_BIDIR_RAD if arc.place in bidir else 0.0))

        for arc in t.outputs:                      # transition -> place
            pc = pn.places[arc.place].position
            tc = t.position
            p0 = _box_edge(tc, pc)
            p1 = _place_edge(pc, tc)
            colour = "#e67e22" if firing_now else "#4a148c"
            lw = 2.0 if firing_now else 1.1
            _arc(ax, p0, p1, colour=colour, lw=lw,
                 label=(f"\u00d7{arc.weight}" if arc.weight > 1 else ""),
                 rad=(_BIDIR_RAD if arc.place in bidir else 0.0))


def _draw_all(ax, pn: PetriNet, firing: Optional[str] = None) -> None:
    _draw_arcs(ax, pn, firing=firing)
    for p in pn.places.values():
        _draw_place(ax, p, pn)
    for tname, t in pn.transitions.items():
        _draw_transition(ax, t, firing=(firing == tname))


def draw_static(pn: PetriNet, out_path: str,
                title: Optional[str] = None,
                subtitle: Optional[str] = None,
                pdf_only: bool = False,
                inc_cell_layer_sim: bool = False) -> None:
    """Render the Petri Net to a single static file.

    `out_path` names the PNG.  With `pdf_only` False (default) the PNG is
    written; with `pdf_only` True a vector PDF with the same basename is
    written *instead* (no PNG).  The two formats are mutually exclusive.

    When ``inc_cell_layer_sim=True``, a realistic synapse cross-section image
    is placed as a narrow side panel to the right of the Petri net diagram.
    """
    fig, ax, bbox = _figure_for(pn)
    _draw_all(ax, pn)
    _draw_legend(ax, bbox[0], bbox[2])
    title = title or pn.name
    fig.suptitle(title, fontsize=17, weight="bold", y=0.985)
    if subtitle or pn.description:
        ax.set_title(subtitle or pn.description,
                     fontsize=10.5, style="italic", color="#333",
                     loc="center", pad=6, wrap=True)

    # tight_layout reflows the main axes, so add the side panel *after* it: the
    # panel is pinned to the axes' final drawn position via live transforms.
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    if inc_cell_layer_sim:
        _add_cell_layer_sidepanel(fig, ax)

    if pdf_only:
        base, _ = os.path.splitext(out_path)
        fig.savefig(base + ".pdf", bbox_inches="tight", facecolor="white")
    else:
        fig.savefig(out_path, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def animate(
    pn: PetriNet,
    trace: List[Tuple[float, str, Dict[str, int]]],
    out_gif: str,
    fps: int = 8,
    n_frames: int = 80,
    title: Optional[str] = None,
    inc_cell_layer_sim: bool = False,
) -> None:
    """Sample the trace into `n_frames` keyframes and write a GIF.

    When ``inc_cell_layer_sim=True``, a realistic synapse cross-section image
    is placed as a narrow side panel to the right of the Petri net animation.
    """
    from matplotlib.animation import FuncAnimation, PillowWriter

    idxs = np.linspace(0, len(trace) - 1, n_frames, dtype=int)
    idxs = sorted(set(idxs.tolist()))

    # Use the same figure size as draw_static so a data-unit maps to the same
    # number of points: the transition-label auto-fit (calibrated for the PNG)
    # then renders identically in the GIF instead of overflowing the boxes.
    fig, ax, bbox = _figure_for(pn)

    # Trim the wide left/right white margins. _figure_for uses a fixed 22x13
    # figure, but with set_aspect("equal") the (narrower) net is centred in it,
    # leaving big L/R whitespace. The static PNG hides this via
    # savefig(bbox_inches="tight"), which GIF writers can't use -- so instead
    # size the figure width to the net's own aspect ratio and pull the axes to
    # the edges, reserving a right strip for the synapse side panel when shown.
    data_aspect = (bbox[1] - bbox[0]) / (bbox[3] - bbox[2])
    left, right = (0.012, 0.88) if inc_cell_layer_sim else (0.012, 0.985)
    bottom, top = 0.05, 0.95
    fig_h = fig.get_figheight()
    fig.set_size_inches(
        data_aspect * fig_h * (top - bottom) / (right - left), fig_h,
        forward=True,
    )
    fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top)

    if title:
        fig.suptitle(title, fontsize=15, weight="bold", y=0.985)

    # Add side panel once (persists across animation frames). Added *after* the
    # resize above so it pins to the net's final drawn edge.
    if inc_cell_layer_sim:
        _add_cell_layer_sidepanel(fig, ax)

    hud = fig.text(0.01, 0.01, "", fontsize=8, family="monospace",
                   color="#222",
                   bbox=dict(facecolor="white", edgecolor="#999",
                             alpha=0.88, pad=4))

    def update(frame_idx):
        ax.clear()
        ax.set_axis_off()
        x_min, x_max = bbox[0], bbox[1]
        ax.set_xlim(bbox[0], bbox[1])
        ax.set_ylim(bbox[2], bbox[3])
        ax.set_aspect("equal")
        _draw_compartments(ax, x_min, x_max)

        i = idxs[frame_idx]
        t_i, tr_i, marking = trace[i]
        for n, v in marking.items():
            pn.places[n].tokens = v

        _draw_all(ax, pn,
                  firing=(tr_i if tr_i not in ("<init>", "<end>") else None))
        _draw_legend(ax, x_min, bbox[2])

        ax.text(x_min + 0.3, bbox[3] - 0.25,
                f"t = {t_i:6.3f} s    fired: {tr_i}",
                fontsize=11, weight="bold", color="#1b1b1b",
                bbox=dict(facecolor="#fff9c4", edgecolor="#666",
                          pad=3, alpha=0.95))

        hud.set_text("  ".join(
            f"{k}={v}" for k, v in sorted(marking.items()) if v > 0
        )[:160])
        return []

    anim = FuncAnimation(fig, update, frames=len(idxs),
                         interval=1000 / fps, blit=False)
    # dpi raised from 110 -> 150 for a sharper GIF (matches the static figure's
    # crispness). Combined with the larger figure above this notably increases
    # GIF resolution and file size.
    anim.save(out_gif, writer=PillowWriter(fps=fps), dpi=150)
    plt.close(fig)


# =========================================================================
# Combined animation: Petri net (left) + dynamic curves (right)
# =========================================================================
def animate_combined(
    pn: PetriNet,
    trace: List[Tuple[float, str, Dict[str, int]]],
    out_gif: str,
    t_max: float,
    fps: int = 8,
    n_frames: int = 80,
    title: Optional[str] = None,
    inc_cell_layer_sim: bool = False,
) -> None:
    """Animated GIF with two synchronised subplots.

    Left  : Petri net token-flow (same rendering as ``animate``).
    Right : Ca²⁺ curves — instantaneous (top) and cumulative (bottom)
            overlaid, growing in lockstep with the Petri net animation.

    When ``inc_cell_layer_sim=True``, a realistic synapse cross-section image
    is placed as a narrow side panel to the right of the Petri net panel.
    """
    from matplotlib.animation import FuncAnimation, PillowWriter

    # -- sample trace for Petri net frames (same strategy as `animate`) --
    idxs = np.linspace(0, len(trace) - 1, n_frames, dtype=int)
    idxs = sorted(set(idxs.tolist()))

    # -- curve data (pre-compute on a dense grid) --
    t_grid = np.linspace(0, t_max, n_frames)
    ca_inst = sample_trace(trace, t_grid, "Ca_in")
    tc, cum_ca = cumulative_ca_influx(trace)
    cum_interp = np.interp(t_grid, tc, cum_ca)

    # -- build the dual-panel figure --
    fig = plt.figure(figsize=(34, 13), dpi=120)

    # Left subplot: Petri net
    ax_pn = fig.add_subplot(1, 2, 1)
    ax_pn.set_axis_off()

    # Right subplot: curves (two stacked axes)
    ax_ca_inst = fig.add_subplot(6, 2, 4)   # top-right  (rows 1-2 of 6)
    ax_ca_cum  = fig.add_subplot(6, 2, 6)   # bottom-right (rows 3-4 of 6)
    # Hide the unused right-side row-5/6 subplot slots
    for idx in range(8, 13, 2):
        ax_hide = fig.add_subplot(6, 2, idx)
        ax_hide.set_visible(False)

    # Titles & styling for curve axes
    for ax in (ax_ca_inst, ax_ca_cum):
        ax.set_xlim(0, t_max)
        ax.grid(True, alpha=0.25)

    ax_ca_inst.set_ylabel("Cytosolic Ca²⁺ tokens", fontsize=10)
    ax_ca_inst.set_title("Instantaneous cytosolic Ca²⁺", fontsize=11, weight="bold")
    ax_ca_inst.set_ylim(0, max(ca_inst.max(), 1) * 1.25 + 0.5)

    ax_ca_cum.set_xlabel("Simulation time (s)", fontsize=10)
    ax_ca_cum.set_ylabel("Σ Ca²⁺ influx events", fontsize=10)
    ax_ca_cum.set_title("Cumulative Ca²⁺ influx", fontsize=11, weight="bold")
    ax_ca_cum.set_ylim(0, max(cum_interp.max(), 1) * 1.25 + 1.0)

    # --- Set up Petri net ---
    xs = [p.position[0] for p in pn.places.values()] + \
         [t.position[0] for t in pn.transitions.values()]
    x_min, x_max = min(xs) - 1.7, max(xs) + 3.6
    y_min, y_max = _vy(-3.7), _vy(11.3)
    ax_pn.set_xlim(x_min, x_max)
    ax_pn.set_ylim(y_min, y_max)
    ax_pn.set_aspect("equal")
    _draw_compartments(ax_pn, x_min, x_max)

    # --- Add cell layer side panel (once, persists across frames) ---
    if inc_cell_layer_sim:
        _add_cell_layer_sidepanel(fig, ax_pn)

    # --- Curve artists (initial empty) ---
    line_inst,      = ax_ca_inst.plot([], [], color="#2e7d32", lw=2.2, label="Ca²⁺ in cytosol")
    line_cum,       = ax_ca_cum.plot([],  [], color="#c62828", lw=2.2, label="Cumulative influx")
    head_inst = ax_ca_inst.scatter([], [], color="#2e7d32", s=60, zorder=5)
    head_cum  = ax_ca_cum.scatter([],  [], color="#c62828", s=60, zorder=5)
    ax_ca_inst.legend(loc="upper left", framealpha=0.9, fontsize=9)
    ax_ca_cum.legend(loc="upper left", framealpha=0.9, fontsize=9)

    # --- Figure-level HUD ---
    if title:
        fig.suptitle(title, fontsize=15, weight="bold", y=0.99)

    hud_pn = fig.text(0.01, 0.01, "", fontsize=7, family="monospace",
                      color="#222",
                      bbox=dict(facecolor="white", edgecolor="#999",
                                alpha=0.88, pad=3))

    hud_curves = fig.text(0.52, 0.01, "", fontsize=8, family="monospace",
                          color="#222",
                          bbox=dict(facecolor="#fff9c4", edgecolor="#888",
                                    alpha=0.92, pad=3))

    # --- Frame update ---
    def update(frame_k):
        # ---- Petri net ----
        ax_pn.clear()
        ax_pn.set_axis_off()
        ax_pn.set_xlim(x_min, x_max)
        ax_pn.set_ylim(y_min, y_max)
        ax_pn.set_aspect("equal")
        _draw_compartments(ax_pn, x_min, x_max)

        i = idxs[frame_k]
        t_i, tr_i, marking = trace[i]
        for n, v in marking.items():
            pn.places[n].tokens = v

        _draw_all(ax_pn, pn,
                  firing=(tr_i if tr_i not in ("<init>", "<end>") else None))
        _draw_legend(ax_pn, x_min, y_min)

        ax_pn.text(x_min + 0.3, y_max - 0.25,
                   f"t = {t_i:6.3f} s    fired: {tr_i}",
                   fontsize=11, weight="bold", color="#1b1b1b",
                   bbox=dict(facecolor="#fff9c4", edgecolor="#666",
                             pad=3, alpha=0.95))

        hud_pn.set_text("  ".join(
            f"{k}={v}" for k, v in sorted(marking.items()) if v > 0
        )[:140])

        # ---- Curves (grow from frame 0 to frame_k) ----
        k = frame_k
        line_inst.set_data(t_grid[:k+1], ca_inst[:k+1])
        line_cum.set_data(t_grid[:k+1], cum_interp[:k+1])

        head_inst.set_offsets([[t_grid[k], ca_inst[k]]])
        head_cum.set_offsets([[t_grid[k], cum_interp[k]]])

        hud_curves.set_text(
            f"t={t_grid[k]:5.2f}s    "
            f"cyt Ca²⁺ = {ca_inst[k]:.0f}    "
            f"Σ influx = {cum_interp[k]:.0f}"
        )

        return []

    anim = FuncAnimation(fig, update, frames=len(idxs),
                         interval=1000 / fps, blit=False)
    anim.save(out_gif, writer=PillowWriter(fps=fps), dpi=120)
    plt.close(fig)


# =========================================================================
# NMDAR comparison: side-by-side Petri nets + cumulative Ca2+ curve
# =========================================================================
def _draw_pn_on_ax(
    ax, pn: PetriNet, firing: Optional[str] = None, *, show_legend: bool = False,
    inc_cell_layer_sim: bool = False,
) -> None:
    """Render a single Petri net onto an existing *Axes* object.

    Reuses the same internal helpers as ``draw_static`` and ``animate``,
    drawing compartments, arcs, places, transitions, and (optionally) the
    legend inside *ax*.

    When ``inc_cell_layer_sim=True``, the side panel is drawn — but only when
    this function is called from the comparison figure (the figure reference is
    accessed via *ax*'s parent figure).
    """
    xs = [p.position[0] for p in pn.places.values()] + \
         [t.position[0] for t in pn.transitions.values()]
    ys = [p.position[1] for p in pn.places.values()] + \
         [t.position[1] for t in pn.transitions.values()]
    x_min, x_max = min(xs) - 1.7, max(xs) + 3.6
    y_min, y_max = _vy(-3.7), _vy(11.3)

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_aspect("equal")
    ax.set_axis_off()

    _draw_compartments(ax, x_min, x_max)
    _draw_arcs(ax, pn, firing=firing)
    for p in pn.places.values():
        _draw_place(ax, p, pn)
    for tname, t in pn.transitions.items():
        _draw_transition(ax, t, firing=(firing == tname))

    if show_legend:
        _draw_legend(ax, x_min, y_min)

    if inc_cell_layer_sim:
        _add_cell_layer_sidepanel(ax.figure, ax)


def draw_nmdar_pn_comparison(
    base: PetriNet,
    pert: PetriNet,
    trace_b: List,
    trace_p: List,
    out_path: str,
    t_max: float,
    pdf_only: bool = False,
    inc_cell_layer_sim: bool = False,
) -> None:
    """Static figure: baseline (left) + perturbation (right) Petri nets with
    cumulative Ca²⁺ influx curve underneath (spanning the full width).

    When ``inc_cell_layer_sim=True``, the realistic synapse cross-section image
    is placed as a narrow side panel to the right of each Petri net panel.
    """
    from .simulation import cumulative_ca_influx

    tc_b, cum_b = cumulative_ca_influx(trace_b)
    tc_p, cum_p = cumulative_ca_influx(trace_p)

    fig = plt.figure(figsize=(26, 16), dpi=130)
    gs = fig.add_gridspec(2, 2, height_ratios=[3, 1.3], hspace=0.22,
                          left=0.01, right=0.99, top=0.93, bottom=0.06)

    ax_base = fig.add_subplot(gs[0, 0])
    ax_pert = fig.add_subplot(gs[0, 1])
    ax_cum  = fig.add_subplot(gs[1, :])

    _draw_pn_on_ax(ax_base, base, show_legend=True,
                   inc_cell_layer_sim=False)    # only one side-panel per figure
    ax_base.set_title("Baseline (healthy NMDAR signalling)", fontsize=13,
                      weight="bold", pad=4)

    _draw_pn_on_ax(ax_pert, pert, show_legend=False,
                   inc_cell_layer_sim=inc_cell_layer_sim)   # rightmost gets the image
    ax_pert.set_title("Perturbation (anti-NMDAR encephalitis)", fontsize=13,
                      weight="bold", pad=4)

    ax_cum.plot(tc_b, cum_b, color="#2e7d32", lw=2.5, label="Baseline (healthy)")
    ax_cum.plot(tc_p, cum_p, color="#c62828", lw=2.5, label="Perturbation (anti-NMDAR)")
    ax_cum.set_xlim(0, t_max)
    ymax_cum = max(cum_b.max(), cum_p.max()) * 1.15 + 1.0
    ax_cum.set_ylim(0, ymax_cum)
    ax_cum.set_xlabel("Simulation time (s)", fontsize=12)
    ax_cum.set_ylabel(r"$\Sigma$ Ca$^{2+}$ influx events", fontsize=12)
    ax_cum.set_title(r"Cumulative $Ca^{2+}$ influx through NMDAR", fontsize=14,
                     weight="bold")
    ax_cum.grid(True, alpha=0.30)
    ax_cum.legend(loc="upper left", framealpha=0.9, fontsize=11)

    if cum_b[-1] > 0:
        ratio = cum_p[-1] / cum_b[-1] * 100
        ax_cum.annotate(
            f"Perturbation / Baseline at t = {t_max:.1f}s : "
            f"{cum_p[-1]:.0f} / {cum_b[-1]:.0f}  =  {ratio:.1f} %",
            xy=(tc_p[-1], cum_p[-1]),
            xytext=(t_max * 0.50, cum_p[-1] + (cum_b[-1] - cum_p[-1]) * 0.5),
            fontsize=11, color="#222", weight="bold",
            arrowprops=dict(arrowstyle="->", color="#666"),
            bbox=dict(facecolor="#fff9c4", edgecolor="#888", pad=5),
        )

    fig.suptitle("Anti-NMDAR Encephalitis \u2014 Petri Net Models and Cumulative Ca\u00b2\u207a Influx",
                 fontsize=16, weight="bold", y=0.97)

    import os
    if pdf_only:
        base_path, _ = os.path.splitext(out_path)
        fig.savefig(base_path + ".pdf", facecolor="white", bbox_inches="tight")
    else:
        fig.savefig(out_path, dpi=150, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def animate_nmdar_pn_comparison(
    base: PetriNet,
    pert: PetriNet,
    trace_b: List,
    trace_p: List,
    out_gif: str,
    t_max: float,
    fps: int = 8,
    n_frames: int = 60,
    title: Optional[str] = None,
    inc_cell_layer_sim: bool = False,
) -> None:
    """Animated GIF: baseline (left) + perturbation (right) Petri nets with
    synchronised cumulative Ca²⁺ influx curve underneath.

    When ``inc_cell_layer_sim=True``, the realistic synapse cross-section image
    is placed as a narrow side panel to the right of each Petri net panel.
    """
    from matplotlib.animation import FuncAnimation, PillowWriter
    from .simulation import cumulative_ca_influx

    import numpy as np

    idxs = np.linspace(0, len(trace_b) - 1, n_frames, dtype=int)
    idxs = sorted(set(idxs.tolist()))

    tc_b, cum_b = cumulative_ca_influx(trace_b)
    tc_p, cum_p = cumulative_ca_influx(trace_p)
    t_grid = np.linspace(0, t_max, n_frames)
    cum_b_interp = np.interp(t_grid, tc_b, cum_b)
    cum_p_interp = np.interp(t_grid, tc_p, cum_p)

    fig = plt.figure(figsize=(26, 16), dpi=120)
    gs = fig.add_gridspec(2, 2, height_ratios=[3, 1.3], hspace=0.22,
                          left=0.01, right=0.99, top=0.93, bottom=0.06)

    ax_base = fig.add_subplot(gs[0, 0])
    ax_pert = fig.add_subplot(gs[0, 1])
    ax_cum  = fig.add_subplot(gs[1, :])

    for ax in (ax_base, ax_pert):
        ax.set_axis_off()

    def _pn_extent(pn):
        xs = [p.position[0] for p in pn.places.values()] + \
             [t.position[0] for t in pn.transitions.values()]
        return min(xs) - 1.7, max(xs) + 3.6, _vy(-3.7), _vy(11.3)

    xmb, xMb, ymb, yMb = _pn_extent(base)
    xmp, xMp, ymp, yMp = _pn_extent(pert)

    for ax, (xm, xM, ym, yM) in zip(
        (ax_base, ax_pert), ((xmb, xMb, ymb, yMb), (xmp, xMp, ymp, yMp))
    ):
        ax.set_xlim(xm, xM)
        ax.set_ylim(ym, yM)
        ax.set_aspect("equal")

    # Add a single side panel to the rightmost (perturbation) Petri net
    if inc_cell_layer_sim:
        _add_cell_layer_sidepanel(fig, ax_pert)

    ymax_cum = max(cum_b.max(), cum_p.max()) * 1.15 + 1.0
    ax_cum.set_xlim(0, t_max)
    ax_cum.set_ylim(0, ymax_cum)
    ax_cum.set_xlabel("Simulation time (s)", fontsize=12)
    ax_cum.set_ylabel(r"$\Sigma$ Ca$^{2+}$ influx events", fontsize=12)
    ax_cum.set_title(r"Cumulative Ca$^{2+}$ influx through NMDAR", fontsize=14,
                     weight="bold")
    ax_cum.grid(True, alpha=0.30)

    line_b,   = ax_cum.plot([], [], color="#2e7d32", lw=2.5, label="Baseline (healthy)")
    line_p,   = ax_cum.plot([], [], color="#c62828", lw=2.5, label="Perturbation (anti-NMDAR)")
    head_b    = ax_cum.scatter([], [], color="#2e7d32", s=70, zorder=5)
    head_p    = ax_cum.scatter([], [], color="#c62828", s=70, zorder=5)
    ax_cum.legend(loc="upper left", framealpha=0.9, fontsize=11)

    hud = fig.text(0.50, 0.01, "", fontsize=9, family="monospace", ha="center",
                   bbox=dict(facecolor="#fff9c4", edgecolor="#888", pad=4))

    if title:
        fig.suptitle(title, fontsize=16, weight="bold", y=0.97)

    def update(frame_k):
        k = frame_k

        ax_base.clear()
        ax_base.set_axis_off()
        ax_base.set_xlim(xmb, xMb)
        ax_base.set_ylim(ymb, yMb)
        ax_base.set_aspect("equal")
        _draw_compartments(ax_base, xmb, xMb)

        i = idxs[min(k, len(idxs) - 1)]
        t_i, tr_i, marking = trace_b[i]
        for n, v in marking.items():
            base.places[n].tokens = v
        _draw_arcs(ax_base, base,
                   firing=(tr_i if tr_i not in ("<init>", "<end>") else None))
        for p in base.places.values():
            _draw_place(ax_base, p, base)
        for tname, t in base.transitions.items():
            _draw_transition(ax_base, t,
                             firing=(tr_i if tr_i not in ("<init>", "<end>") else None))
        ax_base.set_title("Baseline (healthy NMDAR signalling)", fontsize=13,
                          weight="bold", pad=4)

        ax_pert.clear()
        ax_pert.set_axis_off()
        ax_pert.set_xlim(xmp, xMp)
        ax_pert.set_ylim(ymp, yMp)
        ax_pert.set_aspect("equal")
        _draw_compartments(ax_pert, xmp, xMp)

        tp_idx = min(range(len(trace_p)), key=lambda j: abs(trace_p[j][0] - t_i))
        tp_i, trp_i, marking_p = trace_p[tp_idx]
        for n, v in marking_p.items():
            pert.places[n].tokens = v
        _draw_arcs(ax_pert, pert,
                   firing=(trp_i if trp_i not in ("<init>", "<end>") else None))
        for p in pert.places.values():
            _draw_place(ax_pert, p, pert)
        for tname, t in pert.transitions.items():
            _draw_transition(ax_pert, t,
                             firing=(trp_i if trp_i not in ("<init>", "<end>") else None))
        ax_pert.set_title("Perturbation (anti-NMDAR encephalitis)", fontsize=13,
                          weight="bold", pad=4)

        line_b.set_data(t_grid[:k+1], cum_b_interp[:k+1])
        line_p.set_data(t_grid[:k+1], cum_p_interp[:k+1])
        head_b.set_offsets([[t_grid[k], cum_b_interp[k]]])
        head_p.set_offsets([[t_grid[k], cum_p_interp[k]]])

        ratio = (cum_p_interp[k] / cum_b_interp[k]) * 100 if cum_b_interp[k] > 0 else 0.0
        hud.set_text(
            f"t = {t_grid[k]:5.2f}s    "
            f"Baseline cumulative: {cum_b_interp[k]:5.0f}    "
            f"Perturbation cumulative: {cum_p_interp[k]:5.0f}    "
            f"Ratio: {ratio:5.1f}%"
        )
        return []

    anim = FuncAnimation(fig, update, frames=len(idxs),
                         interval=1000 / fps, blit=False)
    anim.save(out_gif, writer=PillowWriter(fps=fps), dpi=120)
    plt.close(fig)
