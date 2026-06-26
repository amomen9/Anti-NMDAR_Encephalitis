"""
Snoopy XML exporters for the anti-NMDAR encephalitis Petri-Net models.

Two file flavours are produced per model so all capabilities of the Python
implementation - and a few extra ones Snoopy can offer - are covered:

* ``.spn``  Stochastic Petri Net.
            - Discrete integer tokens.
            - ``MassAction(k)`` rate functions match the Gillespie SSA in
              ``src/simulation.py`` (rate * combinatorial term over input arcs).
            - Suitable for Snoopy's stochastic simulator, CSL model checking
              and structural analyses (boundedness, P/T-invariants).

* ``.cpn``  Continuous Petri Net.
            - Real-valued markings, mass-action ODE semantics.
            - Lets Snoopy run a deterministic ODE simulation, which the Python
              code does not provide.
            - Same topology, capacities, weights and rate constants as the
              ``.spn``, so trajectories are directly comparable.

Both files use Snoopy's native XML format (``<Snoopy version="2"
revision="2.0">``) reverse-engineered from published example models, with
the full attribute set Snoopy needs to render and simulate (Name, Marking,
ID, Fixed, Logic, Comment for places; Name, FunctionList, ID, Logic,
Comment, Reversible, Node Type for transitions; Multiplicity + Comment on
edges; plus per-attribute and per-node graphics, parent graphic boxes, and
default General / Simulation Properties / Plot metadata).

Run directly with ``python src/snoopy_export.py`` or via ``python run.py``.
"""

from __future__ import annotations

import datetime as _dt
import os
import sys
import xml.etree.ElementTree as ET
from typing import Dict, List, Tuple

# Allow running as `python src/snoopy_export.py` directly (no package context)
# as well as `python -m src.snoopy_export` or `from src.snoopy_export import ...`.
try:
    from .petri_net import Arc, PetriNet, Place, Transition  # type: ignore
except ImportError:                                         # pragma: no cover
    _here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, os.path.dirname(_here))
    from src.petri_net import Arc, PetriNet, Place, Transition  # type: ignore


# ---------------------------------------------------------------------------
# Layout scaling: the Python grid uses unit-spacing integers; Snoopy stores
# pixel coordinates.  The scale is square (X == Y) so the diagram preserves
# the aspect ratio of the Python rendering.  Every pixel-dimension in the PN
# (node w/h, label offsets, annotation-side offsets) is derived from
# ``_GEOM`` so changing ``_SCALE`` rescales the *whole* diagram uniformly --
# nodes grow at the same rate as the layout, exactly like a vector zoom.
# ---------------------------------------------------------------------------
_SCALE = 160.0
_SCALE_X = _SCALE
_SCALE_Y = _SCALE
_OFFSET_X = 80.0
_OFFSET_Y = 80.0
_Y_FLIP = 12.0                # rows above this go up, below it go down

# Pixel-geometry factor: 1.0 reproduces the original Snoopy default look
# (positions at ~60 px / unit, node w/h = 26 px).  At ``_SCALE = 160`` the
# diagram - positions and node sizes alike - is ``160/60`` larger.
_GEOM = _SCALE / 60.0

# Node parent box (the circle / square drawn on the canvas)
_NODE_W = 26.0 * _GEOM
_NODE_H = 26.0 * _GEOM

# Per-attribute label offsets, all scaled with the geometry
_NAME_XOFF      = 4.0  * _GEOM
_NAME_YOFF      = 24.0 * _GEOM
_TNAME_XOFF     = 2.0  * _GEOM
_TNAME_YOFF     = 25.0 * _GEOM
_ID_XOFF        = 20.0 * _GEOM
_PLACE_COMMENT_YOFF = 40.0 * _GEOM
_FUNC_XOFF      = 9.0  * _GEOM
_FUNC_YOFF      = -25.0 * _GEOM


def _snoopy_xy(x: float, y: float) -> Tuple[float, float]:
    return (_OFFSET_X + x * _SCALE_X,
            _OFFSET_Y + (_Y_FLIP - y) * _SCALE_Y)


# Where the small annotation label for a transition should sit relative to
# the transition centre.  Mirrors ``Transition.annotation_side`` in the
# Python model so paired forward/reverse transitions split above/below
# cleanly.  Scaled with the same geometry factor as the nodes themselves.
_ANNOT_OFFSETS = {
    "below": (0.0,            48.0 * _GEOM),
    "above": (0.0,           -48.0 * _GEOM),
    "right": ( 72.0 * _GEOM,   0.0),
    "left":  (-72.0 * _GEOM,   0.0),
}


# ---------------------------------------------------------------------------
# ID allocator: every element in a Snoopy file (node, attribute, graphic,
# edge, point, metadata, ...) gets a unique integer ``id``.  We keep one
# global counter and hand out blocks of consecutive IDs as needed.
# ---------------------------------------------------------------------------
class _IDPool:
    def __init__(self, start: int = 1) -> None:
        self._next = start

    def one(self) -> int:
        i = self._next
        self._next += 1
        return i

    def block(self, n: int) -> List[int]:
        base = self._next
        self._next += n
        return list(range(base, base + n))


# ---------------------------------------------------------------------------
# XML helpers
# ---------------------------------------------------------------------------
def _cdata_attr(parent: ET.Element, name: str, attr_id: int, value) -> ET.Element:
    """Build an ``<attribute name=... id=... net="1">`` whose payload is CDATA.

    ElementTree has no native CDATA support, so we stash the text and rewrite
    it during serialisation (see ``_serialise``).
    """
    a = ET.SubElement(parent, "attribute", {
        "name": name,
        "id":   str(attr_id),
        "net":  "1",
    })
    a.text = f"__CDATA__{value}__/CDATA__"
    return a


def _graphic(parent: ET.Element, gid: int, *, x: float, y: float,
             xoff: float = 0.0, yoff: float = 0.0,
             grparent: int | None = None, show: str = "1",
             w: float | None = None, h: float | None = None,
             pen: str = "0,0,0", brush: str = "255,255,255",
             thickness: str = "1") -> ET.Element:
    attrs = {
        "x":         f"{x:.2f}",
        "y":         f"{y:.2f}",
        "id":        str(gid),
        "net":       "1",
        "show":      show,
        "state":     "1",
        "pen":       pen,
        "brush":     brush,
        "thickness": thickness,
    }
    if xoff:
        attrs["xoff"] = f"{xoff:.2f}"
    if yoff:
        attrs["yoff"] = f"{yoff:.2f}"
    if grparent is not None:
        attrs["grparent"] = str(grparent)
    if w is not None:
        attrs["w"] = f"{w:.2f}"
    if h is not None:
        attrs["h"] = f"{h:.2f}"
    return ET.SubElement(parent, "graphic", attrs)


def _graphics_with_one(parent: ET.Element, **kwargs) -> ET.Element:
    g = ET.SubElement(parent, "graphics", count="1")
    _graphic(g, **kwargs)
    return g


def _graphics_empty(parent: ET.Element) -> ET.Element:
    return ET.SubElement(parent, "graphics", count="0")


# ---------------------------------------------------------------------------
# Per-element writers
# ---------------------------------------------------------------------------
def _write_place(nodeclass_el: ET.Element, ids: _IDPool, idx: int, p: Place,
                 net_type: str) -> Tuple[int, int]:
    """Emit one ``<node>`` for a place. Returns (node_id, parent_graphic_id)."""
    node_type = "Stochastic Place" if net_type == "stoch" else "Continuous Place"

    # 13 consecutive ids: node, Name+gfx, NodeType, ID+gfx, Marking+gfx,
    # Fixed, Logic, Comment+gfx, parent gfx.
    (node_id, name_a, name_g, ntype_a, id_a, id_g, mark_a, mark_g,
     fixed_a, logic_a, comment_a, comment_g, parent_gfx_id) = ids.block(13)

    x, y = _snoopy_xy(*p.position)

    node = ET.SubElement(nodeclass_el, "node", id=str(node_id), net="1")

    # --- Name (Snoopy identifier; must be [A-Za-z][A-Za-z0-9_]*) ---
    a = _cdata_attr(node, "Name", name_a, p.name)
    _graphics_with_one(a, gid=name_g, x=x, y=y, xoff=_NAME_XOFF, yoff=_NAME_YOFF,
                       grparent=parent_gfx_id, show="1")

    # --- Node Type ---
    a = _cdata_attr(node, "Node Type", ntype_a, node_type)
    _graphics_empty(a)

    # --- ID (Snoopy-local 0-based index within nodeclass) ---
    a = _cdata_attr(node, "ID", id_a, idx)
    _graphics_with_one(a, gid=id_g, x=x + _ID_XOFF, y=y, xoff=_ID_XOFF,
                       grparent=parent_gfx_id, show="0")

    # --- Marking ---
    marking = p.tokens if net_type == "stoch" else float(p.tokens)
    a = _cdata_attr(node, "Marking", mark_a, marking)
    _graphics_with_one(a, gid=mark_g, x=x, y=y, grparent=parent_gfx_id, show="1")

    # --- Fixed (continuous) / always 0 ---
    a = _cdata_attr(node, "Fixed", fixed_a, 0)
    _graphics_empty(a)

    # --- Logic (place is not a logic / read place) ---
    a = _cdata_attr(node, "Logic", logic_a, 0)
    _graphics_empty(a)

    # --- Comment (short annotation only, like the Python figure's caption) ---
    comment_text = p.annotation or (p.label.replace("\n", " ") if p.label else "")
    a = _cdata_attr(node, "Comment", comment_a, comment_text)
    _graphics_with_one(a, gid=comment_g, x=x, y=y + _PLACE_COMMENT_YOFF,
                       yoff=_PLACE_COMMENT_YOFF,
                       grparent=parent_gfx_id, show="1")

    # --- Parent graphic (the place circle on the canvas) ---
    g = ET.SubElement(node, "graphics", count="1")
    _graphic(g, gid=parent_gfx_id, x=x, y=y, w=_NODE_W, h=_NODE_H,
             pen="128,128,128", brush="255,255,255", thickness="3")

    return node_id, parent_gfx_id


def _write_transition(nodeclass_el: ET.Element, ids: _IDPool, idx: int,
                      t: Transition, net_type: str) -> Tuple[int, int]:
    """Emit one ``<node>`` for a transition. Returns (node_id, parent_graphic_id)."""
    node_type = "Stochastic Transition" if net_type == "stoch" else "Continuous Transition"

    # 13 consecutive ids: node, Name+gfx, FunctionList+gfx, ID+gfx,
    # Logic, Comment+gfx, Reversible, NodeType, parent gfx.
    (node_id, name_a, name_g, func_a, func_g, id_a, id_g, logic_a,
     comment_a, comment_g, rev_a, ntype_a, parent_gfx_id) = ids.block(13)

    x, y = _snoopy_xy(*t.position)

    node = ET.SubElement(nodeclass_el, "node", id=str(node_id), net="1")

    # --- Name (Snoopy identifier; must be [A-Za-z][A-Za-z0-9_]*) ---
    a = _cdata_attr(node, "Name", name_a, t.name)
    _graphics_with_one(a, gid=name_g, x=x, y=y, xoff=_TNAME_XOFF, yoff=_TNAME_YOFF,
                       grparent=parent_gfx_id, show="1")

    # --- FunctionList = MassAction(rate) ---
    a = ET.SubElement(node, "attribute", {
        "name": "FunctionList", "type": "ColList",
        "id":   str(func_a), "net": "1",
    })
    col = ET.SubElement(a, "colList", {
        "row_count":  "1", "col_count": "2",
        "active_row": "0", "active_col": "0",
    })
    head = ET.SubElement(col, "colList_head")
    h1 = ET.SubElement(head, "colList_colLabel")
    h1.text = "__CDATA__Function set__/CDATA__"
    h2 = ET.SubElement(head, "colList_colLabel")
    h2.text = "__CDATA__Function__/CDATA__"
    body = ET.SubElement(col, "colList_body")
    row = ET.SubElement(body, "colList_row", nr="0")
    c0 = ET.SubElement(row, "colList_col", nr="0")
    c0.text = "__CDATA__Main__/CDATA__"
    c1 = ET.SubElement(row, "colList_col", nr="1")
    c1.text = f"__CDATA__MassAction({t.rate})__/CDATA__"
    # FunctionList graphic hidden by default: the MassAction(k) value stays
    # in the colList so simulation works, but 30+ rate labels overlapping on
    # the canvas would just be noise. Toggle visible via Snoopy's View menu.
    _graphics_with_one(a, gid=func_g, x=x + _FUNC_XOFF, y=y + _FUNC_YOFF,
                       xoff=_FUNC_XOFF, yoff=_FUNC_YOFF,
                       grparent=parent_gfx_id, show="0")

    # --- ID ---
    a = _cdata_attr(node, "ID", id_a, idx)
    _graphics_with_one(a, gid=id_g, x=x + _ID_XOFF, y=y, xoff=_ID_XOFF,
                       grparent=parent_gfx_id, show="0")

    # --- Logic ---
    a = _cdata_attr(node, "Logic", logic_a, 0)
    _graphics_empty(a)

    # --- Comment (rate label + citation, placed per annotation_side so paired
    #     forward/reverse transitions don't collide on top of each other) ---
    parts = []
    if t.annotation:
        parts.append(t.annotation)
    if t.reference and (not t.annotation or t.reference not in t.annotation):
        parts.append(f"[{t.reference}]")
    comment_text = " ".join(parts)
    cx_off, cy_off = _ANNOT_OFFSETS.get(t.annotation_side, _ANNOT_OFFSETS["below"])
    a = _cdata_attr(node, "Comment", comment_a, comment_text)
    _graphics_with_one(a, gid=comment_g,
                       x=x + cx_off, y=y + cy_off,
                       xoff=cx_off, yoff=cy_off,
                       grparent=parent_gfx_id, show="1")

    # --- Reversible ---
    a = _cdata_attr(node, "Reversible", rev_a, 0)
    _graphics_empty(a)

    # --- Node Type ---
    a = _cdata_attr(node, "Node Type", ntype_a, node_type)
    _graphics_empty(a)

    # --- Parent graphic (transition box) ---
    g = ET.SubElement(node, "graphics", count="1")
    _graphic(g, gid=parent_gfx_id, x=x, y=y, w=_NODE_W, h=_NODE_H,
             pen="128,128,128", brush="255,255,255", thickness="3")

    return node_id, parent_gfx_id


def _write_edge(edgeclass_el: ET.Element, ids: _IDPool,
                src_node_id: int, tgt_node_id: int,
                src_gfx_id: int, tgt_gfx_id: int,
                src_pos: Tuple[float, float], tgt_pos: Tuple[float, float],
                weight: int, comment: str = "") -> None:
    (edge_id, mult_a, mult_g, comment_a, comment_g, edge_gfx_id) = ids.block(6)

    edge = ET.SubElement(edgeclass_el, "edge", {
        "source": str(src_node_id),
        "target": str(tgt_node_id),
        "id":     str(edge_id),
        "net":    "1",
    })

    # Multiplicity label sits at the midpoint of the line.
    mx = 0.5 * (src_pos[0] + tgt_pos[0])
    my = 0.5 * (src_pos[1] + tgt_pos[1])

    a = _cdata_attr(edge, "Multiplicity", mult_a, weight)
    _graphics_with_one(a, gid=mult_g, x=mx, y=my, show="1")

    a = _cdata_attr(edge, "Comment", comment_a, comment)
    _graphics_with_one(a, gid=comment_g, x=mx + 30.0, y=my + 10.0, show="1")

    # Edge line graphic with two endpoints
    g = ET.SubElement(edge, "graphics", count="1")
    line = ET.SubElement(g, "graphic", {
        "id":              str(edge_gfx_id),
        "net":             "1",
        "source":          str(src_gfx_id),
        "target":          str(tgt_gfx_id),
        "state":           "1",
        "show":            "1",
        "pen":             "0,0,0",
        "brush":           "0,0,0",
        "edge_designtype": "3",
        "thickness":       "1",
    })
    points = ET.SubElement(line, "points", count="2")
    ET.SubElement(points, "point",
                  x=f"{src_pos[0]:.2f}", y=f"{src_pos[1]:.2f}")
    ET.SubElement(points, "point",
                  x=f"{tgt_pos[0]:.2f}", y=f"{tgt_pos[1]:.2f}")


# ---------------------------------------------------------------------------
# Default Snoopy metadata: General + Simulation Properties + Plot.
# Snoopy will accept a file without these but populating them avoids the
# "incomplete file" warnings and lets the simulator run with sensible defaults.
# ---------------------------------------------------------------------------
def _write_general_metadata(parent: ET.Element, ids: _IDPool,
                            pn: PetriNet) -> None:
    (meta_id, name_a, name_g, created_a, created_g, authors_a, authors_g,
     keys_a, keys_g, desc_a, desc_g, refs_a, refs_g, parent_g) = ids.block(14)

    m = ET.SubElement(parent, "metadata", id=str(meta_id), net="1")

    now = _dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Long-text fields (Description, References) get show="0" so they don't
    # render as one wide overflowing line above the diagram. Short fields
    # (Name, Authors, Keywords) remain visible.
    for nm, attr_id, gfx_id, value, show in [
        ("Name",        name_a,    name_g,    pn.name,                                             "1"),
        ("Created",     created_a, created_g, now,                                                 "0"),
        ("Authors",     authors_a, authors_g, "",                                                  "1"),
        ("Keywords",    keys_a,    keys_g,    "NMDAR; anti-NMDAR encephalitis; Petri net",         "0"),
        ("Description", desc_a,    desc_g,    pn.description,                                      "0"),
        ("References",  refs_a,    refs_g,    "See README.md",                                     "0"),
    ]:
        a = _cdata_attr(m, nm, attr_id, value)
        _graphics_with_one(a, gid=gfx_id, x=20.0, y=20.0, xoff=25.0, yoff=20.0,
                           grparent=parent_g, show=show)

    g = ET.SubElement(m, "graphics", count="1")
    _graphic(g, gid=parent_g, x=20.0, y=20.0, w=1.0, h=1.0,
             pen="255,255,255", brush="255,255,255")


def _write_simulation_properties(parent: ET.Element, ids: _IDPool) -> None:
    (meta_id, start_a, end_a, step_a, sim_a, sem_a, props_a) = ids.block(7)
    m = ET.SubElement(parent, "metadata", id=str(meta_id), net="1")

    for nm, attr_id, val in [
        ("interval start", start_a, 0),
        ("interval end",   end_a,   8),
        ("output step",    step_a,  100),
        ("simulator",      sim_a,   0),
        ("simulator Semantics", sem_a, 0),
    ]:
        a = _cdata_attr(m, nm, attr_id, val)
        _graphics_empty(a)

    a = ET.SubElement(m, "attribute", {
        "name": "simulator properties", "type": "ColList",
        "id":   str(props_a), "net": "1",
    })
    col = ET.SubElement(a, "colList", {
        "row_count": "10", "col_count": "2",
        "active_row": "0", "active_col": "0",
    })
    head = ET.SubElement(col, "colList_head")
    ET.SubElement(head, "colList_colLabel").text = "__CDATA__Property__/CDATA__"
    ET.SubElement(head, "colList_colLabel").text = "__CDATA__Value__/CDATA__"
    body = ET.SubElement(col, "colList_body")
    rows = [
        ("InBetweenVisualization", "1"),
        ("Refreshrate",            "5000"),
        ("IntegratorStepSize",     "0.01"),
        ("LinearSolverType",       "CVDense"),
        ("CalculateStepSize",      "0"),
        ("CheckValues",            "0"),
        ("OutputNoiseValues",      "0"),
        ("UseODEsReduction",       "1"),
        ("RelativeTol",            "1.0e-5"),
        ("AbsTol",                 "1.0e-10"),
    ]
    for i, (k, v) in enumerate(rows):
        row = ET.SubElement(body, "colList_row", nr=str(i))
        c0 = ET.SubElement(row, "colList_col", nr="0")
        c0.text = f"__CDATA__{k}__/CDATA__"
        c1 = ET.SubElement(row, "colList_col", nr="1")
        c1.text = f"__CDATA__{v}__/CDATA__"
    _graphics_empty(a)
    _graphics_empty(m)


def _write_plot_metadata(parent: ET.Element, ids: _IDPool) -> None:
    """Bare minimum 'Default View' plot metadata so Snoopy opens cleanly."""
    plain_attrs = [
        ("Name",                "Default View"),
        ("ID",                  "0"),
        ("Comment",             ""),
        ("Nodeclass",           "Place"),
        ("ViewerType",          "xyPlot"),
        ("Results",             "Marking"),
        ("NodeColour",          "0"),
        ("RegEx",               "."),
        ("RegExInvert",         "0"),
        ("IsCurrent",           "1"),
        ("ViewTitle",           "Default View"),
        ("ChartWidth",          "600"),
        ("ChartHeight",         "600"),
        ("WindowWidth",         "600"),
        ("WindowHeight",        "600"),
        ("XAxisTitle",          "Time"),
        ("YAxisTitle",          ""),
        ("DefaultLineWidth",    "2"),
        ("DefaultLineStyle",    "100"),
        ("ShowLegend",          "1"),
        ("LegendHorizontalPosition", "32"),
        ("LegendVerticalPosition",   "1"),
        ("ShowLines",           "1"),
        ("ShowSymbols",         "0"),
        ("BarWidth",            "20"),
        ("IntervalWidth",       "0"),
        ("XAxisVariable",       "Time"),
        ("XAxisVariableName",   "Simulation Time"),
        ("SumOfNodes",          "0"),
        ("FixedXAdjustment",    "0"),
        ("FixedYAdjustment",    "0"),
        ("X_Axis_Max",          "1"),
        ("X_Axis_Min",          "0"),
        ("Y_Axis_Max",          "1"),
        ("Y_Axis_Min",          "0"),
        ("X_ZOOM",              "1"),
        ("Y_ZOOM",              "1"),
        ("SortCurvesBy",        "0"),
        ("HistogramFrequencyType", "0"),
    ]
    list_attrs = [("PlaceList", 1), ("CurveInfo", 7)]

    block = ids.block(1 + len(plain_attrs) + 2 * len(list_attrs))
    meta_id = block[0]
    plain_ids = block[1:1 + len(plain_attrs)]
    list_ids = block[1 + len(plain_attrs):]

    m = ET.SubElement(parent, "metadata", id=str(meta_id), net="1")
    for (nm, val), aid in zip(plain_attrs, plain_ids):
        a = _cdata_attr(m, nm, aid, val)
        _graphics_empty(a)

    for (nm, col_count), aid in zip(list_attrs, [list_ids[0], list_ids[1]]):
        a = ET.SubElement(m, "attribute", {
            "name": nm, "type": "ColList", "id": str(aid), "net": "1",
        })
        col = ET.SubElement(a, "colList", {
            "row_count": "0", "col_count": str(col_count),
            "active_row": "0", "active_col": "0",
        })
        head = ET.SubElement(col, "colList_head")
        for _ in range(col_count):
            ET.SubElement(head, "colList_colLabel").text = "__CDATA____/CDATA__"
        ET.SubElement(col, "colList_body")
        _graphics_empty(a)
    _graphics_empty(m)


# ---------------------------------------------------------------------------
# Final XML serialisation: ElementTree has no CDATA; replace placeholders.
# ---------------------------------------------------------------------------
def _serialise(root: ET.Element, path: str) -> None:
    ET.indent(root, space="  ")
    raw = ET.tostring(root, encoding="utf-8").decode("utf-8")
    # Splice <![CDATA[...]]> in wherever we tagged it
    raw = raw.replace("__CDATA__", "<![CDATA[").replace("__/CDATA__", "]]>")
    # ElementTree generates self-closing tags for empty attributes; replace
    # the empty-attr pattern <attribute ... /> -- not strictly required by
    # Snoopy, but safer to leave them as-is. CDATA replacement above already
    # injected the marker text only into elements where we set .text.
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        f.write('<?xml-stylesheet type="text/xsl" href="/xsl/spped2svg.xsl"?>\n')
        f.write(raw)
        if not raw.endswith("\n"):
            f.write("\n")


# ---------------------------------------------------------------------------
# Top-level exporter
# ---------------------------------------------------------------------------
def export(pn: PetriNet, path: str, net_type: str = "stoch") -> None:
    """Write a Snoopy XML file for ``pn`` at ``path``.

    ``net_type`` is ``"stoch"`` for .spn (Stochastic Petri Net) or
    ``"continuous"`` for .cpn (Continuous Petri Net).
    """
    if net_type not in ("stoch", "continuous"):
        raise ValueError(f"Unknown Snoopy net type: {net_type!r}")

    netclass_name = "Stochastic Petri Net" if net_type == "stoch" else "Continuous Petri Net"
    place_nc_name = "Place" if net_type == "stoch" else "Place, Continuous"
    trans_nc_name = "Transition" if net_type == "stoch" else "Transition, Continuous"

    root = ET.Element("Snoopy", version="2", revision="2.0")
    ET.SubElement(root, "netclass", name=netclass_name)

    ids = _IDPool(start=1)

    # ---- nodeclasses ----
    nodeclasses = ET.SubElement(root, "nodeclasses", count="4")

    place_nc = ET.SubElement(nodeclasses, "nodeclass",
                             count=str(len(pn.places)), name=place_nc_name)

    place_node_ids: Dict[str, int]   = {}    # place name -> node id (edges' source/target)
    place_gfx_ids:  Dict[str, int]   = {}    # place name -> parent graphic id (edge graphic source/target)
    place_positions: Dict[str, Tuple[float, float]] = {}

    for idx, p in enumerate(pn.places.values()):
        nid, gid = _write_place(place_nc, ids, idx, p, net_type)
        place_node_ids[p.name] = nid
        place_gfx_ids[p.name]  = gid
        place_positions[p.name] = _snoopy_xy(*p.position)

    trans_nc = ET.SubElement(nodeclasses, "nodeclass",
                             count=str(len(pn.transitions)), name=trans_nc_name)

    trans_node_ids: Dict[str, int] = {}
    trans_gfx_ids:  Dict[str, int] = {}
    trans_positions: Dict[str, Tuple[float, float]] = {}

    for idx, t in enumerate(pn.transitions.values()):
        nid, gid = _write_transition(trans_nc, ids, idx, t, net_type)
        trans_node_ids[t.name] = nid
        trans_gfx_ids[t.name]  = gid
        trans_positions[t.name] = _snoopy_xy(*t.position)

    ET.SubElement(nodeclasses, "nodeclass", count="0", name="Coarse Place")
    ET.SubElement(nodeclasses, "nodeclass", count="0", name="Coarse Transition")

    # ---- edgeclasses ----
    edgeclasses = ET.SubElement(root, "edgeclasses", count="4")

    # Count input + output edges
    n_edges = sum(len(t.inputs) + len(t.outputs) for t in pn.transitions.values())
    edge_class = ET.SubElement(edgeclasses, "edgeclass",
                               count=str(n_edges), name="Edge")

    for t in pn.transitions.values():
        for arc in t.inputs:                               # Place -> Transition
            _write_edge(
                edge_class, ids,
                src_node_id=place_node_ids[arc.place], tgt_node_id=trans_node_ids[t.name],
                src_gfx_id=place_gfx_ids[arc.place],   tgt_gfx_id=trans_gfx_ids[t.name],
                src_pos=place_positions[arc.place],    tgt_pos=trans_positions[t.name],
                weight=arc.weight, comment=arc.annotation,
            )
        for arc in t.outputs:                              # Transition -> Place
            _write_edge(
                edge_class, ids,
                src_node_id=trans_node_ids[t.name], tgt_node_id=place_node_ids[arc.place],
                src_gfx_id=trans_gfx_ids[t.name],   tgt_gfx_id=place_gfx_ids[arc.place],
                src_pos=trans_positions[t.name],    tgt_pos=place_positions[arc.place],
                weight=arc.weight, comment=arc.annotation,
            )

    ET.SubElement(edgeclasses, "edgeclass", count="0", name="Inhibitor Edge")
    ET.SubElement(edgeclasses, "edgeclass", count="0", name="Read Edge")
    ET.SubElement(edgeclasses, "edgeclass", count="0", name="Modifier Edge")

    # ---- metadataclasses ----
    metadataclasses = ET.SubElement(root, "metadataclasses", count="10")

    general = ET.SubElement(metadataclasses, "metadataclass", count="1", name="General")
    _write_general_metadata(general, ids, pn)

    ET.SubElement(metadataclasses, "metadataclass", count="0", name="Comment")
    ET.SubElement(metadataclasses, "metadataclass", count="0", name="Image")

    plot = ET.SubElement(metadataclasses, "metadataclass", count="1", name="Plot")
    _write_plot_metadata(plot, ids)

    sim_props = ET.SubElement(metadataclasses, "metadataclass",
                              count="1", name="Simulation Properties")
    _write_simulation_properties(sim_props, ids)

    ET.SubElement(metadataclasses, "metadataclass", count="0", name="Simulation Results")
    ET.SubElement(metadataclasses, "metadataclass", count="0", name="Constant Class")
    ET.SubElement(metadataclasses, "metadataclass", count="0", name="Function Class")
    ET.SubElement(metadataclasses, "metadataclass", count="0", name="Oberserver Class")
    ET.SubElement(metadataclasses, "metadataclass", count="0", name="Constant Class1")

    _serialise(root, path)


def export_spn(pn: PetriNet, path: str) -> None:
    """Write a Snoopy Stochastic Petri Net (.spn) file."""
    export(pn, path, net_type="stoch")


def export_cpn(pn: PetriNet, path: str) -> None:
    """Write a Snoopy Continuous Petri Net (.cpn) file."""
    export(pn, path, net_type="continuous")


def export_pair(pn: PetriNet, base_path: str) -> None:
    """Write both .spn and .cpn next to each other.

    ``base_path`` should be a path *without* extension, e.g.
    ``output/snoopy/baseline``.
    """
    export_spn(pn, base_path + ".spn")
    export_cpn(pn, base_path + ".cpn")


# ---------------------------------------------------------------------------
# Direct-invocation entry point: `python src/snoopy_export.py`
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    from src.models import build_baseline, build_perturbation

    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "output", "snoopy")
    os.makedirs(out, exist_ok=True)
    export_pair(build_baseline(),     os.path.join(out, "baseline"))
    export_pair(build_perturbation(), os.path.join(out, "perturbation"))
    print(f"Wrote Snoopy files to {out}")
