r"""
Charlie / APNN exporter for the Petri-Net models in this project.

Format
======
This exporter mirrors the structure of a real Charlie-exported APNN file
(``cs1.apnn``, the Cottbus production-cell benchmark) so that the emitted nets
load cleanly into Charlie (Heiner, Schwarick & Wegener; BTU Cottbus DSSZ).

The proven element shapes are:

  * net      ``\beginnet{<id>}`` ... ``\endnet``
  * place    ``\place{<id>}{\name{<n>}\capacity{<c>}\init{<m>}}``
             (sub-tokens packed, NO spaces between them; ``\capacity`` omitted
             when the place is unbounded)
  * trans.   ``\transition{<id>}{\name{<n>}}``
  * arc      ``\arc{<id>}{\from{<src>} \to{<tgt>} \weight{<w>} \type{ordinary}}``
             (sub-tokens SPACE-separated)

Line endings are CRLF, and single blank lines separate the place / transition /
arc blocks (two blank lines before ``\endnet``), matching ``cs1.apnn``.

History note (why earlier exports failed)
-----------------------------------------
An earlier version dropped ``\weight`` from arcs, believing Charlie's APNN
reader discarded any arc carrying a ``\weight`` token, and encoded weight>1 as
parallel arcs.  The reference file ``cs1.apnn`` disproves that: every arc there
carries ``\weight{1} \type{ordinary}``.  The real requirement is the
``\type{ordinary}`` token *next to* ``\weight``.  With both present, a weight>1
arc (the bivalent IgG crosslink, ``2x NMDAR_Ab -> T_crosslink``) is written as a
single ``\weight{2}`` arc -- no parallel-arc trick, and the arc count is the
honest logical count (42 baseline / 58 perturbation).

A second earlier issue -- ``PN = null`` -- came from ``\name`` labels that
contained spaces / punctuation (e.g. "Ca2+ (extracellular)").  ``cs1.apnn``
shows ``\name`` is fine as long as it is a clean identifier token, so this
exporter keeps ``\name`` and sanitises it to ``[A-Za-z0-9_]``.

Charlie consumes net *structure* only (rates / kinetics are ignored), so no
rate information is written.
"""

from __future__ import annotations

import os
import re
import sys
from typing import Dict, List

try:
    from .petri_net import PetriNet  # type: ignore
except ImportError:                                         # pragma: no cover
    _here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, os.path.dirname(_here))
    from src.petri_net import PetriNet  # type: ignore


_WS = re.compile(r"\s+")
_NON_ID = re.compile(r"[^A-Za-z0-9_]")
_CRLF = "\r\n"


def _ident(name: str, fallback: str) -> str:
    """Collapse free text to a single ``[A-Za-z0-9_]+`` identifier token."""
    slug = _NON_ID.sub("_", _WS.sub("_", str(name).strip()))
    slug = re.sub(r"_+", "_", slug).strip("_")
    return slug or fallback


def export_apnn(pn: PetriNet, path: str, *, weights: bool = True,
                with_capacity: bool = True) -> Dict[str, int]:
    r"""Write a Charlie-readable APNN file for ``pn`` at ``path``.

    The output follows the ``cs1.apnn`` reference format: ``\name`` labels,
    optional ``\capacity``, explicit ``\init`` markings, and arcs carrying both
    ``\weight{w}`` and ``\type{ordinary}``.  A weight>1 arc is a single arc
    (no parallel-arc encoding).

    Parameters
    ----------
    weights
        Retained for backward compatibility with existing call sites.  It no
        longer changes correctness: arcs always carry an explicit ``\weight``
        plus ``\type{ordinary}`` because that is the form Charlie reads.
    with_capacity
        Emit ``\capacity{c}`` for places whose model capacity is set.  Places
        with no capacity are written without the token (APNN default =
        unbounded), which is correct for this project's open biological nets.

    Returns
    -------
    dict with the element counts actually written ({"places", "transitions",
    "arcs"}), so callers can assert the figures Charlie should report.
    """
    # Stable, readable identifiers: the model's clean .name doubles as the APNN
    # id and the \name label (proven to load for places/transitions already).
    place_id: Dict[str, str] = {}
    used: set = set()
    for i, (key, p) in enumerate(pn.places.items()):
        pid = _ident(p.name or key, f"P_{i}")
        if pid in used:                       # guarantee uniqueness
            pid = f"{pid}_{i}"
        used.add(pid)
        place_id[key] = pid

    trans_id: Dict[str, str] = {}
    for i, (key, t) in enumerate(pn.transitions.items()):
        tid = _ident(t.name or key, f"T_{i}")
        if tid in used:
            tid = f"{tid}_{i}"
        used.add(tid)
        trans_id[key] = tid

    lines: List[str] = [f"\\beginnet{{{_ident(pn.name, 'net')}}}", ""]

    # ---- places ----
    for key, p in pn.places.items():
        body = f"\\name{{{place_id[key]}}}"
        cap = getattr(p, "capacity", None)
        if with_capacity and cap is not None:
            body += f"\\capacity{{{int(cap)}}}"
        body += f"\\init{{{int(p.tokens)}}}"
        lines.append(f"\\place{{{place_id[key]}}}{{{body}}}")
    lines.append("")

    # ---- transitions ----
    for key, t in pn.transitions.items():
        lines.append(f"\\transition{{{trans_id[key]}}}{{\\name{{{trans_id[key]}}}}}")
    lines.append("")

    # ---- arcs (single arc per logical arc; weight + type carried explicitly) ----
    aidx = 0

    def arc(src: str, tgt: str, w: int) -> str:
        nonlocal aidx
        s = (f"\\arc{{A_{aidx}}}{{\\from{{{src}}} \\to{{{tgt}}} "
             f"\\weight{{{int(w)}}} \\type{{ordinary}}}}")
        aidx += 1
        return s

    for key, t in pn.transitions.items():
        tid = trans_id[key]
        for a in t.inputs:                     # place -> transition
            lines.append(arc(place_id[a.place], tid, a.weight))
        for a in t.outputs:                    # transition -> place
            lines.append(arc(tid, place_id[a.place], a.weight))

    lines += ["", "", "\\endnet"]

    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write(_CRLF.join(lines) + _CRLF)

    return {"places": len(pn.places), "transitions": len(pn.transitions),
            "arcs": aidx}


def export_legend(pn: PetriNet, path: str) -> None:
    """Write a TSV mapping APNN ids <-> model names (decode raw Charlie output)."""
    rows = ["id\tkind\tname\tinit"]
    for i, (key, p) in enumerate(pn.places.items()):
        rows.append(f"{_ident(p.name or key, f'P_{i}')}\tplace\t{p.name}\t{int(p.tokens)}")
    for i, (key, t) in enumerate(pn.transitions.items()):
        rows.append(f"{_ident(t.name or key, f'T_{i}')}\ttransition\t{t.name}\t")
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write(_CRLF.join(rows) + _CRLF)


if __name__ == "__main__":
    from src.models import build_baseline, build_perturbation

    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "output", "NMDAR", "charlie")
    for tag, builder in [("baseline", build_baseline), ("perturbation", build_perturbation)]:
        pn = builder()
        counts = export_apnn(pn, os.path.join(out, f"NMDAR_{tag}.apnn"))
        export_legend(pn, os.path.join(out, f"NMDAR_{tag}_legend.tsv"))
        print(f"NMDAR_{tag}.apnn  ->  {counts}")


# ===========================================================================
# ANDL exporter (Abstract Net Description Language; Snoopy/Marcie/Charlie)
# ===========================================================================
# Unlike APNN, ANDL is Charlie's actively-maintained native textual format and
# carries arc weights *natively*: a transition's marking change is written as
#   <name> : <guard> : [P + w] & [P - w] & ... : <rate function> ;
# where [P + w] is an output arc of weight w and [P - w] an input arc of
# weight w.  A weight>1 arc (e.g. the bivalent crosslink, 2 NMDAR_Ab consumed)
# is therefore [NMDAR_Ab - 2] -- a genuine weight Charlie reads, so the
# perturbation correctly registers as non-ordinary.  Format mirrors a real
# Snoopy export (baseline.andl): `spn` header, `constants:` / `places:` /
# `transitions:` sections, CRLF endings, outputs listed before inputs.

def _andl_name(name: str) -> str:
    """Net-name slug: replace each non-word char with '_' (no collapsing),
    matching Snoopy (e.g. 'Baseline: Healthy ...' -> 'Baseline__Healthy_...')."""
    return _NON_ID.sub("_", str(name).strip())


def _andl_rate(r: float) -> str:
    s = f"{float(r):.10g}"
    if "." not in s and "e" not in s and "E" not in s:
        s += ".0"
    return s


def export_andl(pn: PetriNet, path: str, *, net_class: str = "spn") -> Dict[str, int]:
    r"""Write a Charlie-readable ANDL file for ``pn`` at ``path``.

    Arc weights are carried natively via the ``[Place +/- w]`` terms, so this is
    the format to use when any arc has weight > 1 (Charlie's APNN reader does
    not honour such weights and reports the net as ordinary).

    Returns element counts ({"places", "transitions", "arcs"}).
    """
    L: List[str] = [f"{net_class}  [{_andl_name(pn.name)}]", "{", "constants:", "places:"]

    for p in pn.places.values():
        L.append(f"  {p.name} = {int(p.tokens)};")
    L += ["", "transitions:"]

    narcs = 0
    for t in pn.transitions.values():
        # Snoopy convention: outputs (+) first, then inputs (-), in model order.
        terms = [f"[{a.place} + {int(a.weight)}]" for a in t.outputs]
        terms += [f"[{a.place} - {int(a.weight)}]" for a in t.inputs]
        narcs += len(t.outputs) + len(t.inputs)
        effect = " & ".join(terms)
        L.append(f"  {t.name}")
        L.append("    : ")                                   # guard (empty)
        L.append(f"    : {effect}")                          # marking change
        L.append(f"    : MassAction({_andl_rate(t.rate)})")  # rate function
        L.append("    ;")
    L += ["", "}"]

    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write(_CRLF.join(L))                               # no trailing newline
    return {"places": len(pn.places), "transitions": len(pn.transitions), "arcs": narcs}
