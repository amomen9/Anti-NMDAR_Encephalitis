"""
Stochastic (Gillespie SSA) simulation for a `PetriNet`.

Returned trace gives a sequence of (time, transition, marking-snapshot) tuples
suitable for animation and quantitative comparison of Ca2+ influx.
"""

from __future__ import annotations

import math
from typing import List, Tuple, Dict

import numpy as np

from .petri_net import PetriNet


def gillespie(
    pn: PetriNet,
    t_max: float,
    rng: np.random.Generator,
    max_steps: int = 100_000,
) -> List[Tuple[float, str, Dict[str, int]]]:
    """Standard direct-method SSA.

    Records the marking *after* each firing.  When all propensities vanish the
    simulation halts early and pads the trace with a final 'no-op' record so
    that callers can still slice by `t_max`.
    """
    pn.reset()
    t = 0.0
    trace: List[Tuple[float, str, Dict[str, int]]] = [
        (t, "<init>", pn.marking())
    ]

    for _ in range(max_steps):
        props = pn.propensities()
        a0 = sum(props.values())
        if a0 <= 0:
            break

        r1, r2 = rng.random(), rng.random()
        dt = -math.log(r1) / a0
        t += dt
        if t > t_max:
            break

        # choose transition
        threshold = r2 * a0
        running = 0.0
        chosen = None
        for name, a in props.items():
            running += a
            if running >= threshold:
                chosen = name
                break
        if chosen is None:                          # numerical drift safety
            chosen = max(props, key=props.get)

        pn.fire(chosen, time=t)
        trace.append((t, chosen, pn.marking()))

    trace.append((t_max, "<end>", pn.marking()))
    return trace


def sample_trace(
    trace: List[Tuple[float, str, Dict[str, int]]],
    times: np.ndarray,
    place: str,
) -> np.ndarray:
    """Left-continuous step interpolation of a place's token count."""
    out = np.zeros_like(times, dtype=float)
    j = 0
    last = trace[0][2][place]
    for i, ti in enumerate(times):
        while j < len(trace) and trace[j][0] <= ti:
            last = trace[j][2][place]
            j += 1
        out[i] = last
    return out


def cumulative_ca_influx(
    trace: List[Tuple[float, str, Dict[str, int]]],
    influx_transition: str = "T_Ca_in",
) -> Tuple[np.ndarray, np.ndarray]:
    """Cumulative count of Ca2+ ions that have passed through open NMDAR."""
    ts, cum = [0.0], [0.0]
    total = 0
    for time, tr, _ in trace[1:]:
        if tr == influx_transition:
            total += 1
        ts.append(time)
        cum.append(total)
    return np.asarray(ts), np.asarray(cum)
