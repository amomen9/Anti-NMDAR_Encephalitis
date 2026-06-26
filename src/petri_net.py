"""
Lightweight Petri Net core used by the anti-NMDAR encephalitis model.

Supports:
  - Discrete tokens with integer multiplicities
  - Stoichiometric arcs (input + output weights)
  - Mass-action propensities for Gillespie SSA
  - Compartment metadata for layered visualisation
  - Rich annotations (rate constants, citations) carried on transitions/arcs

The class is intentionally self-contained so the visualisation and simulation
layers can introspect every element. SNAKES is referenced in README.md as the
formal modelling counterpart; we keep our own implementation here for fine-grain
control over the matplotlib renderer and the GIF animator.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


# Visualisation palette grouped by compartment
COMPARTMENT_COLORS = {
    "extracellular": "#bde0fe",   # light blue (synaptic cleft)
    "membrane":      "#cdb4db",   # mauve (plasma membrane)
    "scaffold":      "#ffd6a5",   # peach (post-synaptic density)
    "cytosol":       "#caffbf",   # light green (cytoplasm)
    "endosome":      "#ffadad",   # salmon (endo/lysosomal)
    "regulator":     "#fdffb6",   # light yellow
}


@dataclass
class Place:
    name: str
    label: str
    tokens: int
    position: Tuple[float, float]
    compartment: str = "cytosol"
    capacity: Optional[int] = None
    annotation: str = ""
    label_side: str = "below"   # below | above | left | right
    # Hierarchical link: when set, the renderer draws a Snoopy-style "coarse
    # node" block beside this place, signalling that a named subnet is logically
    # attached here but is rendered in a separate diagram (not expanded inline).
    subnet_link: str = ""
    subnet_offset: Tuple[float, float] = (0.0, -2.3)  # block centre relative to place

    def __post_init__(self) -> None:
        self.initial_tokens = self.tokens


@dataclass
class Arc:
    place: str
    weight: int = 1
    annotation: str = ""        # short on-arc text (e.g. rate constant)


@dataclass
class Transition:
    name: str
    label: str
    inputs: List[Arc]
    outputs: List[Arc]
    rate: float                 # mass-action coefficient (s^-1 in stoichiometric units)
    position: Tuple[float, float]
    kind: str = "mass-action"   # mass-action | michaelis | source | sink
    annotation: str = ""
    reference: str = ""         # citation tag, expanded in README.md
    annotation_side: str = "below"  # below | above | left | right


class PetriNet:
    def __init__(self, name: str, description: str = "") -> None:
        self.name = name
        self.description = description
        self.places: Dict[str, Place] = {}
        self.transitions: Dict[str, Transition] = {}
        self.history: List[Tuple[float, str]] = []   # (time, transition)

    # ------------------------------------------------------------------ build
    def add_place(self, place: Place) -> Place:
        self.places[place.name] = place
        return place

    def add_transition(self, transition: Transition) -> Transition:
        for arc in transition.inputs + transition.outputs:
            if arc.place not in self.places:
                raise KeyError(f"Arc references missing place '{arc.place}'")
        self.transitions[transition.name] = transition
        return transition

    # ------------------------------------------------------------------ state
    def marking(self) -> Dict[str, int]:
        return {n: p.tokens for n, p in self.places.items()}

    def reset(self) -> None:
        for p in self.places.values():
            p.tokens = p.initial_tokens
        self.history.clear()

    def is_enabled(self, transition_name: str) -> bool:
        t = self.transitions[transition_name]
        for arc in t.inputs:
            if self.places[arc.place].tokens < arc.weight:
                return False
        for arc in t.outputs:                       # capacity check
            cap = self.places[arc.place].capacity
            if cap is not None and self.places[arc.place].tokens + arc.weight > cap:
                return False
        return True

    def fire(self, transition_name: str, time: Optional[float] = None) -> bool:
        if not self.is_enabled(transition_name):
            return False
        t = self.transitions[transition_name]
        for arc in t.inputs:
            self.places[arc.place].tokens -= arc.weight
        for arc in t.outputs:
            self.places[arc.place].tokens += arc.weight
        if time is not None:
            self.history.append((time, transition_name))
        return True

    # ----------------------------------------------------------- propensities
    def propensity(self, transition_name: str) -> float:
        """Stochastic mass-action propensity: rate * \\prod binom(#tokens, weight)."""
        t = self.transitions[transition_name]
        a = float(t.rate)
        for arc in t.inputs:
            n = self.places[arc.place].tokens
            if n < arc.weight:
                return 0.0
            # combinatorial term for higher-order reactions
            for i in range(arc.weight):
                a *= (n - i) / (i + 1)
        return max(a, 0.0)

    def propensities(self) -> Dict[str, float]:
        return {n: self.propensity(n) for n in self.transitions}
