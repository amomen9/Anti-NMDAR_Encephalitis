"""
Downstream Petri Net models: Calcium -> apoptosis / neurodegeneration.

This module is a standalone extension of the NMDAR ensephalitis model.  It
defines two functions that build self-contained PetriNet objects (one for
the baseline (healthy neuron) and one for the perturbed (anti-NMDAR
encephalitis) state) covering the pathways from cytosolic Calcium through to
CASP3-mediated apoptosis and neurodegeneration output.

Perturbation versus baseline
-----------------------------
Baseline   : Calcium starts at steady resting level (2 tokens); downstream
             pathways are activated but at low basal rates.
Perturbation: Calcium input is reduced (reflects depleted NMDAR pool), but
             because Calcium is chronically perturbed, ROS and ER-stress
             transitions are up-regulated to model the hypo-Calcium 
             maladaptive response (Hughes 2010; Moscato 2014).

Layout strategy
---------------
The diagram is placed on a new coordinate grid that is designed to sit
*below* or *beside* the existing NMDAR model when the two are merged.

  X axis: 5 columns, spacing 2.4
  Y axis:
       3.0  Cytosol interface  (Calcium input comes in here from upstream)
       1.0  Mitochondrial / ROS branch
      -1.0  ER stress branch
      -3.0  CASP cascade row
      -5.0  Output row (Apoptosis / Neurodegeneration)
"""

from __future__ import annotations

from .petri_net import Arc, PetriNet, Place, Transition  # adjust import for your package layout

# ---------------------------------------------------------------------------
# Layout grid
# ---------------------------------------------------------------------------
_DX = 2.4
X = {i: 0.8 + (i - 1) * _DX for i in range(1, 7)}   # 6 columns
H = {i + 0.5: 0.8 + (i - 0.5) * _DX for i in range(1, 7)}

# Y rows
Y_CYTO  =  3.0   # cytosol interface: Ca_in lives here
Y_MITO  =  1.0   # mitochondrial / ROS branch
Y_ER    = -1.0   # ER stress / CHOP branch
Y_CASP  = -3.0   # CASP9 / CASP3 row
Y_OUT   = -5.0   # apoptosis / neurodegeneration output

# Transition lanes
YT_CYTO_MITO  =  2.0
YT_CYTO_ER    =  1.0   # shared lane; transitions offset in X
YT_MITO_CASP  =  0.0   # mito -> CASP9
YT_ER_CASP    = -2.0   # CHOP -> CASP3
YT_CASP_OUT   = -4.0   # CASP3 -> Apoptosis

# ---------------------------------------------------------------------------
# Rate constants
# Reference: Hardy et al., 2023; Gleichmann & Bhatt 2011; scaled for 8 s sim
# ---------------------------------------------------------------------------

# --- Baseline (resting / low-activity) ---
KB_CA_MITO      = 0.30   # Ca2+ uptake by mitochondria
KB_CYTO_C       = 0.12   # cytochrome-c release on loss
KB_CASP9_ON     = 0.20   # CASP9 activation by cyt-c
KB_CASP9_OFF    = 0.15   # CASP9 inactivation
KB_CASP3_ON     = 0.18   # CASP3 activation by CASP9
KB_CASP3_OFF    = 0.10   # CASP3 basal degradation
KB_ROS_ON       = 0.25   # ROS production (NADPH ox / mito ROS)
KB_ROS_OFF      = 0.40   # ROS scavenging (SOD/catalase)
KB_ROS_CASP3    = 0.08   # ROS-driven CASP3 activation
KB_ER_STRESS    = 0.20   # ER stress from Ca2+ perturbation
KB_CHOP_ON      = 0.15   # CHOP induction by ER stress
KB_CHOP_OFF     = 0.12   # CHOP degradation
KB_CHOP_CASP3   = 0.10   # CHOP-driven CASP3 activation
KB_APOPTOSIS    = 0.22   # CASP3 -> apoptosis commitment
KB_NEURODEGEN   = 0.06   # apoptosis -> neurodegeneration accumulation

# --- Perturbation adjustments (anti-NMDAR, hypo-Ca2+ chronic state) ---
#  NMDAR depletion = Ca2+ influx is reduced but ER stress and ROS are elevated 
KP_CA_MITO      = 0.15   # reduced uptake due to low Ca2+
KP_CYTO_C       = 0.18   # slightly elevated: mito already stressed
KP_CASP9_ON     = 0.20
KP_CASP9_OFF    = 0.12   # slower recovery
KP_CASP3_ON     = 0.22   # elevated
KP_CASP3_OFF    = 0.08
KP_ROS_ON       = 0.55   # markedly elevated: maladaptive ROS burst
KP_ROS_OFF      = 0.28   # impaired scavenging
KP_ROS_CASP3    = 0.20   # elevated
KP_ER_STRESS    = 0.50   # elevated: chronic ER Ca2+ depletion
KP_CHOP_ON      = 0.40   # elevated CHOP
KP_CHOP_OFF     = 0.10
KP_CHOP_CASP3   = 0.25   # elevated CHOP-driven apoptosis
KP_APOPTOSIS    = 0.35   # faster commitment
KP_NEURODEGEN   = 0.18   # faster neurodegeneration accumulation


# ===========================================================================
# Shared builder
# ===========================================================================
def _build(name: str, description: str,
           ca_tokens: int,
           k_ca_mito, k_cytoc, k_casp9_on, k_casp9_off,
           k_casp3_on, k_casp3_off,
           k_ros_on, k_ros_off, k_ros_casp3,
           k_er_stress, k_chop_on, k_chop_off, k_chop_casp3,
           k_apop, k_neurodegen) -> PetriNet:

    pn = PetriNet(name=name, description=description)

    # ================================================================
    # PLACES
    # ================================================================

    # --- Cytosol interface ---
    pn.add_place(Place("Ca_in", "Ca²⁺\n(cytosol)", tokens=ca_tokens,
                       position=(X[4], Y_CYTO), compartment="cytosol",
                       annotation="Input from NMDAR model",
                       subnet_link="Block:\nNMDAR subnet",
                       subnet_offset=(0.0, 2.3)))

    # --- Mitochondrial branch ---
    pn.add_place(Place("Ca_mito", "Ca²⁺\n(mitochondria)", tokens=0,
                       position=(X[2], Y_MITO), compartment="endosome",
                       annotation="Mito Ca²⁺ overload"))
    pn.add_place(Place("DeltaPsi",  "DeltaPsi\n(intact)", tokens=6,
                       position=(X[3], Y_MITO), compartment="endosome",
                       annotation="Mito membrane potential",
                       capacity=6))
    pn.add_place(Place("CytoC", "Cytochrome-c\n(released)", tokens=0,
                       position=(X[4], Y_MITO), compartment="cytosol",
                       annotation="Apoptosome trigger"))

    # --- ROS branch ---
    pn.add_place(Place("ROS", "ROS", tokens=0,
                       position=(X[5], Y_MITO), compartment="cytosol",
                       annotation="Reactive oxygen species"))

    # --- ER stress branch ---
    pn.add_place(Place("ER_stress", "ER stress\n(PERK/eIF2alpha)", tokens=0,
                       position=(X[2], Y_ER), compartment="cytosol",
                       annotation="Unfolded protein response"))
    pn.add_place(Place("CHOP", "CHOP\n(DDIT3)", tokens=0,
                       position=(X[3], Y_ER), compartment="cytosol",
                       annotation="Pro-apoptotic TF"))

    # --- CASP cascade ---
    pn.add_place(Place("CASP9", "CASP9\n(active)", tokens=0,
                       position=(X[3], Y_CASP), compartment="cytosol",
                       annotation="Initiator caspase"))
    pn.add_place(Place("CASP3", "CASP3\n(active)", tokens=0,
                       position=(X[4], Y_CASP), compartment="cytosol",
                       annotation="Executioner caspase"))

    # --- Outputs ---
    pn.add_place(Place("Apoptosis", "Apoptosis\n(committed)", tokens=0,
                       position=(X[3], Y_OUT), compartment="endosome",
                       annotation="Cell death signal"))
    pn.add_place(Place("Neurodegeneration", "Neuro-\ndegeneration", tokens=0,
                       position=(X[4], Y_OUT), compartment="endosome",
                       annotation="Cumulative damage marker"))

    # ================================================================
    # TRANSITIONS
    # ================================================================

    # --- Mitochondrial branch ---
    # Ca2+ uptake by mitochondria is the key initiating event that triggers both the cytochrome-c release and ROS production branches
    # so it is placed first in the diagram
    pn.add_transition(Transition(
        "T_Ca_mito", "Ca2+\nmito uptake",
        inputs=[Arc("Ca_in", 1)],
        outputs=[Arc("Ca_mito", 1)],
        rate=k_ca_mito,
        position=(H[2.5], YT_CYTO_MITO),
        annotation=f"k={k_ca_mito:.2f} [REF13]",
        annotation_side="above",
        reference="REF13",
    ))

    # Ca2+ mito overload leads to loss of membrane potential (DeltaPsi) and release of cytochrome-c
    # which are modelled as two separate transitions that both consume DeltaPsi   
    pn.add_transition(Transition(
        "T_DeltaPsi_loss", "DeltaPsi\nloss",
        inputs=[Arc("Ca_mito", 1), Arc("DeltaPsi", 1)],
        outputs=[Arc("Ca_mito", 1)],   # Ca_mito persists, DeltaPsi consumed
        rate=k_cytoc,
        position=(H[2.5], Y_MITO),
        annotation=f"k={k_cytoc:.2f} [REF14]",
        annotation_side="below",
        reference="REF14",
    ))

    # Cytochrome-c release is modelled as a separate transition that also consumes Ca_mito (reflecting the overload trigger) 
    # but does not consume DeltaPsi (reflecting the fact that cytochrome-c release can occur even with partial loss of membrane potential)
    pn.add_transition(Transition(
        "T_CytoC_release", "Cyt-c\nrelease",
        inputs=[Arc("Ca_mito", 1)],
        outputs=[Arc("CytoC", 1)],
        rate=k_cytoc,
        position=(H[3.5], YT_CYTO_MITO),
        annotation=f"k={k_cytoc:.2f} [REF14]",
        annotation_side="above",
        reference="REF14",
    ))

    # --- ROS branch ---
    # ROS production is modelled as a single transition driven by cytosolic Ca2+  
    pn.add_transition(Transition(
        "T_ROS_on", "ROS\nproduction",
        inputs=[Arc("Ca_in", 1)],
        outputs=[Arc("Ca_in", 1), Arc("ROS", 1)],  # Ca_in catalytic
        rate=k_ros_on,
        position=(H[4.5], YT_CYTO_MITO),
        annotation=f"k={k_ros_on:.2f} [REF15]",
        annotation_side="above",
        reference="REF15",
    ))

    # ROS scavenging is modelled as a single transition that consumes ROS
    pn.add_transition(Transition(
        "T_ROS_off", "ROS\nscavenging",
        inputs=[Arc("ROS", 1)],
        outputs=[],
        rate=k_ros_off,
        position=(X[5], YT_CYTO_MITO),
        annotation=f"k={k_ros_off:.2f} [REF15]",
        annotation_side="right",
        reference="REF15",
    ))

    # --- ER stress branch ---
    # ER stress induction is modelled as a single transition driven by cytosolic Ca2+ perturbation that produces ER stress as an output
    pn.add_transition(Transition(
        "T_ER_stress", "ER stress\ninduction",
        inputs=[Arc("Ca_in", 1)],
        outputs=[Arc("Ca_in", 1), Arc("ER_stress", 1)],  # Ca_in catalytic
        rate=k_er_stress,
        position=(H[1.5], YT_CYTO_ER),
        annotation=f"k={k_er_stress:.2f} [REF16]",
        annotation_side="left",
        reference="REF16",
    ))

    # CHOP induction is modelled as a single transition that consumes ER stress and produces CHOP
    pn.add_transition(Transition(
        "T_CHOP_on", "CHOP\ninduction",
        inputs=[Arc("ER_stress", 1)],
        outputs=[Arc("CHOP", 1)],
        rate=k_chop_on,
        position=(H[2.5], Y_ER),
        annotation=f"k={k_chop_on:.2f} [REF16]",
        annotation_side="below",
        reference="REF16",
    ))

    # CHOP degradation is modelled as a single transition that consumes CHOP
    pn.add_transition(Transition(
        "T_CHOP_off", "CHOP\ndegradation",
        inputs=[Arc("CHOP", 1)],
        outputs=[],
        rate=k_chop_off,
        position=(H[2.5], YT_ER_CASP),
        annotation=f"k={k_chop_off:.2f} [REF16]",
        annotation_side="left",
        reference="REF16",
    ))

    # --- CASP9 (intrinsic pathway, cyt-c driven) ---
    # CASP9 activation is modelled as a single transition driven by cytosolic cytochrome-c that produces active CASP9 as an output
    pn.add_transition(Transition(
        "T_CASP9_on", "CASP9\nactivation",
        inputs=[Arc("CytoC", 1)],
        outputs=[Arc("CASP9", 1)],
        rate=k_casp9_on,
        position=(H[3.5], YT_MITO_CASP),
        annotation=f"k={k_casp9_on:.2f} [REF17]",
        annotation_side="right",
        reference="REF17",
    ))

    # CASP9 inactivation is modelled as a single transition that consumes active CASP9
    pn.add_transition(Transition(
        "T_CASP9_off", "CASP9\ninactivation",
        inputs=[Arc("CASP9", 1)],
        outputs=[],
        rate=k_casp9_off,
        position=(H[3.5], YT_MITO_CASP - 0.6),
        annotation=f"k={k_casp9_off:.2f} [REF17]",
        annotation_side="left",
        reference="REF17",
    ))

    # --- CASP3: three activation sources ---
    # CASP3 activation is modelled as three separate transitions driven by 
    # (1) active CASP9
    # (2) ROS
    # (3) CHOP
    # reflecting the fact that all three can independently drive CASP3 activation and apoptosis in neurons
    pn.add_transition(Transition(
        "T_CASP3_from_CASP9", "CASP3 <-\nCASP9",
        inputs=[Arc("CASP9", 1)],
        outputs=[Arc("CASP9", 1), Arc("CASP3", 1)],  # CASP9 catalytic
        rate=k_casp3_on,
        position=(H[3.5], YT_ER_CASP),
        annotation=f"k={k_casp3_on:.2f} [REF17]",
        annotation_side="right",
        reference="REF17",
    ))

    # CASP3 activation from CHOP 
    pn.add_transition(Transition(
        "T_CASP3_from_CHOP", "CASP3 <-\nCHOP",
        inputs=[Arc("CHOP", 1)],
        outputs=[Arc("CHOP", 1), Arc("CASP3", 1)],   # CHOP catalytic
        rate=k_chop_casp3,
        position=(H[2.5], YT_ER_CASP),
        annotation=f"k={k_chop_casp3:.2f} [REF16]",
        annotation_side="left",
        reference="REF16",
    ))

    # CASP3 activation from ROS
    pn.add_transition(Transition(
        "T_CASP3_from_ROS", "CASP3 <-\nROS",
        inputs=[Arc("ROS", 1)],
        outputs=[Arc("ROS", 1), Arc("CASP3", 1)],    # ROS catalytic
        rate=k_ros_casp3,
        position=(H[4.5], YT_ER_CASP),
        annotation=f"k={k_ros_casp3:.2f} [REF15]",
        annotation_side="right",
        reference="REF15",
    ))

    # CASP3 basal degradation is modelled as a single transition that consumes active CASP3
    pn.add_transition(Transition(
        "T_CASP3_off", "CASP3\nclearance",
        inputs=[Arc("CASP3", 1)],
        outputs=[],
        rate=k_casp3_off,
        position=(X[4], YT_ER_CASP - 0.4),
        annotation=f"k={k_casp3_off:.2f}",
        annotation_side="right",
    ))

    # --- Output transitions ---
    pn.add_transition(Transition(
        "T_Apoptosis", "Apoptosis\ncommitment",
        inputs=[Arc("CASP3", 1)],
        outputs=[Arc("Apoptosis", 1)],
        rate=k_apop,
        position=(H[3.5], YT_CASP_OUT),
        annotation=f"k={k_apop:.2f} [REF17]",
        annotation_side="left",
        reference="REF17",
    ))

    pn.add_transition(Transition(
        "T_Neurodegeneration", "Neuro-\ndegeneration",
        inputs=[Arc("Apoptosis", 1)],
        outputs=[Arc("Apoptosis", 1), Arc("Neurodegeneration", 1)],  # cumulative
        rate=k_neurodegen,
        position=(H[4.5], YT_CASP_OUT),
        annotation=f"k={k_neurodegen:.2f} [REF18]",
        annotation_side="right",
        reference="REF18",
    ))

    return pn


# ===========================================================================
# Public builders
# ===========================================================================

def build_downstream_baseline() -> PetriNet:
    """Healthy neuron downstream signalling from resting Ca2+ levels."""
    return _build(
        name="Downstream Baseline: Ca2+ -> CASP3 / Neurodegeneration",
        description=(
            "Downstream consequences of resting cytosolic Ca2+ in a healthy "
            "glutamatergic neuron.  Three parallel pro-apoptotic inputs converge "
            "on CASP3: (i) mitochondrial pathway via cytochrome-c and CASP9, "
            "(ii) ROS production, (iii) ER stress / CHOP.  At baseline Ca2+ "
            "levels all three branches are activated at low rates consistent "
            "with a live, non-apoptotic cell."
        ),
        ca_tokens=35,
        k_ca_mito=KB_CA_MITO,    k_cytoc=KB_CYTO_C,
        k_casp9_on=KB_CASP9_ON,  k_casp9_off=KB_CASP9_OFF,
        k_casp3_on=KB_CASP3_ON,  k_casp3_off=KB_CASP3_OFF,
        k_ros_on=KB_ROS_ON,      k_ros_off=KB_ROS_OFF,
        k_ros_casp3=KB_ROS_CASP3,
        k_er_stress=KB_ER_STRESS,
        k_chop_on=KB_CHOP_ON,    k_chop_off=KB_CHOP_OFF,
        k_chop_casp3=KB_CHOP_CASP3,
        k_apop=KB_APOPTOSIS,     k_neurodegen=KB_NEURODEGEN,
    )


def build_downstream_perturbation() -> PetriNet:
    """Anti-NMDAR encephalitis downstream: reduced Ca2+ but elevated stress."""
    return _build(
        name="Downstream Perturbation: anti-NMDAR Ca2+ -> CASP3 / Neurodegeneration",
        description=(
            "Downstream consequences in an anti-NMDAR encephalitis neuron.  "
            "NMDAR depletion reduces Ca2+ influx (ca_tokens=1) but chronically "
            "elevates ROS production and ER stress through maladaptive "
            "compensatory mechanisms leading to "
            "disproportionate CASP3 accumulation and neurodegeneration output "
            "relative to the reduced Ca2+ signal."
        ),
        ca_tokens=15,             # reduced pool (perturbed NMDAR model output)
        k_ca_mito=KP_CA_MITO,   k_cytoc=KP_CYTO_C,
        k_casp9_on=KP_CASP9_ON, k_casp9_off=KP_CASP9_OFF,
        k_casp3_on=KP_CASP3_ON, k_casp3_off=KP_CASP3_OFF,
        k_ros_on=KP_ROS_ON,     k_ros_off=KP_ROS_OFF,
        k_ros_casp3=KP_ROS_CASP3,
        k_er_stress=KP_ER_STRESS,
        k_chop_on=KP_CHOP_ON,   k_chop_off=KP_CHOP_OFF,
        k_chop_casp3=KP_CHOP_CASP3,
        k_apop=KP_APOPTOSIS,    k_neurodegen=KP_NEURODEGEN,
    )


# ===========================================================================
# Quick self-test  (python downstream_models.py)
# ===========================================================================
if __name__ == "__main__":
    for builder, label in [
        (build_downstream_baseline,     "Baseline"),
        (build_downstream_perturbation, "Perturbation"),
    ]:
        pn = builder()
        print(f"\n{'='*60}")
        print(f"  {label}: {pn.name}")
        print(f"  Places      : {len(pn.places)}")
        print(f"  Transitions : {len(pn.transitions)}")
        print(f"  Initial marking:")
        for n, p in pn.places.items():
            if p.tokens > 0:
                print(f"    {n:25s}: {p.tokens}")
        print(f"  Enabled transitions at t=0:")
        for t in pn.transitions:
            if pn.is_enabled(t):
                print(f"    {t}")
