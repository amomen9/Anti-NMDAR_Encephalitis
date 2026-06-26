"""
Two Petri-Net models around the NMDA receptor:

  * Baseline      - healthy NMDAR signalling at a glutamatergic synapse
  * Perturbation  - anti-NMDAR encephalitis (IgG against GluN1)

Layout strategy
---------------
The diagram is laid out on a strict grid so that no place ever occupies the
same coordinates as a transition.

  * Compartment bands  : Y values   9 (pre-syn) -> -1 (endosomal), step 2
  * Transition lanes   : odd Y values 8, 6, 4, 2, 0 sit *between* bands
  * Place columns      : X[1..8] integer columns, spacing 2.6
  * Half columns       : X[i.5]   sit between place columns and host
                         "in-row" transitions (NMDAR gating chain)
  * Annotation_side    : each transition tells the renderer where its rate
                         label should sit (above / below / left / right)

Only species/processes immediately around the NMDAR are included.  Molecules
mentioned in the brief but not relevant to anti-NMDAR encephalitis are *not*
modelled - their omission, with rationale, is documented in README.md.
"""

from __future__ import annotations

from .petri_net import Arc, PetriNet, Place, Transition

# ---------------------------------------------------------------------------
# Layout grid
# ---------------------------------------------------------------------------
# Place columns
_X_STEP = 2.6
X = {i: 0.8 + (i - 1) * _X_STEP for i in range(1, 9)}
# Half columns (between integer columns)
H = {i + 0.5: 0.8 + (i - 0.5) * _X_STEP for i in range(1, 9)}

# Place rows (anatomical compartments) - spaced 2.5 apart so place radii
# (0.5) and transition boxes (h~0.95) never touch their neighbours.
Y_VESICLE = 10.0
Y_EXTRA   = 7.5
Y_MEMB    = 5.0
Y_PSD     = 2.5
Y_CYTO    = 0.0
Y_ENDO    = -2.5

# Transition lanes (mid-points between place rows)
YT_PRE_EX  = 8.75   # pre-syn <-> extra
YT_EX_MEM  = 6.25   # extra   <-> memb
YT_MEM_PSD = 3.75   # memb    <-> psd
YT_PSD_CY  = 1.25   # psd     <-> cyto
YT_CY_END  = -1.25  # cyto    <-> endo


# ---------------------------------------------------------------------------
# Mass-action rate constants (token-units / s).  Literature values are cited
# in README.md; the values here are scaled so the simulation runs on an 8-s
# window with synaptic events in the sub-second range and antibody-driven
# internalisation visible.
# ---------------------------------------------------------------------------
K_GLU_ON      = 4.0
K_GLU_OFF     = 4.0
K_OPEN        = 9.0
K_CLOSE       = 3.0
K_MG_BLOCK    = 8.0
K_MG_UNBLOCK  = 6.0
K_CA_INFLUX   = 6.0
K_CA_EFFLUX   = 1.2
K_PSD95_BIND  = 0.8
K_EPHB2_BIND  = 0.6
K_CAMKII_ON   = 2.0
K_CAMKII_OFF  = 0.4
K_LTP         = 0.3
K_GLU_RELEASE = 0.7
K_GLU_REUPTAKE = 1.5

K_AB_BIND       = 4.0
K_AB_UNBIND     = 0.02
K_CROSSLINK     = 2.5
K_EPHB2_DISRUPT = 2.5
K_ENDO          = 2.0
K_LYSO          = 0.6


# ===========================================================================
# Baseline (healthy NMDAR signalling)
# ===========================================================================
def build_baseline(include_psd95_bind: bool = False) -> PetriNet:
    pn = PetriNet(
        name="Baseline: Healthy NMDAR signalling",
        description=(
            "Glutamatergic synapse around a single NMDAR pool. "
            "Co-agonists (glycine / D-serine), voltage-dependent Mg2+ block, "
            "Ca2+ influx, PSD-95 scaffolding, EphB2 stabilisation and "
            "downstream CaMKII -> LTP cascade."
        ),
    )

    # ----------- Pre-synaptic terminal -----------
    pn.add_place(Place("Vesicle", "Pre-syn\nvesicle", tokens=8,
                       position=(X[1], Y_VESICLE), compartment="regulator",
                       annotation="Pre-synaptic Glu store"))

    # ----------- Extracellular reservoirs -----------
    pn.add_place(Place("Glu_ext", "Glutamate\n(cleft)", tokens=4,
                       position=(X[2], Y_EXTRA), compartment="extracellular",
                       annotation="NMDAR agonist"))
    pn.add_place(Place("Gly_ext", "Glycine /\nD-Serine", tokens=10,
                       position=(X[4], Y_EXTRA), compartment="extracellular",
                       annotation="GluN1 co-agonist"))
    pn.add_place(Place("Mg_ext", "Mg2+\n(cleft)", tokens=14,
                       position=(X[6], Y_EXTRA), compartment="extracellular",
                       annotation="Voltage-dep. blocker"))
    pn.add_place(Place("Ca_ext", "Ca2+\n(extracellular)", tokens=80,
                       position=(X[7], Y_EXTRA), compartment="extracellular",
                       annotation="~1.8 mM reservoir"))

    # ----------- Plasma membrane: NMDAR conformational states -----------
    pn.add_place(Place("NMDAR", "NMDAR\n(closed)", tokens=10,
                       position=(X[2], Y_MEMB), compartment="membrane",
                       annotation="GluN1/GluN2 tetramer"))
    pn.add_place(Place("NMDAR_GluGly", "NMDAR\n(Glu+Gly\nbound)", tokens=0,
                       position=(X[4], Y_MEMB), compartment="membrane",
                       annotation="Doubly liganded, still closed"))
    pn.add_place(Place("NMDAR_open", "NMDAR\n(open)", tokens=0,
                       position=(X[6], Y_MEMB), compartment="membrane",
                       annotation="Conducting state"))
    pn.add_place(Place("NMDAR_Mg", "NMDAR\n(Mg2+\nblocked)", tokens=0,
                       position=(X[7], Y_MEMB), compartment="membrane",
                       annotation="Pore blocked"))

    # ----------- PSD scaffold -----------
    if include_psd95_bind:
        pn.add_place(Place("PSD95", "PSD-95\n(free)", tokens=6,
                           position=(X[1], Y_PSD), compartment="scaffold",
                           annotation="PDZ scaffold"))
        pn.add_place(Place("NMDAR_PSD95", "NMDAR -\nPSD-95", tokens=0,
                           position=(X[2], Y_PSD), compartment="scaffold",
                           annotation="Synaptic anchor"))
    pn.add_place(Place("EphB2", "EphB2\n(surface)", tokens=5,
                       position=(X[3], Y_PSD), compartment="scaffold",
                       annotation="ECD cross-talk partner"))
    pn.add_place(Place("NMDAR_EphB2", "NMDAR -\nEphB2", tokens=0,
                       position=(X[4], Y_PSD), compartment="scaffold",
                       annotation="Lateral anchor (Mikasova 2012)"))

    # ----------- Cytosol -----------
    pn.add_place(Place("LTP", "LTP\nmarker", tokens=0,
                       position=(X[1], Y_CYTO), compartment="cytosol",
                       annotation="AMPAR insertion"))
    pn.add_place(Place("CaMKII_a", "CaMKII\n(active)", tokens=0,
                       position=(X[3], Y_CYTO), compartment="cytosol",
                       annotation="Ca/CaM-bound"))
    pn.add_place(Place("CaMKII_i", "CaMKII\n(inactive)", tokens=8,
                       position=(X[5], Y_CYTO), compartment="cytosol",
                       annotation="Holoenzyme, no CaM"))
    pn.add_place(Place("Ca_in", "Ca2+\n(cytosol)", tokens=2,
                       position=(X[7], Y_CYTO), compartment="cytosol",
                       annotation="Resting ~100 nM",
                       # label to the left so the area directly below stays
                       # clear for the attached-subnet block + connector.
                       label_side="left",
                       subnet_link="Block:\nNeurodegeneration\nPathway PN"))

    # ====================== TRANSITIONS ======================

    # --- Pre-syn / extracellular boundary (transition lane YT_PRE_EX = 8) ---
    pn.add_transition(Transition("T_release", "Vesicle\nrelease",
        inputs=[Arc("Vesicle", 1)], outputs=[Arc("Glu_ext", 1)],
        rate=K_GLU_RELEASE,
        position=(H[1.5] - 0.55, YT_PRE_EX),
        annotation="k ~ 0.7 /s [REF11]", annotation_side="below",
        reference="REF11"))

    pn.add_transition(Transition("T_reuptake", "EAAT\nreuptake",
        inputs=[Arc("Glu_ext", 1)], outputs=[Arc("Vesicle", 1)],
        rate=K_GLU_REUPTAKE,
        position=(H[1.5] + 0.55, YT_PRE_EX),
        annotation="EAAT clearance [REF12]", annotation_side="above",
        reference="REF12"))

    # --- NMDAR gating: forward in MEMB row, reverse in lane above (YT_EX_MEM) ---
    pn.add_transition(Transition("T_bind", "Glu+Gly\nbind",
        inputs=[Arc("NMDAR", 1), Arc("Glu_ext", 1), Arc("Gly_ext", 1)],
        outputs=[Arc("NMDAR_GluGly", 1)],
        rate=K_GLU_ON,
        position=(H[2.5], Y_MEMB),
        annotation="k_on ~ 5e6 /M/s [REF1]", annotation_side="below",
        reference="REF1"))

    pn.add_transition(Transition("T_unbind", "Glu+Gly\nunbind",
        inputs=[Arc("NMDAR_GluGly", 1)],
        outputs=[Arc("NMDAR", 1), Arc("Glu_ext", 1), Arc("Gly_ext", 1)],
        rate=K_GLU_OFF,
        position=(H[2.5], YT_EX_MEM),
        annotation="k_off ~ 4 /s [REF2]", annotation_side="above",
        reference="REF2"))

    pn.add_transition(Transition("T_open", "Channel\nopen",
        inputs=[Arc("NMDAR_GluGly", 1)],
        outputs=[Arc("NMDAR_open", 1)],
        rate=K_OPEN,
        position=(H[4.5], Y_MEMB),
        annotation="beta ~ 90 /s [REF2]", annotation_side="below",
        reference="REF2"))

    pn.add_transition(Transition("T_close", "Channel\nclose",
        inputs=[Arc("NMDAR_open", 1)],
        outputs=[Arc("NMDAR_GluGly", 1)],
        rate=K_CLOSE,
        position=(H[4.5], YT_EX_MEM),
        annotation="alpha ~ 30 /s [REF2]", annotation_side="above",
        reference="REF2"))

    pn.add_transition(Transition("T_Mg_block", "Mg2+\nblock",
        inputs=[Arc("NMDAR_open", 1), Arc("Mg_ext", 1)],
        outputs=[Arc("NMDAR_Mg", 1)],
        rate=K_MG_BLOCK,
        position=(H[6.5], Y_MEMB),
        annotation="tau < 200 us [REF3]", annotation_side="below",
        reference="REF3"))

    pn.add_transition(Transition("T_Mg_unblock", "Mg2+\nunblock",
        inputs=[Arc("NMDAR_Mg", 1)],
        outputs=[Arc("NMDAR_open", 1), Arc("Mg_ext", 1)],
        rate=K_MG_UNBLOCK,
        position=(H[6.5], YT_EX_MEM),
        annotation="voltage-gated [REF3]", annotation_side="above",
        reference="REF3"))

    # --- Membrane <-> PSD scaffolding (lane YT_MEM_PSD = 4) ---
    if include_psd95_bind:
        pn.add_transition(Transition("T_PSD95", "PSD-95\nbind",
            inputs=[Arc("NMDAR", 1), Arc("PSD95", 1)],
            outputs=[Arc("NMDAR_PSD95", 1)],
            rate=K_PSD95_BIND,
            position=(H[1.5], YT_MEM_PSD),
            annotation="PDZ anchor [REF5]", annotation_side="left",
            reference="REF5"))

    pn.add_transition(Transition("T_EphB2", "EphB2\nbind",
        inputs=[Arc("NMDAR", 1), Arc("EphB2", 1)],
        outputs=[Arc("NMDAR_EphB2", 1)],
        rate=K_EPHB2_BIND,
        position=(H[3.5], YT_MEM_PSD),
        annotation="ECD complex [REF6]", annotation_side="right",
        reference="REF6"))

    # --- Ca2+ influx & clearance (lane YT_PSD_CY = 2) ---
    pn.add_transition(Transition("T_Ca_in", "Ca2+\ninflux",
        inputs=[Arc("NMDAR_open", 1), Arc("Ca_ext", 1)],
        outputs=[Arc("NMDAR_open", 1), Arc("Ca_in", 1)],
        rate=K_CA_INFLUX,
        position=(H[6.5], YT_PSD_CY),
        annotation="g ~ 55 pS [REF4]", annotation_side="left",
        reference="REF4"))

    pn.add_transition(Transition("T_Ca_out", "Ca2+ efflux\n(PMCA/NCX)",
        inputs=[Arc("Ca_in", 1)], outputs=[Arc("Ca_ext", 1)],
        rate=K_CA_EFFLUX,
        position=(H[7.5] - 0.6, YT_PSD_CY),
        annotation="basal [REF4]", annotation_side="right",
        reference="REF4"))

    # --- CaMKII chain in cytosol row (forward and reverse side-by-side) ---
    pn.add_transition(Transition("T_CaMKII_on", "CaMKII\nactivate",
        inputs=[Arc("CaMKII_i", 1), Arc("Ca_in", 1)],
        outputs=[Arc("CaMKII_a", 1)],
        rate=K_CAMKII_ON,
        position=(X[4] - 0.55, Y_CYTO),
        annotation="tau_on ~ 0.1 s [REF7]", annotation_side="below",
        reference="REF7"))

    pn.add_transition(Transition("T_CaMKII_off", "CaMKII\ndeactivate",
        inputs=[Arc("CaMKII_a", 1)], outputs=[Arc("CaMKII_i", 1)],
        rate=K_CAMKII_OFF,
        position=(X[4] + 0.55, Y_CYTO),
        annotation="tau_off ~ 3 s [REF7]", annotation_side="above",
        reference="REF7"))

    pn.add_transition(Transition("T_LTP", "LTP\ninduction",
        inputs=[Arc("CaMKII_a", 1)],
        outputs=[Arc("LTP", 1), Arc("CaMKII_i", 1)],
        rate=K_LTP,
        position=(H[1.5], Y_CYTO),
        annotation="AMPAR insertion [REF7]", annotation_side="above",
        reference="REF7"))

    return pn


# ===========================================================================
# Perturbation (anti-NMDAR encephalitis)
# ===========================================================================
def build_perturbation(include_psd95_bind: bool = False) -> PetriNet:
    pn = build_baseline(include_psd95_bind=include_psd95_bind)
    pn.name = "Perturbation: anti-NMDAR encephalitis"
    pn.description = (
        "Patient IgG against the GluN1 amino-terminal domain binds surface NMDARs, "
        "is bivalently crosslinked, displaces the EphB2-NMDAR interaction "
        "(Mikasova 2012), and drives clathrin-mediated endocytosis -> lysosomal "
        "degradation (Hughes 2010, Moscato 2014). Net effect: progressive loss "
        "of synaptic NMDARs and a fall in Ca2+ influx."
    )

    # --- Disease-specific places, all on column X[8] (right lane) ---
    pn.add_place(Place("Ab", "Anti-GluN1\nIgG", tokens=18,
                       position=(X[8], Y_EXTRA), compartment="regulator",
                       annotation="Patient autoantibody"))
    pn.add_place(Place("NMDAR_Ab", "NMDAR -\nIgG", tokens=0,
                       position=(X[8], Y_MEMB), compartment="membrane",
                       annotation="Single Ab bound",
                       label_side="right"))   # below-label would hide under T_crosslink
    pn.add_place(Place("NMDAR_Ab2", "Crosslinked\nNMDAR-IgG2", tokens=0,
                       position=(X[8], Y_PSD), compartment="membrane",
                       annotation="Bivalent cap (Hughes 2010)",
                       label_side="right"))   # below-label would hide under T_endocytosis
    pn.add_place(Place("Endosome", "Internalised\n(endosome)", tokens=0,
                       position=(X[8], Y_CYTO), compartment="endosome",
                       annotation="Clathrin endocytosis",
                       label_side="right"))   # below-label would hide under T_lyso
    pn.add_place(Place("Lysosome", "Degraded\n(lysosome)", tokens=0,
                       position=(X[8], Y_ENDO), compartment="endosome",
                       annotation="Receptor pool lost"))

    # --- Antibody binding (membrane row, between NMDAR_Mg and NMDAR_Ab) ---
    pn.add_transition(Transition("T_Ab_bind", "IgG\nbinding",
        inputs=[Arc("NMDAR", 1), Arc("Ab", 1)],
        outputs=[Arc("NMDAR_Ab", 1)],
        rate=K_AB_BIND,
        position=(H[7.5], Y_MEMB),
        annotation="k_on ~ 1e5 /M/s [REF8]", annotation_side="below",
        reference="REF8"))

    pn.add_transition(Transition("T_Ab_unbind", "IgG\nunbind",
        inputs=[Arc("NMDAR_Ab", 1)],
        outputs=[Arc("NMDAR", 1), Arc("Ab", 1)],
        rate=K_AB_UNBIND,
        position=(H[7.5], YT_EX_MEM),
        annotation="slow [REF8]", annotation_side="above",
        reference="REF8"))

    # --- Bivalent crosslink (memb-PSD lane, right side) ---
    pn.add_transition(Transition("T_crosslink", "Bivalent\ncrosslink",
        inputs=[Arc("NMDAR_Ab", 2)],
        outputs=[Arc("NMDAR_Ab2", 1)],
        rate=K_CROSSLINK,
        position=(X[8], YT_MEM_PSD),
        annotation="2x NMDAR-IgG -> cap [REF9]", annotation_side="right",
        reference="REF9"))

    # --- EphB2 disruption (lane YT_MEM_PSD, between EphB2 stuff and Ab column) ---
    pn.add_transition(Transition("T_EphB2_disrupt", "EphB2\ndisruption",
        inputs=[Arc("NMDAR_EphB2", 1), Arc("Ab", 1)],
        outputs=[Arc("NMDAR_Ab", 1), Arc("EphB2", 1)],
        rate=K_EPHB2_DISRUPT,
        position=(H[5.5], YT_MEM_PSD),
        annotation="surface cross-talk loss [REF6]", annotation_side="below",
        reference="REF6"))

    # --- Endocytosis (PSD-cyto lane, right column) ---
    pn.add_transition(Transition("T_endocytosis", "Clathrin\nendocytosis",
        inputs=[Arc("NMDAR_Ab2", 1)], outputs=[Arc("Endosome", 1)],
        rate=K_ENDO,
        position=(X[8], YT_PSD_CY),
        annotation="t1/2 ~ 5 h [REF9]", annotation_side="right",
        reference="REF9"))

    # --- Lysosomal degradation (cyto-endo lane) ---
    pn.add_transition(Transition("T_lyso", "Lysosomal\ndegradation",
        inputs=[Arc("Endosome", 1)], outputs=[Arc("Lysosome", 1)],
        rate=K_LYSO,
        position=(X[8], YT_CY_END),
        annotation="receptor lost [REF10]", annotation_side="right",
        reference="REF10"))

    return pn
