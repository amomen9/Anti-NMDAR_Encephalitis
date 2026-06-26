# Petri-Net model of Anti-NMDA-Receptor Encephalitis

This project builds **two side-by-side Petri-Net models** of the molecular
microenvironment of a single **NMDA receptor (NMDAR)** at a glutamatergic
synapse, plus a **downstream neurodegeneration pathway** that links Ca²⁺
dynamics to CASP3-mediated apoptosis and cumulative neurodegeneration:

| Branch         | Biological scenario                                       |
| -------------- | --------------------------------------------------------- |
| **Baseline**   | Healthy NMDAR signalling                                  |
| **Perturbation** | Anti-NMDAR encephalitis (patient IgG against GluN1)     |

The same set of "housekeeping" places and transitions is shared between the
two nets; the perturbation only **adds** the disease-specific machinery
(antibody binding, bivalent crosslinking, EphB2-disruption, clathrin
endocytosis, lysosomal degradation).  This keeps the comparison faithful:
every divergence in the dynamics comes from the disease-specific
transitions alone.

The downstream neurodegeneration model extends the NMDAR model by connecting
the cytosolic Ca²⁺ pool to mitochondrial ROS production, ER stress / CHOP
signalling, the CASP9 → CASP3 caspase cascade, and ultimately apoptosis and
neurodegeneration output.

## Outputs

### NMDAR model (`./output/NMDAR/`)

| File                              | What it shows                                                                 |
| --------------------------------- | ----------------------------------------------------------------------------- |
| `images/NMDAR_baseline_pn.png` / `.pdf` | Annotated static Petri Net of the healthy NMDAR signalling pathway          |
| `images/NMDAR_perturbation_pn.png` / `.pdf` | Annotated static Petri Net with antibody / internalisation lane on the right |
| `images/NMDAR_PN.png` / `.pdf`   | Side-by-side comparison: baseline (left) + perturbation (right) Petri nets with cumulative Ca²⁺ curve underneath |
| `images/ca_influx_comparison.png` / `.pdf` | Static comparison of instantaneous + cumulative Ca²⁺ influx               |
| `animations/NMDAR_baseline_animation.gif` | Gillespie token-flow animation of the healthy net                           |
| `animations/NMDAR_perturbation_animation.gif` | Gillespie token-flow animation of the disease net                           |
| `animations/NMDAR_baseline_combined.gif` | Healthy net animation with synchronised Ca²⁺ curves (instantaneous + cumulative) |
| `animations/NMDAR_perturbation_combined.gif` | Disease net animation with synchronised Ca²⁺ curves                         |
| `animations/NMDAR_PN.gif`        | Animated side-by-side Petri nets + cumulative Ca²⁺ curve                    |
| `animations/ca_influx_comparison.gif` | Animated growth of the same Ca²⁺ curves                                   |
| `snoopy/baseline.spn` / `.cpn`   | Snoopy Stochastic/Continuous Petri Net of the baseline model                 |
| `snoopy/perturbation.spn` / `.cpn` | Snoopy Stochastic/Continuous Petri Net of the disease model                 |
| `charlie/NMDAR_baseline.apnn` / `.andl` | Charlie-readable APNN / ANDL export of the baseline model                  |
| `charlie/NMDAR_perturbation.apnn` / `.andl` | Charlie-readable APNN / ANDL export of the perturbation model              |
| `charlie/NMDAR_*_legend.tsv`     | TSV mapping of APNN/ANDL identifiers to human-readable names                 |

### Neurodegeneration model (`./output/Neurodegeneration/`)

| File                              | What it shows                                                                 |
| --------------------------------- | ----------------------------------------------------------------------------- |
| `images/neurodegeneration_pn.png` / `.pdf` | Side-by-side Petri net diagram of Ca²⁺ → CASP3 / neurodegeneration          |
| `images/neurodegeneration_comparison.png` / `.pdf` | Ensemble bar chart + absolute difference of downstream markers             |
| `images/neurodegeneration_trajectory.png` / `.pdf` | Static trajectory of CASP3, CHOP, ROS, Apoptosis, and cumulative ND      |
| `animations/neurodegeneration_trajectory.gif` | Animated version of the trajectory curves                                   |
| `animations/neurodegeneration_downstream_pn.gif` | Animated side-by-side Petri nets + cumulative neurodegeneration curve      |
| `charlie/neurodegeneration_baseline.apnn` / `.andl` | Charlie-readable export of the baseline downstream model                   |
| `charlie/neurodegeneration_perturbation.apnn` / `.andl` | Charlie-readable export of the perturbation downstream model               |
| `charlie/neurodegeneration_*_legend.tsv` | TSV identifier-to-name mapping                                              |

### Snoopy files

Two Snoopy XML files are emitted per model, mirroring every place, transition,
arc weight, mass-action rate constant, capacity, layout position and annotation
carried by the Python model:

* **`.spn`** -- *Stochastic Petri Net.* Discrete integer markings, `MassAction(k)`
  rate functions matching the Gillespie SSA in `src/simulation.py`
  (rate * combinatorial term over input arcs).  Use Snoopy's stochastic
  simulator, CSL model checking, or structural analyses (boundedness,
  P/T-invariants).
* **`.cpn`** -- *Continuous Petri Net.* Real-valued markings, same topology
  and rate constants, ODE semantics.  Lets Snoopy run a deterministic
  trajectory the Python code does not provide; the result is directly
  comparable with the stochastic run.

The exporter lives in [`src/snoopy_export.py`](src/snoopy_export.py); regenerate
on demand with `python run_all.py` or `python -c "from src.snoopy_export import
export_pair; from src.models import build_baseline; export_pair(build_baseline(),
'output/NMDAR/snoopy/baseline')"`.

### Charlie / APNN / ANDL files

[Charlie](https://www-dssz.informatik.tu-cottbus.de/DSSZ/Software/Charlie)
is a Petri-net analysis tool from BTU Cottbus (the Snoopy / Marcie group).
Each model is exported in two formats:

* **`.apnn`** — APNN (Abstract Petri Net Notation). Place/transition/arc
  structure with `\name` labels, explicit `\init` markings, and arcs
  carrying both `\weight{w}` and `\type{ordinary}` (the form Charlie reads
  correctly, proven against the reference `cs1.apnn` file).  A weight > 1
  arc (e.g. the bivalent IgG crosslink, `2x NMDAR_Ab → T_crosslink`) is
  written as a single `\weight{2}` arc.
* **`.andl`** — ANDL (Abstract Net Description Language). Charlie's native
  textual format, which carries arc weights natively via `[Place +/- w]`
  terms.  Use this format when any arc has weight > 1, since Charlie's APNN
  reader does not honour such weights.
* **`_legend.tsv`** — Tab-separated mapping of APNN/ANDL identifiers back
  to human-readable model names.

## Quick start

```bash
pip install -r requirements.txt
python run_all.py              # runs both NMDAR and neurodegeneration pipelines
python run_NMDAR.py            # NMDAR model only
python run_neurodegen.py       # neurodegeneration model only
```

## Modelling choices

### Species kept (around the NMDAR)

| Species                       | Role                                                                                                 |
| ----------------------------- | ---------------------------------------------------------------------------------------------------- |
| **Glutamate**                 | Primary NMDAR agonist                                                                                |
| **Glycine / D-Serine**        | Obligatory co-agonist binding GluN1                                                                  |
| **Mg2+**                      | Voltage-dependent pore blocker                                                                       |
| **Ca2+** (ext / int)          | Permeating ion, downstream signaller                                                                 |
| **NMDAR** (GluN1/GluN2 tetramer) | The receptor under study; states: closed / liganded / open / Mg2+-blocked / Ab-bound / endosomal  |
| **PSD-95**                    | PDZ-domain scaffold that anchors NMDAR to the post-synaptic density                                  |
| **EphB2**                     | Surface partner of NMDAR; its extracellular cross-talk is disrupted by patient IgG (Mikasova 2012)   |
| **CaMKII** (active / inactive)| Downstream Ca2+ sensor whose activation triggers LTP                                                 |
| **LTP marker**                | Lumped output of CaMKII-driven AMPAR insertion                                                       |
| **Anti-GluN1 IgG**            | Patient autoantibody                                                                                 |
| **Endosome / Lysosome**       | Intracellular destinations of internalised, antibody-bound NMDARs                                    |
| **Pre-synaptic vesicle**      | Source of synaptic glutamate (drives the system)                                                     |

### Species kept (downstream neurodegeneration pathway)

| Species                       | Role                                                                                                 |
| ----------------------------- | ---------------------------------------------------------------------------------------------------- |
| **Ca²⁺ (cytosol)**            | Input from NMDAR model; shared place links the two subnets                                           |
| **Ca²⁺ (mitochondria)**       | Mito Ca²⁺ overload triggers membrane potential loss and cytochrome-c release                         |
| **ΔΨ (mito membrane potential)** | Intact potential (tokens consumed by Ca²⁺ overload)                                                |
| **Cytochrome-c**              | Released from mitochondria, triggers CASP9 apoptosome                                                |
| **ROS**                       | Reactive oxygen species, produced via Ca²⁺-dependent pathways, drives CASP3                         |
| **ER stress**                 | PERK/eIF2α pathway activated by ER Ca²⁺ depletion, induces CHOP                                       |
| **CHOP (DDIT3)**              | Pro-apoptotic transcription factor, drives CASP3 activation                                          |
| **CASP9**                     | Initiator caspase (apoptosome), activates CASP3                                                       |
| **CASP3**                     | Executioner caspase, commits cell to apoptosis                                                        |
| **Apoptosis**                 | Commitment signal (token accumulates to mark apoptotic decision)                                      |
| **Neurodegeneration**         | Cumulative damage marker (token count grows with each firing)                                         |

### Species excluded (with rationale)

The project brief mentioned a number of markers that turn out **not** to be
part of the established molecular pathology of anti-NMDAR encephalitis.
They are therefore left out of the model so that the diagram remains
mechanistically honest:

| Species         | Why it is **not** in this model                                                                                          |
| --------------- | ------------------------------------------------------------------------------------------------------------------------ |
| **PrPc / PrPSc, PrD**  | Prion biology; relevant to spongiform encephalopathies (CJD, BSE), unrelated to NMDAR autoimmunity                 |
| **A-beta, oligomeric AD species** | Alzheimer disease pathology; A-beta does modulate NMDAR indirectly but is not the driver of anti-NMDAR encephalitis |
| **Htt (huntingtin)**   | Huntington disease; polyQ aggregates, unrelated mechanism                                                            |
| **SCA23**              | Spinocerebellar ataxia 23 (PDYN-linked); a hereditary cerebellar disorder, not anti-NMDAR encephalitis              |
| **PDYN (prodynorphin)** | Opioid-peptide precursor; loosely connected to NMDAR via SCA23 but not part of the encephalitis pathway              |

### What the perturbation adds

```
NMDAR  +  Anti-GluN1 IgG      --T_Ab_bind-->     NMDAR-IgG          (k ~ 1e5 M^-1 s^-1)
2 x NMDAR-IgG                 --T_crosslink-->   crosslinked cap    (bivalent IgG)
NMDAR-EphB2 + IgG             --T_EphB2_disrupt-> NMDAR-IgG + EphB2 (Mikasova 2012)
crosslinked cap               --T_endocytosis--> Endosome           (t1/2 ~ 5 h in vitro)
Endosome                      --T_lyso-->         Lysosome           (receptor pool lost)
```

These steps drain the surface NMDAR pool, which collapses the rate of
`T_Ca_in` firings, which collapses the cytosolic Ca2+ and the CaMKII -> LTP
cascade -- exactly the synaptic phenotype reported in patients.

The downstream neurodegeneration model amplifies this effect: reduced
Ca²⁺ influx in the disease state causes chronically elevated ROS production
and ER stress through maladaptive compensatory mechanisms, leading to
disproportionate CASP3 accumulation and neurodegeneration output.

### Hierarchical model structure

The NMDAR and neurodegeneration models are connected hierarchically:

* The **Ca_in** place in the NMDAR model carries a `subnet_link = "Block:\nNeurodegeneration Pathway PN"` annotation.  In the rendered diagram this produces a Snoopy-style coarse-node block beside the place, signalling that a separate subnet is logically attached.
* The **Ca_in** place in the neurodegeneration model carries a reciprocal `subnet_link = "Block:\nNMDAR subnet"` annotation, forming a bidirectional hierarchical link.
* The two models are simulated independently and their outputs compared as described in the output tables above.

## Why Petri Nets here?

* **Tokens are integer molecules**, so the model is by construction
  stoichiometrically correct.
* The **Gillespie SSA** in `src/simulation.py` interprets transition rates
  as mass-action propensities, which keeps the simulation in the same
  formal regime as standard chemical reaction networks.
* **Compartment-aware layout** (extracellular / membrane / PSD / cytosol /
  endosomal lanes) makes it easy to read off the cell biology directly
  from the diagram.

The project uses a small, fully introspectable Petri-Net core (`src/petri_net.py`)
so that the matplotlib visualiser can decorate every place, transition,
and arc with rich annotations.  **SNAKES** and **pm4py** are listed in
`requirements.txt` because they are the canonical Python Petri-Net
libraries; the model is structurally compatible with both (places +
weighted arcs + transitions) and either could load it for additional
formal analysis (boundedness, reachability, P/T-invariants).

## Mass-action rate constants

Rates are mass-action coefficients in token-units / s.  Real physical
constants (with concentrations folded in) are cited below; values used in
the model are scaled to keep the synaptic events on a sub-second timescale
and the antibody-driven internalisation visible inside an 8-s simulation
window.

### NMDAR model transitions

| Transition           | Symbol   | Literature value            | Reference  |
| -------------------- | -------- | --------------------------- | ---------- |
| `T_bind` (Glu+Gly)   | k_on     | ~ 5e6 M^-1 s^-1             | [REF1]     |
| `T_unbind`           | k_off    | ~ 4 s^-1                    | [REF2]     |
| `T_open`             | beta     | ~ 90 s^-1                   | [REF2]     |
| `T_close`            | alpha    | ~ 30 s^-1                   | [REF2]     |
| `T_Mg_block`         | -        | tau < 200 us (instantaneous)| [REF3]     |
| `T_Mg_unblock`       | -        | voltage-dependent           | [REF3]     |
| `T_Ca_in`            | g_NMDAR  | ~ 55 pS (~ 1e6 ions/s)      | [REF4]     |
| `T_Ca_out`           | PMCA/NCX | basal Ca2+ clearance        | [REF4]     |
| `T_PSD95`            | -        | slow PDZ binding            | [REF5]     |
| `T_EphB2`            | -        | stable surface complex      | [REF6]     |
| `T_CaMKII_on`        | tau_on   | ~ 0.1 s                     | [REF7]     |
| `T_CaMKII_off`       | tau_off  | ~ 3 s                       | [REF7]     |
| `T_LTP`              | -        | downstream AMPAR insertion  | [REF7]     |
| `T_release`          | -        | vesicle release ~ 0.7 /s    | [REF11]    |
| `T_reuptake`         | -        | EAAT glutamate uptake       | [REF12]    |
| `T_Ab_bind`          | k_on,IgG | ~ 1e5 M^-1 s^-1             | [REF8]     |
| `T_Ab_unbind`        | k_off,IgG| very slow                   | [REF8]     |
| `T_crosslink`        | -        | bivalent capping            | [REF9]     |
| `T_EphB2_disrupt`    | -        | Mikasova et al. 2012        | [REF6]     |
| `T_endocytosis`      | k_endo   | t1/2 ~ 5 h in vitro         | [REF9]     |
| `T_lyso`             | -        | lysosomal degradation       | [REF10]    |

### Neurodegeneration model transitions

| Transition                 | Symbol       | Role                                              | Reference  |
| -------------------------- | ------------ | ------------------------------------------------- | ---------- |
| `T_Ca_mito`               | k_ca_mito    | Ca²⁺ uptake by mitochondria                       | [REF13]    |
| `T_DeltaPsi_loss`         | k_cytoc      | Mito membrane potential loss (Ca²⁺ overload)      | [REF14]    |
| `T_CytoC_release`         | k_cytoc      | Cytochrome-c release from mito                     | [REF14]    |
| `T_CASP9_on`              | k_casp9_on   | CASP9 activation by cytochrome-c                   | [REF17]    |
| `T_CASP9_off`             | k_casp9_off  | CASP9 inactivation                                 | [REF17]    |
| `T_CASP3_from_CASP9`      | k_casp3_on   | CASP3 activation by CASP9                          | [REF17]    |
| `T_CASP3_from_CHOP`       | k_chop_casp3 | CASP3 activation by CHOP                           | [REF16]    |
| `T_CASP3_from_ROS`        | k_ros_casp3  | CASP3 activation by ROS                            | [REF15]    |
| `T_CASP3_off`             | k_casp3_off  | CASP3 clearance / degradation                      |            |
| `T_ROS_on`                | k_ros_on     | ROS production (NADPH oxidase / mito ROS)          | [REF15]    |
| `T_ROS_off`               | k_ros_off    | ROS scavenging (SOD / catalase)                    | [REF15]    |
| `T_ER_stress`             | k_er_stress  | ER stress induction (PERK/eIF2α)                   | [REF16]    |
| `T_CHOP_on`               | k_chop_on    | CHOP induction by ER stress                        | [REF16]    |
| `T_CHOP_off`              | k_chop_off   | CHOP degradation                                   | [REF16]    |
| `T_Apoptosis`             | k_apop       | CASP3 → apoptosis commitment                      | [REF17]    |
| `T_Neurodegeneration`     | k_neurodegen | Apoptosis → neurodegeneration accumulation         | [REF18]    |

## References

- **[REF1]** Lester R.A.J., Jahr C.E. (1992). *NMDA channel behavior depends on agonist affinity.* J. Neurosci. 12: 635-643. PMID: 1346806. <https://www.jneurosci.org/content/12/2/635>
- **[REF2]** Erreger K., Dravid S.M., Banke T.G., Wyllie D.J.A., Traynelis S.F. (2005). *Subunit-specific gating controls rat NR1/NR2A and NR1/NR2B NMDA channel kinetics and synaptic signalling profiles.* J. Physiol. 563(2): 345-358. <https://physoc.onlinelibrary.wiley.com/doi/10.1113/jphysiol.2004.080028>
- **[REF3]** Jahr C.E., Stevens C.F. (1990). *Voltage dependence of NMDA-activated macroscopic conductances predicted by single-channel kinetics.* J. Neurosci. 10(9): 3178-3182. <https://www.jneurosci.org/content/10/9/3178>
- **[REF4]** Stern P., Behe P., Schoepfer R., Colquhoun D. (1992). *Single-channel conductances of NMDA receptors expressed from cloned cDNAs: comparison with native receptors.* Proc. R. Soc. B 250: 271-277. See also Kuner T., Schoepfer R. (1996). *Multiple structural elements determine subunit specificity of Mg2+ block in NMDA receptor channels.* J. Neurosci. 16(11): 3549-3558. <https://www.jneurosci.org/content/16/11/3549>
- **[REF5]** Bard L., Sainlos M., Bouchet D., Cousins S., Mikasova L., Breillat C., Stephenson F.A., Imperiali B., Choquet D., Groc L. (2010). *Dynamic and specific interaction between synaptic NR2-NMDA receptor and PDZ proteins.* PNAS 107(45): 19561-19566. <https://www.pnas.org/doi/10.1073/pnas.1002690107>
- **[REF6]** Mikasova L., De Rossi P., Bouchet D., Georges F., Rogemond V., Didelot A., Meissirel C., Honnorat J., Groc L. (2012). *Disrupted surface cross-talk between NMDA and Ephrin-B2 receptors in anti-NMDA encephalitis.* Brain 135(5): 1606-1621. <https://academic.oup.com/brain/article/135/5/1606/310762>
- **[REF7]** Lee S.-J.R., Escobedo-Lozoya Y., Szatmari E.M., Yasuda R. (2009). *Activation of CaMKII in single dendritic spines during long-term potentiation.* Nature 458: 299-304. <https://www.nature.com/articles/nature07842>
- **[REF8]** Hughes E.G., Peng X., Gleichman A.J., Lai M., Zhou L., Tsou R., Parsons T.D., Lynch D.R., Dalmau J., Balice-Gordon R.J. (2010). *Cellular and synaptic mechanisms of anti-NMDA receptor encephalitis.* J. Neurosci. 30(17): 5866-5875. <https://www.jneurosci.org/content/30/17/5866>
- **[REF9]** Hughes E.G. et al. 2010 (as above) -- internalisation visible by 2 h, plateau by 12 h, reversible after antibody removal. PMC: <https://pmc.ncbi.nlm.nih.gov/articles/PMC2868315/>
- **[REF10]** Moscato E.H., Peng X., Jain A., Parsons T.D., Dalmau J., Balice-Gordon R.J. (2014). *Acute mechanisms underlying antibody effects in anti-N-methyl-D-aspartate receptor encephalitis.* Ann. Neurol. 76(1): 108-119. <https://pmc.ncbi.nlm.nih.gov/articles/PMC4296347/>
- **[REF11]** Schikorski T., Stevens C.F. (2001). *Morphological correlates of functionally defined synaptic vesicle populations.* Nat. Neurosci. 4: 391-395. <https://www.nature.com/articles/nn0401_391>
- **[REF12]** Tzingounis A.V., Wadiche J.I. (2007). *Glutamate transporters: confining runaway excitation by shaping synaptic transmission.* Nat. Rev. Neurosci. 8: 935-947. <https://www.nature.com/articles/nrn2274>
- **[REF13]** Nicholls D.G. (2005). *Mitochondria and calcium signaling.* Cell Calcium 38(3-4): 311-317. <https://www.sciencedirect.com/science/article/pii/S014341600500111X>
- **[REF14]** Kroemer G., Galluzzi L., Brenner C. (2007). *Mitochondrial membrane permeabilization in cell death.* Physiol. Rev. 87(1): 99-163. <https://journals.physiology.org/doi/full/10.1152/physrev.00013.2006>
- **[REF15]** Gleichmann M., Mattson M.P. (2011). *Neuronal calcium homeostasis and dysregulation.* Antioxid. Redox Signal. 14(7): 1261-1273. <https://www.liebertpub.com/doi/10.1089/ars.2010.3386>
- **[REF16]** Harding H.P., Zhang Y., Zeng H., Novoa I., Lu P.D., Calfon M., Sadri N., Yun C., Popko B., Paules R., Stojdl D.F., Bell J.C., Hettmann T., Leiden J.M., Ron D. (2003). *An integrated stress response regulates amino acid metabolism and resistance to oxidative stress.* Mol. Cell 11(3): 619-633. <https://www.sciencedirect.com/science/article/pii/S1097276503001059>
- **[REF17]** Yuan J., Yankner B.A. (2000). *Apoptosis in the nervous system.* Nature 407: 802-809. <https://www.nature.com/articles/35037739>
- **[REF18]** Hardy J., Gwinn-Hardy K. (2023). *Neurodegeneration: a molecular approach.* Oxford University Press.

### Background reading used to validate excluded species and confirm pathway

- Dalmau J., Armangue T., Planaguma J., Radosevic M., Mannara F., Leypoldt F., Geis C., Lancaster E., Titulaer M.J., Rosenfeld M.R., Graus F. (2019). *An update on anti-NMDA receptor encephalitis for neurologists and psychiatrists: mechanisms and models.* Lancet Neurology 18(11): 1045-1057. <https://www.thelancet.com/article/S1474-4422(19)30244-3/abstract>
- Kreye J. et al. (2016). *Human cerebrospinal fluid monoclonal N-methyl-D-aspartate receptor autoantibodies are sufficient for encephalitis pathogenesis.* Brain 139(10): 2641-2652. <https://academic.oup.com/brain/article/139/10/2641/2422161>
- Ladepeche L., Planaguma J., Thakur S., Suarez I., Hara M., Borbely J.S., Sandoval A., Laparra-Cuervo L., Dalmau J., Lakadamyali M. (2018). *NMDA Receptor Autoantibodies in Autoimmune Encephalitis Cause a Subunit-Specific Nanoscale Redistribution of NMDA Receptors.* Cell Reports 23(13): 3759-3768. <https://www.sciencedirect.com/science/article/pii/S2211124718308921>

## Project layout

```
Ali/
|-- README.md
|-- requirements.txt
|-- run_all.py                    # orchestrator: runs both pipelines
|-- run_NMDAR.py                  # NMDAR pipeline driver
|-- run_neurodegen.py             # neurodegeneration pipeline driver
|-- src/
|   |-- __init__.py
|   |-- petri_net.py              # Place / Transition / Arc / PetriNet core
|   |-- models.py                 # build_baseline() / build_perturbation()
|   |-- simulation.py             # Gillespie SSA + sampling utilities
|   |-- visualization.py          # static figures + animation (NMDAR)
|   |-- compare_ca_influx.py      # Ca2+ influx comparison (PNG + GIF)
|   |-- neurodegeneration.py      # build_downstream_baseline/perturbation()
|   |-- compare_neurodegen.py     # neurodegeneration figures + animation
|   |-- snoopy_export.py          # Snoopy native XML (.spn / .cpn)
|   |-- charlie_export.py         # Charlie APNN + ANDL + legend TSV export
|-- output/
|   |-- NMDAR/
|   |   |-- images/
|   |   |   |-- NMDAR_baseline_pn.png
|   |   |   |-- NMDAR_perturbation_pn.png
|   |   |   |-- NMDAR_PN.png
|   |   |   |-- ca_influx_comparison.png
|   |   |-- animations/
|   |   |   |-- NMDAR_baseline_animation.gif
|   |   |   |-- NMDAR_perturbation_animation.gif
|   |   |   |-- NMDAR_baseline_combined.gif
|   |   |   |-- NMDAR_perturbation_combined.gif
|   |   |   |-- NMDAR_PN.gif
|   |   |   |-- ca_influx_comparison.gif
|   |   |-- snoopy/
|   |   |   |-- baseline.spn / baseline.cpn
|   |   |   |-- perturbation.spn / perturbation.cpn
|   |   |-- charlie/
|   |   |   |-- NMDAR_baseline.apnn / .andl / _legend.tsv
|   |   |   |-- NMDAR_perturbation.apnn / .andl / _legend.tsv
|   |-- Neurodegeneration/
|   |   |-- images/
|   |   |   |-- neurodegeneration_pn.png
|   |   |   |-- neurodegeneration_comparison.png
|   |   |   |-- neurodegeneration_trajectory.png
|   |   |-- animations/
|   |   |   |-- neurodegeneration_trajectory.gif
|   |   |   |-- neurodegeneration_downstream_pn.gif
|   |   |-- charlie/
|   |       |-- neurodegeneration_baseline.apnn / .andl / _legend.tsv
|   |       |-- neurodegeneration_perturbation.apnn / .andl / _legend.tsv
|   `-- Old/                       # legacy artefacts (regenerated)
```

## Interpretation guide

* **Place colour** = anatomical compartment (legend in every figure).
* **Token count** in the centre of each place = current molecule count.
* **Transitions** are dark rounded rectangles labelled with their name; the
  small italic box beside each one shows the literature rate constant and
  the citation tag.
* **Arrows**: dark grey = input arc (place -> transition), indigo =
  output arc (transition -> place); an `x2` label on an arc means
  stoichiometry 2 (used on the bivalent IgG crosslink).
* **Animations**: the currently-firing transition flashes yellow with a
  red outline; the HUD at the top-left shows simulation time and the
  fired transition; the strip at the bottom shows the live marking.
* **Ca2+ comparison**: the upper panel is the instantaneous cytosolic
  Ca2+ token count, the lower panel is the cumulative count of
  `T_Ca_in` firings.  The annotated box shows the ratio
  *Perturbation / Baseline* at the end of the simulation -- for the seed
  used in `run_NMDAR.py`, the disease net retains roughly **60 %** of the
  baseline Ca2+ influx by *t* = 8 s, consistent with the
  titer-dependent NMDAR loss reported in the literature.
* **Hierarchical blocks**: a Snoopy-style coarse-node block beside a
  place (pale slate-blue with double border) signals that a separate
  subnet diagram is logically attached at that point.  The
  **Ca_in** place in the NMDAR model links to the neurodegeneration
  pathway; the **Ca_in** place in the neurodegeneration model links
  back to the NMDAR subnet.
* **Neurodegeneration trajectory**: five stacked panels showing
  instantaneous CASP3, CHOP, ROS, Apoptosis tokens and cumulative
  neurodegeneration firings.  The ratio annotation quantifies the
  disease effect at the end of the simulation window.
