---
title: "Simulation and theory: suggested additions for the OB–PA manuscript"
---

*Draft additions covering the computational simulations and the charge-regulated
adsorption model. Part A slots into the main text (Methods 2.2.10 and Results
3.3.4, with a pointer for 3.3.3); Part B is the matching Supplementary
Information. Symbols are defined at first use; numbers come from the analysis in
`zeta/charge_regulation.py`, `mc/`, and `nmr/`.*

---

# Part A — Suggested main text

## 2.2.10. Computational simulations

We modelled phytate binding to an N-terminal oleosin brush with coarse-grained
Monte Carlo simulations. The N-terminal arm of *Cannabis sativa* oleosin (UniProt
A0A803P977) was modelled one bead per residue with the CALVADOS 3 force field
(Tesei et al., 2021; Cao et al., 2024) and grafted as a fifty-chain brush to one
wall of a slit, representing the cationic protein layer at the oil-body surface.
A single phytate molecule, represented as a rigid seven-bead unit (one inositol
centre, six phosphates), titrated all twelve dissociable protons together with the
peptide in the semi-grand-canonical ensemble, so its net charge $z$ tracks pH. The
Faunus engine (Stenqvist et al., 2013) sampled pH 4.0, 6.0 and 10.0 at 5–150 mM NaCl
and returned, against the distance $r$ from the grafted wall: the potential of mean
force $w(r)$ (bulk-referenced), the electric potential $\varphi(r)$, the phytate
valence $z(r)$, the wall alignment $S_\text{wall}$, and the relative rotational
mobility $D(r)/D_0$ of an adsorbed phytate.

A charge-regulated Poisson–Boltzmann (CR-PB) model then extended this molecular
picture to finite phytate concentration and to the measured $\zeta$-potential. The
simulation supplied the local binding-plane potential and phytate valence; the
measured bare-surface $\zeta$ anchored the diffuse-layer surface charge; and one
intrinsic affinity $K$ per pH captured the non-electrostatic contact term. The full
model, parameters and validity window are in Supplementary Section S-2.

## 3.3.4. Molecular characterization of oleosin-PA binding

The simulations identify a local binding-plane potential that reconciles two
observations: phytate binds near pH 6 even when the net oil-body charge is small, and
phytate can overcharge the surface at pH 4. Figure 10 collects the low-salt
simulation output — the potential of mean force, the phytate charge, the electric
potential, and the rotational mobility, by pH.

![Figure 10. Low-salt (5 mM) simulation composite by pH: (A) potential of mean force $w(r)$, (B) relative rotational mobility $D(r)/D_0$, (C) electric potential with simulated (open) and measured (filled) $\zeta$, (D) phytate charge, versus distance $r$ from the grafted wall.](../mc/composite_I0.005.png)

The ionizable oleosin termini raise the local electric potential at the brush
periphery above the surface average: the potential reaches $+52$ mV at the pH-4
binding minimum and remains positive ($\approx +5$ mV) at pH 6, where the net oil-body
charge is near zero. Phytate can therefore bind near the isoelectric point because the
local cationic field attracts it even when the net charge does not. In the CR-PB model,
the measured bare-surface $\zeta$ sets the net charge that governs the diffuse layer,
whereas the simulated binding-plane potential sets the local phytate enrichment.

Surface coverage reflects a balance between attraction to the binding-plane potential
and repulsion between adsorbed phytates. This lateral repulsion scales as $z^2$ and is
screened by salt, so pH and salt set both the attraction to the brush and the maximum
adsorbed charge. At pH 6, the higher phytate valence limits coverage relative to the
bare attraction and therefore sets how far the surface can be overcharged.

At pH 4, this balance tips far enough that phytate adsorbs past neutralization:
each molecule carries $z \approx -7$ at the surface (Figure 10D; the bulk value
is $\approx -6$), so once the adsorbed charge exceeds the bare surface charge the sign
of the surface charge inverts. This charge inversion
(overcharging) is the mechanism behind the measured $\zeta$-potential crossing from
$+38$ mV to $-4$ mV as phytate is added (Figure 6A). It is self-limiting: as the
surface charges negative, the local potential opposes further anion binding, so the
inversion saturates. The CR-PB model reproduces the crossing (Supplementary Section
S-3); a non-regulating mean-field surface cannot, because it stops binding at
neutralization.
This charge-regulated overcharging explains the high-PA $\zeta$ reversal and is
consistent with the pH-4 destabilization, which begins earlier through charge
reduction and PA-mediated inter-droplet bridging.

The mechanism makes three experimental predictions: phytate should shift
$\zeta$-potential in a salt-dependent manner, slow phytate reorientation without
immobilizing it, and produce more heat than a self-regulating monolayer can explain.

The electrophoretic data test the first prediction. The model reproduces the pH-4
reversal and the low-salt pH-6 phytate series, where added phytate makes $\zeta$ more
negative. The phytate-induced shift collapses as salt is added in both the experiment
and the simulation: the pH-6 PMF well shrinks from $-15$ to $-2.6\ k_BT$ between 5 and
150 mM salt (Supplementary Section S-4), so screening weakens adsorption. At pH 10 the
simulated interaction is repulsive,
matching the negligible direct binding; the high-phytate destabilization there is
double-layer compression, not adsorption.

The NMR data provide an independent mobility and bound-fraction check. The PMF-derived
bound fraction matches the fraction required by NMR without
fitting an NMR-specific parameter. At 0.15 M, the PMF gives a surface partition length
$K_H$; the oil-body specific surface area and adsorption isotherm then give the bound
fraction $f_\text{bound}$ (Supplementary Section S-5). The measured slow-motion contrast
$\Delta R_{1\rho}=R_{1\rho}-R_1$ (Figure 9G) gives the fraction required by a two-state
fast-exchange model. The simulated and required fractions agree at every pH (Table 1),
so the free-energy profile accounts quantitatively for the NMR magnitude at the
measured salt.

Table 1. Bound fraction from the MC PMF at 0.15 M, 7.6 mM PA ($K_H$ the surface partition
length; $f_\text{req}$ the fraction the two-site model requires, using the fully-bound
rate $\Delta R_\text{bound}=3\times10^{4}$ s⁻¹).

| pH | $K_H$ / Å | $f_\text{bound}$ (sim) / % | $f_\text{req}$ / % |
|----|-----------|----------------------------|--------------------|
| 4  | 925  | 1.30 | 1.33 |
| 6  | 1209 | 1.05 | 1.01 |
| 10 | 17   | 0.11 | 0.16 |

The simulation also explains why bound phytate still relaxes like a mobile species.
Its rotational mobility relative to the free molecule, $D(r)/D_0$ (Figure 10B),
is $0.5$–$0.7$ near the surface, so it reorients at roughly half its bulk rate but is
not frozen. The raised $\Delta R_{1\rho}$ at pH 4 and 6 then reflects slow reorientation as
phytate co-tumbles with the micron-scale droplet, not rigid immobilization —
consistent with the absent cross-polarization signal. A bare-hydrodynamic estimate
supports a single free-tumbling reference $\Delta R_\text{free}$: free phytate tumbles
near-isotropically, with rotational correlation time $\tau_r \approx 134$–$147$ ps and
little variation across
protonation states (Supplementary Section S-5).

The calorimetry tests whether the measured heat can be assigned to monolayer
adsorption. The one-site model and the charge-regulated model both fit the pH-4
thermogram (Supplementary Section S-6), but both require $\approx 0.38$ phytate per
oleosin to account for the heat. The surface charge implied by $\zeta$ and the
NMR-scale bound fraction caps the surface at $\approx 0.03$ phytate per oleosin,
leaving a ~13-fold gap. Most of the heat must therefore arise off the surface,
plausibly from the inter-droplet
cross-linking already suggested by the calorimetric shoulder at $\sim 0.3$ PA:oleosin
(Section 3.3.2). The one-site $K_a$, $n$ and $\Delta H$ should therefore be presented as
apparent thermodynamic descriptors of a multiphasic pH-4 process, not as direct
surface-binding constants. A non-aggregating ITC experiment could use the same isotherm
framework to quantify surface binding (Supplementary Section S-6).

The mean-field treatment has one important limitation. It captures the regulation
feedback — binding weakening as the surface charges up — but smears the adsorbed
layer, so it omits the discrete phytate–phytate correlations that matter at low salt,
where the double layer is thick relative to the inter-phytate spacing. There the
mean-field description is least reliable, and the simulation provides the cross-check
(Supplementary Section S-2).

## CONCLUSIONS

The simulations and the charge-regulated model give the electrostatic mechanism behind
the observed behavior. Phytate binding reflects a pH- and salt-dependent balance
between attraction to the binding-plane potential and repulsion within the adsorbed
phytate layer. At pH 4 this balance can push adsorption past neutralization, producing
overcharging and the measured $\zeta$ reversal at high PA; the earlier
destabilization is consistent with charge reduction and PA-mediated inter-droplet
bridging. At pH 6 the weaker, salt-screened binding adds negative charge that
stabilizes the emulsion. At pH 10 direct binding is negligible, and the double layer
governs behavior.

The three molecular probes are consistent with this mechanism. The CR-PB model
reproduces the $\zeta$ reversal, and the NMR bound fraction follows quantitatively from
the simulated free energy with no parameter fitted to the NMR data. The same framework
shows why calorimetry resists a physical surface-binding fit: the surface coverage is
far below what the measured heat requires, so most of the heat must arise off the
surface, plausibly through inter-droplet cross-linking. At low salt, where mean field is
least reliable, the simulation provides the reference.

---

# Part B — Suggested Supplementary Information

## S-1. Monte Carlo simulations

### S-1.1. System and grafting

The modelled fragment is the 14-residue N-terminal arm of
*Cannabis sativa* oleosin (UniProt A0A803P977, sequence MKSFLPESGPSAQQ),
coarse-grained one bead per residue with CALVADOS 3 (Tesei et al., 2021;
Cao et al., 2024). Consecutive beads are joined by harmonic bonds (force constant
80.33 kJ mol⁻¹ Å⁻², equilibrium length 3.8 Å). Fifty chains are grafted at
$3\times10^{-3}$ Å⁻² onto the bottom wall of a $130\times130\times200$ Å³ slit
(periodic in $x$ and $y$; hard walls at $z=\pm100$ Å) through a harmonic tether on the
terminal anchor bead (0.5 kJ mol⁻¹ Å⁻²). Water is implicit, entering through the
temperature-dependent dielectric permittivity at 298.15 K.

### S-1.2. Peptide titration

Four groups protonate through atom-type swap moves against an
implicit proton bath ($a_{\text{H}^+}=10^{-\text{pH}}$): Glu (p$K_a$ 4.24), Lys
(10.4), the N-terminus (8.0) and the C-terminus (3.67) (Grimsley et al., 2009). The
mean chain charge runs from $\approx +1e$ at pH 4 to $\approx -1.3e$ at pH 10, i.e.
the brush passes through its isoelectric region near
pH 6, where the local cationic patch still supports binding even though the net-charge
picture does not (Section 3.3.4).

### S-1.3. Phytate titration

One phytate
(*myo*-inositol hexakisphosphate) is a rigid seven-bead unit: one inositol centre
(INO, diameter 6.2 Å) and six phosphates (PH, 5.8 Å) placed at the crystallographic
phosphorus positions of the chair conformation ($\approx 4.5$ Å from the ring centre,
one axial and five equatorial); bead diameters are scaled from CALVADOS 3 as
$M^{1/3}$. Phytate carries twelve dissociable protons — two per phosphate — giving
thirteen charge states from neutral H₁₂PA to PA¹²⁻. Rather than fixing the charge,
all twelve dissociations are sampled as molecular swap moves in the
semi-grand-canonical ensemble, with macroscopic p$K_a$ values from potentiometry in
NaCl (De Stefano et al., 2003): a low-p$K_a$ sextet 1.1, 1.5, 1.7, 2.1, 2.5, 2.9 (the
first proton on each phosphate) and a higher, spread-out sextet 5.72, 6.28, 6.81,
7.60, 9.94, 11.84 (the second). This bimodal ladder makes the phytate valence climb
steeply with pH ($z\approx-6$ at pH 4 to $\approx-10.6$ at pH 10). It also straddles
the experimental window, so phytate is itself charge-regulated as it approaches the
brush. The net charge is spread uniformly over the six phosphate beads. All thirteen
states share the same seven-bead geometry, so a swap move changes only the charge and
preserves the molecule's position and orientation through gyration-tensor alignment. A
single phytate can therefore explore its full protonation manifold without the
geometric bias introduced by separate fixed-charge models.

### S-1.4. Interactions

Nonbonded pairs combine the Ashbaugh–Hatch potential (short-range,
hydrophobicity-scaled; Ashbaugh & Hatch, 2008) with a screened-Coulomb (Debye–Hückel)
term set by the salt. Phytate–peptide cross-pairs use $\lambda=0$ (purely repulsive
WCA) to suppress a spurious hydrophobic attraction that arithmetic $\lambda$-mixing
between hydrophobic amino acids and the non-hydrophobic phytate would otherwise
introduce. Bonded neighbors and all intramolecular phytate pairs are excluded.

### S-1.5. Monte Carlo moves

Per sweep, Faunus (Stenqvist et al., 2013) applies: peptide
translation in $xy$ ($\delta_\text{max}$ 10 Å) and phytate translation along $z$
(5 Å, 20 repeats); peptide pivot rotation (0.5 rad) and phytate rigid rotation
(1.0 rad); and speciation (molecular swap) for the twelve phytate and four peptide
reactions (10 repeats per draw).

### S-1.6. Potential of mean force

The PMF along the inositol $z$-coordinate is obtained by
multi-walker umbrella sampling with 40 hard-wall windows (6 Å wide, $\approx 2.4$ Å
apart, overlapping) spanning the grafted wall to the slit centre. A harmonic
drive (5.0 kJ mol⁻¹ Å⁻²) steers phytate into each window before the bias is removed
for unbiased production; adjacent windows are stitched by the overlap ratio
$\Delta F_{i\to i+1}=-k_BT\ln(f_i/f_{i+1})$, with $f_i$ the fraction of window-$i$
samples falling in the $i/i{+}1$ overlap. PMFs were computed at pH 4, 6 and 10 and
$I=0.005$, 0.03 and 0.15 M.

### S-1.7. Position-resolved rotational mobility

$D(r)/D_0$ is read from the Monte Carlo
rotations themselves: the mean accepted squared angular step per attempted move,
$\langle\Delta\theta^2\rangle(r)$, normalized to its bulk value near the box centre.
Monte Carlo carries no time scale, so only this relative profile is meaningful — it
measures how the brush hinders phytate tumbling on approach. The
disk-shaped phytate also takes only a weak orientational bias at the wall
($S_\text{wall}\approx-0.005$, the second-rank order parameter: $0$ isotropic,
$-\tfrac12$ face-on, $+1$ edge-on).

### S-1.8. Electric potential

The laterally averaged slab charge $\sigma(z)$ gives the
screened-sum potential

$$\varphi(z)=\frac{2\pi l_B}{\kappa}\sum_{z'}\sigma(z')\,e^{-\kappa|z-z'|},$$

with $l_B$ the Bjerrum length and $1/\kappa$ the Debye length, referenced to zero in
bulk. The simulated $\zeta$ is $\varphi$ at the shear plane, taken as the outermost
$r$ at which the peptide charge $|\sigma(z)|$ has not yet fallen below 2% of its peak
(the open markers in the composite panel C).

## S-2. The charge-regulated Poisson–Boltzmann model

The oil-body surface is a planar Gouy–Chapman wall in a 1:1 electrolyte whose Debye
screening $\kappa$ is fixed by the background salt. The model separates the far-field
monopole from the local binding field. The net surface charge density $\sigma_0$ is
obtained from the measured bare $\zeta$ through the Grahame relation, while the local
binding-plane potential $\psi_\beta$ is read from the simulation at the brush periphery
where phytate binds.

The Grahame relation ties the diffuse-layer (surface) potential $\psi_d$ to the
effective surface charge $\sigma_\text{eff}$,

$$\sigma_\text{eff} = \frac{\kappa\,\varepsilon\varepsilon_0 k_B T}{e}\,
\sinh\!\Big(\frac{e\psi_d}{2k_BT}\Big),$$

with $e$ the elementary charge, $k_BT$ the thermal energy, and
$\varepsilon\varepsilon_0$ the permittivity. Phytate of valence $z$ and bulk
concentration $c_L$ occupies surface sites (areal density $\Gamma$, intrinsic
association constant $K$) at the Boltzmann-enhanced contact concentration
$c_L^{s}=c_L\,e^{-z e\psi_\beta/k_BT}$, giving a Langmuir occupancy $\theta$ and a
regulated charge,

$$\theta=\frac{K c_L^{s}}{1+K c_L^{s}},\qquad
\sigma_\text{eff}=\sigma_0 + z e\,\Gamma\,\theta .$$

Because the adsorbed phytate itself shifts $\psi_\beta$, the two expressions for
$\sigma_\text{eff}$ are solved self-consistently (one-dimensional root finding) at
each $c_L$; the $\zeta$-potential is read as $\psi_d$. This is a charge-regulation
boundary condition on Gouy–Chapman (Grahame, 1947; Trefalt et al., 2016). The
binding-plane potential $\psi_\beta$ carries the cationic-patch dipole from the
simulation, so it is positive at pH 4 and 6 even where $\sigma_0$ is near zero; the
model therefore binds phytate near the isoelectric point while $\sigma_0$ still sets
the diffuse layer.

### S-2.1. Overcharging

The surface charge inverts once the adsorbed charge exceeds the bare
charge, $|z| e\,\Gamma\theta > |\sigma_0|$. With $|z|\approx 7$ each phytate flips a
large charge, so inversion is efficient and self-limiting: as $\psi_d$ crosses zero
the Boltzmann factor suppresses further binding (Grosberg et al., 2002).

### S-2.2. Two roles, and the affinity $K$

The model serves two purposes: it computes the
$\zeta$-potential by combining the simulated local binding field with the measured
bare-surface charge anchor, and it provides an isotherm to fit the calorimetry. $K$ is
the intrinsic (non-electrostatic) affinity — the electrostatics are already carried by
$\psi_\beta$ and $\sigma_\text{eff}$, so $K$ is *not* the simulated Henry constant
(using that would double-count the electrostatics). It reflects the pH-set protonation
of the oleosin and phytate groups and is nearly salt-independent. For the
$\zeta$-potential it is fit once per pH to the $\zeta$-versus-[PA] series and held
across salt (Table S-3a). For the calorimetry it is fit freely, together with the
surface area, as an independent test of whether a surface-adsorption isotherm can
describe the heat at all (Section S-6). The NMR does not use this model: its bound
fraction comes from the MC PMF and a two-dimensional adsorption isotherm (Section S-5).
The saturation coverage $\Gamma_\text{max}=1/(11\ \text{Å})^2$ is a geometric
close-packing prior. The
$\zeta$ sample is diluted 1000-fold, so its coverage is the isotherm re-equilibrated at
$[\text{PA}]/1000$.

### S-2.3. Scope and validity

The mean field captures the regulation feedback but smears the
adsorbed layer, omitting discrete phytate–phytate correlations (Grosberg et al.,
2002). Their relative importance is set by the screened
two-dimensional coupling $\Gamma_\text{eff}\approx z^2 l_B\,e^{-\kappa a}/a$, with
$l_B\approx 7$ Å the Bjerrum length and $a\approx\sqrt{|z|e/|\sigma_0|}$ the
inter-phytate spacing set by the surface charge $\sigma_0$: mean field is adequate at
$\sim 100$ mM
($\Gamma_\text{eff}\lesssim 0.2$) and correlations grow at low salt
($\Gamma_\text{eff}\sim 2$ at 10 mM). The pH-6 high-salt residual (below) sits where
the model is expected to hold; the low-salt corner is where the simulation is the
reference.

## S-3. $\zeta$-potential

Figure S (below) compares the model with the measured $\zeta$-potential
versus phytate (A) and versus salt (B). Panel A shows the pH-4 charge reversal
reproduced; the model uses one intrinsic $K$ per pH (Table S-3a), fit to the
phytate series and held across salt.

![Measured $\zeta$-potential (× with error bars) and the charge-regulated model (lines) coloured by pH (blue 4, orange 6, green 10): (A) versus [PA] at low salt; (B) versus [NaCl] at 1.5 mM PA. The model is an absolute prediction (not anchored to the bare salt curve).](../zeta/fig_crpb.png)

Table S-3a. Intrinsic affinities.

| pH | $K$ / M⁻¹ | binding |
|----|----------------|---------|
| 4  | 64             | overcharges (charge inversion) |
| 6  | $2.9\times10^{4}$ | binds, adds negative charge |
| 10 | 0              | repulsive, no direct binding |

The pH-6 salt series is the model's clearest limitation, and the simulation explains
it. The measured phytate-induced shift screens away almost completely by 20 mM salt;
the model, holding a large salt-independent $K$, retains too much shift at high salt.
The simulated potential of mean force, in contrast, screens correctly (Table S-3b),
so the discrepancy is the reduced fixed-$K$ description, not the underlying physics.

Table S-3b. pH-6 potential-of-mean-force well depth versus salt (simulation).

| $I$ / M | well depth / $k_BT$ |
|---------|---------------------|
| 0.005   | $-15.1$ |
| 0.03    | $-8.4$  |
| 0.15    | $-2.6$  |

## S-4. Simulation composites versus salt

Figures S (below) reproduce the low-salt composite at 0.03 and 0.15 M. The
potential-of-mean-force wells shrink from $\sim -15$ to $\sim -2.6\ k_BT$ and the
electric potential is attenuated as salt increases, screening the attraction that
drives adsorption.

![Simulation composite at $I = 0.03$ M (panels as in the low-salt composite).](../mc/composite_I0.03.png)

![Simulation composite at $I = 0.15$ M. The measured $\zeta$ marker is the highest available salt (100 mM added NaCl); the series does not extend to 150 mM, so it overstates $|\zeta|$ relative to the simulation, which over-screens at high salt.](../mc/composite_I0.15.png)

## S-5. NMR bound fraction and free-phytate hydrodynamics

### S-5.1. Bound fraction from the free energy

The NMR comparison uses no parameter fitted to
the NMR data. Integrating the simulated potential of mean force $w(r)$ over the
adsorption region gives the surface partition length (a Henry constant with units of
length),

$$K_H=\int\big(e^{-w(r)/k_BT}-1\big)\,dr ,$$

evaluated at $I=0.15$ M. The PMF used here includes the same fixed short-range contact
term used throughout the adsorption analysis (a salt-robust $2.5\,k_BT$ well at the
brush-contact position); this term was set before comparison with the NMR relaxation
data. The areal coverage $\Gamma$ follows from a two-dimensional
generalized-van-der-Waals isotherm (Nordholm et al., 2018). In the dilute (Henry) limit
$\Gamma\to K_H\,n_B$
($n_B$ the bulk phytate number density), but at the measured 7.6 mM phytate the adsorbed
anions crowd and repel, so the coverage fraction $\theta=\Gamma d^2$ (phytate hard-disk
diameter $d$) solves the implicit isotherm

$$\theta=\phi_0\,(1-\theta)\,\exp\!\Big(\alpha\,\theta-\frac{\theta}{1-\theta}\Big),
\qquad \phi_0=K_H\,n_B\,d^2 ,$$

with $\phi_0$ the Henry-limit coverage and $\alpha<0$ the dimensionless lateral
coefficient carrying the like-charge phytate–phytate repulsion. That repulsion is the
screened Coulomb interaction between adsorbed multivalent anions, which diverges at
contact in mean field. Forsman's (2004) correlation correction removes the
divergence: below a correlation-hole radius $R_c$ (the nearest-neighbour spacing
$\sim 1/\sqrt{\pi\Gamma}$, floored at the hard core $d$ and capped at $\sim$two Debye
lengths) the pair repulsion is replaced by its tangent, so the layer packs without the
singularity. Because $R_c$ depends on the coverage, $\alpha(\Gamma)$
and $\theta$ are solved together by fixed-point iteration. The bound fraction is then
$f_\text{bound}=\Gamma a/n_B$, where $a=6\phi_\text{OB}/d_{43}$ is the oil-body specific
surface area ($\phi_\text{OB}$ the droplet volume fraction, $d_{43}$ the reported
volume-weighted mean diameter); it reduces to $a K_H$ only as $\theta\to0$.
In the fast-exchange limit, the observed slow-motion contrast is the population-weighted
average of the free and bound states,

$$\Delta R_{1\rho}^{\text{obs}}
=f_\text{bound}\Delta R_\text{bound}
+f_\text{free}\Delta R_\text{free},\qquad
f_\text{free}=1-f_\text{bound}.$$

The bound fraction required by the NMR data is therefore

$$f_\text{req}
=\frac{\Delta R_{1\rho}^{\text{obs}}-\Delta R_\text{free}}
{\Delta R_\text{bound}-\Delta R_\text{free}}.$$

The simulated $f_\text{bound}$ enters this relation with no parameter fitted to the NMR
data, and matches $f_\text{req}$ at every pH (Table 1). Because $K_H$ is read at the
measured salt, nothing is transferred across ionic strength.

### S-5.2. Free-tumbling reference

The remainder of this section details the reference
$\Delta R_\text{free}$.

A bare-hydrodynamic estimate (RDKit conformer ensembles, Landrum et al., with
Rotne–Prager–Yamakawa bead mobilities — the GRPY method of Zuk et al., 2018, via the
PyGRPY package, Waszkiewicz) fixes the reference. The
rotational correlation time $\tau_r$ and its anisotropy vary little across protonation
states (Table S-5), so different phytate conformers tumble alike and a single
$\Delta R_\text{free}$ is well defined. This is a purely geometric measure: the
mobilities are computed from the bare molecular shape with stick boundary conditions
and no ionic-atmosphere friction, so electrostatics are neglected — net charge
enters only through the shape (removal of peripheral protons), not through any
electrostatic friction. It is therefore a rough estimate, adequate only for
establishing that the free reference is conformer-insensitive.

Table S-5. Free-phytate rotation by net charge (geometric estimate, electrostatics
neglected).

| charge | $\tau_r$ / ps | anisotropy | $a_H$ / Å |
|--------|---------------|------------|-----------|
| 0   | 143 | 1.21 | 5.28 |
| $-3$ | 134 | 1.19 | 5.18 |
| $-6$ | 137 | 1.19 | 5.21 |
| $-9$ | 147 | 1.17 | 5.33 |
| $-12$| 164 | 1.14 | 5.55 |

## S-6. Calorimetry

Both models fit the pH-4 thermogram well. The one-site (independent-sites) model returns
$n=0.38$ phytate per oleosin, and the charge-regulated model — used here in its second
role, as an adsorption isotherm fit *freely* to the heat (affinity $K$, surface area and
enthalpy all free) — reproduces it about as well (RMS 6.6 versus 5.9 µJ, Figure S). Yet
both imply the same physically incompatible monolayer coverage.

![pH-4 ITC thermogram: experiment (replicate mean $\pm$ SD), the charge-regulated model (solid), and the one-site model (dashed).](../itc/fig_itc_crpb.png)

The heat demands $\approx 0.38$ phytate bound per oleosin — the one-site $n$, and the
amount the charge-regulated fit also needs to match the enthalpy. But the surface charge
allowed by the $\zeta$ reversal and the NMR-scale bound fraction caps the surface at
$\approx 0.03$ phytate per oleosin: a self-regulating monolayer cannot hold the ~13× more
phytate the heat requires.
The two models absorb this gap in different ways (Table S-6). The
charge-regulated fit inflates the peptide-to-oleosin area to $\approx 7\times10^{3}$ Å²,
21× the physical grafting density (338 Å² per oleosin, the same value in the
simulated brush and in the measured $\sim 1$ µmol oleosin g⁻¹), and the affinity and
area are strongly anti-correlated ($-0.93$), so the fit cannot separate binding
strength from capacity. The one-site model, which has no surface constraint, represents
the discrepancy as an apparent stoichiometry. Either way the excess heat cannot be
monolayer binding;
the calorimetry's own multiphasic shoulder (Section 3.3.2) points to inter-droplet cross-linking
as the most likely source. The one-site thermodynamic parameters should therefore be
interpreted as apparent descriptors of a multiphasic process. A non-aggregating ITC
experiment could use the CR-PB isotherm framework to quantify surface binding.

Table S-6. Free fits of the two models to the pH-4 thermogram (uncertainties are $\pm1$
standard error from the fit covariance; the physical grafting area 338 Å² is *not* a
fit).

|                       | one-site (independent) | charge-regulated (CR-PB) |
|-----------------------|------------------------|--------------------------|
| adsorption            | independent sites      | charge-regulated surface |
| binding $K$           | $K_a=(2.5\pm0.8)\times10^{5}$ M⁻¹ | $K=4$ M⁻¹ (fit) |
| capacity              | $n=0.38\pm0.03$        | area $\approx 7\times10^{3}$ Å² (fit) — 21× the physical 338 Å² |
| $\Delta H$ / kJ mol⁻¹ | $-47\pm5$         | $-38$ |
| largest $|$parameter correlation$|$ | 0.88 ($K_a$–$\Delta H$) | 0.93 ($K$–area) |
| RMS / µJ              | 5.9                    | 6.6 |

---

## Supporting Information Available

All simulation input files, modelling and analysis scripts are available at
<https://github.com/mlund/si-hempoil-phytate>.

## References

- Ashbaugh, H. S., & Hatch, H. W. (2008). Natively unfolded protein stability as a
  coil-to-globule transition in charge/hydropathy space. *J. Am. Chem. Soc.* **130**,
  9536–9542. https://doi.org/10.1021/ja802124e
- Cao, F., von Bülow, S., Tesei, G., & Lindorff-Larsen, K. (2024). A coarse-grained
  model for disordered and multi-domain proteins. *Protein Sci.* **33**, e5172.
  https://doi.org/10.1002/pro.5172
- De Stefano, C., Milea, D., & Sammartano, S. (2003). Speciation of phytate ion in
  aqueous solution: protonation constants in tetraethylammonium iodide and sodium
  chloride. *J. Chem. Eng. Data* **48**, 114–119. https://doi.org/10.1021/je020124m
- Forsman, J. (2004). A simple correlation-corrected Poisson–Boltzmann theory.
  *J. Phys. Chem. B* **108**, 9236–9245. https://doi.org/10.1021/jp049571u
- Grahame, D. C. (1947). The electrical double layer and the theory of
  electrocapillarity. *Chem. Rev.* **41**, 441–501.
  https://doi.org/10.1021/cr60130a002
- Grimsley, G. R., Scholtz, J. M., & Pace, C. N. (2009). A summary of the measured
  pK values of the ionizable groups in folded proteins. *Protein Sci.* **18**,
  247–251. https://doi.org/10.1002/pro.19
- Grosberg, A. Yu., Nguyen, T. T., & Shklovskii, B. I. (2002). Colloquium: The physics
  of charge inversion in chemical and biological systems. *Rev. Mod. Phys.* **74**,
  329–345. https://doi.org/10.1103/RevModPhys.74.329
- Landrum, G., et al. RDKit: Open-source cheminformatics.
  https://doi.org/10.5281/zenodo.591637
- Nordholm, S., Forsman, J., Woodward, C., Freasier, B., & Abbas, Z. (2018).
  *Generalized van der Waals Theory of Molecular Fluids in Bulk and at Surfaces*.
  Elsevier. ISBN 978-0-12-811136-9.
- Stenqvist, B., Thuresson, A., Kurut, A., Vácha, R., & Lund, M. (2013). Faunus — a
  flexible framework for Monte Carlo simulation. *Mol. Simul.* **39**, 1233–1239.
  https://doi.org/10.1080/08927022.2013.828207
- Tesei, G., Schulze, T. K., Crehuet, R., & Lindorff-Larsen, K. (2021). Accurate model
  of liquid–liquid phase behavior of intrinsically disordered proteins from
  optimization of single-chain properties. *Proc. Natl. Acad. Sci. USA* **118**,
  e2111696118. https://doi.org/10.1073/pnas.2111696118
- Trefalt, G., Behrens, S. H., & Borkovec, M. (2016). Charge regulation in the
  electrical double layer: ion adsorption and surface interactions. *Langmuir* **32**,
  380–400. https://doi.org/10.1021/acs.langmuir.5b03611
- Waszkiewicz, R. PyGRPY: a Python implementation of the generalized
  Rotne–Prager–Yamakawa bead mobility. https://github.com/RadostW/PyGRPY
- Zuk, P. J., Cichocki, B., & Szymczak, P. (2018). GRPY: an accurate bead method for
  calculation of hydrodynamic properties of rigid biomacromolecules. *Biophys. J.*
  **115**, 782–800. https://doi.org/10.1016/j.bpj.2018.07.015
