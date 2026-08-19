---
title: "Simulation and theory for the OB–PA manuscript"
---

# Part A — Main text

## 2.2.10. Computational simulations

We used coarse-grained Monte Carlo (MC) simulations to model phytate binding to an
N-terminal oleosin brush. The model represented the N-terminal arm of *Cannabis
sativa* oleosin (UniProt A0A803P977) with one bead per residue and the CALVADOS 3
force field (Tesei et al., 2021; Cao et al., 2024). We grafted 50 chains to one wall
of a slit to represent the cationic protein layer at the oil-body surface. A rigid
seven-bead phytate comprised one inositol centre and six phosphates. We sampled
phytate and peptide protonation in the semi-grand-canonical ensemble, so the signed
phytate valence $Z_{\mathrm{PA}}=q_{\mathrm{PA}}/e$ responded to pH.

We sampled pH 4.0, 6.0, and 10.0 at 5–150 mM NaCl with Faunus 3.0.0-beta.2
(Git revision v3.0.0-beta.2-12-g75193d4; Stenqvist et al., 2013). As functions of
distance $r$ from the grafted wall, we calculated the bulk-referenced potential of
mean force $w(r)$, electric potential $\varphi(r)$, phytate valence
$Z_{\mathrm{PA}}(r)$, wall alignment $S_\text{wall}$, and ensemble orientational
order $S^2_\text{ens}$ of P–O proxy vectors. Because equilibrium MC has no physical
clock, $S^2_\text{ens}$ is not a rotational correlation time.

We then extended the one-phytate simulation to finite surface coverage with a
lateral-correlation adsorption-regulated Poisson–Boltzmann model (LC–AR-PB). The MC
simulation supplied the peptide-only potential and mean phytate charge at the binding
coordinate. The measured bare-surface $\zeta$-potential set the diffuse-layer surface
charge, and one intrinsic affinity $K$ per pH represented the remaining contact free
energy. Nonlinear PB described the collective mean-field electrostatics. A matched
residual free energy described lateral correlations without counting the full screened
pair interaction twice. Supplementary Section S-2 gives the assumptions, equations,
parameters, and validity range. We fitted the same grand-potential isotherm separately
to the pH-4 isothermal titration calorimetry (ITC) thermogram, with $K$, area per
oleosin, and adsorption enthalpy as free parameters (Supplementary Section S-6).

## 3.3.4. Molecular characterization of oleosin-PA binding

The simulations resolve a local binding-plane potential that is not apparent from the
net oil-body charge. This local field provides attraction near pH 6, where the net
charge is small, and permits overcharging at pH 4. Figure 10 shows the low-salt MC
profiles underlying this interpretation.

![Figure 10. MC profiles at 5 mM ionic strength versus distance $r$ from the grafted wall: (A) potential of mean force $w(r)$; (B) ensemble P–O proxy order $S^2_\mathrm{ens}$ (solid) and whole-molecule wall-pose order $S_\mathrm{wall}$ (dashed); (C) electric potential, simulated $\zeta$ at the shear plane (open marker), and measured $\zeta$ (dashed line, shaded band $\pm$ one standard deviation over replicates; drawn across the panel because the measurement has no $r$ coordinate); and (D) phytate charge. The shear plane is the outermost position where the peptide fixed-charge density exceeds 2% of its peak, 36–43 Å from the wall across pH. Panel B reports equilibrium order, not rotational rates.](../mc/composite_I0.005.png)

The ionizable oleosin termini make the brush periphery more positive than the surface
average. The peptide-only binding potential is $+62$ mV at pH 4 and remains positive
at pH 6 ($\approx +19$ mV), where the net oil-body charge is near zero. The local
cationic field can therefore attract phytate near the isoelectric point. LC–AR-PB uses
the measured bare-surface $\zeta$-potential to set the net charge of the diffuse layer
and the simulated peptide-only potential to set local one-particle enrichment. The
external field excludes the probe phytate's contribution.

Four terms determine the surface coverage: attraction to the peptide binding-plane
potential, nonlinear PB repulsion, site exclusion, and discrete lateral correlations.
The PB repulsion increases with phytate charge and weakens with salt screening. The
residual correlation term corrects short-range PB over-repulsion; it does not add a
second full pair interaction. Thus, pH and salt control both attraction to the brush
and the net adsorbed charge. At pH 6, the larger MC phytate valence strengthens the
electrostatic feedback and limits overcharging.

At pH 4, the fitted coverage passes the neutralization point. Each adsorbed phytate
carries $Z_{\mathrm{PA}}\approx-7$ at the surface, compared with $\approx-6$ in bulk
(Figure 10D). The surface charge therefore changes sign when the adsorbed charge exceeds the bare
positive charge. The fitted charge inversion describes the measured $\zeta$-potential
shift from $+38$ to $-4$ mV as phytate is added (Figure 6A). The inversion is
self-limiting because a negative diffuse-layer potential opposes further anion
binding. LC–AR-PB fits the crossing with the lateral correlation correction included
(Grosberg et al., 2002; Supplementary Section S-3). The resulting charge reduction is
also consistent with the earlier pH-4 destabilization and possible PA-mediated
inter-droplet bridging.

The three experiments constrain different parts of this interpretation.
$\zeta$-potential reports the electrokinetic response, NMR relaxation probes the local
phytate environment, and ITC measures the total heat of interaction. A single-field
NMR measurement cannot determine an adsorbed population or rotational correlation
time.

The low-salt fit captures the pH-4 reversal and the pH-6 shift towards more negative
$\zeta$. Without refitting $K$, the model also captures the weakening of this shift
with added salt. Over the same range, the simulated pH-6 PMF well decreases from
$-13.9$ to $-2.4\,k_BT$ between 5 and 150 mM salt (Supplementary Section S-4). Both
results support salt screening of adsorption. At pH 10, the simulated interaction is
repulsive, and adsorption is not identifiable from the $\zeta$ fit. The high-phytate
destabilization at pH 10 is therefore consistent with double-layer compression rather
than direct adsorption.

The NMR relaxation data show a pH-dependent perturbation of phytate. At the single
24 kHz spin-lock field, $R_{1\rho}-R_1$ contains transverse relaxation,
spectral-density, and possible exchange contributions (Quinn & McDermott, 2012). The
OB-induced increment,
$[R_{1\rho}-R_1]_\mathrm{OB-PA}-[R_{1\rho}-R_1]_\mathrm{aq}$, is $397\pm51$,
$302\pm14$, and $46.9\pm2.1$ s⁻¹ at pH 4, 6, and 10, respectively. The relaxation
perturbation is therefore larger at pH 4 and 6 than at pH 10. This pH dependence
qualitatively follows the simulated interaction: the MC PMF is attractive at pH 4 and
6 but repulsive at pH 10. The comparison assigns neither an adsorbed fraction nor a
rotational timescale. Quantitative population analysis would require a bound-state
relaxation function and an exchange model (Supplementary Section S-5).

The Widom calculation provides a separate equilibrium comparison. At the nominal
binding plane and $I=0.15$ M, $S^2_\mathrm{ens}$ is
$3\times10^{-5}$–$7\times10^{-5}$, and the conditional orientational-cage cost is
$0.9$–$1.5\,k_BT$. The one-particle model therefore gives little persistent alignment
and a modest angular cage. These quantities neither predict $R_{1\rho}$ nor exclude
slow isotropic reorientation. Independently, the lack of a detectable
cross-polarization signal argues against a dominant rigid phytate population.

The calorimetry tests whether monolayer adsorption can account for the measured heat.
The one-site instrument model and the unconstrained LC–AR-PB surface fit describe the pH-4
thermogram with root-mean-square (RMS) residuals of 5.2 and 6.1 µJ, respectively
(Supplementary Section S-6). The one-site model gives $0.34$ phytate per oleosin.
LC–AR-PB reaches a similar maximum of $0.41$ only by fitting
$7.40\times10^3$ Å² per oleosin, 22 times the
physical grafting area. At the physical area of 338 Å², LC–AR-PB gives
$0.02$–$0.03$ phytate per oleosin, consistent with the scale inferred from
$\zeta$-potential. Thus, the LC–AR-PB curve fit does not support a physical monolayer
interpretation. Off-surface processes likely contribute to the heat; one candidate is
inter-droplet cross-linking, as suggested by the shoulder near 0.3 PA:oleosin
(Section 3.3.2). The fit parameters are therefore diagnostics of a multiphasic pH-4
process, not direct monolayer constants. A non-aggregating ITC experiment could apply
the same framework with an independently measured, fixed area (Supplementary Section
S-6).

The lateral closure remains a reduced two-dimensional approximation. It includes a
density-dependent correlation hole but holds the single-probe peptide field and
phytate charge fixed as coverage grows. The model therefore omits brush restructuring
and charge regulation caused by neighbouring phytates. Any cooperative kink from the
residual correlation term must be tested against the no-correlation baseline and
alternative closures (Supplementary Section S-2).

## CONCLUSIONS

The simulations and LC–AR-PB support an electrostatic interpretation of the observed
behaviour. Phytate binding depends on attraction to the peptide-only binding-plane
potential and on collective and correlated interactions within the adsorbed layer. pH
and salt control both contributions. At pH 4, the fitted coverage passes
neutralization and describes the high-PA $\zeta$ reversal. The earlier
destabilization is consistent with charge reduction and possible PA-mediated
inter-droplet bridging. At pH 6, weaker and salt-sensitive binding adds negative
charge, consistent with the observed stabilization. At pH 10, direct binding is not
resolved, and double-layer compression is consistent with the observed behaviour.

The three probes constrain different parts of this interpretation. LC–AR-PB fits the
$\zeta$ reversal with a PB-matched lateral correlation correction. NMR spectra and
rotating-frame relaxation show larger changes at pH 4 and 6 than at pH 10, but they do
not determine the adsorbed fraction. Both ITC models fit the thermogram,
yet the LC–AR-PB fit requires 21 times the physical area. The thermogram therefore
cannot be assigned to a physical monolayer within this model, and off-surface processes
likely contribute. The main model uncertainty at low salt is the reduced
correlation-hole closure, which neglects coverage-induced changes in the peptide field
and phytate charge.

---

# Part B — Supplementary Information

## S-1. Monte Carlo simulations

### S-1.1. System and grafting

We represented the 14-residue N-terminal arm of *Cannabis sativa* oleosin (UniProt
A0A803P977, sequence MKSFLPESGPSAQQ) with one bead per residue and CALVADOS 3
(Tesei et al., 2021; Cao et al., 2024). Harmonic bonds joined consecutive beads with
an equilibrium length of 3.8 Å and a force constant of 80.33 kJ mol⁻¹ Å⁻². We
grafted 50 chains at $3\times10^{-3}$ Å⁻² to the bottom wall of a
$130\times130\times200$ Å³ slit. The cell was periodic in $x$ and $y$, with hard
walls at $z=\pm100$ Å. A harmonic potential tethered each terminal anchor bead to the
bottom wall with a force constant of 0.5 kJ mol⁻¹ Å⁻². Water was implicit and entered
through its temperature-dependent dielectric permittivity at 298.15 K.

### S-1.2. Peptide titration

We sampled four protonation reactions by atom-type swaps against an implicit proton
bath, $a_{\mathrm{H}^+}=10^{-\mathrm{pH}}$: Glu (p$K_a=4.24$), Lys (10.4), the
N-terminus (8.0), and the C-terminus (3.67) (Grimsley et al., 2009). The mean chain
charge changed from $\approx+1e$ at pH 4 to $\approx-1.3e$ at pH 10. The brush
therefore passed through its isoelectric region near pH 6, where a local cationic field
remained despite the small net charge (Section 3.3.4).

### S-1.3. Phytate titration

We represented one phytate (*myo*-inositol hexakisphosphate) as a rigid seven-bead
unit. One inositol-centre bead (INO; diameter 6.2 Å) was surrounded by six phosphate
beads (PH; 5.8 Å) in an idealized chair geometry. The phosphate beads were
$\approx4.5$ Å from the ring centre, with one axial and five equatorial. We scaled the
bead diameters from CALVADOS 3 as $M^{1/3}$.

The twelve dissociable protons gave 13 charge states from neutral H₁₂PA to PA¹²⁻.
Molecular swap moves sampled all twelve dissociations in the semi-grand-canonical
ensemble. We used macroscopic p$K_a$ values measured potentiometrically in NaCl: 1.1,
1.5, 1.7, 2.1, 2.5, and 2.9 for the first proton on each phosphate, and 5.72, 6.28,
6.81, 7.60, 9.94, and 11.84 for the second (De Stefano et al., 2003). This p$K_a$
ladder increased the mean phytate valence from $Z_{\mathrm{PA}}\approx-6$ at pH 4 to
$Z_{\mathrm{PA}}\approx-10.6$ at pH 10 and allowed the valence to respond to the
peptide environment. For each integer charge state, charges occupied the maximally
separated subset of phosphate beads. Above $-6e$, every phosphate carried $-1e$, and
the additional charges followed the same placement rule. All 13 states shared one
geometry, so each swap changed only the charge and preserved position and orientation
through gyration-tensor alignment. Thus, one phytate sampled the full protonation
manifold without combining separate fixed-charge simulations.

### S-1.4. Interactions

Nonbonded pairs combined the short-range, hydrophobicity-scaled Ashbaugh–Hatch
potential (Ashbaugh & Hatch, 2008) with a salt-dependent Debye–Hückel term.
Phytate–peptide cross-pairs used $\lambda=0$, corresponding to a repulsive
Weeks–Chandler–Andersen potential. This choice prevented arithmetic
$\lambda$-mixing from introducing a hydrophobic attraction between phytate and
hydrophobic amino acids. Peptide nearest neighbours retained their Coulomb interaction
but were excluded from the short-range nonbonded term. All intramolecular phytate pairs
were excluded because their effective interaction was already represented by the
macroscopic p$K_a$ ladder.

### S-1.5. Monte Carlo moves

Each propagation cycle made 40 stochastic move selections. The available moves were
peptide translation in $xy$ ($\delta_\text{max}=10$ Å), peptide pivot rotation
(0.5 rad), phytate translation along $z$ (5 Å; 20 repeats), and rigid phytate rotation
(1.0 rad). Molecular swap moves sampled the twelve phytate and four peptide
protonation reactions with 10 repeats per selection (Stenqvist et al., 2013).

### S-1.6. Potential of mean force

We obtained the PMF along the inositol $z$-coordinate by multi-walker umbrella
sampling. The calculation used 38 overlapping hard-wall windows, each 6 Å wide and
spaced by 2.5 Å, from the grafted wall to the slit centre. A harmonic drive with a
force constant of 5.0 kJ mol⁻¹ Å⁻² moved phytate into each window. We then removed
the drive and sampled an unbiased distribution conditional on the hard-wall window.
For adjacent windows, we calculated the free-energy offset from
$\Delta F_{i\to i+1}=-k_BT\ln(f_i/f_{i+1})$, where $f_i$ is the fraction of
window-$i$ samples in the shared overlap. We computed PMFs at pH 4, 6, and 10 and at
$I=0.005$, 0.03, and 0.15 M.

### S-1.7. Position-resolved equilibrium orientational diagnostics

In each umbrella window, Widom trial rotations gave the ensemble second-rank order
$S^2_\mathrm{ens}$ of six phosphate-to-inositol-centre vectors, used as coarse P–O
orientation proxies. The estimator first formed one ensemble second-moment tensor and
then evaluated its rotational invariant. Consequently, $S^2_\mathrm{ens}$ approaches
zero when the preferred local axis varies among brush configurations, even if each
configuration orients the molecule locally. The same analysis gave the conditional
orientational-cage free energy $W_\mathrm{rot}$ relative to the best trial orientation
in each configuration.

Both quantities are equilibrium observables and contain no physical time. They cannot
be converted to rotational diffusion or $R_{1\rho}$. The separate whole-molecule
wall-pose order $S_\text{wall}$ equals $0$ for isotropic, $-\tfrac12$ for face-on,
and $+1$ for edge-on orientations. Its mean remained close to zero at the binding
plane. The umbrella PMF already contains the orientational partition function;
therefore, adding $W_\mathrm{rot}$ to $w(r)$ would count orientational entropy twice.

### S-1.8. Electric potential

Let $\sigma(z)$ be the laterally averaged slab charge density in $e$ Å$^{-2}$. Its
screened planar-slab sum gives the reduced potential

$$
y(z)=\frac{e\varphi(z)}{k_BT}
=\frac{2\pi l_B}{\kappa}\sum_{z'}\sigma(z')\,e^{-\kappa|z-z'|}.
$$

Here, $l_B$ is the Bjerrum length and $1/\kappa$ is the Debye length. We converted
$y$ to $\varphi=(k_BT/e)y$, referenced $\varphi$ to zero in bulk, and evaluated the
simulated $\zeta$ at a model shear plane. We defined this plane as the outermost $r$
at which the magnitude of the peptide charge density remained above 2% of its peak
(open markers in composite panel C).

The cell was fixed and non-periodic in $z$, matching the open-boundary axial kernel;
we applied no finite-disk correction. Because the shear plane is phenomenological, we
used the simulated $\zeta$ values as profile diagnostics rather than as a complete
hydrodynamic prediction.

## S-2. Lateral-correlation adsorption-regulated PB model

The single-phytate MC simulation resolves the local peptide field but cannot represent
finite lateral coverage. We therefore coupled it to LC–AR-PB to calculate
$\zeta$-potential and, separately, to test a surface interpretation of the ITC data.
The model represents the oil-body surface as a planar Gouy–Chapman diffuse layer in a
1:1 electrolyte. It separates the local peptide field acting on one phytate at the
brush edge from the collective diffuse-layer response to adsorbed charge.

The molecular coupling is one-way: MC supplies a fixed local peptide field and a
fixed phytate valence to the continuum model, but the continuum model does not update
the MC simulation. Self-consistency refers instead to the surface problem. For each
trial coverage, PB calculates the charging free energy from the corresponding
effective surface charge; that free energy changes the equilibrium coverage. We solve
this feedback by minimizing one surface grand potential globally (Figure S2).

Mathematically, LC–AR-PB is a reduced surface classical density-functional theory
(cDFT): its sole density variable is the laterally uniform $\Gamma$, and PB integrates
out the diffuse ions. It is not a nonlocal cDFT for lateral density variations.

![Figure S2. LC–AR-PB construction. Blue denotes fixed MC input, purple experimental and prescribed surface inputs, grey the bare PB state, and red the coverage-dependent PB response. (A) The one-phytate MC calculation supplies the peptide-only potential at $R_\beta$, the minimum of the sampled PMF. Small fluctuations in the blue schematic curve indicate stochastic MC sampling and do not represent measured uncertainty. The planar PB profiles begin at an effective diffuse/shear plane and illustrate the coverage-dependent continuum response. The model uses the PB boundary-potential difference rather than adding the two radial curves point by point. (B) MC inputs enter the surface grand potential without a return loop. Self-consistency is restricted to the feedback among coverage $\Gamma$, effective charge $\sigma_{\mathrm{eff}}$, and the PB charging term. Global minimization selects the stable solution.](../zeta/fig_lc_arpb_schematic.png)

### S-2.1. MC input and stated one-particle approximation

We defined the binding coordinate as $R_\beta=26$ Å from the grafting wall. The
external field was the peptide-only potential
$\psi_{\mathrm{pep}}^{MC}(R_\beta)$, excluding the probe phytate's contribution. We
averaged the peptide-only profiles across umbrella windows to reduce noise. This
approximation assumes that one dilute phytate perturbs the mean brush profile
negligibly. At the same coordinate, the MC charge profile gave the positive phytate
charge magnitude $z=-q_{\mathrm{PA}}/e$. We interpolated both quantities linearly in
$\log I$ between the three simulated ionic strengths.

This MC-to-continuum handoff includes phytate charge regulation by pH and by the
one-particle peptide environment. It then holds the peptide field and $z(R_\beta)$
fixed as coverage increases. The model therefore omits charge regulation and brush
restructuring induced by neighbouring adsorbed phytates.

### S-2.2. Bare charge and common surface grand potential

We inferred the bare surface charge density $\sigma_0$, in $e$ Å$^{-2}$, from the
measured zero-PA $\zeta$-potential at $I_0=5$ mM and held $\sigma_0$ fixed across
salt concentrations. For the reduced diffuse-layer potential
$y_d=e\psi_d/(k_BT)$, the Grahame relation is (Grahame, 1947)

$$
\sigma=\frac{\kappa}{2\pi l_B}\sinh\!\left(\frac{y_d}{2}\right),
\qquad
y_d(\sigma)=2\,\operatorname{asinh}\!\left(\frac{2\pi l_B\sigma}{\kappa}\right),
\qquad\text{(S1)}
$$

where $l_B$ is the Bjerrum length and $\kappa$ is the inverse Debye length. For an
adsorbed number density $\Gamma$,

$$
\sigma_{\mathrm{eff}}=\sigma_0-z\Gamma,
\qquad
\zeta=\frac{k_BT}{e}\,y_d(\sigma_{\mathrm{eff}}).
\qquad\text{(S2)}
$$

Adsorption thus changes both the effective surface charge and the diffuse-layer
potential. This coupling is the adsorption-regulated boundary condition (Trefalt et
al., 2016). The mapping equates the Gouy–Chapman boundary potential with the measured
$\zeta$-potential; it does not calculate a hydrodynamic offset between the PB plane
and the shear plane.

To distinguish the two electrostatic inputs, we define the incremental diffuse-layer
response and the effective binding-plane potential as

$$
\Delta y_d(\Gamma)=y_d(\sigma_{\mathrm{eff}})-y_d(\sigma_0),
\qquad
y_\beta(\Gamma)=y_{\mathrm{pep}}^{MC}(R_\beta)+\Delta y_d(\Gamma).
\qquad\text{(S3)}
$$

Thus, $\psi_\beta=(k_BT/e)y_\beta$ is not the one-phytate MC potential alone. It is
the fixed local MC field plus the change in the continuum potential caused by
adsorption. The subtraction of the bare PB potential also prevents the measured
bare-surface monopole field from being added twice.

We obtained the coverage by minimizing one surface grand potential, which avoids
applying sequential electrostatic corrections. Let
$y_{\mathrm{pep}}^{MC}=e\psi_{\mathrm{pep}}^{MC}/(k_BT)$ be the reduced peptide
field, $c_L$ the free ligand concentration, and $K$ the intrinsic affinity in M⁻¹.
For surface area $A$, $\beta=1/(k_BT)$, and fractional coverage
$\theta=\Gamma/\Gamma_{\max}$, the grand-potential density is

$$
\frac{\beta\Omega}{A}
=\Gamma_{\max}\left[\theta\ln\theta+(1-\theta)\ln(1-\theta)\right]
-\Gamma\ln(Kc_L)
-z\,y_{\mathrm{pep}}^{MC}(R_\beta)\Gamma
+g_{\mathrm{PB}}(\Gamma)
+\beta f_{\mathrm{corr,res}}(\Gamma).
\qquad\text{(S4)}
$$

The five terms describe Langmuir mixing and site exclusion, intrinsic adsorption,
the one-phytate MC field, nonlinear PB charging, and residual lateral correlations,
respectively. The product $Kc_L$ is dimensionless. We referenced the PB term to its
tangent at the bare state:

$$
g_{\mathrm{PB}}(\Gamma)
=F_{\mathrm{PB}}(\sigma_{\mathrm{eff}})-F_{\mathrm{PB}}(\sigma_0)
-y_d(\sigma_0)(\sigma_{\mathrm{eff}}-\sigma_0),
\qquad
\frac{\partial F_{\mathrm{PB}}}{\partial\sigma}=y_d(\sigma).
\qquad\text{(S5)}
$$

For the planar Grahame relation, the corresponding primitive is

$$
F_{\mathrm{PB}}(\sigma)
=2\left[
\sigma\,\operatorname{asinh}(a\sigma)
-\frac{\sqrt{1+(a\sigma)^2}}{a}
\right],
\qquad a=\frac{2\pi l_B}{\kappa}.
\qquad\text{(S6)}
$$

Differentiation with respect to coverage gives

$$
\frac{\partial g_{\mathrm{PB}}}{\partial\Gamma}
=-z\,[y_d(\sigma_{\mathrm{eff}})-y_d(\sigma_0)].
\qquad\text{(S7)}
$$

Thus, the same nonlinear PB response determines the adsorption equilibrium and the
reported $\zeta$-potential. Subtracting the bare-state tangent prevents the monopole
field inferred from the measured bare $\zeta$ from being added again to the
peptide-only MC binding field.

The stationarity condition makes the feedback explicit. If
$\beta\mu_{\mathrm{corr,res}}=\partial(\beta f_{\mathrm{corr,res}})/\partial\Gamma$,
then an interior equilibrium obeys

$$
\frac{\theta}{1-\theta}
=Kc_L\exp\!\left\{
z\left[y_{\mathrm{pep}}^{MC}(R_\beta)+\Delta y_d(\Gamma)\right]
-\beta\mu_{\mathrm{corr,res}}(\Gamma)
\right\}.
\qquad\text{(S8)}
$$

Equation (S8) could be solved by fixed-point iteration: a trial $\Gamma$ determines
$\sigma_{\mathrm{eff}}$ and $\Delta y_d$, which update the adsorption probability.
The implementation instead minimizes Equation (S4) globally because the lateral
closure can produce more than one stationary branch. This procedure retains the MC
field in every free-energy evaluation while selecting the thermodynamically stable
coverage.

### S-2.3. PB-matched lateral correlation free energy

PB contains the long-wavelength mean-field repulsion between adsorbed charges but
omits their discrete lateral correlation hole. Adding the complete screened pair
interaction would count the mean-field repulsion twice. LC–AR-PB therefore adds only
the residual free-energy density

$$
\beta f_{\mathrm{corr,res}}(\Gamma)
=\frac{\Gamma^2}{2}
\left[J_{\mathrm{F}}(\Gamma)-J_{\mathrm{Y,MF}}\right].
\qquad\text{(S9)}
$$

The hard-core Yukawa mean-field reference is

$$
J_{\mathrm{Y,MF}}=\int_{D_{\mathrm{PA}}}^\infty 2\pi r\,
u_Y(r)\,dr,
\qquad
u_Y(r)=z^2l_B\frac{e^{-\kappa r}}{r}.
\qquad\text{(S10)}
$$

Here, $D_{\mathrm{PA}}=10.6$ Å is the phytate hard-core diameter. To obtain
$J_{\mathrm{F}}$, we used the Forsman tangent below the correlation-hole radius and
the Yukawa interaction above it (Forsman, 2004):

$$
u_{\mathrm{F}}(r;\Gamma)=
\begin{cases}
u_Y(R_c)+u_Y^{\prime}(R_c)(r-R_c), & D_{\mathrm{PA}}\le r<R_c,\\
u_Y(r), & r\ge R_c,
\end{cases}
\qquad
J_{\mathrm{F}}=\int_{D_{\mathrm{PA}}}^{\infty}2\pi r\,
u_{\mathrm{F}}(r;\Gamma)\,dr.
\qquad\text{(S11)}
$$

The cutoff depends on the solved coverage:

$$
R_c(\Gamma)=\max\!\left[
D_{\mathrm{PA}},
\min\!\left(\frac{2}{\kappa},\frac{1}{\sqrt{\pi\Gamma}}\right)
\right].
\qquad\text{(S12)}
$$

The cutoff cannot fall below the phytate core. At low coverage it is capped at two
Debye lengths; at higher coverage it follows the Wigner–Seitz radius
$1/\sqrt{\pi\Gamma}$. The residual vanishes as $\Gamma^2$ at zero coverage and also
when $R_c=D_{\mathrm{PA}}$.

The tangent softens the short-range PB over-repulsion, so
$J_{\mathrm{F}}-J_{\mathrm{Y,MF}}\leq0$. The negative residual is a correction to the
PB repulsion, not an independent attraction. We minimized the complete
density-dependent free energy globally, thereby including the dependence of
$R_c(\Gamma)$ on coverage. Langmuir site exclusion already represents finite
capacity, so we added no second hard-disk packing free energy.

### S-2.4. Parameters, charge inversion and scope

We set the geometric site-capacity prior to
$\Gamma_{\max}=1/(11\ \text{Å})^2$. At each pH, we fitted one $K$ to the low-salt
$\zeta$-versus-[PA] series and held it fixed across salt concentrations. If allowing
adsorption improved $\chi^2$ by less than 1 relative to $K=0$, we reported $K=0$ as
not identifiable; this criterion applied at pH 10. The $\zeta$ sample was diluted
1000-fold. We approximated its re-equilibrated free
concentration by the diluted analytical concentration. The independent ITC fit left
$K$, area per oleosin, and a constant adsorption enthalpy free; it imposed neither the
$\zeta$-fitted $K$ nor the physical area. The ITC calculation included injection
dilution and ligand mass balance explicitly (Section S-6). We interpreted the NMR
relaxation separately and extracted no adsorption coverage from it (Section S-5).

For a positively charged bare surface, charge inversion occurs when
$z\Gamma>\sigma_0$. After inversion, the diffuse-layer potential opposes adsorption of
another anion, making the response self-limiting. The residual correlation term can
make the coverage more cooperative, especially at low salt. Any predicted kink or
branch arises from this correlation-hole closure and must be compared with the matched
$f_{\mathrm{corr,res}}=0$ baseline. The main omitted physics is the coverage
dependence of the peptide field and phytate charge, together with the finite thickness
and permeability of the oleosin brush. Only the diluted $\zeta$ calculation equates
free and analytical ligand concentrations; the ITC fit includes ligand mass balance.

## S-3. $\zeta$-potential

Figure S3 compares LC–AR-PB with the measured $\zeta$-potential as a function of
phytate (A) and salt (B). Panel A includes the low-salt data used to fit one intrinsic
$K$ per pH (Table S-3a), including the pH-4 charge reversal. Panel B transfers those
affinities across salt without refitting. Both calculations use the peptide-only MC
field, the MC phytate charge at $R_\beta$, and the PB-matched residual lateral term.

![Measured $\zeta$-potential (× with error bars) and LC–AR-PB with residual lateral correlations (lines), coloured by pH (blue 4, orange 6, green 10): (A) versus [PA] at low salt; (B) versus [NaCl] at 1.5 mM PA. The zero-PA low-salt value fixes the bare surface charge; the PA-dependent curves and salt transfer follow from the common surface grand potential.](../zeta/fig_lateral_arpb.png)

Table S-3a. Intrinsic affinities.

| pH | $K$ / M⁻¹ | binding |
|----|----------------|---------|
| 4  | 1.29           | overcharges (charge inversion) |
| 6  | 47.3           | binds, adds negative charge |
| 10 | 0 (not identifiable) | direct binding not resolved |

The salt series tests transfer because $K$ was fitted only to the low-salt phytate
series. At pH 6, the transferred model captures the measured weakening of the phytate
response with salt. The simulated PMF wells independently weaken with salt (Table
S-3b), supporting the same screening trend. Deviations at pH 4 and 10 may arise from
the planar, impermeable Gouy–Chapman representation of the finite brush or from holding
the one-particle field and charge fixed as coverage changes.

Table S-3b. pH-6 potential-of-mean-force well depth versus salt (simulation).

| $I$ / M | well depth / $k_BT$ |
|---------|---------------------|
| 0.005   | $-13.9$ |
| 0.03    | $-8.3$  |
| 0.15    | $-2.4$  |

## S-4. Simulation composites versus salt

The composites below extend Figure 10 to $I=0.03$ and 0.15 M. Across the full salt
series, the pH-6 PMF well decreases from $-13.9$ to $-2.4\,k_BT$, and the electric
potential decreases in magnitude. Both changes show progressive electrostatic
screening of adsorption.

![Simulation composite at $I = 0.03$ M (panels as in the low-salt composite).](../mc/composite_I0.03.png)

![Simulation composite at $I = 0.15$ M. The measured $\zeta$ line corresponds to the highest measured salt concentration, 100 mM added NaCl. Because $|\zeta|$ falls with salt, the line overestimates $|\zeta|$ at the simulated ionic strength, although the Debye length differs by only a factor of $1.2$ between the two conditions.](../mc/composite_I0.15.png)

## S-5. NMR interpretation and simulation scope

We compare the $^{31}$P magic-angle-spinning (MAS) spectra and relaxation data with
the simulations only qualitatively. The larger OB-induced line broadening and
rotating-frame relaxation changes at pH 4 and 6 show a larger spectral and relaxation
perturbation than at pH 10. This pH ordering matches the simulated PMFs, which are
attractive at pH 4 and 6 and repulsive at pH 10. The lack of a detectable
cross-polarization (CP) signal argues against a dominant rigid phytate population. A
minor rigid population below the detection limit remains possible.

We extracted no adsorbed fraction from $R_{1\rho}-R_1$. Because only one spin-lock
field was measured, intrinsic transverse relaxation, exchange, and motional
spectral-density contributions cannot be separated (Quinn & McDermott, 2012). The
equilibrium MC simulations provide neither an exchange rate nor a physical rotational
correlation time. We therefore report $S^2_\mathrm{ens}$ and $W_\mathrm{rot}$ only as equilibrium
orientational diagnostics and do not fit them to the NMR rates. A quantitative
connection requires relaxation-dispersion measurements and either a dynamical model
or simulations of the relevant orientational time-correlation functions.

## S-6. Calorimetry: LC–AR-PB versus the one-site instrument model

We fitted the pH-4 thermogram with two models. First, we refitted the one-site model
used in the instrument analysis to the extracted injection heats, including the cell
overflow correction. Both fits excluded the first injection, whose heat is unreliable
because the syringe tip exchanges with the cell during thermal equilibration and
because it was a smaller priming injection (5.1 against 8 µL). The phytate it
delivered was retained in the dilution and mass balance of all later injections. We then fitted LC–AR-PB as a surface isotherm with three free
parameters: intrinsic $K$, area $A_o$ per oleosin, and constant molar adsorption
enthalpy $\Delta H$. This diagnostic LC–AR-PB fit was independent of the
$\zeta$-potential fit and did not constrain $K$ to the $\zeta$-derived value.

For injection $i$, let $M_{t,i}$ and $L_{t,i}$ be the diluted total oleosin and
phytate concentrations. We evaluated the LC–AR-PB coverage at the free phytate
concentration $c_{\mathrm{free},i}$ and solved the ligand balance explicitly:

$$
L_{t,i}=c_{\mathrm{free},i}
+M_{t,i}A_o\Gamma(c_{\mathrm{free},i};K).
$$

The bound concentration was
$B_i=M_{t,i}A_o\Gamma(c_{\mathrm{free},i};K)$. For injection volume $\Delta V_i$ and
cell volume $V_0$, the displacement-corrected heat was

$$
q_i=\Delta H\,V_0\left[
B_i-B_{i-1}
+\frac{\Delta V_i}{V_0}\frac{B_i+B_{i-1}}{2}
\right].
$$

This heat model assigns one constant enthalpy to each net adsorption event. It does
not resolve coverage-dependent electrostatic enthalpy, correlation enthalpy, or heat
from aggregation.

![pH-4 ITC thermogram: extracted injection heats (filled points; the excluded first injection is open), the unconstrained LC–AR-PB surface fit with residual lateral correlations (solid), and the independent one-site instrument model refitted with the same dilution correction (dashed).](../itc/fig_itc_lateral_arpb.png)

Both models described the thermogram with similar residuals (Table S-6). The one-site
model gave
$K_a=(1.9\pm0.5)\times10^5$ M⁻¹, $n=0.336\pm0.03$ phytate per oleosin and
$\Delta H=-56\pm7$ kJ mol⁻¹, with RMS 5.2 µJ, reproducing the instrument fit
($K_a=1.87\times10^5$ M⁻¹, $n=0.322$, $\Delta H=-56.4$ kJ mol⁻¹). LC–AR-PB gave
$K=0.263\pm0.23$ M⁻¹, $A_o=(7.40\pm0.95)\times10^3$ Å² per oleosin and
$\Delta H=-41.2\pm3.1$ kJ mol⁻¹, with RMS 6.1 µJ. The two affinity constants are
not directly comparable. The one-site $K_a$ contains the complete apparent binding
free energy. By contrast, LC–AR-PB represents the peptide field and PB response
explicitly, leaving $K$ to describe the remaining intrinsic activity.

The fitted LC–AR-PB area was 21.9 times the physical grafting area of 338 Å² per
oleosin. Moreover, $K$ and area were strongly anticorrelated ($-0.90$). The inflated
area allowed the fit to reach 0.405 bound phytate per oleosin, close to the one-site
stoichiometry. Imposing the physical area reduced the endpoint to 0.020 phytate per
oleosin, or 0.026 with the pH-4 $\zeta$-fitted $K$. These values match the
$\sim0.03$ scale estimated from the $\zeta$-calibrated surface model. Thus, the LC–AR-PB
curve fit absorbs the missing capacity into an unphysical area and does not establish
a monolayer adsorption enthalpy.

The capacity mismatch and multiphasic calorimetric shoulder (Section 3.3.2) suggest
that inter-droplet cross-linking or aggregation contributes additional heat. A
non-aggregating calorimetric experiment could test the surface model with an
independently fixed area and explicit ligand mass balance.

Table S-6. Unconstrained fits to the pH-4 ITC thermogram. Uncertainties are formal
one-standard-error estimates from the local fit covariance and do not include model
uncertainty.

| quantity | one-site instrument model | LC–AR-PB surface model |
|----------|---------------------------|--------------------------|
| affinity | $K_a=(1.9\pm0.5)\times10^5$ M⁻¹ | $K=0.263\pm0.23$ M⁻¹ |
| capacity parameter | $n=0.336\pm0.03$ PA/oleosin | $A_o=(7.40\pm0.95)\times10^3$ Å²/oleosin |
| maximum fitted bound amount | $\approx0.34$ PA/oleosin | 0.405 PA/oleosin |
| $\Delta H$ | $-56\pm7$ kJ mol⁻¹ | $-41.2\pm3.1$ kJ mol⁻¹ |
| largest $|$parameter correlation$|$ | 0.93 ($n$–$\Delta H$) | 0.90 ($K$–$A_o$; anticorrelated) |
| RMS | 5.2 µJ | 6.1 µJ |
| physical interpretation | apparent multiphasic stoichiometry | area is 21.9× the physical 338 Å² |

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
- Quinn, C. M., & McDermott, A. E. (2012). Quantifying conformational dynamics using
  solid-state $R_{1\rho}$ experiments. *J. Magn. Reson.* **222**, 1–7.
  https://doi.org/10.1016/j.jmr.2012.05.014
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
