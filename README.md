# si-hempoil-phytate

Monte Carlo simulation and analysis of **phytate (IP6)** interacting with an
**oleosin N-terminal peptide brush** grafted to an oil-body surface, plus a
hydrodynamic estimate of free-phytate diffusion. This is the reproduction repository
for the simulation and theory results in the oil-body–phytate manuscript.

The peptide is the N-terminal arm of *Cannabis sativa* oleosin, modeled as a
coarse-grained brush on one wall of a slit. A titrating phytate molecule is pulled
along the slit normal by umbrella sampling to obtain the potential of mean force
(PMF), together with its charge regulation, shape, and the electric potential
profile. The [Faunus](https://github.com/mlund/faunus) MC engine runs the
simulations; the Python here sets up the runs and analyzes the output.

A **charge-regulated Poisson–Boltzmann** model (`zeta/charge_regulation.py`) carries
the simulated PMF to finite concentration and to the measured ζ-potential: it
reproduces the pH-4 ζ charge reversal (overcharging) at the cost of one intrinsic
affinity per pH, and it also fits the ITC thermogram. The ³¹P NMR bound fraction is
read directly from the simulated free energy (`nmr/binding_fraction.py`). Shared
MC-grounded utilities live in `zeta/mc_shared.py`.

## Geometry conventions

The Faunus `!Slit` cell is centered at `z = 0`:

- **decorated wall** (peptide brush) at `z = -Lz/2`
- **neutral wall** (bulk side, implicit salt) at `z = +Lz/2`

Analysis works in **`r = z + Lz/2`**, the distance from the decorated wall
(`r = 0` at the brush, `r → Lz/2` toward bulk). The phytate is sampled over
`z ∈ [-Lz/2, 0]`, so its reaction coordinate spans `r ∈ [0, Lz/2]`.

## Reproduce the paper results

Uses [uv](https://docs.astral.sh/uv/) (Python ≥ 3.12). The small CSV layer is
committed, so the figures and tables reproduce from a clean clone **without** the
large raw MC output:

```bash
uv sync                      # create .venv from uv.lock

uv run make_figures.py       # CSV layer -> all figures + tables (Table N, S-3a, S-6)
uv run writeup/build.py      # embed the figures into the Word/ODF write-up
```

`make_figures.py` regenerates every paper figure (`zeta/fig_crpb.png`,
`itc/fig_itc_crpb.png`, `mc/composite_I*.png`) and prints the tables to stdout.

To rebuild the CSV layer itself from the raw Faunus output — needs `mc/runs/`, which
is large and git-ignored (obtain it separately) — run Stage 1 first:

```bash
uv run extract_csvs.py       # mc/runs/ -> mc/{pmf,window,potential}_profiles.csv, simulated_zeta.csv
```

## Layout

```
mc/                       Faunus MC simulation + extraction
  input.yaml, state.yaml, phytate.xyz   simulation inputs
  generate_runs.py        write per-(pH, salt) run dirs + SLURM scripts
  extract_pmf.py          runs/ PMF -> pmf_profiles.csv (+ pmf.png), screened-Coulomb tail fit
  extract_profiles.py     runs/ -> window_profiles, potential_profiles, simulated_zeta
  plot_composites.py      composite figures (PMF, mobility, potential, charge) per salt
  pmf_profiles.csv        PMF(r) per (pH, ionic)                             [committed]
  window_profiles.csv     phytate/peptide charge, rotational mobility, shape vs r  [committed]
  potential_profiles.csv  electric potential phi(r)                          [committed]
  simulated_zeta.csv      shear-plane zeta per condition                     [committed]
  runs/                   raw Faunus output (large, git-ignored)
zeta/                     charge-regulated model + shared utilities
  mc_shared.py            shared MC-grounded utilities (PMF, valence, Grahame, Debye)
  charge_regulation.py    CR-PB model: zeta vs [PA]/[NaCl] (fig_crpb.png) + ITC fit
  parse_zeta_workbook.py  tidy the measured-zeta workbook -> measured_zeta.csv
  Zeta_potential_tidied.xlsx   measured-zeta workbook (source)
  measured_zeta.csv       measured zeta (mean, SD, n)                        [committed]
nmr/                      NMR bound fraction + free-phytate hydrodynamics
  binding_fraction.py     simulated PMF -> NMR bound fraction (Table N)
  rotational_diffusion.py RDKit conformers + PyGRPY -> free_phytate_diffusion.csv
  free_phytate_diffusion.csv   free-phytate hydrodynamics                    [committed]
itc/                      isothermal titration calorimetry
  extract_thermograms.py  load NanoAnalyze .xls/.nitc exports
  fit_thermogram.py       one-site (Wiseman) fit + shared titration engine
  pH4/ pH6/ pH10/         raw ITC exports
writeup/                  simulation-support text for the manuscript
  sim_for_mauricio.md     Word/ODF draft (main text + SI)
  build.py                sim_for_mauricio.md -> .docx + .odt (pandoc)
extract_csvs.py           Stage 1: mc/runs/ -> the committed CSV layer
make_figures.py           Stage 2: CSV layer -> all figures + tables
pyproject.toml            uv project (dependencies + ruff config)
```

## Simulate (optional)

Running the simulations needs the `faunus` binary on `PATH` (installed separately
from the analysis environment). Generate and (optionally) submit a scan over pH and
ionic strength:

```bash
cd mc
uv run generate_runs.py            # write runs/pH{pH}_I{salt}/ with SLURM scripts
uv run generate_runs.py --submit   # ... and sbatch each job
```

Each job equilibrates, then runs umbrella sampling, producing `pmf.csv` and the
per-window analysis under `runs/pH{pH}_I{salt}/`. Then run `extract_csvs.py` to
rebuild the CSV layer.

## Development

```bash
uv run ruff check .    # lint
uv run ruff format .   # format
```
