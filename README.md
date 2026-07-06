# si-hempoil-phytate

Monte Carlo simulation and analysis of **phytate (IP6)** interacting with an
**oleosin N-terminal peptide brush** grafted to an oil-body surface, plus a
hydrodynamic estimate of free-phytate diffusion.

The peptide is the N-terminal arm of *Cannabis sativa* oleosin, modeled as a
coarse-grained brush on one wall of a slit. A titrating phytate molecule is pulled
along the slit normal by umbrella sampling to obtain the potential of mean force
(PMF), together with its charge regulation, shape, and the electric potential
profile. The [Faunus](https://github.com/mlund/faunus) MC engine runs the
simulations; the Python here sets up scans and analyzes the output.

## Geometry conventions

The Faunus `!Slit` cell is centered at `z = 0`:

- **decorated wall** (peptide brush) at `z = -Lz/2`
- **neutral wall** (bulk side, implicit salt) at `z = +Lz/2`

Analysis works in **`r = z + Lz/2`**, the distance from the decorated wall
(`r = 0` at the brush, `r → Lz/2` toward bulk). The phytate is sampled over
`z ∈ [-Lz/2, 0]`, so its reaction coordinate spans `r ∈ [0, Lz/2]`.

## Layout

```
mc/                     Faunus MC simulation + analysis
  input.yaml            umbrella-sampling input (MiniJinja template: pH, salt, ...)
  setup_scan.py         generate per-(pH, I) run dirs + SLURM scripts
  plot_pmf.py           PMF plots, faceted by pH, with a screened-Coulomb tail fit
  extract_data.py       per-window charge & polymer shape + averaged potential
  peptide_density_z.py  peptide number density along z from an XTC trajectory
  analysis.ipynb        merge the extracted CSVs and plot per condition
  runs/                 simulation output (large files git-ignored)
nmr/                    NMR binding bridge + free-phytate hydrodynamics
  nmr_binding_bridge.py
  rotational_diffusion.py (RDKit conformers + PyGRPY)
latex/                  manuscript
pyproject.toml          uv project (dependencies + ruff config)
```

## Setup

Uses [uv](https://docs.astral.sh/uv/). Requires Python ≥ 3.12.

```bash
uv sync            # create .venv from uv.lock
```

Running the simulations additionally needs the `faunus` binary on `PATH`
(installed separately from the analysis environment).

## Workflow

### 1. Simulate

Generate and (optionally) submit a scan over pH and ionic strength:

```bash
cd mc
uv run setup_scan.py            # write runs/pH{pH}_I{salt}/ with SLURM scripts
uv run setup_scan.py --submit   # ... and sbatch each job
```

Each job equilibrates, then runs umbrella sampling, producing `pmf.csv` and the
per-window analysis under `runs/pH{pH}_I{salt}/`.

### 2. Analyze

```bash
cd mc
uv run plot_pmf.py         # -> all_pmfs.csv, pmf.png (PMF fit, referenced to bulk)
uv run extract_data.py     # -> charge.csv, polymer_shape.csv, potential.csv
uv run jupyter lab analysis.ipynb   # merge CSVs and plot per (pH, I) condition
```

`plot_pmf.py` fits the bulk tail of each PMF to a screened-Coulomb (Yukawa) form
`PMF(r) = a/r · exp(-r/D) + b` with the Debye length `D` fixed by the ionic
strength, and subtracts `b` to reference the PMF to zero in the bulk.

### 3. Free-phytate diffusion

Rotational and translational diffusion of free phytate across protonation states,
from RDKit conformers and PyGRPY bead hydrodynamics:

```bash
uv run nmr/rotational_diffusion.py      # -> nmr/diffusion.csv
```

Gives τ_r ≈ 140–165 ps (near-isotropic tumbler). This is a bare-hydrodynamic
estimate (stick boundary conditions, no ionic-atmosphere friction); see the
notebook and script docstring for caveats.

## Development

```bash
uv run ruff check .    # lint .py and .ipynb
uv run ruff format .   # format
```
