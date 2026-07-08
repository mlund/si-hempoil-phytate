# Supporting Information for Hempoil-Phytate Interactions: Theory and Simulations

Reproducibility repository for the simulation, charge-regulation, NMR,
zeta-potential, and ITC analyses used in the oil-body/phytate manuscript.

The large raw Faunus output is expected under `mc/runs/` and is not needed for the
standard reproduction path. The extracted CSV layer is tracked, so the published
figures and reported tables can be rebuilt from a clean clone.

## Reproduce Published Data and Figures

Requirements: [uv](https://docs.astral.sh/uv/) and Python 3.12 or newer.
`writeup/build.py` also requires `pandoc` on `PATH`.

```bash
uv sync
uv run make_figures.py
uv run writeup/build.py
```

`make_figures.py` is the main reproduction target. It rebuilds:

- `zeta/measured_zeta.csv` from `zeta/Zeta_potential_tidied.xlsx`
- `nmr/free_phytate_diffusion.csv`
- `mc/composite_I0.005.png`, `mc/composite_I0.03.png`, and
  `mc/composite_I0.15.png`
- `zeta/fig_crpb.png`
- `itc/fig_itc_crpb.png`

It also prints the reported tables to stdout: NMR binding fractions, zeta model
parameters, and ITC fit parameters.

To rebuild the tracked CSV layer from raw MC output, obtain `mc/runs/` separately
and run this first:

```bash
uv run extract_csvs.py
uv run make_figures.py
```

`extract_csvs.py` writes `mc/pmf_profiles.csv`, `mc/window_profiles.csv`,
`mc/potential_profiles.csv`, and `mc/simulated_zeta.csv`.

## Optional Simulation Rerun

Running the simulations requires the `faunus` binary on `PATH`; submission uses
SLURM.

```bash
cd mc
uv run generate_runs.py
uv run generate_runs.py --submit
cd ..
uv run extract_csvs.py
```

The MC analysis coordinate is `r = z + Lz/2`, the distance from the decorated
peptide-brush wall.

## Repository Layout

- `mc/`: Faunus inputs, run generation, PMF/profile extraction, tracked MC CSVs,
  and composite MC figures. `mc/runs/` holds large raw outputs and is git-ignored.
- `zeta/`: measured zeta-potential data, workbook parser, shared MC utilities, and
  the charge-regulated Poisson-Boltzmann model.
- `nmr/`: PMF-to-bound-fraction calculation and free-phytate hydrodynamic
  diffusion estimate.
- `itc/`: raw NanoAnalyze ITC exports, thermogram extraction, and one-site
  thermogram fitting.
- `writeup/`: manuscript/SI support text and the `pandoc` build script for Word
  and ODF output.
  `extract_csvs.py` for rebuilding the MC CSV layer, and `pyproject.toml` /
  `uv.lock` for the Python environment.

## Development

```bash
uv run ruff check .
uv run ruff format .
```
