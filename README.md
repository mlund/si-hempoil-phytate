[![Code: MIT](https://img.shields.io/badge/code-MIT-blue.svg)](LICENSE)
[![Data: CC BY 4.0](https://img.shields.io/badge/data-CC%20BY%204.0-lightgrey.svg)](LICENSE-DATA)
<!-- TODO: Zenodo concept-DOI badge (all versions) after the first release -->

# Supporting Information for Hempoil-Phytate Interactions: Theory and Simulations

Reproducibility repository for the simulation, adsorption-regulated PB,
zeta-potential, and ITC analyses in *Phytate-hempseed oil bodies interactions:
effect of pH on colloidal, thermodynamic, and spectroscopic properties*
(Oyarzún, Lund, Bonilla, Fuentealba, Topgaard, Gouseti, Jensen, Risbo; manuscript
in preparation).

The large raw Faunus output is expected under `mc/runs/` and is not needed for the
standard reproduction path. The extracted CSV layer is tracked, so the published
figures and reported tables can be rebuilt from a clean clone.

## Reproduce Published Data and Figures

Requirements: [uv](https://docs.astral.sh/uv/) (installs Python 3.12).

```bash
uv sync
uv run make_figures.py
```

`make_figures.py` is the main reproduction target. It rebuilds:

- `zeta/measured_zeta.csv` from `zeta/Zeta_potential_tidied.xlsx`
- `mc/composite_I0.005.png`, `mc/composite_I0.03.png`, and
  `mc/composite_I0.15.png`
- `zeta/fig_lc_arpb_schematic.png`
- `zeta/fig_lateral_arpb.png` and `zeta/model_zeta.csv`
- `itc/itc_injections.csv` and `itc/fig_itc_lateral_arpb.png`

It also prints the zeta-model and ITC fit parameters. `zeta/model_zeta.csv` holds the
predicted LC--AR-PB zeta for both scans, in long format: `sampling == "curve"` is the
dense grid drawn in the figure, `sampling == "measured"` the rows matching the measured
conditions (Table S-3a).

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
  the LC–AR-PB model.
- `itc/`: raw NanoAnalyze ITC exports, thermogram extraction, and one-site
  thermogram fitting.
- `extract_csvs.py` and `make_figures.py`: top-level reproduction drivers.
- `pyproject.toml` and `uv.lock`: Python environment.

## Data

Text files are UTF-8 CSV with a header row; units are in the column names
(`_kT`, `_mV`, `_mM`, `_e_A2` = e/Å², `_A3` = 1/Å³). `r` is in Å, `ionic_strength`
in mol/L.

| File | Content | Origin |
|---|---|---|
| `mc/pmf_profiles.csv` | Phytate PMF `w(r)` and standard error per pH and ionic strength | `extract_csvs.py` from `mc/runs/` |
| `mc/window_profiles.csv` | Per-window phytate and peptide charge, orientational order, peptide shape | `extract_csvs.py` |
| `mc/potential_profiles.csv` | Mean electric potential and peptide charge density vs. `r` | `extract_csvs.py` |
| `mc/simulated_zeta.csv` | Potential at the shear plane `r_shear` | `extract_csvs.py` |
| `zeta/Zeta_potential_tidied.xlsx` | Measured zeta potentials | experiment |
| `zeta/measured_zeta.csv` | Same, tidy: mean and SD over `n` replicates | `make_figures.py` |
| `itc/pH*/*.nitc` | Raw TA Instruments ITC runs (NanoAnalyze format) | experiment |
| `itc/pH*/*.xls`, `Summary_*.xlsx` | NanoAnalyze exports: integrated heats (µJ) per injection | experiment |
| `itc/itc_injections.csv` | All ITC replicates, one row per injection: volume, heat, concentrations; `analysed`/`fitted` flag what enters the reported fits | `make_figures.py` |

Raw MC output (`mc/runs/`, ~130 MB zipped) is not deposited: `mc/input.yaml`,
`mc/state.yaml`, and `mc/generate_runs.py` regenerate it, and the CSVs above
reproduce every figure.

## Citation and license

Cite via the Zenodo DOI of this repository and the associated article;
`CITATION.cff` holds the metadata. Code is under the MIT license
([LICENSE](LICENSE)); data, figures, and text are under CC BY 4.0
([LICENSE-DATA](LICENSE-DATA)).

## Development

```bash
uv run pytest        # physics invariants and regression tests
uv run ruff check .
uv run ruff format .
uv run clean.py      # remove figures and transient caches
```

The clean command preserves the tracked CSV analysis layer, raw simulation output,
experimental workbooks, and archives. Use `uv run clean.py --dry-run` to inspect its
cleanup set first.

## Release checklist (delete before publishing)

- [ ] Add funding (`grants` in `.zenodo.json`), if any.
- [ ] Confirm the CC BY 4.0 data license with the co-authors (experimental data).
- [ ] After acceptance: article DOI in `.zenodo.json` (`related_identifiers`,
      `isSupplementTo`) and `CITATION.cff` (`preferred-citation`).
- [ ] Make the repository public and enable it on https://zenodo.org (GitHub integration).
- [ ] Create a GitHub release (SemVer) to trigger the deposit; add the concept-DOI badge.
