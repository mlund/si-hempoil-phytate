#!/usr/bin/env python3
"""Plot phytate PMF(s) from umbrella sampling and fit the screened-Coulomb tail.

The Faunus slit is centred at z = 0, with the peptide-decorated wall at z = -Lz/2
and the bulk at the box centre. We work in r = z + Lz/2, the distance from the
decorated wall (r = 0 at the wall, r ≈ Lz/2 in the bulk). The brush is a planar
double layer, so the bulk tail is fit to the planar (Gouy-Chapman) screened form

    PMF(r) = a * exp(-r / D) + b       (not the spherical exp(-r/D)/r)

with the Debye length D fixed by the ionic strength I read from the run name
(e.g. I0.005). Subtracting the fitted b references the PMF to zero in the bulk,
consistent with PMF -> 0 far from the wall (and everywhere at strong screening).
"""

import argparse
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

RUNS_DIR = Path(__file__).parent / "runs"

# Water at the simulation temperature (298.15 K); reproduces Faunus's reported
# debye_length/Å of 42.98 for I = 0.005 M.
BJERRUM_A = 7.1528  # Bjerrum length, Å
N_A = 6.02214076e23  # Avogadro's number, 1/mol


def debye_length(ionic_strength_molar: float) -> float:
    """Debye length in Å from ionic strength in mol/L (1:1 electrolyte, water)."""
    number_density = ionic_strength_molar * N_A / 1e27  # molar -> 1/Å^3
    kappa = math.sqrt(8 * math.pi * BJERRUM_A * number_density)
    return 1.0 / kappa


def ionic_strength_from_name(name: str) -> float:
    """Parse the ionic strength (mol/L) encoded as I<value> in a run name."""
    match = re.search(r"I([0-9.]+)", name)
    if not match:
        raise ValueError(f"no ionic strength (I<value>) found in run name {name!r}")
    return float(match.group(1))


def ph_from_name(name: str) -> float:
    """Parse the pH encoded as pH<value> in a run name."""
    match = re.search(r"pH([0-9.]+)", name)
    if not match:
        raise ValueError(f"no pH (pH<value>) found in run name {name!r}")
    return float(match.group(1))


def box_length_z(output_yaml: Path) -> float:
    """Read Lz, the slit's z-length, from the `cell: !Slit [Lx, Ly, Lz]` block."""
    # Explicit UTF-8: read_text() would otherwise use the locale default (cp1252 on Windows).
    lines = output_yaml.read_text(encoding="utf-8").splitlines()
    start = next(i for i, line in enumerate(lines) if line.strip().startswith("cell: !Slit"))
    # The three dimensions follow as a YAML list ("- <value>"); Lz is the third.
    dims = [
        float(m.group(1))
        for line in lines[start + 1 : start + 4]
        if (m := re.match(r"\s*-\s*([0-9.eE+-]+)", line))
    ]
    return dims[2]


def fit_tail(r, pmf, err, debye, window):
    """Fit PMF(r) = a*exp(-r/D) + b over the `window` Å of bulk (largest r).

    The brush is a planar double layer, so the far-field PMF decays as the planar
    (Gouy-Chapman) screened form exp(-r/D) -- NOT the spherical/point-charge Yukawa
    exp(-r/D)/r. (The two agree at high salt, where the tail is flat within the
    window, but differ by ~0.3-0.8 kT in b at low salt where D is large.)
    With D fixed the model is linear in (a, b): a weighted least-squares fit.
    """
    # The bulk is the box centre, at the largest sampled r; fit the screened tail
    # in the `window` Å approaching it.
    mask = r >= r.max() - window
    rw, pmfw, errw = r[mask], pmf[mask], err[mask]

    basis = np.column_stack([np.exp(-rw / debye), np.ones_like(rw)])
    weights = 1.0 / errw
    coeffs, *_ = np.linalg.lstsq(basis * weights[:, None], pmfw * weights, rcond=None)
    a, b = coeffs
    return a, b, rw


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    # Each run is a subdirectory of runs/ named by conditions, e.g. pH4_I0.005.
    parser.add_argument(
        "runs",
        nargs="*",
        help="run names to plot (default: every run under runs/ with a pmf.csv)",
    )
    parser.add_argument("-o", "--output", default="pmf.png", help="output image path")
    parser.add_argument(
        "--csv",
        default="all_pmfs.csv",
        help="tidy long-format CSV of the referenced PMFs (default: all_pmfs.csv)",
    )
    parser.add_argument(
        "-w",
        "--window",
        type=float,
        default=10.0,
        help="width in Å of the tail fit window (default: 10)",
    )
    parser.add_argument(
        "--max-ionic-strength",
        type=float,
        default=1.0,
        help="skip runs above this I in mol/L (default: 1.0, drops the I1000 sentinels)",
    )
    args = parser.parse_args()

    if args.runs:
        pmf_files = [RUNS_DIR / name / "pmf.csv" for name in args.runs]
    else:
        pmf_files = sorted(RUNS_DIR.glob("*/pmf.csv"))

    runs = [(p.parent.name, p) for p in pmf_files]
    runs = [(n, p) for n, p in runs if ionic_strength_from_name(n) <= args.max_ionic_strength]
    if not runs:
        raise SystemExit(f"no runs (I <= {args.max_ionic_strength}) found under {RUNS_DIR}")

    # One panel per pH; within a panel, colour curves by ionic strength.
    ph_values = sorted({ph_from_name(n) for n, _ in runs})
    ionic_values = sorted({ionic_strength_from_name(n) for n, _ in runs})
    shades = np.linspace(0, 0.85, len(ionic_values))
    colors = {i: plt.cm.viridis(x) for i, x in zip(ionic_values, shades, strict=True)}

    ncols = 2
    nrows = -(-len(ph_values) // ncols)  # ceil division
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(5 * ncols, 3.2 * nrows), sharex=True, sharey=True, squeeze=False
    )
    panels = {ph: axes.flat[k] for k, ph in enumerate(ph_values)}

    records = []  # accumulate one long-format block (pH, I, r, pmf, stderr) per run
    for name, path in runs:
        ph = ph_from_name(name)
        ionic = ionic_strength_from_name(name)
        ax = panels[ph]
        data = np.genfromtxt(path, delimiter=",", names=True)
        z, pmf, err = data["cv"], data["pmf_kT"], data["stderr_kT"]
        # Decorated wall at z = -Lz/2; measure r from it (r = 0 at wall, bulk at r = Lz/2).
        lz = box_length_z(path.parent / "output.yaml")
        r = z + lz / 2

        debye = debye_length(ionic)
        a, b, _ = fit_tail(r, pmf, err, debye, args.window)
        print(f"{name}: D = {debye:.2f} Å, a = {a:.3g} kT·Å, b = {b:.3g} kT")

        # The PMF's zero is arbitrary. Subtract only the constant offset b, leaving
        # PMF(r) - b that approaches the screened-Coulomb tail a/r*exp(-r/D) in the
        # bulk (which itself decays to zero) rather than being pinned to zero.
        pmf = pmf - b

        color = colors[ionic]
        ax.fill_between(r, pmf - err, pmf + err, alpha=0.2, color=color)
        ax.plot(r, pmf, color=color, lw=1.5, label=f"I = {ionic:g} M")

        records.append(np.column_stack([np.full_like(r, ph), np.full_like(r, ionic), r, pmf, err]))

    for ph, ax in panels.items():
        ax.set_title(f"pH {ph:g}")
        ax.axhline(0, color="gray", ls="--", lw=0.5)
        ax.legend(fontsize="small")
    for ax in axes.flat[len(ph_values) :]:  # hide any unused panels
        ax.set_visible(False)

    fig.supxlabel("r = distance from decorated wall (Å)")
    fig.supylabel("PMF (kT)")
    fig.suptitle("Phytate PMF near oleosin brush")
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)
    print(f"Saved {args.output}")

    # Tidy long format: read with pandas.read_csv(...) then groupby([pH, ionic_strength]).
    table = np.vstack(records)
    np.savetxt(
        args.csv,
        table,
        delimiter=",",
        comments="",
        header="pH,ionic_strength,r,pmf_kT,stderr_kT",
        fmt=["%g", "%g", "%.4f", "%.6f", "%.6f"],
    )
    print(f"Saved {args.csv} ({len(table)} rows)")

    plt.show()


if __name__ == "__main__":
    main()
