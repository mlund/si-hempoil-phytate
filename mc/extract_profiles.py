#!/usr/bin/env python3
"""Extract per-window mean charge and polymer shape as functions of r.

Each umbrella window's <run>/umbrella_states/window*_output.yaml holds an
`analysis` block with the Multipole (mean charge, once per selection) and
PolymerShape analyses. The window's position is the mean phytate (INO) z from its
AtomPosition collective variable, converted to r = z + Lz/2 (distance from the
decorated wall at z = -Lz/2; see the extract_pmf.py geometry).

The ElectricPotentialProfile analyses (potential_all.csv, potential_peptide.csv)
are written per window over the full slit; they are averaged over all windows.

The MC RotateMolecule moves report a per-attempt mean-square angular displacement;
pooled over the phytate protonation states this is a rough rotational-mobility
proxy, giving D(r)/D0 ~ msd_rot(r) / msd_rot(bulk).

The potential analysis also samples the charge density sigma(z). We estimate a
zeta potential as the total potential at the hydrodynamic shear plane, taken at
the outer edge of the fixed brush charge (largest r where |sigma_peptide| exceeds
a small fraction of its peak). Magnitude is shear-plane dependent; the sign and
pH-trend are the meaningful output.

Writes three tidy CSVs for pandas.read_csv. The three per-window observables (phytate
and peptide charge, rotational MSD, and the gyration-tensor shape) all live on the same
window-centre r-grid, so they are merged into one wide table:
  window_profiles.csv    pH, ionic_strength, r, q_ino, q_ino_err, q_peptide,
                         q_peptide_err, msd_rot, <shape descriptors...>
  potential_profiles.csv pH, ionic_strength, r, potential_all_mV, potential_peptide_mV,
                         sigma_peptide_e_A2   (full-slit profile, its own r-grid)
  simulated_zeta.csv     pH, ionic_strength, r_shear, zeta_mV   (one row per condition)
"""

import argparse
import csv
import re
from pathlib import Path

import numpy as np
import yaml

RUNS_DIR = Path(__file__).parent / "runs"

# Faunus output.yaml carries custom !Tags (!Slit, !Every, ...); map every tag to
# its plain scalar/sequence/mapping so a SafeLoader can read the file.
_Base = getattr(yaml, "CSafeLoader", yaml.SafeLoader)


class _Loader(_Base):
    pass


def _ignore_tag(loader, suffix, node):
    if isinstance(node, yaml.ScalarNode):
        return loader.construct_scalar(node)
    if isinstance(node, yaml.SequenceNode):
        return loader.construct_sequence(node)
    return loader.construct_mapping(node)


_Loader.add_multi_constructor("!", _ignore_tag)


def load_yaml(path: Path):
    return yaml.load(path.read_text(), Loader=_Loader)


def ph_from_name(name: str) -> float:
    return float(re.search(r"pH([0-9.]+)", name).group(1))


def ionic_strength_from_name(name: str) -> float:
    return float(re.search(r"I([0-9.]+)", name).group(1))


def parse_pm(value: str) -> tuple[float, float]:
    """Split a Faunus 'mean ± error' string into floats."""
    mean, error = value.split("±")
    return float(mean), float(error)


def window_index(path: Path) -> int:
    return int(re.search(r"window(\d+)_output", path.name).group(1))


def rotation_msd(doc) -> float:
    """Pooled mean-square angular displacement per attempted phytate rotation.

    The tag-stripping loader discards the move type, so a phytate RotateMolecule is
    identified as a non-peptide molecule move without a `directions` key (its
    TranslateMolecule has `directions: z`). Summing accepted <Δθ²> over every
    protonation state's rotation move and dividing by the (position-independent)
    total attempts gives a consistent MC rotational-mobility proxy across windows.
    """
    total_sq = 0.0
    total_n = 0
    for collection in doc["propagate"]["collections"]:
        for move in collection.get("moves", []):
            mol = move.get("molecule")
            if mol is None or mol == "peptide" or "directions" in move:
                continue
            msd = move["statistics"]["mean_square_displacement"]
            total_sq += msd["avg"] * msd["n"]
            total_n += msd["n"]
    return total_sq / total_n if total_n else float("nan")


def parse_window(output_yaml: Path):
    """Return (z, [(selection, charge, err), ...], shape_dict, rot_msd) for one window."""
    doc = load_yaml(output_yaml)
    z = None
    charges = []
    shape = None
    for item in doc["analysis"]:
        ((key, val),) = item.items()  # each analysis entry is a single-key mapping
        if key == "collectivevariable" and val.get("property") == "AtomPosition":
            z = val["mean"]  # mean phytate (INO) z for this window
        elif key == "multipole":
            mean, error = parse_pm(val["charge"])
            charges.append((val["selection"], mean, error))
        elif key == "polymershape":
            shape = {k: v for k, v in val.items() if k != "num_samples"}
    return z, charges, shape, rotation_msd(doc)


def mean_profile(profile_csvs, cols) -> tuple:
    """Average the requested data columns over a set of same-grid profiles.

    Columns of potential_*.csv: z/Å (0), slab_charge_density/e·Å⁻² (1),
    charge_density/e·Å⁻³ (2), potential/mV (3), error (4), field (5). The grid is
    identical across windows. Returns (z, mean_col, ...) for each column in `cols`.
    """
    profiles = [
        np.genfromtxt(p, delimiter=",", skip_header=1, usecols=(0, *cols)) for p in profile_csvs
    ]
    avg = np.mean(profiles, axis=0)
    return (avg[:, 0], *(avg[:, i + 1] for i in range(len(cols))))


def shear_zeta(r, sigma_peptide, potential_all, half_lz, threshold=0.02):
    """Zeta = total potential at the fixed-charge edge of the brush.

    The shear plane is the outermost point (brush side, r < Lz/2) where the peptide
    areal charge density exceeds `threshold` of its peak; beyond it the fixed charge
    has vanished and the potential is purely diffuse. Returns (r_shear, zeta_mV).
    """
    sig = np.abs(sigma_peptide)
    brush = (r < half_lz) & (sig > threshold * sig.max())
    r_shear = r[brush].max()
    return r_shear, float(np.interp(r_shear, r, potential_all))


# Faunus charge-selection labels -> the window_profiles column stems they map to.
_CHARGE_COLS = {"atomtype INO": "q_ino", "molecule peptide": "q_peptide"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--window", default="window_profiles.csv", help="output window-profile CSV")
    parser.add_argument(
        "--potential", default="potential_profiles.csv", help="output potential CSV"
    )
    parser.add_argument("--zeta", default="simulated_zeta.csv", help="output shear-plane zeta CSV")
    parser.add_argument(
        "--max-ionic-strength",
        type=float,
        default=1.0,
        help="skip runs above this I in mol/L (default: 1.0, drops the I1000 sentinels)",
    )
    args = parser.parse_args()

    window_rows = []
    shape_fields = None
    potential_rows = []
    zeta_rows = []
    for output_yaml in sorted(RUNS_DIR.glob("*/output.yaml")):
        run = output_yaml.parent
        name = run.name
        ph = ph_from_name(name)
        ionic = ionic_strength_from_name(name)
        if ionic > args.max_ionic_strength:  # skip the I1000 sentinel (matches extract_pmf)
            continue
        half_lz = load_yaml(output_yaml)["cell"][2] / 2  # Lz is the slit z-length

        # Potential profiles span the full slit; average each selection over windows.
        all_csvs = sorted(run.glob("umbrella_states/window*/potential_all.csv"))
        pep_csvs = sorted(run.glob("umbrella_states/window*/potential_peptide.csv"))
        z, pot_all = mean_profile(all_csvs, (3,))
        _, pot_peptide, sigma_peptide = mean_profile(pep_csvs, (3, 1))  # potential + charge density
        r_prof = z + half_lz
        for ri, pa, pp, sp in zip(r_prof, pot_all, pot_peptide, sigma_peptide, strict=True):
            potential_rows.append([ph, ionic, ri, pa, pp, sp])

        # Zeta at the fixed-charge edge of the brush (sigma_peptide -> 0).
        r_shear, zeta = shear_zeta(r_prof, sigma_peptide, pot_all, half_lz)
        zeta_rows.append([ph, ionic, r_shear, zeta])

        # Per-window observables (charge, shape, rotational MSD) share the window-centre
        # r-grid, so each window yields one merged row.
        windows = sorted(run.glob("umbrella_states/window*_output.yaml"), key=window_index)
        for win in windows:
            z, charges, shape, msd = parse_window(win)
            if z is None:
                continue
            r = z + half_lz
            q = dict.fromkeys(_CHARGE_COLS.values(), "")
            for selection, charge, error in charges:
                stem = _CHARGE_COLS.get(selection)
                if stem is not None:
                    q[stem], q[f"{stem}_err"] = charge, error
            clean = {}
            if shape is not None:
                # Sanitise field names (e.g. Re2/Rg2 -> Re2_Rg2) for tidy columns.
                clean = {k.replace("/", "_"): v for k, v in shape.items()}
                if shape_fields is None:
                    shape_fields = list(clean)
            # Keep the shape dict; expand it against the final shape_fields at write time so
            # rows stay aligned even if a shapeless window is seen before the first shaped one.
            window_rows.append(
                [
                    ph,
                    ionic,
                    r,
                    q["q_ino"],
                    q.get("q_ino_err", ""),
                    q["q_peptide"],
                    q.get("q_peptide_err", ""),
                    msd,
                    clean,
                ]
            )

    window_rows.sort(key=lambda row: (row[0], row[1], row[2]))
    potential_rows.sort(key=lambda row: (row[0], row[1], row[2]))
    zeta_rows.sort(key=lambda row: (row[0], row[1]))

    shape_fields = shape_fields or []
    with open(args.window, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "pH",
                "ionic_strength",
                "r",
                "q_ino",
                "q_ino_err",
                "q_peptide",
                "q_peptide_err",
                "msd_rot",
                *shape_fields,
            ]
        )
        for *base, clean in window_rows:
            w.writerow([*base, *(clean.get(fld, "") for fld in shape_fields)])
    print(f"Saved {args.window} ({len(window_rows)} rows)")

    with open(args.potential, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "pH",
                "ionic_strength",
                "r",
                "potential_all_mV",
                "potential_peptide_mV",
                "sigma_peptide_e_A2",
            ]
        )
        w.writerows(potential_rows)
    print(f"Saved {args.potential} ({len(potential_rows)} rows)")

    with open(args.zeta, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pH", "ionic_strength", "r_shear", "zeta_mV"])
        w.writerows(zeta_rows)
    print(f"Saved {args.zeta} ({len(zeta_rows)} rows)")


if __name__ == "__main__":
    main()
