#!/usr/bin/env python3
"""Extract per-window mean charge and polymer shape as functions of r.

Each umbrella window's <run>/umbrella_states/window*/output.yaml holds an
`analysis` block with the Multipole (mean charge, once per selection) and
PolymerShape analyses. The window's position is the mean phytate (INO) z from its
AtomPosition collective variable, converted to r = z + Lz/2 (distance from the
decorated wall at z = -Lz/2; see the extract_pmf.py geometry).

The ElectricPotentialProfile analyses (potential_all.csv, potential_peptide.csv)
are written per window over the full slit; they are averaged over all windows.

The WidomRotation analysis adds step-size-independent equilibrium quantities: the
ensemble second-rank orientational order s2 (0 = isotropic, 1 = perfectly aligned)
and the conditional orientational-cage free energy widom_w [kJ/mol], referenced to
the best trial orientation. These are static diagnostics, not correlation times or
NMR relaxation rates. Older runs lacking the analysis leave these columns blank.

The PolymerShape gyration tensor additionally yields s_wall = (3 Szz/Rg^2 - 1)/2,
the P2 order parameter of the disk-like phytate's mass distribution about the wall
normal z (0 isotropic, -1/2 face-on, +1 edge-on) -- a lab-frame *pose* diagnostic,
distinct from the P--O-proxy ensemble order s2 (see wall_alignment).

The potential analysis also samples the charge density sigma(z). We estimate a
zeta potential as the total potential at the hydrodynamic shear plane, taken at
the outer edge of the fixed brush charge (largest r where |sigma_peptide| exceeds
a small fraction of its peak). Magnitude is shear-plane dependent; the sign and
pH-trend are the meaningful output.

The DensityProfile analysis gives the brush monomer density rho_peptide(r), retained as
a structural profile. It comes from the monomer density rather than sigma_peptide,
which changes sign at the anchor and depends on pH. Older runs lacking the analysis
leave the column blank.

Writes three tidy CSVs for pandas.read_csv. The per-window charge, Widom, and
gyration-tensor observables all live on the same
window-centre r-grid, so they are merged into one wide table:
  window_profiles.csv    pH, ionic_strength, r, q_ino, q_ino_err, q_peptide,
                         q_peptide_err, s2, s2_err, widom_w, s_wall,
                         <shape descriptors...>
  potential_profiles.csv pH, ionic_strength, r, potential_all_mV, potential_peptide_mV,
                         sigma_peptide_e_A2, rho_peptide_A3  (full-slit profile, own r-grid)
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


def parse_pm(value) -> tuple[float, float]:
    """Return (mean, error) from a Faunus multipole moment.

    Current Faunus (PR #fd0422e) emits a {mean, error} mapping; older runs wrote a
    'mean ± error' string. Accept either.
    """
    if isinstance(value, dict):
        return float(value["mean"]), float(value.get("error", 0.0))
    mean, error = value.split("±")
    return float(mean), float(error)


def window_index(path: Path) -> int:
    # Faunus nests per-window files as umbrella_states/window{i}/output.yaml (PR #61ccb6a);
    # older runs used umbrella_states/window{i}_output.yaml. Read the index from either.
    return int(re.search(r"window(\d+)", str(path)).group(1))


def wall_alignment(shape) -> float | str:
    """P2 order parameter of the disk-like phytate's gyration tensor about the wall
    normal z: s_wall = (3 Szz/Rg^2 - 1)/2, with Rg^2 = Sxx + Syy + Szz.

    0 = isotropic; -1/2 = face-on (disk plane parallel to the wall); +1 = edge-on.
    A lab-frame *pose* diagnostic, distinct from the P--O-proxy ensemble s2.
    Returns "" when the gyration-tensor columns are absent (shapeless window).
    """
    try:
        sxx, syy, szz = shape["Sxx"], shape["Syy"], shape["Szz"]
    except KeyError:
        return ""
    rg2 = sxx + syy + szz
    return (3 * szz / rg2 - 1) / 2 if rg2 else ""


def parse_window(output_yaml: Path):
    """Return (z, charges, shape_dict, widom) for one window.

    `charges` is a list of (selection, charge, err); `widom` is the tuple
    (s2, s2_err, widom_w) from the WidomRotation analysis, or (None, None, None)
    for runs predating it.
    """
    doc = load_yaml(output_yaml)
    z = None
    charges = []
    shape = None
    widom = (None, None, None)
    for item in doc["analysis"]:
        ((key, val),) = item.items()  # each analysis entry is a single-key mapping
        # Faunus PR #74 switched output-section names to snake_case (widomrotation ->
        # widom_rotation, ...). Strip underscores so both the pre- and post-#74 layouts
        # parse -- lets the archived runs.zip and freshly generated runs share this code.
        key = key.replace("_", "")
        if key == "collectivevariable" and val.get("property") == "AtomPosition":
            z = val["mean"]  # mean phytate (INO) z for this window
        elif key == "multipole":
            mean, error = parse_pm(val["charge"])
            charges.append((val["selection"], mean, error))
        elif key == "polymershape":
            shape = {k: v for k, v in val.items() if k != "num_samples"}
        elif key == "widomrotation":
            widom = (val["mean_S2"]["mean"], val["mean_S2"]["error"], val["W/kJ/mol"]["mean"])
    return z, charges, shape, widom


def mean_profile(profile_csvs, cols) -> tuple:
    """Average the requested data columns over a set of same-grid profiles.

    Columns of potential_*.csv: z/Å (0), slab_charge_density/e·Å⁻² (1),
    charge_density/e·Å⁻³ (2), potential/mV (3), error (4), field (5).
    Columns of peptide_profile.csv: z/Å (0), density/Å⁻³ (1), its error (2), then
    molarity and mass density with their errors.
    The grid is identical across windows. Returns (z, mean_col, ...) for each column in `cols`.
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

        # Brush monomer density -> the Zwanzig roughness profile. DensityProfile lays out its
        # own slabs, so map it onto the potential grid rather than assuming shared bin centres.
        dens_csvs = sorted(run.glob("umbrella_states/window*/peptide_profile.csv"))
        if dens_csvs:
            z_rho, rho = mean_profile(dens_csvs, (1,))
            rho_peptide = np.interp(z, z_rho, rho)
        else:
            rho_peptide = [""] * len(z)  # older runs predate the DensityProfile analysis

        for ri, pa, pp, sp, rho_i in zip(
            r_prof, pot_all, pot_peptide, sigma_peptide, rho_peptide, strict=True
        ):
            potential_rows.append([ph, ionic, ri, pa, pp, sp, rho_i])

        # Zeta at the fixed-charge edge of the brush (sigma_peptide -> 0).
        r_shear, zeta = shear_zeta(r_prof, sigma_peptide, pot_all, half_lz)
        zeta_rows.append([ph, ionic, r_shear, zeta])

        # Per-window observables share the window-centre r-grid, so each window yields
        # one merged row.
        windows = sorted(
            [
                *run.glob("umbrella_states/window*/output.yaml"),  # nested layout (current)
                *run.glob("umbrella_states/window*_output.yaml"),  # flat layout (archived runs)
            ],
            key=window_index,
        )
        for win in windows:
            z, charges, shape, (s2, s2_err, widom_w) = parse_window(win)
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
                    "" if s2 is None else s2,
                    "" if s2_err is None else s2_err,
                    "" if widom_w is None else widom_w,
                    wall_alignment(clean),
                    clean,
                ]
            )

    window_rows.sort(key=lambda row: (row[0], row[1], row[2]))
    potential_rows.sort(key=lambda row: (row[0], row[1], row[2]))
    zeta_rows.sort(key=lambda row: (row[0], row[1]))

    shape_fields = shape_fields or []
    with open(args.window, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(
            [
                "pH",
                "ionic_strength",
                "r",
                "q_ino",
                "q_ino_err",
                "q_peptide",
                "q_peptide_err",
                "s2",
                "s2_err",
                "widom_w",
                "s_wall",
                *shape_fields,
            ]
        )
        for *base, clean in window_rows:
            w.writerow([*base, *(clean.get(fld, "") for fld in shape_fields)])
    print(f"Saved {args.window} ({len(window_rows)} rows)")

    with open(args.potential, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(
            [
                "pH",
                "ionic_strength",
                "r",
                "potential_all_mV",
                "potential_peptide_mV",
                "sigma_peptide_e_A2",
                "rho_peptide_A3",
            ]
        )
        w.writerows(potential_rows)
    print(f"Saved {args.potential} ({len(potential_rows)} rows)")

    with open(args.zeta, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["pH", "ionic_strength", "r_shear", "zeta_mV"])
        w.writerows(zeta_rows)
    print(f"Saved {args.zeta} ({len(zeta_rows)} rows)")


if __name__ == "__main__":
    main()
