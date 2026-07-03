#!/usr/bin/env python3
"""Hydrodynamic diffusion of free phytate across protonation states, via PyGRPY.

For each net charge, an RDKit conformer ensemble (ETKDGv3 + MMFF94, UFF fallback
for the polyanion) of myo-inositol hexakisphosphate is built into an all-atom bead
model (one bead per atom, Bondi vdW radii). From the 6x6 GRPY grand mobility matrix
we report the ensemble-averaged

    tau_r  rotational correlation time   1 / (6 <D_r>)        [ps]
    D_r    rotational diffusion          isotropic average    [rad^2/s]
    D_t    translational diffusion       isotropic average    [m^2/s]
    a_H    Stokes (hydrodynamic) radius   kT / (6 pi eta D_t)  [A]

and writes a tidy diffusion.csv (one row per charge).

Rough-estimate caveats: PyGRPY uses stick boundary conditions (slip could lower
tau_r by up to ~3x for a sub-nm solute) and models no ionic-atmosphere friction;
net charge therefore enters only through the bare hydrodynamic shape (it removes
peripheral H beads). See the si-hempoil-phytate analysis notes.
"""

import argparse
from pathlib import Path

import numpy as np
from pygrpy.grpy_tensors import conglomerateMobilityMatrixAtCentre
from rdkit import Chem
from rdkit.Chem import AllChem

KB = 1.380649e-23  # J / K
ANG = 1e-10  # m

# Bondi vdW radii [Angstrom]
VDW = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "P": 1.80, "S": 1.80}

PHYTIC_ACID = (
    "[C@@H]1([C@@H]([C@@H]([C@@H]([C@H]([C@@H]1OP(=O)(O)O)"
    "OP(=O)(O)O)OP(=O)(O)O)OP(=O)(O)O)OP(=O)(O)O)OP(=O)(O)O"
)


def phytate_smiles(charge):
    """SMILES for myo-inositol hexakisphosphate carrying net charge -`charge`.

    `charge` is the number of dissociated protons (0 = neutral H12PA, 12 = PA^12-);
    protons come off the six phosphates evenly (the first from each before any
    second), matching the pKa ordering. The exact pattern is hydrodynamically
    immaterial.
    """
    if not 0 <= charge <= 12:
        raise SystemExit("charge must be an integer in 0..12")
    doubly = max(0, charge - 6)
    singly = charge - 2 * doubly
    groups = (
        ["P(=O)([O-])[O-]"] * doubly
        + ["P(=O)([O-])O"] * singly
        + ["P(=O)(O)O"] * (6 - doubly - singly)
    )
    parts = PHYTIC_ACID.split("P(=O)(O)O")  # 7 fragments around the 6 phosphates
    # 7 parts vs 6 groups: the trailing fragment is appended separately, so the
    # length mismatch is intentional (strict=False).
    return "".join(part + grp for part, grp in zip(parts, groups, strict=False)) + parts[-1]


def make_conformers(smiles, n, seed):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.pruneRmsThresh = 0.5
    AllChem.EmbedMultipleConfs(mol, numConfs=n, params=params)
    for cid in range(mol.GetNumConformers()):
        # MMFF94 lacks parameters for phosphate anions (return -1 = no setup);
        # fall back to UFF so charged states get a relaxed geometry.
        if AllChem.MMFFOptimizeMolecule(mol, confId=cid) == -1:
            AllChem.UFFOptimizeMolecule(mol, confId=cid)
    return mol


def coords_and_radii(mol, conf_id, hydration):
    conf = mol.GetConformer(conf_id)
    centres, radii = [], []
    for atom in mol.GetAtoms():
        p = conf.GetAtomPosition(atom.GetIdx())
        centres.append([p.x, p.y, p.z])
        radii.append(VDW[atom.GetSymbol()] + hydration)
    return np.asarray(centres), np.asarray(radii)


def diffusion(centres, radii, eta, T):
    """Isotropic (D_t, D_r) and rotational anisotropy from the GRPY mobility."""
    M = conglomerateMobilityMatrixAtCentre(centres, radii)
    D_t = KB * T * np.linalg.eigvalsh(M[0:3, 0:3] / (eta * ANG))  # m^2/s
    D_r = KB * T * np.linalg.eigvalsh(M[3:6, 3:6] / (eta * ANG**3))  # rad^2/s
    return D_t.mean(), D_r.mean(), D_r.max() / D_r.min()


def gyration_radius(centres):
    return float(np.sqrt(np.mean(np.sum((centres - centres.mean(0)) ** 2, axis=1))))


def scan_charge(charge, args):
    """Ensemble-averaged diffusion observables for one protonation state."""
    mol = make_conformers(phytate_smiles(charge), args.n, args.seed)
    rows = []
    for i in range(mol.GetNumConformers()):
        centres, radii = coords_and_radii(mol, i, args.hydration)
        d_t, d_r, aniso = diffusion(centres, radii, args.viscosity, args.temperature)
        rows.append((gyration_radius(centres), d_t, d_r, aniso))
    rg, d_t, d_r, aniso = np.mean(rows, axis=0)
    a_H = KB * args.temperature / (6 * np.pi * args.viscosity * d_t) / ANG
    return {
        "charge": -charge,
        "n_conf": mol.GetNumConformers(),
        "Rg_A": rg,
        "tau_r_ps": 1e12 / (6 * d_r),
        "D_rot_rad2_s": d_r,
        "D_rot_anisotropy": aniso,
        "D_trans_m2_s": d_t,
        "a_H_A": a_H,
    }


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--charges",
        default="0,3,6,9,12",
        help="comma-separated net-charge magnitudes to scan (0..12)",
    )
    ap.add_argument("--n", type=int, default=20, help="conformers per charge")
    ap.add_argument("--seed", type=int, default=1, help="RDKit seed (avoid 0: embeds 1 conf)")
    ap.add_argument("--temperature", type=float, default=298.15, help="K")
    ap.add_argument("--viscosity", type=float, default=8.9e-4, help="solvent [Pa*s]")
    ap.add_argument(
        "--hydration",
        type=float,
        default=0.0,
        help="shell thickness [A] added to every bead (try 1.4 for a water shell)",
    )
    ap.add_argument("--csv", type=Path, default=Path(__file__).parent / "diffusion.csv")
    args = ap.parse_args()

    charges = [int(c) for c in args.charges.split(",")]
    results = [scan_charge(c, args) for c in charges]

    fields = list(results[0])
    print(
        f"{'q':>4} {'n':>3} {'Rg/A':>6} {'tau_r/ps':>9} {'D_rot':>10} "
        f"{'aniso':>6} {'D_t/m2s':>10} {'a_H/A':>6}"
    )
    for r in results:
        print(
            f"{r['charge']:>4d} {r['n_conf']:>3d} {r['Rg_A']:>6.2f} "
            f"{r['tau_r_ps']:>9.1f} {r['D_rot_rad2_s']:>10.2e} "
            f"{r['D_rot_anisotropy']:>6.2f} {r['D_trans_m2_s']:>10.2e} {r['a_H_A']:>6.2f}"
        )

    import csv

    with open(args.csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(results)
    print(f"\nSaved {args.csv}")


if __name__ == "__main__":
    main()
