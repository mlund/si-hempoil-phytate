#!/usr/bin/env python3
"""Compute peptide bead number density along z from an XTC trajectory.

Usage:
    python peptide_density_z.py [run_dir]

Reads `{run_dir}/traj.psf` + `{run_dir}/traj.xtc` (default: cwd) and writes
`{run_dir}/peptide_density_z.csv` with columns `z_A, density_per_A3`.

Peptide atoms = everything except phytate beads (INO, PH0..PH12) and
ghost H+ particles. Density is averaged over all frames.
"""
import sys
from pathlib import Path

import MDAnalysis as mda
import numpy as np

NON_PEPTIDE = {"INO", "H+"} | {f"PH{i}" for i in range(13)}
BIN_WIDTH = 1.0  # Å


def main(run_dir: Path) -> None:
    u = mda.Universe(str(run_dir / "traj.psf"), str(run_dir / "traj.xtc"))
    is_peptide = ~np.isin(u.atoms.names, list(NON_PEPTIDE))
    peptide = u.atoms[is_peptide]
    if peptide.n_atoms == 0:
        sys.exit("no peptide atoms found")

    Lx, Ly, Lz = u.dimensions[:3]
    # XTC stores positions in [0, L]; z=0 is the grafting wall.
    z_edges = np.arange(0.0, Lz + BIN_WIDTH, BIN_WIDTH)
    counts = np.zeros(len(z_edges) - 1)
    for _ in u.trajectory:
        c, _ = np.histogram(peptide.positions[:, 2], bins=z_edges)
        counts += c
    n_frames = u.trajectory.n_frames

    density = counts / (Lx * Ly * BIN_WIDTH * n_frames)
    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])

    out = run_dir / "peptide_density_z.csv"
    np.savetxt(out, np.c_[z_centers, density],
               delimiter=",", header="z_A,density_per_A3", comments="")
    print(f"wrote {out}  ({n_frames} frames, {peptide.n_atoms} peptide atoms)")


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else Path("."))
