#!/usr/bin/env python3
"""Composite MC figure (PMF, rotational mobility, potential, phytate charge) by pH,
for one ionic strength -- the salt-resolved companion to ``composite_lowsalt.png``.

Reproduces the four-panel low-salt composite from ``analysis.ipynb`` for any simulated
ionic strength, so the SI can show how the well and the double layer screen with salt.
Run: ``uv run python plot_composites.py`` (renders I = 0.005, 0.03, 0.15 M).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
PHS = (4.0, 6.0, 10.0)
PH_COLOR = {4.0: "C0", 6.0: "C1", 10.0: "C2"}  # blue / orange / green, matching the manuscript


def _sel(df: pd.DataFrame, ionic: float, pH: float) -> pd.DataFrame:
    return df[(df.ionic_strength == ionic) & (df.pH == pH)].sort_values("r")


def plot_composite(ionic: float, path: Path) -> None:
    """Four-panel composite at ionic strength ``ionic`` (M)."""
    pmf = pd.read_csv(HERE / "pmf_profiles.csv")
    win = pd.read_csv(HERE / "window_profiles.csv")  # phytate charge + rotational mobility
    pot = pd.read_csv(HERE / "potential_profiles.csv")
    zeta = pd.read_csv(HERE / "simulated_zeta.csv")
    # Measured bare-OB zeta vs added salt (PA = 0), interpolated to this ionic strength.
    exp = pd.read_csv(HERE.parent / "zeta" / "measured_zeta.csv")
    exp = exp[(exp.series == "NaCl") & (exp.PA_mM == 0)]
    nacl = max(0.0, (ionic - 0.003) * 1000)  # ionic ~ 3 mM buffer + added NaCl

    fig, axes = plt.subplots(2, 2, figsize=(8, 6), sharex=True)
    (axA, axB), (axC, axD) = axes
    for pH in PHS:
        c = PH_COLOR[pH]
        p = _sel(pmf, ionic, pH)
        axA.fill_between(p.r, p.pmf_kT - p.stderr_kT, p.pmf_kT + p.stderr_kT, color=c, alpha=0.15)
        axA.plot(p.r, p.pmf_kT, color=c, label=f"pH {pH:g}")

        m = _sel(win, ionic, pH)
        axB.plot(m.r, m.msd_rot / m.msd_rot.iloc[-3:].mean(), color=c)

        v = _sel(pot, ionic, pH)
        axC.plot(v.r, v.potential_all_mV, color=c)
        zc = zeta[(zeta.ionic_strength == ionic) & (zeta.pH == pH)].iloc[0]
        axC.plot(zc.r_shear, zc.zeta_mV, "o", mfc="white", mec=c, zorder=5)  # simulated zeta
        e = exp[exp.pH == pH].sort_values("NaCl_mM")
        if not e.empty:
            # np.interp clamps beyond the measured range (max 100 mM added NaCl), so the
            # 150 mM panel shows the highest-salt measured value (the screened plateau).
            xp = e.NaCl_mM.to_numpy()
            mean = float(np.interp(nacl, xp, e.mean_mV.to_numpy()))
            sd = float(np.interp(nacl, xp, e.sd_mV.to_numpy()))
            # Measured zeta has no intrinsic r; the +4 A offset only separates it from the
            # simulated marker. Read this point by height (zeta), not by its x-position.
            axC.errorbar(zc.r_shear + 4, mean, yerr=sd, fmt="s", color=c, capsize=3, zorder=5)

        q = _sel(win, ionic, pH)
        axD.plot(q.r, q.q_ino, color=c, marker="o", ms=3, lw=1)

    axA.axhline(0, color="gray", ls="--", lw=0.5)
    axA.set_ylabel(r"PMF ($k_\mathrm{B}T$)")
    axB.axhline(1, color="gray", ls="--", lw=0.5)
    axB.set_ylabel(r"$D(r)/D_0$ (rotational)")
    axC.axhline(0, color="gray", ls="--", lw=0.5)
    axC.set_ylabel("electric potential (mV)")
    axC.set_ylim(-90, 60)
    axD.set_ylabel(r"phytate charge ($e$)")
    axC.set_xlabel(r"$r$ (Å)")
    axD.set_xlabel(r"$r$ (Å)")
    for ax, letter in zip(axes.flat, "ABCD", strict=True):
        ax.set_xlim(0, 100)
        ax.grid(axis="x", linestyle="--", alpha=0.5)
        ax.grid(axis="y", linestyle="--", alpha=0.5)
        ax.text(0.03, 0.94, letter, transform=ax.transAxes, fontweight="bold", va="top")
    axA.legend(frameon=False, title=f"I = {ionic:g} M")
    axC.legend(
        handles=[
            Line2D([], [], marker="o", mfc="white", mec="gray", ls="", label="sim ζ"),
            Line2D([], [], marker="s", color="gray", ls="", label="exp ζ"),
        ],
        frameon=False,
        fontsize="small",
        loc="upper right",
    )
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    for ionic in (0.005, 0.03, 0.15):
        out = HERE / f"composite_I{ionic:g}.png"
        plot_composite(ionic, out)
        print(f"saved {out.name}")
