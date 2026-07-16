#!/usr/bin/env python3
"""Draw the LC--AR--PB coupling schematic used in Supplementary Section S-2."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.offsetbox import AnchoredOffsetbox, HPacker, TextArea, VPacker
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

HERE = Path(__file__).resolve().parent
BLUE = "#356aa0"
GREY = "#777777"
RED = "#c44e52"
INPUT = "#73579d"
INPUT_FILL = "#eee8f5"


def _box(ax, xy, width, height, text, *, color="#f3f3f3", fontsize=8.5):
    patch = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle="round,pad=0.02,rounding_size=0.025",
        facecolor=color,
        edgecolor="none",
        linewidth=0,
    )
    ax.add_patch(patch)
    ax.text(
        xy[0] + width / 2, xy[1] + height / 2, text, ha="center", va="center", fontsize=fontsize
    )
    return patch


def _arrow(ax, start, end, *, color="#333333", style="-|>", connection="arc3", linewidth=1.2):
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle=style,
        mutation_scale=11,
        color=color,
        linewidth=linewidth,
        connectionstyle=connection,
    )
    ax.add_patch(arrow)


def _binding_equation(ax) -> None:
    """Add the binding-potential decomposition with mechanism-matched colors."""

    text = {"size": 8.2}
    first_line = HPacker(
        children=[
            TextArea(r"$\psi_\beta(\Gamma)=$", textprops=text),
            TextArea(
                r"$\psi_{\mathrm{pep}}^{\mathrm{MC}}(R_\beta)$", textprops={**text, "color": BLUE}
            ),
        ],
        align="center",
        pad=0,
        sep=1,
    )
    second_line = HPacker(
        children=[
            TextArea(r"$+[$", textprops=text),
            TextArea(r"$\psi_d(\sigma_{\mathrm{eff}})$", textprops={**text, "color": RED}),
            TextArea(r"$-$", textprops=text),
            TextArea(r"$\psi_d(\sigma_0)$", textprops={**text, "color": GREY}),
            TextArea(r"$]$", textprops=text),
        ],
        align="center",
        pad=0,
        sep=1,
    )
    equation = VPacker(children=[first_line, second_line], align="left", pad=0, sep=2)
    ax.add_artist(
        AnchoredOffsetbox(
            loc="lower left",
            child=equation,
            frameon=False,
            bbox_to_anchor=(0.045, 0.085),
            bbox_transform=ax.transAxes,
            borderpad=0,
        )
    )


def _potential_panel(ax) -> None:
    r = np.linspace(0.0, 70.0, 500)
    r_beta = 26.0
    r_shear = 41.0

    # Illustrative shapes only: the calculation samples the MC curve at R_beta and
    # uses the PB boundary-potential difference, rather than adding two radial grids.
    peptide = 29.0 * np.exp(-(((r - 15.0) / 18.0) ** 2)) - 1.5
    rng = np.random.default_rng(7)
    mc_fluctuation = np.convolve(rng.normal(size=r.size), np.ones(9) / 9, mode="same")
    peptide += 0.65 * mc_fluctuation
    diffuse_r = r[r >= r_shear]
    pb_bare = 20.0 * np.exp(-(diffuse_r - r_shear) / 15.0)
    pb_covered = -10.0 * np.exp(-(diffuse_r - r_shear) / 15.0)

    ax.add_patch(Rectangle((0, -20), 23, 56, facecolor="#dbe9f6", edgecolor="none", alpha=0.65))
    ax.spines["left"].set_linewidth(1.4)
    ax.axvline(r_beta, color=BLUE, linestyle=":", linewidth=1.1)
    ax.axvline(r_shear, color=GREY, linestyle=":", linewidth=1.0)
    ax.axhline(0, color=GREY, linewidth=0.6)
    ax.text(
        2,
        1.3,
        r"bulk reference, $\psi=0$",
        ha="left",
        va="bottom",
        fontsize=7.5,
        color="#666666",
    )

    ax.plot(
        r,
        peptide,
        color=BLUE,
        linewidth=2.0,
        label=r"fixed $\psi_{\mathrm{pep}}^{\mathrm{MC}}(r)$",
    )
    ax.plot(
        diffuse_r,
        pb_bare,
        color=GREY,
        linewidth=1.6,
        linestyle="--",
        label=r"$\psi_d(r;\sigma_0)$",
    )
    ax.plot(
        diffuse_r,
        pb_covered,
        color=RED,
        linewidth=1.8,
        linestyle="-.",
        label=r"$\psi_d(r;\sigma_{\mathrm{eff}})$",
    )

    y_beta = np.interp(r_beta, r, peptide)
    ax.plot(r_beta, y_beta, "o", color=BLUE, markersize=4.5)
    ax.annotate(
        "MC potential sampled at\n" r"PMF minimum $R_\beta$",
        xy=(r_beta, y_beta),
        xytext=(33, 28),
        arrowprops={"arrowstyle": "->", "color": BLUE, "lw": 1.0},
        fontsize=8,
        color="#244d75",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9, "pad": 1.5},
    )
    ax.text(8, 33, "peptide brush", ha="center", va="top", fontsize=8, color="#244d75")
    ax.text(
        r_beta + 1.2,
        -18,
        r"$R_\beta=26$ Å",
        rotation=90,
        va="bottom",
        fontsize=7.5,
        color=BLUE,
    )
    ax.text(
        r_shear + 1.2,
        -18,
        "effective PB / shear plane",
        rotation=90,
        va="bottom",
        fontsize=7.5,
        color=GREY,
    )
    _binding_equation(ax)

    ax.set_xlim(0, 70)
    ax.set_ylim(-21, 37)
    ax.set_xlabel(r"distance from grafting wall, $r$")
    ax.set_ylabel("potential (schematic)")
    ax.set_xticks([])
    ax.set_yticks([])
    legend = ax.legend(frameon=False, fontsize=7.6, loc="upper right")
    for label, color in zip(legend.get_texts(), (BLUE, GREY, RED), strict=True):
        label.set_color(color)
    ax.set_title(
        "A  Molecular input and continuum potential", loc="left", fontsize=10, fontweight="bold"
    )
    ax.text(
        0.99,
        0.02,
        "curves not to scale",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=7,
        color="#666666",
    )


def _flow_panel(ax) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    _box(ax, (0.03, 0.72), 0.28, 0.18, "", color="#dbe9f6")
    ax.text(0.17, 0.84, "fixed MC input", ha="center", va="center", fontsize=8.5)
    ax.text(
        0.17,
        0.79,
        r"$\psi_{\rm pep}^{MC}(R_\beta)$ and $z$",
        ha="center",
        va="center",
        fontsize=8.5,
        color=BLUE,
    )
    _box(ax, (0.36, 0.72), 0.28, 0.18, "", color=INPUT_FILL)
    ax.text(0.50, 0.84, "experiment", ha="center", va="center", fontsize=8.5)
    ax.text(
        0.50,
        0.79,
        r"$\zeta_0\rightarrow\sigma_0$",
        ha="center",
        va="center",
        fontsize=8.5,
        color=INPUT,
    )
    _box(ax, (0.69, 0.72), 0.28, 0.18, "", color=INPUT_FILL)
    ax.text(0.83, 0.84, "solution and surface", ha="center", va="center", fontsize=8.5)
    ax.text(
        0.83,
        0.79,
        r"$c_L,K,I,\Gamma_{\max}$",
        ha="center",
        va="center",
        fontsize=8.5,
        color=INPUT,
    )
    _box(ax, (0.30, 0.45), 0.32, 0.16, "", color="#f2f2f2")
    ax.text(0.46, 0.55, r"evaluate $\Omega(\Gamma)$", ha="center", va="center", fontsize=8.5)
    ax.text(0.46, 0.50, r"at trial $\Gamma$", ha="center", va="center", fontsize=8.5)
    _box(ax, (0.69, 0.45), 0.28, 0.16, "", color="#f7dddd")
    ax.text(
        0.83,
        0.55,
        r"$\sigma_{\rm eff}=\sigma_0-z\Gamma$",
        ha="center",
        va="center",
        fontsize=8.5,
    )
    ax.text(
        0.83,
        0.50,
        r"$\psi_d$ and $g_{\rm PB}$",
        ha="center",
        va="center",
        fontsize=8.5,
        color=RED,
    )
    _box(
        ax,
        (0.30, 0.14),
        0.32,
        0.15,
        "minimize $\\Omega(\\Gamma)$\n$\\Gamma^*$ and $\\zeta$",
        color="#f2f2f2",
    )

    _arrow(ax, (0.17, 0.72), (0.36, 0.61), color=BLUE)
    _arrow(ax, (0.50, 0.72), (0.46, 0.61), color=INPUT)
    _arrow(ax, (0.83, 0.72), (0.56, 0.61), color=INPUT)
    _arrow(ax, (0.62, 0.54), (0.69, 0.54), color=RED)
    ax.text(0.655, 0.555, r"$\Gamma$", fontsize=7.5, ha="center", va="bottom")
    _arrow(ax, (0.46, 0.45), (0.46, 0.29))
    _arrow(
        ax,
        (0.91, 0.45),
        (0.62, 0.48),
        color=RED,
        connection="arc3,rad=-0.55",
    )

    ax.set_title(
        "B  Information flow and self-consistency", loc="left", fontsize=10, fontweight="bold"
    )


def main() -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.1), gridspec_kw={"width_ratios": [1.15, 1.0]})
    fig.patch.set_facecolor("white")
    for axis in axes:
        axis.set_facecolor("white")
    _potential_panel(axes[0])
    _flow_panel(axes[1])
    fig.tight_layout(w_pad=1.8)
    fig.savefig(
        HERE / "fig_lc_arpb_schematic.png",
        dpi=220,
        bbox_inches="tight",
        facecolor="white",
        transparent=False,
    )
    plt.close(fig)


if __name__ == "__main__":
    main()
