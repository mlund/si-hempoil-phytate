#!/usr/bin/env python3
"""Generate per-(pH, I) run directories with SLURM scripts for a parameter scan.

For each (pH, I) combination this writes `runs/pH{pH}_I{salt}/` containing:
  - input.yaml     copy of input.yaml with pH and salt MiniJinja vars substituted
  - input_eq.yaml  derived from input.yaml for the short equilibration
                   (repeat=1000, umbrella section disabled)
  - state.yaml     copy of the master state (each job equilibrates its own copy)
  - phytate.xyz    symlink to the shared structure
  - job.slurm      SLURM batch script: short equilibration then umbrella sampling

The script only *creates* files; it never runs faunus itself. Pass --submit to
`sbatch` each job from its own directory (skipping points that already have a
pmf.csv, unless --rerun).

Each SLURM job requests a single node with --cores CPUs and runs, in its dir:
    faunus run     -i input_eq.yaml -s state.yaml       # short equilibration
    faunus umbrella -i input.yaml --state state.yaml \
                    -o pmf.csv -j <cores>               # umbrella sampling -> PMF
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).parent
STRUCTURE = ROOT / "phytate.xyz"

# Kept as strings so the folder name and the YAML value stay byte-identical
# (e.g. "0.030", "1000") rather than risking float reformatting.
PH_RANGE = [4, 6, 10]
SALT_RANGE = ["0.005", "0.030", "0.15"]
MATRIX = [(pH, salt) for pH in PH_RANGE for salt in SALT_RANGE]


def render(template: str, pH: int, salt: str) -> str:
    out = re.sub(r"\{% set pH = [^%]*%\}", f"{{% set pH = {pH}.0 %}}", template, count=1)
    out = re.sub(r"\{% set salt = [^%]*%\}", f"{{% set salt = {salt} %}}", out, count=1)
    return out


def to_equilibration(template: str) -> str:
    """Derive the short-equilibration template from the production input:
    shorten the outer propagate cycle and disable the umbrella section."""
    out = re.sub(r"^  repeat: \d+$", "  repeat: 1000", template, count=1, flags=re.M)
    out = re.sub(r"^umbrella:$", "_umbrella:", out, count=1, flags=re.M)
    return out


def slurm_script(pH: int, salt: str, args: argparse.Namespace) -> str:
    directives = [
        f"#SBATCH --job-name=pH{pH}_I{salt}",
        f"#SBATCH -A {args.project}",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        f"#SBATCH --cpus-per-task={args.cores}",
        f"#SBATCH --time={args.time}",
    ]
    if args.partition:
        directives.append(f"#SBATCH --partition={args.partition}")
    directives.append("#SBATCH --output=slurm-%j.out")
    header = "\n".join(directives)

    return f"""#!/bin/bash
{header}

set -euo pipefail

# input.yaml/state.yaml are per-dir, so resolve paths against the submit dir;
# fall back to the script's own location when run by hand.
cd "${{SLURM_SUBMIT_DIR:-$(dirname "$0")}}"

FAUNUS="{args.faunus}"

# Relax this (pH, I) state first — every umbrella window seeds from it.
$FAUNUS run -i input_eq.yaml -s state.yaml

# One walker per window, one thread per core.
$FAUNUS umbrella -i input.yaml --state state.yaml -o pmf.csv \\
    -j "${{SLURM_CPUS_PER_TASK:-{args.cores}}}"
"""


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--input",
        type=Path,
        default=ROOT / "input.yaml",
        help="template input file (default: input.yaml)",
    )
    p.add_argument(
        "--state",
        type=Path,
        default=ROOT / "state.yaml",
        help="master state file copied into each run (default: state.yaml)",
    )
    p.add_argument(
        "--project",
        default="lu2026-2-71",
        help="SLURM project/account ID, -A (default: lu2026-2-71)",
    )
    p.add_argument(
        "-N",
        "--cores",
        type=int,
        default=40,
        help="CPU cores per job on one node (default: 40, matches num_windows)",
    )
    p.add_argument("--time", default="24:00:00", help="SLURM wall-clock limit (default: 24:00:00)")
    p.add_argument("--partition", default=None, help="SLURM partition/queue (optional)")
    p.add_argument("--faunus", default="faunus", help="faunus executable/command (default: faunus)")
    p.add_argument(
        "--submit", action="store_true", help="sbatch each job from its directory after generating"
    )
    p.add_argument(
        "--rerun",
        action="store_true",
        help="with --submit, resubmit even if pmf.csv exists (clears umbrella_states)",
    )
    args = p.parse_args()

    for f in (args.input, args.state, STRUCTURE):
        if not f.is_file():
            p.error(f"{f} not found")

    template = args.input.read_text()
    template_eq = to_equilibration(template)
    runs = ROOT / "runs"
    runs.mkdir(exist_ok=True)

    for pH, salt in MATRIX:
        d = runs / f"pH{pH}_I{salt}"
        d.mkdir(exist_ok=True)

        (d / "input.yaml").write_text(render(template, pH, salt))
        (d / "input_eq.yaml").write_text(render(template_eq, pH, salt))
        shutil.copy(args.state, d / "state.yaml")
        (d / "job.slurm").write_text(slurm_script(pH, salt, args))

        link = d / "phytate.xyz"
        link.unlink(missing_ok=True)
        link.symlink_to(Path("..") / ".." / STRUCTURE.name)

        print(f"  {d.relative_to(ROOT)}")

    if not args.submit:
        print("\nGenerated. Submit with --submit, or: cd runs/<dir> && sbatch job.slurm")
        return 0

    print("\nSubmitting:")
    for pH, salt in MATRIX:
        d = runs / f"pH{pH}_I{salt}"
        if (d / "pmf.csv").is_file() and not args.rerun:
            print(f"  SKIP {d.name} (pmf.csv exists; use --rerun to force)")
            continue
        if args.rerun:
            shutil.rmtree(d / "umbrella_states", ignore_errors=True)
        try:
            out = subprocess.run(
                ["sbatch", "job.slurm"], cwd=d, check=True, capture_output=True, text=True
            )
            print(f"  {d.name}: {out.stdout.strip()}")
        except FileNotFoundError:
            p.error("sbatch not found — run --submit on the cluster")
        except subprocess.CalledProcessError as e:
            print(f"  FAILED {d.name}: {e.stderr.strip()}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
