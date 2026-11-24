#!/usr/bin/env python3
"""Flatten the tidied zeta-potential workbook into one long CSV.

The workbook has three sheets (vs [PA], vs pH, vs [NaCl]); each holds several
labelled blocks of ``x, n, mean, SD`` rows. We parse them into tidy records

    series, pH, PA_mM, NaCl_mM, n, mean_mV, sd_mV

where ``series`` is 'PA', 'pH' or 'NaCl' (the quantity varied). Downstream the
zeta model reads means *and* SDs so the model-data comparison can be weighted by
the measurement uncertainty (SEM = SD/sqrt(n)).
"""

import re
from pathlib import Path

import pandas as pd

SRC = Path(__file__).parent / "Zeta_potential_tidied.xlsx"
OUT = Path(__file__).parent / "measured_zeta.csv"


def _num(pattern: str, text: str) -> float | None:
    m = re.search(pattern, text)
    return float(m.group(1)) if m else None


def parse() -> pd.DataFrame:
    """Return the long-format table of every (condition, mean, SD) in the workbook."""
    xl = pd.ExcelFile(SRC)
    rows = []
    for sheet in xl.sheet_names:
        raw = xl.parse(sheet, header=None)
        block = None
        for _, row in raw.iterrows():
            c0, c1 = row.iloc[0], row.iloc[1]
            if isinstance(c1, str) and c1.strip() == "n":  # column-header row
                continue
            if isinstance(c0, str) and pd.isna(c1):  # a block label (rest of row empty)
                if not c0.startswith("Effect of"):  # skip the sheet title
                    block = c0.strip()
                continue
            if isinstance(c0, str) or pd.isna(c0) or pd.isna(c1):  # not a data row
                continue

            x, n, mean, sd = float(c0), int(c1), float(row.iloc[2]), float(row.iloc[3])
            if "PA" in sheet and "pH" not in block.lower()[:3]:  # 'Effect PA' -> block 'pH X'
                rows.append(("PA", _num(r"pH\s*(\d+)", block), x, 0.0, n, mean, sd))
            elif "pH" in sheet:  # 'Effect pH' -> block '[PA] = X mM'
                rows.append(("pH", x, _num(r"=\s*([\d.]+)\s*mM", block), 0.0, n, mean, sd))
            else:  # 'Effect ionic strength - NaCl' -> block 'X mM PA - pH Y'
                pa = _num(r"([\d.]+)\s*mM PA", block)
                ph = _num(r"pH\s*(\d+)", block)
                rows.append(("NaCl", ph, pa, x, n, mean, sd))

    return pd.DataFrame(rows, columns=["series", "pH", "PA_mM", "NaCl_mM", "n", "mean_mV", "sd_mV"])


def main() -> None:
    """Write the parsed zeta workbook as a tidy CSV."""
    df = parse()
    df.to_csv(OUT, index=False)
    print(f"wrote {OUT.name}: {len(df)} rows, pH {sorted(df.pH.unique())}")


if __name__ == "__main__":
    main()
