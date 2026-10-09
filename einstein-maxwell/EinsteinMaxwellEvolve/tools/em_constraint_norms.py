#!/usr/bin/env python3
"""Normalize fixed-shell CarpetIOScalar norm2 outputs without new 3-D dumps.

This is a CarpetReduce-weighted discrete RMS, NOT a proper-volume integral.
Restart duplicates are resolved by later output-NNNN files. Require matching
timestamps for both fields and the mask; never silently interpolate residuals.
"""
import argparse
import csv
import math
import re
from pathlib import Path


def read_series(root, variable, reduction):
    files = sorted(Path(root).rglob(f"*.{reduction}.asc"))
    if not files:
        raise ValueError(f"Missing {variable}.{reduction}.asc under {root}")
    values = {}
    for path in files:
        column = None
        with path.open() as stream:
            for line in stream:
                if line.lstrip().startswith("#"):
                    if "data columns:" in line:
                        for number, name in re.findall(r"(\d+):([^\s]+)", line.split("data columns:", 1)[1]):
                            if name.split("::")[-1].lower() == variable.lower():
                                column = int(number) - 1
                    continue
                if not line.strip():
                    continue
                if column is None:
                    if path.name.lower().endswith(f"{variable.lower()}.{reduction}.asc"):
                        column = 2
                    else:
                        continue
                row = [float(value) for value in line.split()]
                if len(row) <= column:
                    raise ValueError(f"Missing declared data column in {path}")
                iteration, time, value = row[0], row[1], row[column]
                values[int(iteration)] = (time, value)
    if not values:
        raise ValueError(f"Missing or empty scalar series for {variable}.{reduction}")
    return values


def normalized_rows(root, region):
    names = [f"em_ce_{region}", f"em_cb_{region}", f"em_mask_{region}"]
    series = [read_series(root, name, "norm2") for name in names[:2]]
    series += [read_series(root, names[2], "average")]
    series += [read_series(root, name, "norm_inf") for name in names[:2]]
    keys = set(series[0])
    if any(set(item) != keys for item in series[1:]):
        raise ValueError("Residual/mask time coverage differs; restrict to matched outputs first")
    for iteration in sorted(keys):
        samples = [item[iteration] for item in series]
        time = samples[0][0]
        if any(not math.isclose(t, time, rel_tol=1e-12, abs_tol=1e-10) for t, _ in samples):
            raise ValueError(f"Timestamp mismatch at iteration {iteration}")
        ce, cb, fraction, ce_max, cb_max = [v for _, v in samples]
        if not all(math.isfinite(v) for v in (time, ce, cb, fraction, ce_max, cb_max)):
            raise ValueError(f"Nonfinite value at iteration {iteration}")
        if not 0 < fraction <= 1 + 1e-12:
            raise ValueError(f"Empty/invalid monitoring mask at iteration {iteration}")
        if min(ce, cb, ce_max, cb_max) < 0:
            raise ValueError(f"Negative norm at iteration {iteration}")
        yield iteration, time, ce / math.sqrt(fraction), cb / math.sqrt(fraction), ce_max, cb_max, fraction


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", help="One simulation directory, including output-NNNN subdirectories")
    parser.add_argument("--region", choices=("near", "far"), required=True)
    parser.add_argument("--output", required=True, help="CSV destination; an existing file is not overwritten")
    args = parser.parse_args()
    rows = list(normalized_rows(args.root, args.region))
    with Path(args.output).open("x", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("iteration", "time", "CE_L2_RMS", "CB_L2_RMS", "CE_max_abs", "CB_max_abs", "mask_average"))
        writer.writerows(rows)


if __name__ == "__main__":
    main()
