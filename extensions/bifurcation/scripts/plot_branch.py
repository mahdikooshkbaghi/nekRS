#!/usr/bin/env python3
"""Plot accepted continuation rows; no CFD assumptions are encoded here."""
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt

if len(sys.argv) != 3:
    raise SystemExit(f"usage: {sys.argv[0]} branch.csv output.png")
rows = list(csv.DictReader(Path(sys.argv[1]).open()))
re = [float(r["Re"]) for r in rows]
sigma = [float(r["sigma"]) for r in rows]
st = [float(r["St"]) for r in rows]
fig, axes = plt.subplots(2, 1, sharex=True, constrained_layout=True)
axes[0].plot(re, sigma, "o-"); axes[0].axhline(0, color="k", lw=.7); axes[0].set_ylabel("growth rate")
axes[1].plot(re, st, "o-"); axes[1].set_xlabel("Re"); axes[1].set_ylabel("St")
fig.savefig(sys.argv[2], dpi=150)
