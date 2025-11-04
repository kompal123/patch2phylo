#!/usr/bin/env python
"""
Compute per-position and sliding-window Shannon entropy from a multi-FASTA alignment.
Writes two CSV outputs (raw and windowed) and a PNG plot.
"""

import os
import math
from collections import Counter
from Bio import AlignIO

import matplotlib
matplotlib.use("Agg")  # headless backend to avoid Qt
import matplotlib.pyplot as plt

# ---- Snakemake I/O ----
msa_path   = snakemake.input["msa"]
raw_csv    = snakemake.output["txt"]
win_csv    = snakemake.output["windowed"]
plot_png   = snakemake.output["plot"]

# Params from config (with sensible defaults)
win  = int(snakemake.config.get("parameters", {}).get("variability", {}).get("window_size", 100))
step = int(snakemake.config.get("parameters", {}).get("variability", {}).get("step_size", win))

# Ensure output dirs exist
for p in (raw_csv, win_csv, plot_png):
    os.makedirs(os.path.dirname(p), exist_ok=True)

def shannon_entropy(col_iterable):
    """Shannon entropy over non-gap characters in a column."""
    counts = Counter(b for b in col_iterable if b not in (b"-", b"."))
    n = sum(counts.values())
    if n == 0:
        return 0.0
    return -sum((c/n) * math.log2(c/n) for c in counts.values())

# Load alignment
aln = AlignIO.read(msa_path, "fasta")
L = aln.get_alignment_length()

# Per-position entropy
raw_ent = []
for i in range(L):
    # aln[:, i] returns a str of column residues; convert to bytes to treat uniformly
    col = aln[:, i]
    col_bytes = col.encode() if isinstance(col, str) else col
    raw_ent.append(shannon_entropy(col_bytes))

# Write raw CSV
with open(raw_csv, "w") as f:
    f.write("position,entropy\n")
    for pos, e in enumerate(raw_ent, 1):
        f.write(f"{pos},{e:.6f}\n")

# Sliding-window means
starts, means = [], []
if L > 0:
    idx = 0
    while idx < L:
        block = raw_ent[idx: idx + win]
        if not block:
            break
        avg = sum(block) / len(block)
        starts.append(idx + 1)
        means.append(avg)
        idx += max(step, 1)

# Write windowed CSV
with open(win_csv, "w") as f:
    f.write("window_start,avg_entropy\n")
    for s, m in zip(starts, means):
        f.write(f"{s},{m:.6f}\n")

# Plot (windowed curve; add raw outline if not too long)
plt.figure(figsize=(12, 4))
if raw_ent and L <= 10000:
    plt.plot(range(1, L + 1), raw_ent, linewidth=0.8, alpha=0.5, label="per-position")
if means:
    plt.plot(starts, means, linewidth=1.6, marker="o", markersize=2.5, label=f"window={win}")
plt.xlabel("Alignment position")
plt.ylabel("Shannon entropy (bits)")
plt.title("Sliding-window Shannon entropy")
plt.legend(loc="upper right")
plt.tight_layout()
plt.savefig(plot_png, dpi=150)

