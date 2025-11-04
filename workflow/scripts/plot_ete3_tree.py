#!/usr/bin/env python3
# plot_ete3_tree.py — handle IQ-TREE "UFboot/SH-aLRT" labels like 100/100:0.01

import os, re, sys
from ete3 import Tree, TreeStyle, NodeStyle
from ete3.parser.newick import NewickError

# Snakemake I/O (CLI fallback)
try:
    in_tree  = snakemake.input["tree"]
    out_svg  = snakemake.output["svg"]
    out_pdf  = snakemake.output["pdf"]
    out_png  = snakemake.output["png"]
except NameError:
    if len(sys.argv) < 5:
        print("Usage: plot_ete3_tree.py <in.newick> <out.svg> <out.pdf> <out.png>", file=sys.stderr)
        sys.exit(2)
    in_tree, out_svg, out_pdf, out_png = sys.argv[1:5]

os.makedirs(os.path.dirname(out_svg), exist_ok=True)

with open(in_tree, "r") as f:
    nw = f.read().strip()

# Replace “)X/Y:” or “)X/Y,” or “)X/Y)” with “)X:” / “)X,” / “)X)”.
# Keep the first support (X) and drop the second (Y).
nw_clean = re.sub(
    r'\)([0-9eE+\-\.]+)\/([0-9eE+\-\.]+)(?=[:\),])',
    r')\1',
    nw
)

# Try to parse with support-aware format first
try:
    t = Tree(nw_clean, format=1)
except NewickError:
    # Fallback: IQ-TREE sometimes writes node labels slightly differently;
    # try the default parser as a last resort
    t = Tree(nw_clean)

# Style (compact, no node circles)
ts = TreeStyle()
ts.show_leaf_name = True
ts.scale = 120
ts.branch_vertical_margin = 12

for n in t.traverse():
    ns = NodeStyle()
    ns["size"] = 0
    n.set_style(ns)
    # Show support as integer if present
    if not n.is_leaf() and getattr(n, "support", None) is not None:
        try:
            n.name = f"{float(n.support):.0f}"
        except Exception:
            pass

t.render(out_svg, tree_style=ts)
t.render(out_pdf, tree_style=ts)
t.render(out_png, tree_style=ts, w=1200)


