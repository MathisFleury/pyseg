"""The same map drawn every way pyseg offers: layout x outline weight.

Synthetic data, with enough regions outlined that the trade-off shows: a thin
outline disappears in a crowd, a thick one swallows the small parcels.

Writes PNG and PDF for each combination into examples/gallery/.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pyseg
from dk_example import synthetic

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "gallery")

LAYOUTS = [("dispersed", {}),                       # ggseg's default row of four
           ("stacked", dict(position="stacked")),   # views x hemispheres
           ("left", dict(hemisphere="left")),
           ("lateral", dict(view="lateral"))]
WIDTHS = [("thin", 0.8), ("default", 2.0), ("bold", 3.5)]
DIMS = [("solid", None), ("dim", 0.3), ("faded", 0.12)]

t, sig = synthetic("dk", seed=3)
os.makedirs(OUT, exist_ok=True)
print(f"{len(sig)}/{len(t)} regions outlined")

# how far to fade what is not significant, at one layout and weight
for weight, d in DIMS:
    fig, _ = pyseg.plot_brain(t, sig=sig, cmap="RdBu_r", vmin=-5, vmax=5,
                              ylabel="t (synthetic)", dim=d,
                              title=f"Synthetic t-map — dim={d}")
    stem = os.path.join(OUT, f"dk_dispersed_{weight}")
    for ext in ("png", "pdf"):
        fig.savefig(f"{stem}.{ext}", dpi=200, bbox_inches="tight")
    print(f"  {os.path.basename(stem):24} dim={d}")

for layout, kw in LAYOUTS:
    for weight, lw in WIDTHS:
        fig, _ = pyseg.plot_brain(
            t, sig=sig, cmap="RdBu_r", vmin=-5, vmax=5, ylabel="t (synthetic)",
            sig_lw=lw, lw=max(0.3, lw / 4),      # region borders scale with it
            title=f"Synthetic t-map — {layout}, sig_lw={lw}", **kw)
        stem = os.path.join(OUT, f"dk_{layout}_{weight}")
        for ext in ("png", "pdf"):
            fig.savefig(f"{stem}.{ext}", dpi=200, bbox_inches="tight")
        print(f"  {os.path.basename(stem):24} sig_lw={lw:<4} "
              f"{os.path.getsize(stem + '.pdf') / 1024:5.0f} kB pdf")
