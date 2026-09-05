"""The surface renderer's version of the gallery: layout x outline weight.

Same synthetic map as examples/gallery.py, so the two can be compared side by
side -- flat ggseg panels against an inflated surface.

Writes PNG and PDF for each combination into examples/gallery/.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dk_surface_render as R          # noqa: E402
from dk_example import synthetic       # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "gallery")
WIDTHS = [("thin", 0.7), ("default", 1.6), ("bold", 3.0)]

values, sig = synthetic("dk", seed=3)
os.makedirs(OUT, exist_ok=True)
print(f"{len(sig)}/{len(values)} regions outlined on {R.SURFACE}, "
      f"INFLATE={R.INFLATE}")

R.VLIM = 5.0
R.cmap, R.norm = R.mpl.colormaps["RdBu_r"], R.mpl.colors.Normalize(-5, 5)
for position in ("dispersed", "stacked"):
    panels = None                      # geometry is identical across widths
    for weight, lw in WIDTHS:
        stem = os.path.join(OUT, f"surface_{position}_{weight}")
        panels = R.main(position, values, sig, stem=stem, lw=lw,
                        label="t (synthetic)",
                        title=f"Synthetic t-map on {R.SURFACE} — "
                              f"{position}, lw={lw}", panels=panels)
        print(f"  {os.path.basename(stem):28} lw={lw:<4} "
              f"{os.path.getsize(stem + '.pdf') / 1024:5.0f} kB pdf")
