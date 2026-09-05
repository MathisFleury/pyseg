"""pyseg on synthetic data: a t-map with a handful of "significant" regions.

Values are random, not a result -- the point is the figure. Swap `synthetic()`
for your own dict or pandas Series keyed by region name; FreeSurfer labels
(lh_bankssts) work as-is.
"""
import numpy as np
import pyseg

CUT = 2.3          # |t| above this counts as significant, for the outline


def synthetic(atlas="dk", seed=0):
    """-> (values, significant). Reproducible nonsense in the shape of a t-map."""
    regions = pyseg.brain_regions(atlas)
    t = np.random.default_rng(seed).normal(0, 1.5, len(regions))
    values = dict(zip(regions, t))
    return values, [r for r, v in values.items() if abs(v) > CUT]


def main():
    for atlas, title in [("dk", "Cortical t-map"), ("aseg", "Subcortical t-map")]:
        values, sig = synthetic(atlas)
        fig, _ = pyseg.plot_brain(
            values, atlas=atlas, sig=sig, cmap="RdBu_r", vmin=-4, vmax=4,
            ylabel="t (synthetic)",
            title=f"{title} — black outline: |t| > {CUT}")
        for ext in ("pdf", "png"):
            fig.savefig(f"examples/{atlas}_example.{ext}", dpi=200, bbox_inches="tight")
        print(f"{atlas}: {len(sig)}/{len(values)} outlined "
              f"-> examples/{atlas}_example.pdf")

    # The same structures as surfaces rather than slices.
    names = pyseg.subcortical_regions()
    t = np.random.default_rng(1).normal(0, 1.5, len(names))
    values = dict(zip(names, t))
    sig = [r for r, v in values.items() if abs(v) > 1.4]
    fig, _ = pyseg.plot_subcortical(
        values, sig=sig, cmap="RdBu_r", vmin=-4, vmax=4, ylabel="t (synthetic)",
        title=f"Subcortical surfaces — black outline: |t| > 1.4")
    for ext in ("pdf", "png"):
        fig.savefig(f"examples/sctx_example.{ext}", dpi=200, bbox_inches="tight")
    print(f"sctx: {len(sig)}/{len(values)} outlined -> examples/sctx_example.pdf")


if __name__ == "__main__":
    main()
