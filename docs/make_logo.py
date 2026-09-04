"""Regenerate docs/logo.png from the atlas itself."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import brainseg

hi = ["superiortemporal_left", "supramarginal_left", "parsopercularis_left"]
data = {r: i for i, r in enumerate(hi)}
fig, ax = brainseg.plot_brain(data, hemisphere="left", view="lateral", sig=hi,
                              cmap="Spectral_r", vmin=-2, vmax=4, sig_lw=6,
                              na_color="0.92", edgecolor="w", lw=2,
                              colorbar=False, figsize=(6, 4))
fig.patch.set_alpha(0)
fig.savefig("docs/logo.png", dpi=110, bbox_inches="tight", transparent=True)
