"""Run: python test_brainseg.py"""
import matplotlib
matplotlib.use("Agg")
from matplotlib.patches import PathPatch
import brainseg

r = brainseg.brain_regions("dk")
assert len(r) == 70 and "frontalpole_left" in r, len(r)   # 35 per hemisphere
assert not any(n.startswith("_") for n in r)              # silhouettes aren't regions
assert brainseg.brain_views("dk") == [
    ("left_lateral", "left", "lateral"), ("left_medial", "left", "medial"),
    ("right_medial", "right", "medial"), ("right_lateral", "right", "lateral")]
assert brainseg.brain_views("aseg")[0] == ("coronal", "", "coronal")   # both hemis, one slice

data = {n: i for i, n in enumerate(r[:20])}
sig = r[:3]
fig, ax = brainseg.plot_brain(data, atlas="dk", sig=sig)

patches = [p for p in ax.patches if isinstance(p, PathPatch)]
fills = [p for p in patches if p.get_facecolor()[3] > 0]
outlines = [p for p in patches if p.get_facecolor()[3] == 0]
sig_patches = [p for p in outlines if p.get_linewidth() == 2.0]

ps = brainseg._atlas("dk")[0]
n_shapes = sum(len(p.paths) for p in ps)        # a region can appear in 2 views
n_walls = sum(len(p.wall) for p in ps)

# The whole point: every significant outline sits above every fill and outline.
assert len(sig) <= len(sig_patches) <= 2 * len(sig), sig_patches
assert min(p.get_zorder() for p in sig_patches) > max(p.get_zorder() for p in fills)
assert min(p.get_zorder() for p in sig_patches) > max(
    p.get_zorder() for p in outlines if p not in sig_patches)
assert len(fills) == n_shapes + n_walls
assert len(outlines) == n_shapes + len(sig_patches)

# Hemisphere prefix/suffix and separators don't matter.
c = brainseg._canon
assert len({c(n) for n in ["lh_bankssts", "bankssts_left", "L-bankssts",
                           "Left bankssts", "bankssts.lh"]}) == 1
assert c("rh_insula") == c("insula_right") != c("insula_left")
assert c("Left-Lateral-Ventricle") == c("lateral ventricle L")
assert c("Left-Thalamus") == c("Left-Thalamus-Proper")     # FreeSurfer 6 vs 7
assert c("Cerebellum-Cortex") != c("Left-Cerebellum-Cortex")  # L/R merge stays explicit

for a in ("dk", "aseg", "jhu", "glasser", "schaefer7_400"):
    names = brainseg.brain_regions(a)                # no two regions collapse onto each other
    assert len({c(n) for n in names}) == len(names), a
    brainseg.plot_brain({names[0]: 1.0}, atlas=a, sig=[names[0]])

import warnings
with warnings.catch_warnings():
    warnings.simplefilter("error")   # would warn if the aliases failed
    brainseg.plot_brain({"lh_bankssts": 1.0}, sig=["rh_insula"])
    brainseg.plot_brain({"Left-Thalamus": 1.0, "Brain-Stem": 2.0}, atlas="aseg")

# Layout options.
def spread(**kw):
    f, a = brainseg.plot_brain(data, sig=sig, **kw)
    v = [p.get_path().get_extents() for p in a.patches]
    return (max(b.x1 for b in v) - min(b.x0 for b in v),
            max(b.y1 for b in v) - min(b.y0 for b in v))

wide, tall = spread(), spread(position="stacked")
assert wide[0] / wide[1] > 4 and 1 < tall[0] / tall[1] < 2, (wide, tall)
assert spread(hemisphere="left")[0] < wide[0] / 2
assert len(brainseg.plot_brain(data, view="lateral")[1].patches) < len(patches)
for bad in [dict(view="nope"), dict(position="sideways")]:
    try:
        brainseg.plot_brain(data, **bad); raise AssertionError(bad)
    except ValueError:
        pass

# Stray-geometry cleanup: ggseg sheds sub-pixel rings and slivers of a parcel in
# views it barely reaches. Outlined, each one is a lone speck or dash.
import numpy as np
from matplotlib.path import Path as MPath

def _ring(x, y, w, h):
    v = [(x, y), (x + w, y), (x + w, y + h), (x, y + h), (x, y)]
    return v, [MPath.MOVETO] + [MPath.LINETO] * 3 + [MPath.CLOSEPOLY]

def _shape(*rings):
    v, c = zip(*rings)
    return MPath(np.concatenate([np.array(a) for a in v]), np.concatenate(c))

body = _ring(0, 0, 100, 100)
hole = _ring(40, 40, 5, 5)          # small, but inside the body
chip = _ring(400, 0, 5, 5)          # same size, off on its own
kept = brainseg._despeckle(_shape(body, hole, chip), tol=0.1)
assert len(brainseg._rings(kept)) == 2, "hole dropped, or fragment kept"
assert brainseg._despeckle(_shape(chip), tol=0.1, ref=10000) is None   # sliver of a big parcel
assert brainseg._despeckle(_shape(chip), tol=0.1) is not None          # a parcel of its own
assert brainseg._despeckle(_shape(_ring(0, 0, 0.01, 0.01)), tol=0.1) is None   # sub-pixel

for a, n in [("dk", 70), ("aseg", 26), ("jhu", 21), ("glasser", 359),
             ("schaefer7_400", 400)]:
    ps = brainseg._atlas(a)[0]
    assert len(brainseg.brain_regions(a)) == n, (a, len(brainseg.brain_regions(a)))
    for p in ps:                       # nothing left that is all sliver
        for k, path in p.paths.items():
            big = max(brainseg._area(v) for q in ps if k in q.paths
                      for v, _ in brainseg._rings(q.paths[k]))
            mine = max(brainseg._area(v) for v, _ in brainseg._rings(path))
            assert mine >= 0.02 * big or mine == big, (a, p.name, k)

# Corner cutting: coarse traces (aseg) look spiky once outlined.
def _perim(v):
    return float(np.hypot(*np.diff(v, axis=0).T).sum())

saw = [(0, 0)] + [(4 * i, 3 * (i % 2)) for i in range(1, 25)] + [(100, 100), (0, 100), (0, 0)]
zig = MPath(np.array(saw, float),
            [MPath.MOVETO] + [MPath.LINETO] * (len(saw) - 2) + [MPath.CLOSEPOLY])
sm = brainseg._smooth(zig, 2)
assert brainseg._smooth(zig, 0) is zig
assert _perim(sm.vertices) < _perim(zig.vertices)          # the sawtooth rounds off
assert 0.97 < brainseg._area(sm.vertices) / brainseg._area(zig.vertices) <= 1.0
assert sm.vertices.min() >= -1e-9 and sm.vertices.max() <= 100 + 1e-9   # never grows

# Cutting a distance set by each ring, not a fraction of each edge: sparsely
# traced parcels would otherwise shrink badly (glasser lost a third of its area).
for a in ("dk", "aseg", "glasser", "schaefer7_400"):
    kept = [brainseg._area(q.paths[k].vertices) / brainseg._area(p.paths[k].vertices)
            for p, q in zip(brainseg._atlas(a, 0)[0], brainseg._atlas(a, 2)[0])
            for k in p.paths]
    assert min(kept) > 0.97, (a, min(kept))

raw = next(p for p in brainseg._atlas("aseg", 0)[0] if "brainstem" in p.paths)
soft = next(p for p in brainseg._atlas("aseg", 2)[0] if "brainstem" in p.paths)
assert len(soft.paths["brainstem"].vertices) > len(raw.paths["brainstem"].vertices)
curved = [pl for p in brainseg._atlas("jhu")[0] for pl in p.paths.values()
          if MPath.CURVE4 in pl.codes]
assert curved and all(MPath.CURVE4 in pl.codes for pl in curved)   # left alone

# ggplot layer
df = brainseg.as_brain_df("dk")
assert set(df.columns) == {"region", "hemi", "side", "panel", "ring", "x", "y"}
assert df.region.isna().any() and df.hemi.isin(["left", "right"]).all()
p = brainseg.geom_brain(data, sig=sig, position="stacked")
assert len(p.layers) == 2 and p.facet.__class__.__name__ == "facet_grid"
assert len(brainseg.geom_brain(data).layers) == 1   # no sig, no outline layer

# Connectome: one arc per suprathreshold edge, nodes at region centroids.
import pandas as pd
from matplotlib.patches import FancyArrowPatch
left = [x for x in brainseg.brain_regions("dk") if x.endswith("_left")][:6]
w = np.arange(36, dtype=float).reshape(6, 6); w = (w + w.T) / 2
np.fill_diagonal(w, 0)
edges = pd.DataFrame(w, index=left, columns=left)
_, cax = brainseg.plot_connectome(edges, hemisphere="left", threshold=40)
arcs = [p for p in cax.patches if isinstance(p, FancyArrowPatch)]
assert len(arcs) == int((np.triu(np.abs(w), 1) > 40).sum()), len(arcs)
assert len(arcs) < 15, "threshold ignored"
assert brainseg.brain_atlases()["dk"] == 70
assert "lh_bankssts" in brainseg.brain_labels("dk")

print("ok")
