"""Run: python test_pyseg.py"""
import matplotlib
matplotlib.use("Agg")
from matplotlib.patches import PathPatch
import pyseg

r = pyseg.brain_regions("dk")
assert len(r) == 70 and "frontalpole_left" in r, len(r)   # 35 per hemisphere
assert not any(n.startswith("_") for n in r)              # silhouettes aren't regions
assert pyseg.brain_views("dk") == [
    ("left_lateral", "left", "lateral"), ("left_medial", "left", "medial"),
    ("right_medial", "right", "medial"), ("right_lateral", "right", "lateral")]
assert pyseg.brain_views("aseg")[0] == ("coronal", "", "coronal")   # both hemis, one slice

data = {n: i for i, n in enumerate(r[:20])}
sig = r[:3]
fig, ax = pyseg.plot_brain(data, atlas="dk", sig=sig)

patches = [p for p in ax.patches if isinstance(p, PathPatch)]
fills = [p for p in patches if p.get_facecolor()[3] > 0]
outlines = [p for p in patches if p.get_facecolor()[3] == 0]
sig_patches = [p for p in outlines if p.get_linewidth() == 2.0]

ps = pyseg._atlas("dk")[0]
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
c = pyseg._canon
assert len({c(n) for n in ["lh_bankssts", "bankssts_left", "L-bankssts",
                           "Left bankssts", "bankssts.lh"]}) == 1
assert c("rh_insula") == c("insula_right") != c("insula_left")
assert c("Left-Lateral-Ventricle") == c("lateral ventricle L")
assert c("Left-Thalamus") == c("Left-Thalamus-Proper")     # FreeSurfer 6 vs 7
assert c("Cerebellum-Cortex") != c("Left-Cerebellum-Cortex")  # L/R merge stays explicit

for a in ("dk", "aseg", "jhu", "glasser", "schaefer7_400"):
    names = pyseg.brain_regions(a)                # no two regions collapse onto each other
    assert len({c(n) for n in names}) == len(names), a

import warnings
with warnings.catch_warnings():
    warnings.simplefilter("error")   # would warn if the aliases failed
    pyseg.plot_brain({"lh_bankssts": 1.0}, sig=["rh_insula"])
    pyseg.plot_brain({"Left-Thalamus": 1.0, "Brain-Stem": 2.0}, atlas="aseg")

# Layout options.
def spread(**kw):
    f, a = pyseg.plot_brain(data, sig=sig, **kw)
    v = [p.get_path().get_extents() for p in a.patches]
    return (max(b.x1 for b in v) - min(b.x0 for b in v),
            max(b.y1 for b in v) - min(b.y0 for b in v))

wide, tall = spread(), spread(position="stacked")
assert wide[0] / wide[1] > 4 and 1 < tall[0] / tall[1] < 2, (wide, tall)
assert spread(hemisphere="left")[0] < wide[0] / 2
assert len(pyseg.plot_brain(data, view="lateral")[1].patches) < len(patches)
for bad in [dict(view="nope"), dict(position="sideways")]:
    try:
        pyseg.plot_brain(data, **bad); raise AssertionError(bad)
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
kept = pyseg._despeckle(_shape(body, hole, chip), tol=0.1)
assert len(pyseg._rings(kept)) == 2, "hole dropped, or fragment kept"
assert pyseg._despeckle(_shape(chip), tol=0.1, ref=10000) is None   # sliver of a big parcel
assert pyseg._despeckle(_shape(chip), tol=0.1) is not None          # a parcel of its own
assert pyseg._despeckle(_shape(_ring(0, 0, 0.01, 0.01)), tol=0.1) is None   # sub-pixel

for a, n in [("dk", 70), ("aseg", 26), ("jhu", 21), ("glasser", 359),
             ("schaefer7_400", 400)]:
    ps = pyseg._atlas(a)[0]
    assert len(pyseg.brain_regions(a)) == n, (a, len(pyseg.brain_regions(a)))
    for p in ps:                       # nothing left that is all sliver
        for k, path in p.paths.items():
            big = max(pyseg._area(v) for q in ps if k in q.paths
                      for v, _ in pyseg._rings(q.paths[k]))
            mine = max(pyseg._area(v) for v, _ in pyseg._rings(path))
            assert mine >= 0.02 * big or mine == big, (a, p.name, k)

# Corner cutting: coarse traces (aseg) look spiky once outlined.
def _perim(v):
    return float(np.hypot(*np.diff(v, axis=0).T).sum())

saw = [(0, 0)] + [(4 * i, 3 * (i % 2)) for i in range(1, 25)] + [(100, 100), (0, 100), (0, 0)]
zig = MPath(np.array(saw, float),
            [MPath.MOVETO] + [MPath.LINETO] * (len(saw) - 2) + [MPath.CLOSEPOLY])
sm = pyseg._smooth(zig, 2)
assert pyseg._smooth(zig, 0) is zig
assert _perim(sm.vertices) < _perim(zig.vertices)          # the sawtooth rounds off
assert 0.97 < pyseg._area(sm.vertices) / pyseg._area(zig.vertices) <= 1.0
assert sm.vertices.min() >= -1e-9 and sm.vertices.max() <= 100 + 1e-9   # never grows

# Cutting a distance set by each ring, not a fraction of each edge: sparsely
# traced parcels would otherwise shrink badly (glasser lost a third of its area).
for a in ("dk", "aseg", "glasser", "schaefer7_400"):
    kept = [pyseg._area(q.paths[k].vertices) / pyseg._area(p.paths[k].vertices)
            for p, q in zip(pyseg._atlas(a, 0)[0], pyseg._atlas(a, 2)[0])
            for k in p.paths]
    assert min(kept) > 0.97, (a, min(kept))

raw = next(p for p in pyseg._atlas("aseg", 0)[0] if "brainstem" in p.paths)
soft = next(p for p in pyseg._atlas("aseg", 2)[0] if "brainstem" in p.paths)
assert len(soft.paths["brainstem"].vertices) > len(raw.paths["brainstem"].vertices)
curved = [pl for p in pyseg._atlas("jhu")[0] for pl in p.paths.values()
          if MPath.CURVE4 in pl.codes]
assert curved and all(MPath.CURVE4 in pl.codes for pl in curved)   # left alone

# Panels are a partition. ggseg's atlases were digitised region by region, so
# neighbours overlapped -- 2.3% of dk's left lateral panel, drawn as doubled,
# uneven borders. tools/export_ggseg_atlas.R cuts each region out of what is
# already placed; this is what stops that regressing.
def overlap(atlas, panel, n=300):
    p = next(q for q in pyseg._atlas(atlas, 0)[0] if q.name == panel)
    v = np.vstack([x.vertices for x in p.paths.values()])
    (x0, y0), (x1, y1) = v.min(0), v.max(0)
    gx, gy = np.meshgrid(np.linspace(x0, x1, n), np.linspace(y0, y1, n))
    pts = np.column_stack([gx.ravel(), gy.ravel()])
    hits = sum(q.contains_points(pts).astype(int) for q in p.paths.values())
    return (hits >= 2).sum() / max((hits >= 1).sum(), 1)

for a in ("dk", "aseg", "glasser", "schaefer7_400"):     # jhu is hand-vendored,
    for pan, _, _ in pyseg.brain_views(a):               # not exported, not a partition
        f = overlap(a, pan)
        assert f < 0.002, f"{a}/{pan} regions overlap over {f:.2%} of the panel"

# Every shipped atlas loads, has regions and panels, and plots.
cat = pyseg.brain_atlases()
assert len(cat) >= 20 and cat["dk"] == 70 and cat["schaefer7_400"] == 400, cat
import matplotlib.pyplot as plt
for a, n in cat.items():
    names, views = pyseg.brain_regions(a), pyseg.brain_views(a)
    assert len(names) == n > 0 and views, a
    assert len({c(x) for x in names}) == len(names), f"{a}: names collide"
# one of each shape rather than all 23: two hemispheres, slices, one panel, dense
for a in ("dk", "aseg", "jhu", "glasser", "schaefer17_1000"):
    names = pyseg.brain_regions(a)
    fig, _ = pyseg.plot_brain({names[0]: 1.0}, atlas=a, sig=[names[0]])
    plt.close(fig)

# dim fades everything that is not significant, so the result carries the figure
fig, dax = pyseg.plot_brain(data, sig=sig, dim=0.3)
alpha = {round(float(p.get_facecolor()[3]), 2) for p in dax.patches
         if p.get_facecolor()[3] > 0}
assert alpha == {0.3, 1.0}, alpha              # faded fills and solid ones, nothing else
plt.close(fig)
fig, uax = pyseg.plot_brain(data, sig=sig)     # unchanged without it
assert {round(float(p.get_facecolor()[3]), 2) for p in uax.patches
        if p.get_facecolor()[3] > 0} == {1.0}
plt.close(fig)
assert len(pyseg.geom_brain(data, sig=sig, dim=0.3).layers) == 2

# Subcortical surfaces. Structures are separate closed meshes, so a structure
# is marked out by its silhouette, not by a boundary shared with a neighbour.
from matplotlib.collections import LineCollection
sr = pyseg.subcortical_regions()
assert len(sr) == 16 and "hippocampus_left" in sr, sr
assert len({c(x) for x in sr}) == 16                    # no two names collapse
assert c("Left-Hippocampus") == c("hippocampus_left")   # FreeSurfer names land
assert c("Left-Thalamus") == c("thalamusproper_left")

vals = {r: i - 8 for i, r in enumerate(sr)}
fig, sax = pyseg.plot_subcortical(vals, sig=["hippocampus_left"])
lines = [x for x in sax.collections if isinstance(x, LineCollection)]
assert lines and sum(len(x.get_segments()) for x in lines) > 20, "no outline drawn"
plt.close(fig)
fig, fade = pyseg.plot_subcortical(vals, sig=["hippocampus_left"], dim=0.2)
fc = np.vstack([x.get_facecolor() for x in fade.collections
                if not isinstance(x, LineCollection)])
assert set(np.round(np.unique(fc[:, 3]), 2)) == {0.2, 1.0}, np.unique(fc[:, 3])
plt.close(fig)
fig, bare = pyseg.plot_subcortical(vals)                # nothing significant
assert not [x for x in bare.collections if isinstance(x, LineCollection)]
plt.close(fig)
for kw in (dict(hemisphere="left"), dict(view="medial"), dict(position="stacked")):
    fig, _ = pyseg.plot_subcortical(vals, sig=sr[:2], **kw)
    plt.close(fig)
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    pyseg.plot_subcortical({"not-a-structure": 1.0})
    assert w and "not-a-structure" in str(w[0].message)
plt.close("all")

# ggplot layer
df = pyseg.as_brain_df("dk")
assert set(df.columns) == {"region", "hemi", "side", "panel", "ring", "x", "y"}
assert df.region.isna().any() and df.hemi.isin(["left", "right"]).all()
p = pyseg.geom_brain(data, sig=sig, position="stacked")
assert len(p.layers) == 2 and p.facet.__class__.__name__ == "facet_grid"
assert len(pyseg.geom_brain(data).layers) == 1   # no sig, no outline layer

# Connectome: one arc per suprathreshold edge, nodes at region centroids.
import pandas as pd
from matplotlib.patches import FancyArrowPatch
left = [x for x in pyseg.brain_regions("dk") if x.endswith("_left")][:6]
w = np.arange(36, dtype=float).reshape(6, 6); w = (w + w.T) / 2
np.fill_diagonal(w, 0)
edges = pd.DataFrame(w, index=left, columns=left)
_, cax = pyseg.plot_connectome(edges, hemisphere="left", threshold=40)
arcs = [p for p in cax.patches if isinstance(p, FancyArrowPatch)]
assert len(arcs) == int((np.triu(np.abs(w), 1) > 40).sum()), len(arcs)
assert len(arcs) < 15, "threshold ignored"
assert pyseg.brain_atlases()["dk"] == 70
assert "lh_bankssts" in pyseg.brain_labels("dk")

print("ok")
