"""Systematic check of the surface renderer: every region, every view.

Outlines all 70 DK regions -- not just the handful the dummy data calls
significant -- and asserts that what gets drawn is what should be. Thresholds
are set from the measured distribution, not guessed. Run:

    python examples/dk_surface_check.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dk_surface_render as R          # noqa: E402

ALL = set(R.REGION.values())
report, panels = {}, {}                # (region, view) -> (front, drawn, arcs)
for hemi, view, azim in R.PANELS:
    rep = {}
    panels[hemi, view] = R.render(hemi, azim, {}, ALL, report=rep)
    for r, v in rep.items():
        report[r, f"{hemi[0]}h_{view}"] = v

drawn = {k: v for k, v in report.items() if v[1]}
best = {}
for (r, p), (f, d, a) in report.items():
    if f > best.get(r, (0,))[0]:
        best[r] = (f, d, p)

# 1. Every region is outlined somewhere. A parcel that never appears is a
#    label-mapping or filtering bug -- this is what caught temporalpole_left
#    being deleted by an absolute speckle threshold.
seen = {r for (r, _) in drawn}
assert seen == ALL, f"never drawn: {sorted(ALL - seen)}"

# 2. In the view where a region mostly lives, the depth test must keep nearly
#    all of its front-facing contour. Measured range is 0.94-1.00; at
#    INFLATE=80 insula_right sat at 0.48, with half its outline behind the
#    operculum, which is the failure this guards.
ratio = {r: d / f for r, (f, d, _) in best.items() if f}
low = sorted(ratio.items(), key=lambda kv: kv[1])[:5]
assert min(ratio.values()) >= 0.90, f"outline eaten in its own view: {low}"

# 3. No half-drawn views. A region either shows in a view or it doesn't; a
#    handful of surviving faces means far-side leak -- the dots in the middle
#    of nowhere. Measured: 0 views fall between.
partial = {k: v for k, v in report.items() if v[0] and 0 < v[1] < 0.5 * v[0]}
assert not partial, f"partly drawn: {sorted(partial)[:5]}"

# 4. Continuity: a drawn outline is one arc, or a few where it wraps past the
#    silhouette -- never confetti. Measured: 85 single-arc, 8 two, 1 three.
assert max(v[2] for v in drawn.values()) <= 3, \
    f"fragmented: {sorted(((v[2], k) for k, v in drawn.items()), reverse=True)[:3]}"

# 5. Nothing drawn that faces away.
assert all(d <= f for f, d, _ in report.values())

# 6. DK is bilateral: a region drawn on the left must be drawn on the right.
side = lambda p: p[:2]                                          # noqa: E731
stem = lambda r: r.rsplit("_", 1)[0]                            # noqa: E731
lh = {stem(r) for (r, p) in drawn if side(p) == "lh"}
rh = {stem(r) for (r, p) in drawn if side(p) == "rh"}
assert lh == rh, f"asymmetric: {sorted(lh ^ rh)}"

# 7. The mesh itself is sane in every panel.
for key, (tri, face, seg) in panels.items():
    assert len(tri) == len(face) > 5000, (key, len(tri))
    assert np.isfinite(seg).all() and np.isfinite(tri).all(), key

print(f"{len(ALL)} regions x {len(R.PANELS)} views on {R.SURFACE}, INFLATE={R.INFLATE}")
print(f"  {len(drawn)} region-views carry outline, "
      f"{sum(v[1] for v in drawn.values())} contour faces drawn")
print(f"  retention in own view: {min(ratio.values()):.2f}-{max(ratio.values()):.2f} "
      f"(worst {low[0][0]})")
print(f"  arcs per drawn view: "
      f"{ {n: sum(1 for v in drawn.values() if v[2] == n) for n in sorted({v[2] for v in drawn.values()})} }")
for hemi, view, _ in R.PANELS:
    got = sum(1 for (r, p), v in drawn.items() if p == f"{hemi[0]}h_{view}")
    print(f"  {hemi[0]}h_{view:8} {got:2} regions outlined")
print("ok")
