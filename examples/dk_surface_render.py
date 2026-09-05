"""DK on a real cortical surface, laid out the way pyseg lays out its panels.

pyseg's panels are traced 2D polygons. This renders the same four views from a
real mesh instead, to see whether the layout machinery carries over. It does: a
view is just a camera angle, and panels are still tiled by their own bounding
boxes.

No 3D toolkit and no new dependency -- project the vertices, cull back faces,
sort the rest back to front, and hand matplotlib one PolyCollection per view.
Significant regions keep pyseg's trick: their outlines are a separate pass on
top, so a neighbouring patch cannot paint over them.

Surfaces and the DK label vector come from a local ENIGMA toolbox checkout.
Run examples/dk_surface_check.py to verify every region actually comes out.
"""
import sys
from functools import lru_cache

import matplotlib as mpl
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection
from scipy import sparse
from scipy.sparse.csgraph import connected_components

ENIGMA = "/Users/mfleury/POSTDOC/LIBRAIRY/ENIGMA/enigmatoolbox/datasets/"
SURFACE = "conte69"             # conte69 = 32k vertices a side, fsa5 = 10k
INFLATE = 200                   # smoothing rounds; 0 renders the folded surface,
                                # but then half of every outline hides in a sulcus.
                                # 80 is not enough either: the sylvian fissure stays
                                # shut and the insula loses half its outline behind
                                # the operculum. By 200 every parcel is fully out.
RELAX = 4                       # rounds of smoothing on the outline itself
DEPTH_TOL = 12.0                # mm; see the depth test below
MIN_RUN = 8                     # shorter visible runs are strays, not outline
OUTLINE_LW = 1.6                # significance outline; pyseg's flat plot_brain
SEAM_LW = 0.25                  # takes the same as sig_lw / lw
CMAP, VLIM = "RdBu_r", 4.0

MESH = {"fsa5": ("surfaces/fsa5_{h}h.surf.gii", "parcellations/aparc_fsa5.csv", 10242),
        "conte69": ("surfaces/conte69_32k_{h}h.gii", "parcellations/aparc_conte69.csv", 32492)}

# FreeSurfer's aparc order: labels 1-35 are the left hemisphere, 36-70 the right.
LUT = """bankssts caudalanteriorcingulate caudalmiddlefrontal corpuscallosum cuneus
entorhinal fusiform inferiorparietal inferiortemporal isthmuscingulate lateraloccipital
lateralorbitofrontal lingual medialorbitofrontal middletemporal parahippocampal paracentral
parsopercularis parsorbitalis parstriangularis pericalcarine postcentral posteriorcingulate
precentral precuneus rostralanteriorcingulate rostralmiddlefrontal superiorfrontal
superiorparietal superiortemporal supramarginal frontalpole temporalpole transversetemporal
insula""".split()
REGION = {i + 1: f"{n}_left" for i, n in enumerate(LUT)}
REGION.update({i + 36: f"{n}_right" for i, n in enumerate(LUT)})
CODE = {v: k for k, v in REGION.items()}

# hemisphere, view, azimuth. Lateral means looking at the hemisphere from outside.
PANELS = [("left", "lateral", 180), ("left", "medial", 0),
          ("right", "medial", 180), ("right", "lateral", 0)]

_lab = np.loadtxt(ENIGMA + MESH[SURFACE][1], dtype=int)
_n = MESH[SURFACE][2]
labels = {"left": _lab[:_n], "right": _lab[_n:]}


def camera(azim, elev=0.0):
    """-> (2x3 screen basis, forward vector pointing at the camera)."""
    a, e = np.radians(azim), np.radians(elev)
    fwd = np.array([np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), np.sin(e)])
    right = np.array([-np.sin(a), np.cos(a), 0.0])
    return np.stack([right, np.cross(fwd, right)]), fwd


def inflate(v, f, rounds, lam=0.6):
    """Laplacian smoothing, rescaled back to size -- a poor man's ?h.inflated.
    Region boundaries stay on the outside where you can see them."""
    if not rounds:
        return v
    e = np.vstack([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]])
    e = np.vstack([e, e[:, ::-1]])
    A = sparse.coo_matrix((np.ones(len(e)), (e[:, 0], e[:, 1])),
                          shape=(len(v), len(v))).tocsr()
    deg = np.asarray(A.sum(1))
    v, r0 = v.copy(), np.linalg.norm(v - v.mean(0), axis=1).mean()
    for _ in range(rounds):
        v += lam * (A @ v / deg - v)
    c = v.mean(0)
    return c + (v - c) * (r0 / np.linalg.norm(v - c, axis=1).mean())


@lru_cache(maxsize=None)      # inflation is the slow part; views reuse it
def surface(hemi):
    v, f = nib.load(ENIGMA + MESH[SURFACE][0].format(h=hemi[0])).agg_data()
    f = f.astype(int)
    v = inflate(v.astype(float), f, INFLATE)
    n = np.cross(v[f[:, 1]] - v[f[:, 0]], v[f[:, 2]] - v[f[:, 0]])
    n /= np.linalg.norm(n, axis=1, keepdims=True) + 1e-12
    if (n * (v[f].mean(1) - v.mean(0))).sum() < 0:      # normals point outwards
        n, f = -n, f[:, ::-1]
    return v, f, n


def contour(uv, f, inside):
    """Marching triangles: where the boundary of `inside` crosses a face, the
    segment joining the midpoints of the two edges it cuts. Neighbouring faces
    share those midpoints, so the result is a continuous line lying exactly on
    the label boundary -- unlike walking mesh edges, which zig-zags and depends
    on how each face was assigned to a region.

    The line still steps from one edge midpoint to the next, which reads as a
    staircase at mesh resolution, so relax it along its own graph afterwards.

    -> (segments (m,2,2), the face each belongs to, the endpoint ids (m,2)).
    """
    m = inside[f]
    k = m.sum(1)
    hit = np.flatnonzero((k == 1) | (k == 2))
    if not len(hit):
        return np.zeros((0, 2, 2)), hit, np.zeros((0, 2), int)
    mm, ff = m[hit], f[hit]
    odd = np.where(mm.sum(1) == 1, mm.argmax(1), (~mm).argmax(1))   # the lone vertex
    r = np.arange(len(ff))
    i = [ff[r, odd], ff[r, (odd + 1) % 3], ff[r, (odd + 2) % 3]]
    seg = np.stack([(uv[i[0]] + uv[i[1]]) / 2, (uv[i[0]] + uv[i[2]]) / 2], 1)

    # Name each contour point by the mesh edge it cuts, so neighbouring faces
    # agree on it and the segments form one graph.
    n = len(uv)
    key = np.stack([np.minimum(i[0], j).astype(np.int64) * n + np.maximum(i[0], j)
                    for j in i[1:]], 1)
    _, ids = np.unique(key, return_inverse=True)
    ids = ids.reshape(key.shape)

    pos = np.zeros((ids.max() + 1, 2))
    pos[ids.reshape(-1)] = seg.reshape(-1, 2)
    a, b = ids[:, 0], ids[:, 1]
    A = sparse.coo_matrix((np.ones(2 * len(a)), (np.r_[a, b], np.r_[b, a])),
                          shape=(len(pos),) * 2).tocsr()

    # Vertex label maps speckle: a handful of stray vertices inside another
    # parcel make their own little closed contour. Same nuisance pyseg drops
    # from the flat atlas, and the same rule applies -- judge a loop against the
    # region's own biggest loop, never against the hemisphere. A fixed cutoff
    # deletes small parcels outright: temporalpole_left is 95 vertices and sits
    # just under 6% of the hemisphere, speckle by any absolute measure.
    _, comp = connected_components(A, directed=False)
    lab_seg = comp[ids[:, 0]]
    c, cnt = np.unique(lab_seg, return_counts=True)
    top = cnt.max()
    keep = np.isin(lab_seg, c[(cnt >= max(12, 0.15 * top)) | (cnt == top)])

    deg = np.maximum(np.asarray(A.sum(1)), 1)
    for _ in range(RELAX):
        pos += 0.5 * (A @ pos / deg - pos)
    return pos[ids][keep], hit[keep], ids[keep]


def render(hemi, azim, values, sig, report=None):
    """One view -> (triangles, face colours, outline segments).

    `report`, if given, collects per-region (front-facing, drawn) contour-face
    counts -- what examples/dk_surface_check.py asserts on.
    """
    v, f, n = surface(hemi)
    P, fwd = camera(azim)
    uv, dep, lab = v @ P.T, v @ fwd, labels[hemi]

    front = (n @ fwd) > 0
    tri = uv[f][front]
    order = np.argsort(dep[f].mean(1)[front])                # painter's algorithm

    # Colour a face by the label at least two of its vertices agree on, so the
    # fills line up with the contours below.
    lf = lab[f]
    maj = np.where(lf[:, 0] == lf[:, 1], lf[:, 0],
                   np.where(lf[:, 0] == lf[:, 2], lf[:, 0], lf[:, 1]))
    light = fwd + 0.45 * P[1] - 0.25 * P[0]
    shade = 0.55 + 0.45 * np.clip(n[front] @ (light / np.linalg.norm(light)), 0, 1)
    face = np.array([cmap(norm(values[r])) if r in values else (.85, .85, .85, 1.)
                     for r in (REGION.get(i, "") for i in maj[front])])
    face[:, :3] *= shade[:, None]

    # An outline behind the surface *should* be hidden here, unlike on the flat
    # panels. Culling back faces is not enough even when inflated: on a rounded
    # hemisphere plenty of lateral faces still tilt towards a medial camera, and
    # their outlines print straight through the medial wall. So depth-test every
    # contour face against what is actually painted, at the mesh's own spacing --
    # a finer grid is mostly empty cells and tests nothing.
    res = 100
    cen, cdep = uv[f[front]].mean(1), dep[f[front]].mean(1)
    lo, span = cen.min(0), (cen.max(0) - cen.min(0)).max()
    buf = np.full((res, res), -np.inf)
    g = ((cen - lo) / span * (res - 1)).astype(int).clip(0, res - 1)
    np.maximum.at(buf, (g[:, 1], g[:, 0]), cdep)

    segs = []
    for r in sorted(sig):
        code = CODE.get(r)
        if code is None or not (lab == code).any():
            continue
        seg, fid, pid = contour(uv, f, lab == code)
        if not len(fid):
            continue
        c, d = uv[f[fid]].mean(1), dep[f[fid]].mean(1)
        q = ((c - lo) / span * (res - 1)).astype(int).clip(0, res - 1)
        near = buf[q[:, 1], q[:, 0]]
        # Generous on purpose: a contour face on the same surface sits a few mm
        # behind the local front at worst, while one leaking from the far side
        # is tens of mm back. Tight tolerances just chew holes in the outline.
        vis = front[fid] & np.isfinite(near) & (d >= near - DEPTH_TOL)

        # Cluster what is left *after* culling, not the whole loop: a contour
        # that wraps from lateral to medial keeps one long visible arc, and
        # grouping by the original loop lets its stray far-side faces ride along
        # as dots in the middle of nowhere. Short runs are those strays.
        if vis.any():
            e = pid[vis]
            G = sparse.coo_matrix((np.ones(2 * len(e)),
                                   (np.r_[e[:, 0], e[:, 1]], np.r_[e[:, 1], e[:, 0]])),
                                  shape=(pid.max() + 1,) * 2).tocsr()
            _, cc = connected_components(G, directed=False)
            u, cnt = np.unique(cc[e[:, 0]], return_counts=True)
            runs = u[cnt >= MIN_RUN]
            vis[np.flatnonzero(vis)] = np.isin(cc[e[:, 0]], runs)
        if report is not None:                  # front-facing, drawn, arcs drawn
            report[r] = (int(front[fid].sum()), int(vis.sum()),
                         int(len(runs)) if vis.any() else 0)
        segs.append(seg[vis])
    return tri[order], face[order], np.concatenate(segs) if segs else np.zeros((0, 2, 2))


def tile(panels, position):
    """Place the panels by their own bounding boxes, exactly as pyseg does."""
    box = [np.array([t.reshape(-1, 2).min(0), t.reshape(-1, 2).max(0)])
           for t, _, _ in panels]
    cw = max(b[1, 0] - b[0, 0] for b in box)
    ch = max(b[1, 1] - b[0, 1] for b in box)
    pad = 0.04 * cw
    cells = ([(0, i) for i in range(len(panels))] if position == "dispersed"
             else [(0, 0), (1, 0), (1, 1), (0, 1)])
    off = [np.array([c * (cw + pad) + (cw - (b[1, 0] - b[0, 0])) / 2 - b[0, 0],
                     -r * (ch + pad) + (ch - (b[1, 1] - b[0, 1])) / 2 - b[0, 1]])
           for (r, c), b in zip(cells, box)]
    return off, cw + pad, ch + pad


cmap, norm = mpl.colormaps[CMAP], mpl.colors.Normalize(-VLIM, VLIM)


def dummy_data(seed=0, cut=2.6):
    rng = np.random.default_rng(seed)
    values = {r: v for r, v in zip(sorted(REGION.values()),
                                   rng.normal(0, 2, len(REGION)))}
    return values, {r for r, v in values.items() if abs(v) > cut}


def main(position="dispersed", values=None, sig=None, stem=None, title=None,
         label="t (dummy)", lw=None, panels=None):
    if values is None:
        values, sig = dummy_data()
        print(f"{len(sig)} of {len(values)} regions 'significant': {sorted(sig)[:4]} ...")
    stem = stem or f"examples/dk_surface_{position}"
    lw = OUTLINE_LW if lw is None else lw
    panels = panels or [render(h, a, values, sig) for h, _, a in PANELS]
    off, cw, ch = tile(panels, position)

    rows = 1 if position == "dispersed" else 2
    fig, ax = plt.subplots(figsize=(16, 16 * rows * ch / (len(panels) / rows * cw)))
    for (tri, face, seg), d in zip(panels, off):
        mesh = PolyCollection(tri + d, facecolors=face, edgecolors=face,
                              linewidths=SEAM_LW, zorder=1)       # seal the seams
        mesh.set_rasterized(True)   # 260k triangles stay an image inside the PDF;
        ax.add_collection(mesh)     # the outlines and text stay vector
        if len(seg):
            ax.add_collection(LineCollection(seg + d, colors="k", lw=lw,
                                             capstyle="round", zorder=3))
    ax.autoscale_view()
    ax.set(aspect=1, title=title or
           f"DK on {SURFACE} — black outline: |t| > 2.6   ({position})")
    ax.axis("off")
    fig.colorbar(mpl.cm.ScalarMappable(norm, cmap), ax=ax, fraction=0.02, pad=0.02
                 ).set_label(label)
    for ext in ("pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return panels


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "dispersed")
