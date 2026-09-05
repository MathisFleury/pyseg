"""ggseg-style brain plots, with significance outlines that stay on top.

The one thing this fixes: ggseg draws each region as a single patch carrying
both its fill and its edge, so a neighbour drawn later paints over the black
outline of a significant region. Here fills and outlines are separate passes.
"""
import gzip
import os
import re
import warnings
from collections import namedtuple
from functools import lru_cache

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, PathPatch
from matplotlib.path import Path
from matplotlib.transforms import Affine2D

__version__ = "0.3"
__all__ = ["plot_brain", "plot_dk", "plot_aseg", "plot_tracts", "plot_connectome",
           "geom_brain", "as_brain_df", "brain_atlases", "brain_regions",
           "brain_labels", "brain_views"]

_DATA = os.path.join(os.path.dirname(__file__), "data")
_CMD = {"M": Path.MOVETO, "L": Path.LINETO, "Q": Path.CURVE3,
        "C": Path.CURVE4, "Z": Path.CLOSEPOLY}
_META = ("_aliases", "_panels")
# jhu is hand-vendored and flat; it names its own silhouette. Exported atlases
# are one directory per panel, with the silhouette in "_wall".
_FLAT_WALL = {"jhu": ["NA"]}
_HEMI = {"lh": "left", "l": "left", "left": "left", "lft": "left",
         "rh": "right", "r": "right", "right": "right", "rgt": "right"}

Panel = namedtuple("Panel", "name hemi side paths wall")


# ponytail: SVG subset the atlases use (M/L/Q/C/Z, absolute).
def _parse(svg):
    codes, verts, start = [], [], np.zeros(2)
    tokens = re.split("([A-Za-z])", svg)[1:]
    for cmd, values in zip(tokens[::2], tokens[1::2]):
        if cmd.upper() == "Z":
            points = start.reshape(1, 2)   # CLOSEPOLY repeats the subpath start,
        else:                              # so no vertex sits at a fake origin
            points = np.reshape([list(map(float, e.split(","))) for e in values.split() if e],
                                (-1, 2))
            if cmd.upper() == "M" and len(points):
                start = points[0]
        codes += [_CMD[cmd.upper()]] * len(points)
        verts.append(points)
    return np.array(codes), np.concatenate(verts)


def _canon(name):
    """Region name -> spelling-independent 'stem|hemi' key: 'lh_bankssts',
    'Left bankssts' and 'bankssts_left' all collapse to the same thing."""
    parts = re.split(r"[ _.\-]+", str(name).strip().lower())
    hemi = ""
    if len(parts) > 1 and parts[0] in _HEMI:
        hemi = _HEMI[parts.pop(0)]
    elif len(parts) > 1 and parts[-1] in _HEMI:
        hemi = _HEMI[parts.pop()]
    # FreeSurfer 7 dropped the "-Proper" suffix ggseg still carries on the thalamus.
    return "".join(parts).removesuffix("proper") + "|" + hemi


def _hemi(name):
    return _canon(name).rsplit("|", 1)[1] if name else ""


def _rings(path):
    """-> [(vertices, codes)], one per subpath."""
    at = [i for i, c in enumerate(path.codes) if c == Path.MOVETO] + [len(path.codes)]
    return [(path.vertices[a:b], path.codes[a:b]) for a, b in zip(at, at[1:])]


def _area(v):
    x, y = v[:, 0], v[:, 1]
    return abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1))) / 2


def _despeckle(path, tol, ref=0.0, frac=0.02):
    """ggseg's traced geometry sheds strays: sub-pixel rings, and slivers of a
    parcel left behind in a view where the parcel barely reaches. Invisible when
    filled, but each still takes an outline stroke and shows up as a lone speck
    or dash.

    Drops a ring under `tol` across, or under `frac` of `ref` -- the area of the
    biggest ring that parcel has anywhere in the atlas. A small ring *inside* the
    shape's own main ring is a hole, and holes stay. None if nothing survives.
    """
    parts = _rings(path)
    area = [_area(v) for v, _ in parts]
    main = int(np.argmax(area))
    if area[main] < frac * ref:          # the parcel barely reaches this view
        return None
    inside = Path(*parts[main]).contains_point
    keep = [(v, c) for i, ((v, c), a) in enumerate(zip(parts, area))
            if (v.max(0) - v.min(0)).max() >= tol
            and (i == main or a >= frac * area[main] or inside(v.mean(0)))]
    if not keep:
        return None
    return Path(np.concatenate([v for v, _ in keep]), np.concatenate([c for _, c in keep]))


def _smooth(path, n, frac=0.006):
    """Corner cutting (Chaikin), cutting a distance set by the ring's own size
    rather than a fraction of each edge -- ggseg's traces are coarse enough to
    look spiky once you stroke them, but a fractional cut would visibly shrink
    the sparsely-traced parcels. Short zigzag edges round right off, long ones
    barely move. Curves are left alone."""
    if n <= 0:
        return path
    verts, codes = [], []
    for v, c in _rings(path):
        if not set(c) <= {Path.MOVETO, Path.LINETO, Path.CLOSEPOLY}:
            verts.append(v); codes.append(c); continue
        p = v[:-1]                                    # drop the CLOSEPOLY vertex
        if len(p) > 2 and np.allclose(p[0], p[-1]):
            p = p[:-1]                                # ...and the repeated first
        d = frac * (p.max(0) - p.min(0)).max()        # scaled to this ring
        for _ in range(n):
            edge = np.roll(p, -1, axis=0) - p
            t = np.minimum(0.25, d / np.maximum(np.hypot(*edge.T), 1e-12))[:, None]
            cut = np.empty((2 * len(p), 2))
            cut[0::2], cut[1::2] = p + t * edge, p + (1 - t) * edge
            p = cut
        verts.append(np.vstack([p, p[0]]))
        codes.append([Path.MOVETO] + [Path.LINETO] * (len(p) - 1) + [Path.CLOSEPOLY])
    return Path(np.concatenate(verts), np.concatenate(codes))


def _load(atlas):
    """Read an atlas -> ([(panel, hemi, side, {region: Path}, [wall])], aliases).

    Built-ins are one gzipped table each (see tools/export_ggseg_atlas.R); a
    directory of one file per region still works, for atlases you trace yourself.
    """
    packed = atlas if os.path.isfile(atlas) else os.path.join(_DATA, f"{atlas}.tsv.gz")
    if os.path.isfile(packed):
        panels, alias, paths, wall = [], {}, {}, {}
        with gzip.open(packed, "rt") as fh:
            for line in fh:
                kind, *rest = line.rstrip("\n").split("\t")
                if kind == "P":
                    panels.append(tuple(rest[:3]))
                elif kind == "A":
                    alias[rest[0]] = rest[1]
                elif kind in "RW":
                    pan, name, svg = rest
                    d = wall if kind == "W" else paths
                    d.setdefault(pan, {})[name] = Path(*_parse(svg)[::-1])
        return ([(p, h, s, paths.get(p, {}), list(wall.get(p, {}).values()))
                 for p, h, s in panels], alias)

    wd = atlas if os.path.isdir(atlas) else os.path.join(_DATA, atlas)
    if not os.path.isdir(wd):
        raise ValueError(f"unknown atlas {atlas!r}; built-in: {sorted(brain_atlases())}")
    name = os.path.basename(wd.rstrip("/"))
    spec = os.path.join(wd, "_panels")
    rows = ([ln.rstrip("\n").split("\t") for ln in open(spec) if ln.strip()]
            if os.path.isfile(spec) else [("", "", "")])
    wall_names = _FLAT_WALL.get(name, [])
    raw = []
    for pname, hemi, side in rows:
        pdir = os.path.join(wd, pname)
        paths, wall = {}, []
        for f in sorted(os.listdir(pdir)):
            if f in _META or os.path.isdir(os.path.join(pdir, f)):
                continue
            with open(os.path.join(pdir, f)) as fh:
                path = Path(*_parse(fh.read())[::-1])
            (wall.append(path) if f.startswith("_") or f in wall_names
             else paths.__setitem__(f, path))
        raw.append((pname or name, hemi, side, paths, wall))
    alias = {}
    if os.path.isfile(os.path.join(wd, "_aliases")):
        with open(os.path.join(wd, "_aliases")) as fh:
            alias = dict(ln.rstrip("\n").split("\t") for ln in fh if "\t" in ln)
    return raw, alias


@lru_cache(maxsize=None)
def _atlas(atlas, smooth=2):
    """-> (tuple of Panel, dict of ggseg label -> region name)."""
    raw, alias = _load(atlas)

    ref = {}                    # biggest each parcel gets anywhere in the atlas
    for r in raw:
        for k, pl in r[3].items():
            ref[k] = max(ref.get(k, 0.0), max(_area(v) for v, _ in _rings(pl)))

    out = []
    for pname, hemi, side, paths, wall in raw:
        # 0.003 of the panel: ggseg's strays are all under 0.002, real rings
        # start around 0.005.
        span = Path.make_compound_path(*(list(paths.values()) + wall)).get_extents()
        tol = 0.003 * max(span.width, span.height)
        out.append(Panel(pname, hemi, side,
                         {k: _smooth(c, smooth) for k, v in paths.items()
                          if (c := _despeckle(v, tol, ref[k])) is not None},
                         [_smooth(c, smooth) for w in wall
                          if (c := _despeckle(w, tol)) is not None]))

    # A parcel small enough to be all speck is still a parcel: put it back whole,
    # in whichever panel shows most of it.
    kept = {k for p in out for k in p.paths}
    for k in {k for r in raw for k in r[3]} - kept:
        i = max((i for i, r in enumerate(raw) if k in r[3]),
                key=lambda i: raw[i][3][k].get_extents().size.max())
        out[i].paths[k] = _smooth(raw[i][3][k], smooth)
    return tuple(out), alias


def brain_atlases():
    """Every built-in atlas, with its region count."""
    return {a: len(brain_regions(a)) for a in
            sorted(f[:-7] for f in os.listdir(_DATA) if f.endswith(".tsv.gz"))}


def brain_regions(atlas="dk"):
    """Region names this atlas draws."""
    return sorted({k for p in _atlas(atlas)[0] for k in p.paths})


def brain_labels(atlas="dk"):
    """The atlas' own labels (FreeSurfer names, where ggseg carries them)."""
    return sorted(_atlas(atlas)[1])


def brain_views(atlas="dk"):
    """(panel, hemisphere, view) of each panel, in ggseg's left-to-right order."""
    return [(p.name, p.hemi, p.side) for p in _atlas(atlas)[0]]


def _resolve(atlas, data, sig):
    """Map user region names onto atlas ones, warning about what didn't land."""
    ps, alias = _atlas(atlas)
    names = {k for p in ps for k in p.paths}
    lut = {_canon(a): k for a, k in alias.items() if k in names}
    lut.update({_canon(k): k for k in names})       # a region's own name wins
    data = {lut.get(_canon(k), k): v for k, v in data.items()}
    sig = [lut.get(_canon(k), k) for k in sig]
    unknown = [k for k in list(data) + sig if k not in names]
    if unknown:
        warnings.warn(f"not in atlas {atlas!r}, ignored: {unknown}")
    return data, sig


def _select(ps, hemisphere, view):
    sel = [p for p in ps if view in (None, p.side)
           and (hemisphere is None or p.hemi in (hemisphere, ""))]
    if not sel:
        raise ValueError(f"no panel matches hemisphere={hemisphere!r} view={view!r}; "
                         f"have {[(p.hemi, p.side) for p in ps]}")
    return sel


def _grid(sel, position):
    """-> per-panel (dx, dy), tiling the panels ggseg-style."""
    box = [Path.make_compound_path(*(list(p.paths.values()) + p.wall)).get_extents()
           for p in sel]
    cw, ch = max(b.width for b in box), max(b.height for b in box)
    pad = 0.03 * cw
    if position == "stacked":
        rows = list(dict.fromkeys(p.side for p in sel))
        cols = list(dict.fromkeys(p.hemi for p in sel))
        cells = [(rows.index(p.side), cols.index(p.hemi)) for p in sel]
    elif position == "dispersed":
        cells = [(0, i) for i in range(len(sel))]
    else:
        raise ValueError(f"position must be 'dispersed' or 'stacked', not {position!r}")
    return [(c * (cw + pad) + (cw - b.width) / 2 - b.x0,
             r * (ch + pad) + (ch - b.height) / 2 - b.y0)
            for (r, c), b in zip(cells, box)]


def _lay_out(atlas, hemisphere, view, position, smooth):
    """-> ([(region, Path)], [wall Path]) with the panels tiled into place."""
    sel = _select(_atlas(atlas, smooth)[0], hemisphere, view)
    fills, walls = [], []
    for panel, (dx, dy) in zip(sel, _grid(sel, position)):
        move = Affine2D().translate(dx, dy).transform_path
        fills += [(k, move(v)) for k, v in panel.paths.items()
                  if not hemisphere or _hemi(k) in ("", hemisphere)]
        walls += [move(w) for w in panel.wall]
    return fills, walls


def plot_brain(data=None, atlas="dk", sig=(), hemisphere=None, view=None,
               position="dispersed", smooth=2, cmap="Spectral", vmin=None, vmax=None,
               na_color="0.85", edgecolor="w", lw=0.5, sig_color="k", sig_lw=2.0,
               background="w", figsize=None, title="", ylabel="", colorbar=True,
               ax=None):
    """Plot values on a brain atlas.

    data : dict / pandas Series mapping region name -> value (`brain_regions()`).
    sig  : region names to outline in `sig_color`, drawn above every fill and
           every other outline -- no more fixing the stacking in Illustrator.
           Accepts a list/set of names, or a dict of name -> bool.
    hemisphere : "left" / "right", or None for both.
    view : one of `brain_views(atlas)`'s views ("lateral", "medial", ...), or None.
    position : "dispersed" (one row, as ggseg) or "stacked" (views x hemispheres).
    smooth : rounds of corner cutting on the outlines; 0 for ggseg's raw traces.
    Returns (fig, ax); nothing is shown or saved for you.
    """
    if isinstance(sig, dict):
        sig = [k for k, v in sig.items() if v]
    data, sig = _resolve(atlas, dict(data.items()) if data is not None else {}, sig)

    fills, walls = _lay_out(atlas, hemisphere, view, position, smooth)

    if ax is None:
        verts = np.concatenate([p.vertices for _, p in fills] + [w.vertices for w in walls])
        (x0, y0), (x1, y1) = verts.min(0) - 1, verts.max(0) + 1
        figsize = figsize or (14, 14 * (y1 - y0) / (x1 - x0))
        fig, ax = plt.subplots(figsize=figsize, facecolor=background)
        ax.set(xlim=(x0, x1), ylim=(y1, y0), aspect=1)  # y flipped: SVG coords
        ax.set_facecolor(background)
        ax.axis("off")
    fig = ax.figure
    ax.set_title(title)

    values = [v for k, v in data.items() if k in dict(fills)]
    cmap = cmap if isinstance(cmap, mpl.colors.Colormap) else mpl.colormaps[cmap]
    norm = mpl.colors.Normalize(min(values, default=0) if vmin is None else vmin,
                                max(values, default=1) if vmax is None else vmax)

    for w in walls:      # silhouette / medial wall: the "no data here" fill
        ax.add_patch(PathPatch(w, facecolor=na_color, edgecolor="none", zorder=0))
    for name, p in fills:
        c = cmap(norm(data[name])) if name in data else na_color
        ax.add_patch(PathPatch(p, facecolor=c, edgecolor="none", zorder=1))
    for _, p in fills:                                   # every outline, then
        ax.add_patch(PathPatch(p, facecolor="none", edgecolor=edgecolor,
                               lw=lw, zorder=2))
    for name, p in fills:                                # the significant ones on top
        if name in sig:
            ax.add_patch(PathPatch(p, facecolor="none", edgecolor=sig_color,
                                   lw=sig_lw, zorder=3))

    if colorbar and values:
        cb = fig.colorbar(mpl.cm.ScalarMappable(norm, cmap), ax=ax,
                          fraction=0.025, pad=0.02)
        cb.set_label(ylabel)
    return fig, ax


def plot_dk(data=None, **kw):
    """Cortex, Desikan-Killiany."""
    return plot_brain(data, atlas="dk", **kw)


def plot_aseg(data=None, **kw):
    """Subcortical structures, FreeSurfer aseg."""
    return plot_brain(data, atlas="aseg", **kw)


def plot_tracts(data=None, **kw):
    """White matter tracts, JHU."""
    return plot_brain(data, atlas="jhu", **kw)


def plot_connectome(edges, atlas="dk", nodes=None, threshold=0.0, hemisphere=None,
                    view=None, position="dispersed", smooth=2, cmap="RdBu_r",
                    vmin=None, vmax=None, node_color="k", node_size=12,
                    lw=(0.4, 3.0), curve=0.2, na_color="0.9", background="w",
                    figsize=None, title="", ylabel="", colorbar=True, ax=None):
    """Draw a connectivity matrix as arcs between region centroids.

    edges : square pandas DataFrame indexed by region name, or a square array
            with `nodes` giving the region name of each row.
    threshold : edges with |weight| at or below this are not drawn.
    lw : (thinnest, thickest) line width, mapped from |weight|.
    Regions live in separate panels, so an edge between two views is a straight
    hop across the figure -- pass `hemisphere` or `view` to keep it in one.
    Returns (fig, ax).
    """
    import numpy as np
    w = np.asarray(getattr(edges, "values", edges), dtype=float)
    nodes = list(getattr(edges, "index", nodes) if nodes is None else nodes)
    if w.shape != (len(nodes), len(nodes)):
        raise ValueError(f"edges is {w.shape}, but there are {len(nodes)} nodes")
    nodes = _resolve(atlas, {}, nodes)[1]

    fig, ax = plot_brain(atlas=atlas, hemisphere=hemisphere, view=view,
                         position=position, smooth=smooth, na_color=na_color,
                         background=background, figsize=figsize, title=title,
                         colorbar=False, ax=ax)

    # A region can show up in two views; put its node in the bigger one.
    seen = {}
    for k, p in _lay_out(atlas, hemisphere, view, position, smooth)[0]:
        e = p.get_extents()
        if e.size.prod() > seen.get(k, (0, None))[0]:
            seen[k] = (e.size.prod(), np.array([e.x0 + e.width / 2, e.y0 + e.height / 2]))
    xy = {k: v[1] for k, v in seen.items()}

    i, j = np.triu_indices(len(nodes), k=1)
    hit = np.abs(w[i, j]) > threshold
    i, j, val = i[hit], j[hit], w[i, j][hit]
    cmap = cmap if isinstance(cmap, mpl.colors.Colormap) else mpl.colormaps[cmap]
    lim = max(abs(val).max(), 1e-12) if len(val) else 1.0
    norm = mpl.colors.Normalize(-lim if vmin is None else vmin,
                                lim if vmax is None else vmax)
    span = max(abs(val).max(), 1e-12) if len(val) else 1.0

    for a, b, v in zip(i, j, val):
        if nodes[a] not in xy or nodes[b] not in xy:
            continue
        ax.add_patch(FancyArrowPatch(xy[nodes[a]], xy[nodes[b]], arrowstyle="-",
                                     connectionstyle=f"arc3,rad={curve}",
                                     color=cmap(norm(v)), zorder=4,
                                     lw=lw[0] + (lw[1] - lw[0]) * abs(v) / span))
    if node_size:
        pts = np.array([xy[n] for n in nodes if n in xy])
        if len(pts):
            ax.scatter(*pts.T, s=node_size, color=node_color, zorder=5)
    if colorbar and len(val):
        cb = fig.colorbar(mpl.cm.ScalarMappable(norm, cmap), ax=ax,
                          fraction=0.025, pad=0.02)
        cb.set_label(ylabel)
    return fig, ax


def as_brain_df(atlas="dk"):
    """The atlas as tidy polygons, like ggseg's as.data.frame(): one row per
    vertex, columns region / hemi / side / panel / ring / x / y. Silhouette rows
    have region = None. Feed it to plotnine, or to anything else."""
    import pandas as pd
    out = []
    for p in _atlas(atlas)[0]:
        for name, path in list(p.paths.items()) + [(None, w) for w in p.wall]:
            for j, poly in enumerate(path.to_polygons()):
                out.append(pd.DataFrame({
                    "region": name, "hemi": _hemi(name) or p.hemi, "side": p.side,
                    "panel": p.name, "ring": f"{p.name}/{name}/{j}",
                    "x": poly[:, 0], "y": -poly[:, 1]}))   # ggplot y grows up
    return pd.concat(out, ignore_index=True)


def geom_brain(data=None, atlas="dk", sig=(), hemisphere=None, view=None,
            position="dispersed", edgecolor="w", lw=0.3, sig_color="k",
            sig_lw=0.9):
    """The same plot as a plotnine object, so ggplot2's grammar does the rest:
    add your own scale_fill_*, theme_*, labs, or replace the facetting.

    >>> geom_brain(t_values, sig=fdr_hits) + scale_fill_gradient2() + labs(fill="t")
    """
    import pandas as pd
    from plotnine import (aes, coord_equal, facet_grid, facet_wrap, geom_polygon,
                          ggplot, theme_void)

    if isinstance(sig, dict):
        sig = [k for k, v in sig.items() if v]
    data, sig = _resolve(atlas, dict(data.items()) if data is not None else {}, sig)

    keep = [p.name for p in _select(_atlas(atlas)[0], hemisphere, view)]
    df = as_brain_df(atlas)
    df = df[df.panel.isin(keep)]
    if hemisphere:
        df = df[df.hemi.isin(["", hemisphere])]
    df["panel"] = pd.Categorical(df.panel, categories=keep, ordered=True)
    df["value"] = df.region.map(data)

    # The fix, as ggplot layers: outlines are a later layer, so they cannot be
    # painted over by a neighbour's fill.
    p = (ggplot(df, aes("x", "y", group="ring"))
         + geom_polygon(aes(fill="value"), colour=edgecolor, size=lw)
         + coord_equal() + theme_void())
    if len(hits := df[df.region.isin(sig)]):
        p += geom_polygon(hits, fill="none", colour=sig_color, size=sig_lw)
    return p + (facet_grid("side ~ hemi") if position == "stacked"
                else facet_wrap("~panel", nrow=1))
