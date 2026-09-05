# Atlas data format

One file per atlas, `<atlas>.tsv.gz`: a gzipped, tab-separated table whose first
column is the record type. Nothing but `gzip` is needed to read it.

```sh
gzcat dk.tsv.gz | cut -f1-3 | head
```

| type | columns | meaning |
|---|---|---|
| `P` | panel, hemisphere, view | one panel of the figure, in left-to-right order |
| `A` | label, region | an alias: ggseg's own label for a region |
| `R` | panel, region, path | a region's outline in that panel |
| `W` | panel, `_wall`, path | the "no data here" shape: silhouette or medial wall |

```
P	left_lateral	left	lateral
A	lh_bankssts	bankssts_left
R	left_lateral	bankssts_left	M 217.46,99.91 L 220.62,100.97 222.88,103.53 …
W	left_lateral	_wall	M 4.02,101.28 L 4.71,97.24 …
```

A region appears once per panel it is visible in, so an atlas has more `R`
records than regions — dk has 86 for 70 regions, because some show on both the
lateral and the medial view.

`hemisphere` is `left`, `right`, or empty when a panel holds both (aseg's
coronal slice). A region name carries its own side (`bankssts_left`) unless the
atlas already encodes it (`7Networks_LH_Vis_6`) or the structure is midline
(`brainstem`).

## Paths

The `path` column is a subset of SVG path data: absolute `M`, `L`, `Q`, `C` and
`Z`, with `x,y` pairs separated by spaces. Multiple `M…Z` runs in one record are
the rings of one region — an outer ring plus any holes.

Coordinates are **panel-local**: each panel's own bounding box starts at the
origin, so panels can be tiled into any layout. They follow the SVG convention
of **y increasing downward**, which is why `pyseg` flips the axis when it draws
and negates `y` in `as_brain_df()`.

## Reading one without pyseg

```python
import gzip
regions = {}
for line in gzip.open("dk.tsv.gz", "rt"):
    kind, *rest = line.rstrip("\n").split("\t")
    if kind == "R":
        panel, name, path = rest
        regions.setdefault(panel, {})[name] = path
```

## Where it comes from

`tools/export_ggseg_atlas.R` generates these from the R `ggseg` packages, and is
a maintainer's tool — using `pyseg` never requires R. The exporter also makes
each panel a true partition, since ggseg's atlases were digitised region by
region and neighbours genuinely overlap.

Geometry is derived from [ggseg](https://github.com/ggseg/ggseg) (MIT) and
[python-ggseg](https://github.com/ggseg/python-ggseg).

## Why this format, and what it should probably become

It is not a standard. It exists because the obvious layout — a directory per
atlas, a file per region — came to 10813 files for 23 atlases: 3 MB of path data
sitting in 46 MB of 4 kB blocks, which is not something to put in a wheel.

The format these polygons *should* be in is GeoJSON. ggseg already stores them
as `sf`, R's implementation of OGC Simple Features, so GeoJSON is closer to the
source than this is, and it would open the atlases to geopandas, shapely, QGIS
and D3 with no custom parser. Gzipped it costs about the same: 17 kB against
16 kB for dk. The blocker is small — GeoJSON has no curves, so `jhu`'s Bézier
segments would need flattening. If you are reading this before that migration
happened, it is still the right move.
