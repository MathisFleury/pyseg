#!/usr/bin/env Rscript
# Export ggseg atlases to pyseg's format. Run from the repo root:
#   Rscript tools/export_ggseg_atlas.R dk aseg glasser schaefer7_400
#
#   data/<atlas>/_aliases          ggseg's own labels -> region file names
#   data/<atlas>/_panels           panel dir, hemi, side (left to right)
#   data/<atlas>/<panel>/<region>  SVG path data, panel-local coordinates
#   data/<atlas>/<panel>/_wall     the "no data here" shape for that panel
suppressMessages(library(sf))
for (p in c("ggseg", "ggsegSchaefer", "ggsegGlasser", "ggsegExtra"))
  suppressWarnings(suppressMessages(try(library(p, character.only = TRUE), silent = TRUE)))

out_root <- "pyseg/data"   # run from the repo root

# ggseg y grows upward, pyseg reads SVG coordinates (y down), hence oy - Y.
as_path <- function(geom, ox, oy) {
  cs <- st_coordinates(geom)
  ring <- do.call(paste, c(as.data.frame(cs[, setdiff(colnames(cs), c("X", "Y")), drop = FALSE]), sep = "-"))
  parts <- split(sprintf("%.2f,%.2f", cs[, "X"] - ox, oy - cs[, "Y"]), ring)
  paste(vapply(parts, function(p) paste0("M ", p[1], " L ", paste(p[-1], collapse = " "), " Z"), ""),
        collapse = " ")
}
clean <- function(x) gsub("[^A-Za-z0-9_.+-]", "", gsub(" ", "", x))

# ggseg's atlases were digitised region by region, not built as a partition, so
# neighbours genuinely overlap -- in dk's left lateral panel, 43 of 210 pairs do,
# fusiform and inferior temporal by 15.6% of the smaller one. Drawn as separate
# patches that shows up as doubled, uneven borders. Cutting each region out of
# what has already been placed makes the panel a true partition, so every border
# is one line. Smallest first: a contested strip costs a small parcel a large
# share of itself and a big one almost nothing.
TOPOLOGY <- TRUE

polys <- function(g) {
  g <- suppressWarnings(st_collection_extract(g, "POLYGON"))
  if (length(g) == 0) return(g)
  p <- suppressWarnings(st_cast(g, "POLYGON"))
  a <- as.numeric(st_area(p))
  keep <- a >= 0.02 * sum(a)          # st_difference sheds slivers; drop them
  if (!any(keep)) keep[which.max(a)] <- TRUE
  st_union(p[keep])
}

partition <- function(geo) {
  ord <- order(vapply(geo, function(g) as.numeric(sum(st_area(g))), 0))
  acc <- NULL
  for (k in names(geo)[ord]) {
    g <- geo[[k]]
    if (!is.null(acc)) {
      cut <- polys(suppressWarnings(st_difference(g, acc)))
      if (length(cut) && !all(st_is_empty(cut))) g <- cut
    }
    geo[[k]] <- g
    acc <- if (is.null(acc)) g else suppressWarnings(st_union(acc, g))
  }
  geo
}

for (name in commandArgs(TRUE)) {
  d <- get(name)$data
  d$side <- if ("side" %in% names(d)) as.character(d$side) else ""
  d$hemi <- as.character(d$hemi)

  # Regions present in one hemisphere only keep their bare name (Schaefer already
  # encodes the side, aseg's brainstem has none); the rest get _left / _right.
  reg <- !is.na(d$region)
  bilateral <- any(tapply(d$hemi[reg], d$region[reg], function(h) length(unique(h))) > 1)
  paired <- bilateral & d$hemi %in% c("left", "right")
  d$key <- ifelse(reg, ifelse(paired, paste0(clean(d$region), "_", d$hemi), clean(d$region)), NA)

  # Panels are (hemi, side) groups, merged where their x-ranges touch: aseg's two
  # coronal halves are one picture, dk's four views are four separate ones.
  d$pane <- paste(d$hemi, d$side)
  xr <- t(vapply(split(seq_len(nrow(d)), d$pane),
                 function(i) as.numeric(st_bbox(d[i, ])[c("xmin", "xmax")]), numeric(2)))
  xr <- xr[order(xr[, 1]), , drop = FALSE]
  tol <- 0.02 * diff(range(xr))   # aseg's coronal halves abut, they don't overlap
  grp <- rep(1L, nrow(xr)); reach <- xr[1, 2]
  for (i in seq_len(nrow(xr))[-1]) {
    disjoint <- xr[i, 1] - reach > tol
    grp[i] <- grp[i - 1] + disjoint
    reach <- if (disjoint) xr[i, 2] else max(reach, xr[i, 2])
  }
  d$grp <- grp[match(d$pane, rownames(xr))]

  out <- character()
  panels <- character()
  for (g in sort(unique(d$grp))) {
    sub <- d[d$grp == g, ]
    hemis <- unique(sub$hemi); side <- unique(sub$side)[1]
    hemi <- if (length(hemis) == 1) hemis else ""       # "" = the panel holds both
    pdir <- clean(if (nchar(hemi)) paste(hemi, side, sep = "_") else side)
    b <- st_bbox(sub)
    keys <- unique(na.omit(sub$key))
    geo <- lapply(keys, function(k)
      st_union(st_make_valid(st_geometry(sub)[which(sub$key == k)])))
    names(geo) <- keys
    if (TOPOLOGY) geo <- partition(geo)
    for (k in keys)
      out <- c(out, paste("R", pdir, k, as_path(geo[[k]], b["xmin"], b["ymax"]), sep = "\t"))
    if (any(is.na(sub$key)))
      out <- c(out, paste("W", pdir, "_wall",
                          as_path(st_geometry(sub)[which(is.na(sub$key))],
                                  b["xmin"], b["ymax"]), sep = "\t"))
    panels <- c(panels, paste("P", pdir, hemi, side, sep = "\t"))
  }
  if ("label" %in% names(d)) {
    a <- unique(data.frame(label = as.character(d$label), key = d$key))
    a <- a[!is.na(a$label) & !is.na(a$key), ]
    panels <- c(panels, paste("A", a$label, a$key, sep = "\t"))
  }
  dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
  con <- gzfile(file.path(out_root, paste0(name, ".tsv.gz")), "w")
  writeLines(c(panels, out), con)
  close(con)
  np <- sum(startsWith(panels, "P"))
  cat(sprintf("%-16s %4d regions, %d panels, %5.0f kB\n", name,
              length(unique(na.omit(d$key))), np,
              file.size(file.path(out_root, paste0(name, ".tsv.gz"))) / 1024))
}
