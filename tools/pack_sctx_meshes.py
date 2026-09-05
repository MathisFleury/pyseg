"""Pack ENIGMA's subcortical surfaces into pyseg/data/sctx.npz.

A maintainer's tool, like tools/export_ggseg_atlas.R. Point it at an ENIGMA
toolbox checkout:

    ENIGMA_DATA=~/ENIGMA/enigmatoolbox/datasets python3 tools/pack_sctx_meshes.py

The sctx surface is 16 structures in fixed contiguous vertex blocks (ENIGMA's
own ordering); names below follow pyseg's aseg so the same data plots on the
flat atlas and on the surfaces without renaming.
"""
import os

import nibabel as nib
import numpy as np

E = os.environ.get("ENIGMA_DATA",
                   os.path.expanduser("~/POSTDOC/LIBRAIRY/ENIGMA/enigmatoolbox/datasets"))
STRUCTS = ["accumbensarea", "amygdala", "caudate", "hippocampus",
           "pallidum", "putamen", "thalamusproper", "lateralventricle"]
BLOCKS = {"left": [867, 1419, 3012, 3784, 1446, 4003, 3726, 7653],
          "right": [838, 1457, 3208, 3742, 1373, 3871, 3699, 7180]}

out = {"structs": np.array(STRUCTS)}
for hemi in ("left", "right"):
    v, f = nib.load(f"{E}/surfaces/sctx_{hemi[0]}h.gii").agg_data()
    v, f = np.asarray(v, np.float32), np.asarray(f, np.int32)
    lab = np.repeat(np.arange(len(STRUCTS), dtype=np.int8), BLOCKS[hemi])
    assert len(lab) == len(v), f"{hemi}: {len(lab)} labels for {len(v)} vertices"
    assert (lab[f].max(1) == lab[f].min(1)).all(), "a face straddles two structures"
    out[f"{hemi}_v"], out[f"{hemi}_f"], out[f"{hemi}_lab"] = v, f, lab
    print(f"{hemi}: {len(v)} vertices, {len(f)} faces")

dst = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "pyseg", "data", "sctx.npz")
np.savez_compressed(dst, **out)
print(f"wrote {dst}  {os.path.getsize(dst) / 1024:.0f} kB")
