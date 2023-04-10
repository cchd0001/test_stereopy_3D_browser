# test_stereopy_3D_browser

## loading data
```
import anndata as ad
adata = ad.read_h5ad("D:/L3_b.h5ad")
```
## start a 3D atlas server without mesh
```
from stereopy_3D_browser import launch
launch(adata,meshes={},cluster_label='annotation',spatial_label='spatial')
```

## start a 3D atlas server with meshes
```
from stereopy_3D_browser import launch
launch(adata,meshes={'shell':'D:/shell.obj','midgut':'D:/midgut.obj'},cluster_label='annotation',spatial_label='spatial')
```
