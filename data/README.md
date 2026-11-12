# Private inputs: location and meaning

This folder intentionally contains documentation only. No project grassland
mask or grid is supplied. Use `prepare_inputs.py` to derive the following from
your own valid grassland polygons, or supply compatible files yourself:

Alternatively, use `moweo init-aoi` to create `aoi.gpkg` and a full-square tile
grid without masks, then `moweo masks` to rasterize grasslands on that same
grid later. These AOI-first masks are restricted to the AOI and retain stable
IDs such as `eo_32633_380000_5090000`. See `docs/aoi_workflow.md`.

```text
data/
  inputs/grasslands.gpkg              # user-provided grassland polygons
  EOPatch_PET.gpkg                    # full 5 km square grid
  grass_detection/
    grass_masks_manifest.csv
    <PATCH_ID>_grass.tif
  grass_stats/
    grass_masks_manifest.csv
    <PATCH_ID>_grass.tif
```

Paths can be changed in local `project_config.json`. The input polygons define
where grassland exists; they contain no mowing labels. Their date, land-use
definition, accuracy and license are the user's responsibility. Keep a record
of the source and vintage, especially when applying one mask to multiple years.

The detection mask is rasterized from polygons buffered outward by 5 m.
The statistics mask is rasterized from polygons buffered inward by 7.5 m.
Buffers are applied to source features in the selected metric CRS before tile
clipping. Pixel-centre rasterization (`all_touched=False`) is used. Shared feature
boundaries are also eroded; dissolve adjacent polygons beforehand if that is
not appropriate for your grassland definition. Small/narrow polygons may have
no inner pixels. The inner mask must be a subset of the detection mask.

Every nonempty mask is a one-band uint8 GeoTIFF with exactly 500 rows and 500
columns, north-up 10 m pixels, no rotation, and the configured projected CRS.
`1` means grassland included in the respective analysis; `0` means excluded
(also the file's NoData value). Do not encode exclusion as 255. The mask extent
must equal its 5 km tile, whose west/south coordinates are multiples of 5000 m.
The downloaded NDVI pixel centres must match the mask exactly. There is no
silent mask resampling.

Both manifests must contain one row for every configured tile, including empty
tiles. Required columns and an illustrative row (not distributed data):

```csv
PATCH_ID,filename,has_grass,grass_pixels
tile_500000_5100000,tile_500000_5100000_grass.tif,True,1234
tile_505000_5100000,,False,0
```

`grass_pixels` is the actual number of pixels equal to 1. For empty tiles use
an empty filename, `False`, zero pixels, and no TIFF. Both manifests use exactly
the same IDs as the grid. Legacy manifests may instead have `name`, `split_x`
and `split_y`; these are joined with underscores to construct `PATCH_ID`.
IDs may contain only letters, digits, underscores and hyphens.

The grid GeoPackage must have a CRS and unique `PATCH_ID` values. Legacy grid
fields are supported. `PET`, `BIOGEO_R`, and `AREA_PCT` are optional grouping
metadata for summary reports; the generic helper uses `all`, `all`, and 100.
These placeholder groups do not claim regional ecological classifications.

Acquisition needs the tile grid only. Both mask sets are needed for the full
local detection, finalization and QA pipeline. Do not place private inputs in
public Git history; all files under this folder except this README are ignored.
