# AOI-first workflow and country boundaries

MowEO can initialize the acquisition grid independently of grassland data.
The original grassland-first `moweo init` command remains available.

## 1. Supply an area of interest

Provide valid Polygon/MultiPolygon features in GeoPackage, GeoJSON or another
GeoPandas-readable vector format with a declared CRS. For a country, use your
chosen administrative boundary dataset. MowEO does not bundle or automatically
download country borders; their precision, vintage and license remain explicit.

```sh
moweo --project my_project init-aoi --aoi /absolute/path/aoi.gpkg --epsg 32633 --year 2026
```

To select a country from a multi-country file:

```sh
moweo --project my_project init-aoi --aoi /absolute/path/countries.gpkg --country Slovenia --country-field NAME --epsg 32633 --year 2026
```

Change `NAME` and the name/code to match the source dataset. Matching is exact,
case-insensitive; it is not geocoding. Use `--layer` for a specific source layer.
Use a suitable **projected CRS in metres**, typically the local UTM zone.
EPSG:32633 is an example for the original study area, not a global default.
For extensive countries spanning several projection zones, choose a suitable
regional projection or partition the area into separate projects. Web Mercator
is not appropriate for reliable 5 km ground-distance interpretation. All grid
sizes/buffers are measured in the selected projection's metre units.

## 2. Grid construction and identifiers

The AOI is reprojected, its polygons combined and intersecting 5,000 m square
tiles generated. Tiles with only point/line contact are excluded. Holes and
disconnected polygons are respected. Borders are **not clipped from the tile
geometries**: the downloader needs full 500-by-500 NDVI and 250-by-250 auxiliary
grids. The mask later excludes pixels outside the AOI using pixel-centre rules.

Identifiers are deterministic:

```text
eo_32633_380000_5090000
   EPSG  west   south
```

Coordinates are lower-left tile coordinates, snapped to multiples of 5,000 m.
The identifier is independent of year and country name and unique for a tile
in that EPSG grid. Overlapping AOIs in the same grid share tile IDs. Changing
the CRS changes IDs; do not reuse downloaded cubes across different grids.

The project gets `data/aoi.gpkg`, `data/EOPatch_PET.gpkg` and
`project_config.json`. `AREA_PCT` records the AOI's share of each full tile.
The command refuses to overwrite an existing project grid.

For a small pilot, append `--limit 10` to `init-aoi`. It persists the first
ten tiles in west/south coordinate order as the whole project grid. Unlike
`acquire --limit 10`, later processing is restricted to those same tiles too.

## 3. Acquire before adding masks

```sh
moweo --project my_project prepare
moweo --project my_project acquire
moweo --project my_project acquire --run --parallel 2
```

The second command is a dry run. Actual acquisition requires CDSE login and
credits and computes two cloud-masked annual index cubes per selected tile.
No grassland classification or mask is needed at this stage. The full tile is
downloaded, including any area outside the AOI; this increases cost at borders.
Inspect the configured tile count before submitting a whole-country job.

## 4. Add the grassland definition later

```sh
moweo --project my_project masks --grasslands /absolute/path/private_grasslands.gpkg
moweo --project my_project validate --require-stats
moweo --project my_project process --workers 4
```

The grassland input must contain only the grassland polygons to analyse.
Source features are buffered by +5 m for detection and -7.5 m for statistics
**before clipping to the AOI**, then rasterized by pixel centres on the fixed
10 m grid. This avoids artificial erosion along tile seams. The two manifests
include every configured tile, even when it has no grass pixels. Empty tiles
have no mask TIFF. The command refuses to overwrite existing masks.

This command takes vector grasslands. If users already have binary rasters,
they can supply the per-tile aligned GeoTIFFs and manifests described in
`data/README.md`; arbitrary full-country raster reprojection/tiling is not
implemented by this command. The boundary of a country alone is not a grassland
mask. Keep a record of the mask's definition, source, date and rights.

After masks exist, `moweo run --run` can connect acquisition and local processing.
Without masks, use acquisition alone; the full processing stage intentionally
fails rather than treating all AOI pixels as grassland. The AOI-first workflow
does not modify existing projects or the active private 2026 ten-tile pilot.
