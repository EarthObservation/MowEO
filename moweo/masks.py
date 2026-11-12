"""Add grassland masks to an existing fixed tile grid, without changing IDs."""
import argparse
import re
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import rasterize
from rasterio.transform import from_origin
from shapely.ops import unary_union

from . import config as cfg
from .area import polygons, metric_crs
from .acquisition_prepare import _make_patch_id, _snap_tile_from_centroid


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--grasslands', required=True, type=Path, help='Only grassland Polygon/MultiPolygon features; not mixed land cover.')
    p.add_argument('--layer')
    args = p.parse_args()
    outputs = [(cfg.GRASS_ROOT, cfg.GRASS_MANIFEST, 5),
               (cfg.STATS_GRASS_ROOT, cfg.STATS_GRASS_MANIFEST, -7.5)]
    for root, manifest, _ in outputs:
        if manifest.exists() or (root.exists() and any(root.glob('*_grass.tif'))):
            p.error('Masks already exist. Use a separate project to change the mask definition.')
    crs = metric_crs(cfg.TARGET_EPSG)
    grid = gpd.read_file(cfg.PATCH_GPKG, layer=cfg.PATCH_LAYER).to_crs(crs)
    grid['PATCH_ID'] = _make_patch_id(grid)
    if len(grid) != cfg.EXPECTED_PATCHES or grid.PATCH_ID.duplicated().any():
        p.error('Grid count/IDs do not match the configured project.')
    if not grid.PATCH_ID.map(lambda value: bool(re.fullmatch(r'[A-Za-z0-9_-]+', value))).all():
        p.error('PATCH_ID may contain only letters, digits, underscores and hyphens.')
    source = polygons(args.grasslands, args.layer).to_crs(crs)
    aoi_path = cfg.SETTINGS.get('aoi_file')
    if aoi_path:
        path = Path(aoi_path)
        path = path if path.is_absolute() else cfg.PROJECT_DIR / path
        area = unary_union(polygons(path).to_crs(crs).geometry)
    else:
        area = unary_union(grid.geometry)
    # Buffer complete features before clipping at the AOI or tile boundaries.
    detection = source.geometry.buffer(5).intersection(area)
    inner = source.geometry.buffer(-7.5).intersection(area)
    all_rows = [[], []]
    for row in grid.itertuples():
        west, south, east, north, _ = _snap_tile_from_centroid(row.geometry, row.PATCH_ID)
        transform = from_origin(west, north, 10, 10)
        arrays = []
        for shapes in [detection, inner]:
            indexes = shapes.sindex.query(row.geometry, predicate='intersects')
            selected = [geom for geom in shapes.iloc[indexes] if not geom.is_empty]
            arr = rasterize(((geom, 1) for geom in selected), out_shape=(500, 500),
                            transform=transform, fill=0, dtype='uint8', all_touched=False) if selected else np.zeros((500, 500), dtype='uint8')
            arrays.append(arr)
        if np.any((arrays[1] == 1) & (arrays[0] != 1)):
            raise ValueError(f'{row.PATCH_ID}: inner mask is not a detection-mask subset.')
        for index, ((root, manifest, _), arr) in enumerate(zip(outputs, arrays)):
            root.mkdir(parents=True, exist_ok=True)
            count = int(arr.sum())
            filename = f'{row.PATCH_ID}_grass.tif' if count else ''
            if count:
                with rasterio.open(root / filename, 'w', driver='GTiff', count=1,
                                   width=500, height=500, dtype='uint8', nodata=0,
                                   crs=crs, transform=transform, compress='deflate') as dst:
                    dst.write(arr, 1)
            all_rows[index].append(dict(PATCH_ID=row.PATCH_ID, filename=filename,
                                       has_grass=bool(count), grass_pixels=count))
    if not any(row['grass_pixels'] for row in all_rows[0]):
        print('WARNING: no detection pixels inside the selected tiles and AOI. Check your grassland data/CRS.')
    for (_, manifest, _), rows in zip(outputs, all_rows):
        manifest.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(manifest, index=False)
    print(f'Created both grass-mask manifests for {len(grid)} tiles; grid and IDs unchanged.')
    print('Next: moweo validate --require-stats, then moweo process --workers 4')


if __name__ == '__main__':
    main()
