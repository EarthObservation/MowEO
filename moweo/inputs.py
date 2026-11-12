"""Build portable 5 km tiles and private masks from a user's grassland polygons."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from pyproj import CRS
from rasterio.features import rasterize
from rasterio.transform import from_origin
from shapely.geometry import box

HERE = Path.cwd().resolve()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--grasslands', required=True, type=Path,
                   help='Private polygon vector file: ALL features must represent grassland.')
    p.add_argument('--layer', default=None)
    p.add_argument('--epsg', required=True, type=int, help='Local metric projected CRS, normally the local UTM zone.')
    p.add_argument('--year', required=True, type=int)
    args = p.parse_args()
    if not 2017 <= args.year <= 2100:
        p.error('Choose a Sentinel-2 L2A year from 2017 onward; verify backend coverage.')
    crs = CRS.from_epsg(args.epsg)
    if not crs.is_projected or any(a.unit_name != 'metre' for a in crs.axis_info):
        p.error('--epsg must be a projected CRS with metre units.')
    config_path = HERE / 'project_config.json'
    data = HERE / 'data'
    if config_path.exists() or (data / 'EOPatch_PET.gpkg').exists():
        p.error('Existing project found. Use a separate repository copy for a new area.')
    grass = gpd.read_file(args.grasslands, layer=args.layer)
    if grass.crs is None or grass.empty:
        p.error('Input must have a CRS and at least one grassland polygon.')
    if not grass.geom_type.isin(['Polygon', 'MultiPolygon']).all() or not grass.is_valid.all():
        p.error('Input must contain valid Polygon/MultiPolygon geometries. Repair the source first.')
    grass = grass.to_crs(crs)
    # Buffer complete source features before clipping to tiles to avoid tile-edge seams.
    detection = grass.geometry.buffer(5)
    stats = grass.geometry.buffer(-7.5)
    minx, miny, maxx, maxy = detection.total_bounds
    tiles, manifests = [], {'grass_detection': [], 'grass_stats': []}
    for west in range(math.floor(minx / 5000) * 5000, math.ceil(maxx / 5000) * 5000, 5000):
        for south in range(math.floor(miny / 5000) * 5000, math.ceil(maxy / 5000) * 5000, 5000):
            tile = box(west, south, west + 5000, south + 5000)
            nearby = detection.sindex.query(tile, predicate='intersects')
            if not len(nearby):
                continue
            pid = f'tile_{west}_{south}'
            tiles.append({'PATCH_ID': pid, 'PET': 'all', 'BIOGEO_R': 'all',
                          'AREA_PCT': 100.0, 'geometry': tile})
            transform = from_origin(west, south + 5000, 10, 10)
            for folder, shapes in [('grass_detection', detection), ('grass_stats', stats)]:
                geometries = [geom for geom in shapes.iloc[nearby] if not geom.is_empty]
                arr = rasterize(((geom, 1) for geom in geometries), out_shape=(500, 500),
                                transform=transform, fill=0, dtype='uint8', all_touched=False) if geometries else np.zeros((500, 500), dtype='uint8')
                count = int(arr.sum())
                filename = f'{pid}_grass.tif' if count else ''
                out = data / folder
                out.mkdir(parents=True, exist_ok=True)
                if count:
                    with rasterio.open(out / filename, 'w', driver='GTiff', width=500, height=500,
                                       count=1, dtype='uint8', crs=crs, transform=transform,
                                       nodata=0, compress='deflate') as dst:
                        dst.write(arr, 1)
                manifests[folder].append({'PATCH_ID': pid, 'filename': filename,
                                          'has_grass': bool(count), 'grass_pixels': count})
    if not tiles:
        p.error('No intersecting tiles found.')
    gpd.GeoDataFrame(tiles, crs=crs).to_file(data / 'EOPatch_PET.gpkg', driver='GPKG')
    for folder, rows in manifests.items():
        pd.DataFrame(rows).to_csv(data / folder / 'grass_masks_manifest.csv', index=False)
    settings = json.loads(Path(__file__).with_name('project_config.example.json').read_text(encoding='utf-8'))
    settings.update(year=args.year, expected_patches=len(tiles), target_epsg=args.epsg,
                    study_area_name='Study area', output_prefix='AREA')
    config_path.write_text(json.dumps(settings, indent=2) + '\n', encoding='utf-8')
    print(f'Created {len(tiles)} tiles, private masks, manifests and {config_path.name}.')
    print('Next: python 01_prepare_acquisition.py')


if __name__ == '__main__':
    main()
