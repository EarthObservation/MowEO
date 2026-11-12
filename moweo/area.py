"""Initialize a project from an AOI or a country selected from a boundary file."""
import argparse
import json
import math
from pathlib import Path

import geopandas as gpd
from pyproj import CRS
from shapely.geometry import box
from shapely.ops import unary_union


def metric_crs(epsg):
    crs = CRS.from_epsg(epsg)
    if not crs.is_projected or any(axis.unit_name != 'metre' for axis in crs.axis_info):
        raise ValueError('Choose a projected CRS with metre units; not EPSG:4326.')
    return crs


def polygons(path, layer=None):
    frame = gpd.read_file(path, layer=layer)
    if frame.empty or frame.crs is None:
        raise ValueError('The polygon file must be nonempty and have a CRS.')
    if frame.geometry.isna().any() or frame.geometry.is_empty.any():
        raise ValueError('Remove missing or empty geometries from the source.')
    if not frame.geom_type.isin(['Polygon', 'MultiPolygon']).all() or not frame.is_valid.all():
        raise ValueError('Use valid Polygon/MultiPolygon features; repair the source first.')
    return frame


def build_tiles(area, epsg, limit=None):
    """Full square tiles with positive-area overlap; stable CRS/coordinate IDs."""
    geometries = list(area.geoms) if area.geom_type == 'MultiPolygon' else [area]
    candidates = set()
    for geom in geometries:
        minx, miny, maxx, maxy = geom.bounds
        for west in range(math.floor(minx / 5000) * 5000, math.ceil(maxx / 5000) * 5000, 5000):
            for south in range(math.floor(miny / 5000) * 5000, math.ceil(maxy / 5000) * 5000, 5000):
                if (west, south) in candidates:
                    continue
                tile = box(west, south, west + 5000, south + 5000)
                if geom.intersection(tile).area > 0:
                    candidates.add((west, south))
    rows = []
    for west, south in sorted(candidates):
        tile = box(west, south, west + 5000, south + 5000)
        rows.append(dict(PATCH_ID=f'eo_{epsg}_{west}_{south}', PET='all', BIOGEO_R='all',
                         AREA_PCT=100 * area.intersection(tile).area / tile.area, geometry=tile))
        if limit and len(rows) >= limit:
            break
    return rows, len(candidates)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--aoi', required=True, type=Path, help='Polygon AOI or boundary dataset (GeoPackage, GeoJSON, Shapefile).')
    p.add_argument('--layer')
    p.add_argument('--country', help='Exact country name/code from a column in the supplied boundary file.')
    p.add_argument('--country-field', default='NAME', help='Column to filter when --country is given (default NAME).')
    p.add_argument('--epsg', required=True, type=int, help='Suitable metric projected CRS for the selected area.')
    p.add_argument('--year', required=True, type=int)
    p.add_argument('--limit', type=int, help='Keep only the first N tiles in coordinate order, including for later processing.')
    args = p.parse_args()
    project = Path.cwd()
    if (project / 'project_config.json').exists() or (project / 'data/EOPatch_PET.gpkg').exists():
        p.error('Existing project found; use a new project directory for a different AOI.')
    if not 2017 <= args.year <= 2100 or (args.limit is not None and args.limit < 1):
        p.error('Use year >= 2017 and a positive --limit.')
    crs = metric_crs(args.epsg)
    frame = polygons(args.aoi, args.layer)
    if args.country is not None:
        if args.country_field not in frame:
            p.error(f'Country column {args.country_field!r} not found; available: {list(frame.columns)}')
        frame = frame[frame[args.country_field].astype(str).str.casefold() == args.country.casefold()]
        if frame.empty:
            p.error('No country matched; use the name/code as stored in the boundary file.')
    area = unary_union(frame.to_crs(crs).geometry)
    rows, total = build_tiles(area, args.epsg, args.limit)
    if not rows:
        p.error('The AOI produced no tiles.')
    data = project / 'data'
    data.mkdir(parents=True, exist_ok=True)
    gpd.GeoDataFrame(geometry=[area], crs=crs).to_file(data / 'aoi.gpkg', layer='aoi', driver='GPKG')
    gpd.GeoDataFrame(rows, crs=crs).to_file(data / 'EOPatch_PET.gpkg', driver='GPKG')
    settings = json.loads(Path(__file__).with_name('project_config.example.json').read_text(encoding='utf-8'))
    settings.update(year=args.year, target_epsg=args.epsg, expected_patches=len(rows),
                    study_area_name=args.country or 'User AOI', output_prefix='AREA',
                    aoi_file='data/aoi.gpkg', aoi_source=str(args.aoi.resolve()),
                    country_selection=args.country, country_field=args.country_field if args.country else None)
    (project / 'project_config.json').write_text(json.dumps(settings, indent=2) + '\n', encoding='utf-8')
    print(f'AOI intersects {total} tiles; configured {len(rows)} full 5 km tiles.')
    print('No grassland masks needed for acquisition. Next: moweo prepare, then moweo acquire.')
    print('Before local detection: moweo masks --grasslands /path/to/grasslands.gpkg')


if __name__ == '__main__':
    main()
