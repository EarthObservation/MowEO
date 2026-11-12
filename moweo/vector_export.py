"""Export tile-level QA summaries as a GeoPackage, not event polygons."""
import geopandas as gpd
import pandas as pd
from . import config as cfg


def main():
    grid = gpd.read_file(cfg.PATCH_GPKG, layer=cfg.PATCH_LAYER)
    from .processing_prepare import _patch_id
    grid['PATCH_ID'] = _patch_id(grid)
    qa = pd.read_csv(cfg.NATIONAL_QA, dtype={'PATCH_ID': str})
    summary = grid[['PATCH_ID', 'geometry']].merge(qa, on='PATCH_ID', validate='one_to_one')
    output = cfg.ROOT / 'national_maps' / f'{cfg.OUTPUT_PREFIX}_TILE_QA_{cfg.YEAR}.gpkg'
    output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_file(output, layer='tile_qa', driver='GPKG')
    print('Tile-level summaries:', output)


if __name__ == '__main__':
    main()
