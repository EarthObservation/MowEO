"""AOI/country -> fixed grid -> mask-independent download -> later private masks."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from shapely.geometry import box, Point


class AOIWorkflow(unittest.TestCase):
    def test_country_grid_and_masks_added_later(self):
        source = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as tmp:
            project = Path(tmp) / 'project'
            ignored = ['build', '*.egg-info', '__pycache__', 'workspace', 'project_config.json']
            if os.environ.get('MOWEO_TEST_INSTALLED') == '1':
                ignored.append('moweo')
            shutil.copytree(source, project, ignore=shutil.ignore_patterns(*ignored))

            def run(*args, success=True):
                result = subprocess.run([sys.executable, '-m', 'moweo', *args], cwd=project,
                                        capture_output=True, text=True)
                self.assertEqual(result.returncode == 0, success, (result.stdout + result.stderr)[-5000:])

            countries = project / 'countries.gpkg'
            aoi = box(500100, 5100100, 509900, 5100300)
            gpd.GeoDataFrame({'NAME': ['Testland', 'Otherland']},
                             geometry=[aoi, box(520000, 5100000, 525000, 5105000)],
                             crs=32633).to_file(countries, driver='GPKG')
            run('init-aoi', '--aoi', str(countries), '--country', 'testland',
                '--epsg', '32633', '--year', '2026')
            grid = gpd.read_file(project / 'data/EOPatch_PET.gpkg')
            before = grid.PATCH_ID.tolist()
            self.assertEqual(set(before), {'eo_32633_500000_5100000', 'eo_32633_505000_5100000'})
            self.assertTrue(np.allclose(grid.geometry.area, 25000000))
            self.assertTrue(grid.AREA_PCT.between(0, 100).all())
            run('prepare')
            run('acquire')  # no masks, authentication or credits needed for dry run
            run('validate', '--require-stats', success=False)
            # Test persistent subset selection, with the same stable first tile ID.
            limited = project / 'limited'
            run('--project', str(limited), 'init-aoi', '--aoi', str(countries),
                '--country', 'Testland', '--epsg', '32633', '--year', '2026', '--limit', '1')
            limited_grid = gpd.read_file(limited / 'data/EOPatch_PET.gpkg')
            self.assertEqual(limited_grid.PATCH_ID.tolist(), [before[0]])
            self.assertEqual(json.loads((limited / 'project_config.json').read_text())['expected_patches'], 1)
            grass = project / 'private_grass.gpkg'
            gpd.GeoDataFrame(geometry=[box(499900, 5100000, 510100, 5100400)],
                             crs=32633).to_file(grass, driver='GPKG')
            run('masks', '--grasslands', str(grass))
            run('validate', '--require-stats')
            self.assertEqual(gpd.read_file(project / 'data/EOPatch_PET.gpkg').PATCH_ID.tolist(), before)
            for pid in before:
                with rasterio.open(project / f'data/grass_detection/{pid}_grass.tif') as det:
                    arr = det.read(1)
                    self.assertEqual(set(np.unique(arr)), {0, 1})
                    rows, cols = np.where(arr == 1)
                    xs, ys = rasterio.transform.xy(det.transform, rows, cols)
                    self.assertTrue(all(aoi.covers(Point(x, y)) for x, y in zip(xs, ys)))
                with rasterio.open(project / f'data/grass_stats/{pid}_grass.tif') as stat:
                    self.assertFalse(np.any((stat.read(1) == 1) & (arr != 1)))
            manifest = pd.read_csv(project / 'data/grass_detection/grass_masks_manifest.csv')
            self.assertEqual(len(manifest), 2)
            # Never silently overwrite masks already used by a detector run.
            run('masks', '--grasslands', str(grass), success=False)


if __name__ == '__main__':
    unittest.main()
