"""Offline integration test: private inputs -> synthetic cubes -> final raster.

Run with: python -m unittest discover -s tests -v
No CDSE calls, credentials, credits or real grassland data are used.
"""
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
from shapely.geometry import box
import xarray as xr


class OfflineWorkflow(unittest.TestCase):
    def test_inputs_cubes_detection_and_mosaic(self):
        source = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as tmp:
            project = Path(tmp) / 'project'
            excluded = ['workspace', 'project_config.json', '__pycache__', 'build', '*.egg-info']
            if os.environ.get('MOWEO_TEST_INSTALLED') == '1':
                excluded.append('moweo')
            shutil.copytree(source, project, ignore=shutil.ignore_patterns(*excluded))

            def run(*args, success=True):
                commands = {'prepare_inputs.py': 'init', '01_prepare_acquisition.py': 'prepare',
                            '09_prepare_processing.py': 'validate', '03_acquire.py': 'acquire',
                            '14_run_processing.py': 'process', '13_mosaic.py': 'export'}
                result = subprocess.run([sys.executable, '-m', 'moweo', commands[args[0]], *args[1:]], cwd=project,
                                        capture_output=True, text=True)
                if success:
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                else:
                    self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
                return result

            # Two distant grasslands leave a gap in the final mosaic.
            # EPSG:32632 exercises a region outside the original zone 33 setup.
            vector = project / 'private.gpkg'
            gpd.GeoDataFrame(geometry=[box(500100, 5100100, 500140, 5100140),
                                      box(510100, 5100100, 510140, 5100140)],
                             crs=32632).to_file(vector, driver='GPKG')
            run('prepare_inputs.py', '--grasslands', str(vector), '--epsg', '32632', '--year', '2021')
            run('01_prepare_acquisition.py')
            run('09_prepare_processing.py', '--require-stats')
            settings = json.loads((project / 'project_config.json').read_text())
            self.assertEqual(settings['expected_patches'], 2)
            table = pd.read_csv(project / 'workspace/2021/acquisition/national_patches_2021_acquisition.csv')
            times = pd.date_range('2021-03-15', '2021-10-25', freq='10D')
            for tile_index, row in enumerate(table.itertuples()):
                dest = project / 'workspace/2021/acquisition/process_eopatch' / row.PATCH_ID
                dest.mkdir(parents=True)
                for res, bands, suffix in [(10, {'NDVI': 0.75}, 'NDVI'),
                                           (20, {'NDII': 0.4, 'SWIRMean': 1500}, 'NDII_SWIRMean_explicit20m')]:
                    size = 5000 // res
                    coords = {'t': times, 'y': row.NORTH - res / 2 - np.arange(size) * res,
                              'x': row.WEST + res / 2 + np.arange(size) * res}
                    ds = xr.Dataset({name: (('t', 'y', 'x'), np.full((len(times), size, size),
                                                                   value, dtype='float32'))
                                     for name, value in bands.items()}, coords=coords)
                    if tile_index == 1:
                        # One pronounced drop with slow regrowth in the second tile.
                        series = np.full(len(times), bands.get('NDVI', 0.4), dtype='float32')
                        if res == 10:
                            series[:8] = [0.55, 0.65, 0.8, 0.7, 0.35, 0.5, 0.65, 0.75]
                            ds['NDVI'].values[:] = series[:, None, None]
                        else:
                            ds['NDII'].values[4:6] = 0.1
                            ds['SWIRMean'].values[4:6] = 2200
                    ds.to_netcdf(dest / f'{row.PATCH_ID}_{suffix}_2021.nc', engine='h5netcdf')
            run('03_acquire.py')  # dry run must not need account authentication
            run('14_run_processing.py', '--workers', '1')
            event_tables = list((project / 'workspace/2021/processing/process_eopatch').glob('*/*FINAL_events*.csv'))
            self.assertEqual(len(event_tables), 2)
            frames = [pd.read_csv(path) for path in event_tables]
            events = pd.concat([frame for frame in frames if len(frame)], ignore_index=True)
            spectral_score = events.spectral_code.map({0: 0, 1: 50, 2: 50, 3: 100})
            expected_score = (0.30 * events.signal_quality_score + 0.25 * events.observation_density_score
                              + 0.20 * events.child_support + 0.25 * spectral_score)
            self.assertTrue(np.allclose(events.event_confidence, expected_score))
            checks = pd.read_csv(project / 'workspace/2021/processing/national_analysis_2021/tables/event_scan_checks_2021.csv')
            self.assertEqual(set(checks.STATUS), {'OK'})
            self.assertTrue((checks.EXPECTED == checks.SCANNED).all())
            mosaic = project / 'workspace/2021/processing/national_maps/AREA_MOWING_COUNT_2021.tif'
            self.assertTrue(mosaic.with_name('AREA_TILE_QA_2021.gpkg').exists())
            self.assertTrue(list((project / 'workspace/2021/provenance').glob('*.json')))
            with rasterio.open(mosaic) as ds:
                values = ds.read(1)
                self.assertEqual(ds.crs.to_epsg(), 32632)
                self.assertEqual(ds.nodata, 255)
                self.assertEqual(set(np.unique(values)), {0, 1, 255})
                self.assertTrue((values[:, 500:1000] == 255).all())
            manifest = json.loads((mosaic.parent / 'products_manifest.json').read_text())
            self.assertEqual(len(manifest), 12)
            for product in manifest:
                with rasterio.open(mosaic.parent / product['file']) as src:
                    self.assertEqual(src.dtypes[0], product['dtype'])
                    self.assertEqual(src.nodata, product['nodata'])
                    self.assertTrue((src.read(1)[:, 500:1000] == product['nodata']).all())
            with rasterio.open(mosaic.with_name('AREA_MOWING_DOY_1_2021.tif')) as src:
                self.assertEqual(set(np.unique(src.read(1))), {0, times[4].dayofyear, 65535})
            with rasterio.open(mosaic.with_name('AREA_SCL_MASKED_PERCENT_2021.tif')) as src:
                self.assertEqual(set(np.unique(src.read(1))), {-9999, 0})
            # A missing final product must fail rather than publish a partial mosaic.
            next((project / 'workspace/2021/processing/process_eopatch').glob('*/*MOWING_FINAL_v*.tif')).unlink()
            run('13_mosaic.py', success=False)


if __name__ == '__main__':
    unittest.main()
