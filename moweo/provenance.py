"""Record effective settings and installed versions without credentials."""
from datetime import datetime, timezone
import hashlib
from importlib.metadata import version, PackageNotFoundError
import json
import platform


def record(stage):
    from . import config as cfg
    libraries = {}
    for name in ('moweo', 'openeo', 'numpy', 'pandas', 'scipy', 'rasterio', 'xarray',
                 'geopandas', 'h5netcdf', 'h5py', 'pyproj', 'shapely'):
        try:
            libraries[name] = version(name)
        except PackageNotFoundError:
            libraries[name] = 'not installed'
    effective = dict(cfg.SETTINGS, year=cfg.YEAR)
    digest = hashlib.sha256(json.dumps(effective, sort_keys=True).encode()).hexdigest()
    out = cfg.YEAR_ROOT / 'provenance'
    out.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    data = dict(stage=stage, timestamp_utc=stamp, python=platform.python_version(),
                platform=platform.platform(), versions=libraries, settings=effective,
                settings_sha256=digest, backend=cfg.OPENEO_BACKEND_URL,
                collection='SENTINEL2_L2A', cloud_process='to_scl_dilation_mask',
                cloud_parameters=dict(kernel1_size=17, kernel2_size=77,
                                      mask1_values=[2, 4, 5, 6, 7], mask2_values=[3, 8, 9, 10, 11],
                                      erosion_kernel_size=3), max_cloud_cover=95)
    (out / f'{stamp}_{stage}.json').write_text(json.dumps(data, indent=2) + '\n', encoding='utf-8')
