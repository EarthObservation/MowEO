import importlib

from .config import (
    EXPECTED_PATCHES,
    GRASS_ROOT,
    GRASS_MANIFEST,
    STATS_GRASS_ROOT,
    STATS_GRASS_MANIFEST,
    PATCH_GPKG,
    PATCH_LAYER,
    PROCESS_ROOT,
    ROOT,
    YEAR,
)

PACKAGES = [
    "numpy", "pandas", "scipy", "rasterio", "xarray", "h5netcdf", "h5py",
    "geopandas", "pyogrio", "shapely", "pyproj", "tqdm", "matplotlib",
    "attrs", "openeo",
]

print(f"YEAR: {YEAR}")
print(f"EXPECTED_PATCHES: {EXPECTED_PATCHES}")
print(f"ROOT: {ROOT}")
print(f"ROOT exists: {ROOT.exists()}")
print(f"PROCESS_ROOT: {PROCESS_ROOT}")
print(f"PROCESS_ROOT exists: {PROCESS_ROOT.exists()}")
print(f"GRASS_ROOT: {GRASS_ROOT}")
print(f"GRASS_ROOT exists: {GRASS_ROOT.exists()}")
print(f"Grass manifest: {GRASS_MANIFEST}")
print(f"Detection manifest exists: {GRASS_MANIFEST.exists()}")
print(f"STATS_GRASS_ROOT: {STATS_GRASS_ROOT}")
print(f"STATS_GRASS_ROOT exists: {STATS_GRASS_ROOT.exists()}")
print(f"Statistics manifest: {STATS_GRASS_MANIFEST}")
print(f"Statistics manifest exists: {STATS_GRASS_MANIFEST.exists()}")
print(f"National GPKG: {PATCH_GPKG}")
print(f"National GPKG exists: {PATCH_GPKG.exists()}")
print(f"Layer: {PATCH_LAYER or '<first layer>'}")
print("\nPackages:")

failed = []
for name in PACKAGES:
    try:
        module = importlib.import_module(name)
        version = getattr(module, "__version__", "unknown")
        print(f"  OK  {name:12s} {version}")
    except Exception as exc:
        failed.append((name, str(exc)))
        print(f"  ERR {name:12s} {exc}")

if failed:
    raise SystemExit("\nEnvironment is NOT ready. Fix packages marked ERR.")

print("\nEnvironment imports OK.")
