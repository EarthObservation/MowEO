# Release notes

## 0.2.0a2 — prerelease

- Installable Python package with the `moweo` command and modular acquisition,
  detection, finalization, analysis and export stages.
- Connected acquisition-to-processing execution, with separate stages retained.
- AOI and country selection from supplied vector boundaries.
- Full 5 km tiles with deterministic CRS/coordinate IDs and persistent tile limits.
- Mask-independent acquisition and later grassland rasterization on the fixed grid.
- Private detection/inner statistics masks with binary, alignment and subset checks.
- GeoTIFF mosaics, event CSVs and tile-level GeoPackage summaries.
- Runtime configuration/version records and saved openEO process graphs.
- Complete CSV schema for tiles without detected events, optional cost metadata,
  text report groups and explicit rejection of incomplete count mosaics.
- English documentation, command-line output and plot labels.

## Verification

Both offline tests passed from source and against an installed 0.2.0a2 wheel
on Windows/Python 3.13.5. The installed-wheel run completed in 70.82 seconds.
They cover AOI/country preparation and later masks, plus a complete synthetic
local detection/report/export example. They use no CDSE credentials or credits.

Verified dependencies included openEO 0.52.0, NumPy 2.5.4, pandas 2.3.3,
rasterio 1.4.4, xarray 2025.1.1 and GeoPandas 1.2.0. Rasterio emitted a NumPy 2.5
deprecation warning; the tests passed. Later editorial changes translate the
remaining plot labels and reorganize documentation without changing detector
logic. Python syntax, local documentation links and archive contents are checked
for the publication package.

The supplied Conda environment specifies Python 3.11 but has not been executed
here. CI targets Python 3.11 on Windows and Linux and has not yet run on GitHub.
macOS has not been tested. Offline tests establish software behavior, not
real-world mowing accuracy. No completed real-data validation of this exact
prerelease is claimed. External backend availability is not tested by CI.

## Distribution

Private AOIs, masks, satellite cubes, credentials, reference observations and
results are not distributed. License selection, authorship/citation metadata,
public GitHub publication and DOI archival remain pending. See
[release status](docs/release_status.md).
