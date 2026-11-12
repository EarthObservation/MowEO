# Release notes

## 0.3.0a1 — delivery-reference prerelease

- Ported the delivered v0.3 NDVI/event-quality core, with `grow_r=1.3` and
  event-local signal/density evidence in confidence. This changes results from
  MowEO 0.2.x and requires local recalculation, not relabeling old products.
- Matched the delivery archive's twelve SI raster products: raw, SEIVE and
  spatially harmonized counts; DOY 1–6; rounded UInt8 confidence; valid
  observations; unavailable/SCL-masked percentage.
- Kept original DOY dates distinct from the filtered count. Retained the source
  per-tile, multivalued-date SEIVE and confidence-aware harmonization rules.
- Added profile/schema migration guards, a separate processing subdirectory
  for reuse of downloaded cubes, output staging, NoData-gap initialization and
  a machine-readable twelve-product inventory.
- Fixed analysis to select event tables by configured version and stop if the
  event scan is incomplete or inconsistent with QA.
- Source provenance records imported software version and core/config/export
  hashes, avoiding confusion with another wheel installed in the environment.
- Tests cover event evidence weights, postprocessing behavior and every product's
  data type, NoData gaps and event-date/observation semantics.

The reference core and export/config files were matched across two identical
local source copies and the 2017–2025 delivery archive inventory/metadata. The
detector, quality, spectral, SEIVE and harmonization functions were compared
structurally with the reference implementation. See [products](docs/products.md).

### Verification of this update

All three offline tests passed from source (92.10 s) and from an installed
0.3.0a1 wheel (84.00 s) on Windows/Python 3.13.5. Checks include all twelve
product schemas/NoData gaps, event-local evidence weights, SEIVE small-region
behavior, harmonization of low-confidence pixels and consistency of scanned
event tables with QA. Dependencies match the verification environment below.

The existing 2026 ten-tile pilot was recalculated locally with v0.3 using its
twenty already downloaded index cubes and original masks; no additional cloud
jobs were submitted. All ten tile QA rows and event-scan checks were OK and
all twelve products were written in a separate processing directory. SEIVE and
harmonization rasters were compared pixel-for-pixel with independently executed
functions from the reference export script on the same final tile rasters.
Rounded confidence was also compared pixel-for-pixel. These comparisons passed.
This establishes implementation/output consistency, not independent scientific
accuracy or equivalence of detections between the old and new profiles.

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
