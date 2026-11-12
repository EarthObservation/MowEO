# Delivery-reference detector and raster products

MowEO 0.3.0a1 adopts the `delivery_v0_3` detector profile and the twelve raster
products identified in the latest local national delivery archive. This is a
methodological update from MowEO 0.2.x, not just additional filenames.

## Reference audit

The delivery archive inventory and manifest identify count, SEIVE, harmonized
count, six event-date layers, confidence, valid observations and masked/unavailable
percentage. The delivered 2025 raster metadata identify detector version v0.3.
The two available local v0.3 source copies were byte-identical for core,
configuration and mosaicking. Reference file SHA-256 digests:

| Reference file | SHA-256 |
|---|---|
| `mowing_core.py` | `fb6af459c5f8fc87a7e74ac0a7c4d122ce0192efcbfa4a69e211a104da01715f` |
| `config.py` | `8fe5b93ee252e01a57aa459c5c51a5f08a6eeebde41889ac6f6a64b8490a492e` |
| `06_mosaic_national_products.py` | `63043e72d13a209c444bf710da38c006e4bacdc3c5af4cc1ae1ea3b8c697d854` |

The AST representations of the ported NDVI detector, event-context/quality
functions, temporal postprocessing, spectral response, majority rule, spatial
harmonization and SEIVE functions were checked against these references.
Adaptations cover packaging, configurable paths/prefixes, input guards, staged
export, explicit NoData gaps and the delivery's rounded confidence encoding.
Source data and client rasters are not distributed with the software.

The reference notes call v0.3 a sensitivity version. Its use in a delivery does
not establish universal accuracy or calibrated probabilities. Reuse still
requires regional validation and confirmed software/data rights.

## Detector changes from MowEO 0.2.x

- `grow_r=1.3` instead of 1.4; the default recovery threshold is `0.093*1.3`.
- Event-local observation gaps/density, recovery and pre-event stability are
  recorded. Pre-event stability is evidence, not a hard candidate rejection.
- Pre-event defaults: 30 days, at least two observations, maximum NDVI range
  0.25 and start-to-median tolerance 0.06. The local density window is +/-20 days.
- Event evidence score is `0.30*signal_quality + 0.25*observation_density +
  0.20*child_support + 0.25*spectral_score`, with weights renormalized over finite
  components. Spectral scores are 0, 50 or 100 for neither, one or both directions.
- Baselines use `v1_3`; final/diagnostic rasters use `v0_3`. The final raster
  retains 17 bands and the companion diagnostics retain 44 bands. The baseline
  adds event-local quality bands; event CSVs include the additional metrics.

## Twelve raster products

Files are written under `workspace/<YEAR>/<processing_subdir>/national_maps/`.
`<PREFIX>` defaults to `AREA` and is configurable. All products use the same
10 m grid and the detection mask. `products_manifest.json` records the inventory.

| File pattern | Type / NoData | Meaning |
|---|---|---|
| `<PREFIX>_MOWING_COUNT_<YEAR>.tif` | UInt8 / 255 | Original event count, 0–6. |
| `<PREFIX>_MOWING_COUNT_SEIVE_<YEAR>.tif` | UInt8 / 255 | Count recomputed after filtering each event-date layer. |
| `<PREFIX>_MOWING_COUNT_SPATIAL_HARMONIZED_<YEAR>.tif` | UInt8 / 255 | Confidence-aware spatial harmonization of the original count. |
| `<PREFIX>_MOWING_DOY_1_<YEAR>.tif` through `_DOY_6_` | UInt16 / 65535 | Six original event-date slots; 1–366 is day of year, 0 means no event in the slot. |
| `<PREFIX>_MOWING_CONFIDENCE_<YEAR>.tif` | UInt8 / 255 | Mean event evidence score, rounded to an integer 0–100. Undefined means, including pixels without events, are NoData. |
| `<PREFIX>_VALID_OBSERVATIONS_<YEAR>.tif` | UInt16 / 65535 | Valid daily NDVI observations per pixel. |
| `<PREFIX>_SCL_MASKED_PERCENT_<YEAR>.tif` | Float32 / -9999 | `100*(1-VALID_OBS/N_UNIQUE_NDVI_DATES)`, clipped to 0–100. |

Three count rasters + six date rasters + confidence + observations + percentage
make twelve products. Event CSVs, multi-band tile rasters, analysis figures,
tile QA CSV and tile-summary GeoPackage remain available as additional outputs.

## SEIVE semantics

The historical spelling `SEIVE` is retained for compatibility with the delivery.
Each original `MOWING_DOY_1..6` integer date layer is processed with Rasterio/GDAL
`sieve(size=25, connectivity=8, mask=detection_mask)`. This is **multivalued date
region** filtering, not binary event filtering and not filtering the count raster.
Small regions may be reassigned to neighboring values. The count is then the
number of filtered date slots with values in 1–366. At 10 m, 25 pixels represent
0.25 ha in the projected grid.

Filtering occurs independently within each tile, matching the reference code;
connectivity does not cross tile edges. The exported DOY 1–6 rasters remain
**original, unfiltered dates**, as in the delivery. Thus raw DOY 1 and SEIVE count
are not guaranteed to describe the same filtered event sequence. Filtered date
rasters or recomputed post-SEIVE confidence are not claimed as delivery products.

## Harmonization and quality semantics

Only original-count pixels with missing confidence or confidence below 60 can
change. A 3-by-3 local majority needs at least five agreeing valid pixels. Then
low-confidence pixels can be assigned the dominant count of their 8-connected
grass component within the tile. Components are not field/parcel boundaries.
The original count is retained separately. Harmonization is not SEIVE and does
not imply independently verified accuracy.

`SCL_MASKED_PERCENT` measures unavailable dates after preprocessing, not pure
cloud cover: geometric coverage, invalid data and SCL/dilation rejection can
all contribute. Confidence is an evidence score, not a probability.
The comparison archive's aligned GRAME/CLMS data are external products and
are not generated or redistributed by this workflow.

## Migrating a 0.2.x project

Retain the existing acquisition cubes and masks. Create a separate project
configuration using the current example and update its paths, year and tile
count. Select `detector_profile=delivery_v0_3`, `grow_r=1.3`,
`baseline_version=v1_3`, `final_version=v0_3` and a new `processing_subdir`,
for example `processing_v0_3`. `workspace_base` can point to the old workspace
to reuse its unchanged acquisition inputs. Recalculate with `moweo process`.

Do not merely rename old products or change their version tags. MowEO rejects
old profile/version configurations and exports require v0.3 final raster tags.
Parameter tuning still requires intentional recalculation; not every configurable
threshold is fingerprinted in resumable products. Keep the original results
and provenance for comparison. No repeated cloud acquisition is needed when
the original cubes have the same grid, year and temporal coverage.
