# Parameters and processing settings

`project_config.json` is private project configuration. Relative paths resolve
against the project directory. `project_config.example.json` supplies defaults;
its example year and EPSG are not requirements for all projects.

## Project and data

| Setting | Default | Meaning |
|---|---|---|
| `year` | 2021 in the template | Selected acquisition/detection year; initialization accepts the user's year. Verify backend coverage. |
| `expected_patches` | 1 in the template | Must equal the configured grid and both complete mask manifests; initialization sets this automatically. |
| `target_epsg` | 32633 | User-selected projected CRS in metres. |
| `patch_size_m` | 5000 | Fixed supported tile size. |
| `workspace_base` | `workspace` | Year-separated cubes, products, logs and provenance. |
| `processing_subdir` | `processing` | Output subdirectory within the year workspace; use a different name for detector migration while reusing acquisition. |
| `detector_profile` | `delivery_v0_3` | Required current detector profile; old projects must be migrated. |
| `patch_gpkg` | `data/EOPatch_PET.gpkg` | Full-square tile grid with unique IDs. |
| `patch_layer` | null | First GeoPackage layer unless specified. |
| `patch_selection_file` | `data/patch_ids.txt` | Optional acquisition-only subset; use `init-aoi --limit` for a persistent full-project subset. |
| `detection_grass_root`, `detection_grass_manifest` | `data/grass_detection`, its manifest CSV | Detection masks and complete tile inventory. |
| `stats_grass_root`, `stats_grass_manifest` | `data/grass_stats`, its manifest CSV | Inner statistics masks and complete tile inventory. |
| `study_area_name`, `output_prefix` | `Study area`, `AREA` | Report label and filename prefix. |
| `openeo_backend_url` | `https://openeo.dataspace.copernicus.eu` | Backend must provide the collection and SCL dilation process used by the code. |
| `openeo_parallel_jobs` | 2 | Maximum remote concurrency accepted by the local scheduler configuration. |

AOI initialization also records `aoi_file`, `aoi_source`, and optional country
selection metadata. These remain private, including any machine-specific paths.

## Seasonal windows

| Setting | Default |
|---|---|
| `load_start_mmdd`, `load_end_mmdd` | `03-15`, `10-31` |
| `mowing_start_mmdd`, `mowing_end_mmdd` | `04-01`, `10-31` |

Dates are combined with `year`. Review local phenology and backend endpoint
semantics. The configured within-year windows do not implement a cross-year
growing-season interface. Changing acquisition windows requires a fresh
workspace; existing products are not invalidated by a parameter fingerprint.

## Detector and evidence settings

| Setting | Default | Role |
|---|---|---|
| `drop_thr` | 0.093 | NDVI drop threshold. |
| `grow_r` | 1.3 | Multiplier for the regrowth threshold in the delivery reference. |
| `min_ndvi_valid` | null | Optional lower bound on valid NDVI observations. |
| `min_event_spacing_days` | 15 | Spacing rule for retained candidate events. |
| `rapid_rebound_days` | 5 | Rapid-rebound rejection window. |
| `rapid_rebound_threshold` | 0.15 | NDVI rebound threshold for that rule. |
| `date_tolerance_days` | 5 | Date tolerance in neighboring-event support. |
| `pre_start_days`, `pre_end_days` | 30, 5 | Spectral pre-event observation window. |
| `event_before_days`, `event_after_days` | 4, 10 | Spectral event-response window around the inferred date. |
| `child_date_tolerance_days` | 5 | Date tolerance for 10 m child-pixel support. |
| `max_events` | 6 | Fixed supported event-slot count, not a freely variable raster schema. |
| `baseline_version`, `final_version` | `v1_3`, `v0_3` | Required schema labels for the current profile; not a migration by renaming. |
| `event_confirm_pre_window_days`, `event_confirm_pre_min_obs` | 30, 2 | Pre-event support window and minimum observations. |
| `event_confirm_pre_range_max`, `event_confirm_start_median_tol` | 0.25, 0.06 | Pre-event stability evidence; not a hard rejection rule. |
| `local_obs_window_days` | 20 | Event-local observation-density window half-width. |

These defaults are project-calibrated parameters, not a universal mowing model.
Detector changes require deliberate recalculation and independent validation.
Resumable outputs are not automatically invalidated for every parameter change.

## Fixed cloud preprocessing

The collection is `SENTINEL2_L2A`, maximum scene cloud cover is 95, and the
`to_scl_dilation_mask` process uses `kernel1_size=17`, `kernel2_size=77`,
`mask1_values=[2,4,5,6,7]`, `mask2_values=[3,8,9,10,11]` and
`erosion_kernel_size=3`. These are implemented settings, not JSON options.
NDVI uses a 10 m grid; NDII and mean SWIR use an explicit 20 m grid.

Missing observations can reflect SCL rejection, dilation, NoData or acquisition
coverage and should not be interpreted solely as meteorological cloud cover.
The CLI records implemented cloud settings and actual submitted process graphs.

## Command-line controls

`init-aoi --limit N` fixes the project grid to at most N tiles. By contrast,
`acquire --limit N` limits only that acquisition invocation. `--parallel`
controls remote concurrency; `process --workers` controls local worker count.
`acquire` without `--run` and `process --dry-run` inspect work without submission
or processing respectively. `--skip-analysis` omits summary figures but retains
final products and exports. Keep one acquisition client per project workspace.
