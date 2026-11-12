from __future__ import annotations

import json
import os
from pathlib import Path

PROJECT_DIR = Path(os.environ.get("MOWEO_PROJECT_DIR", Path.cwd())).resolve()
SETTINGS_FILE = PROJECT_DIR / "project_config.json"
SETTINGS_EXAMPLE_FILE = PROJECT_DIR / "project_config.example.json"


def _load_settings() -> dict:
    if not SETTINGS_FILE.exists():
        raise RuntimeError(
            f"Project is not configured: {SETTINGS_FILE}\n"
            "Run: python 00_setup_project.py"
        )
    data = json.loads(SETTINGS_FILE.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise RuntimeError("project_config.json must contain a JSON object.")
    return data


def _path(value: str | Path) -> Path:
    p = Path(value)
    return p if p.is_absolute() else (PROJECT_DIR / p).resolve()


SETTINGS = _load_settings()

YEAR = int(os.environ.get("MOWING_YEAR", SETTINGS["year"]))
EXPECTED_PATCHES = int(SETTINGS.get("expected_patches", 1))
STUDY_AREA_NAME = str(SETTINGS.get("study_area_name", "Study area"))
OUTPUT_PREFIX = str(SETTINGS.get("output_prefix", "AREA"))
TARGET_EPSG = int(SETTINGS.get("target_epsg", 32633))
PATCH_SIZE_M = float(SETTINGS.get("patch_size_m", 5000.0))
NDVI_RES_M = float(SETTINGS.get("ndvi_resolution_m", 10.0))
INDICES20_RES_M = float(SETTINGS.get("indices20_resolution_m", 20.0))
NDVI_WIDTH = int(round(PATCH_SIZE_M / NDVI_RES_M))
NDVI_HEIGHT = NDVI_WIDTH
INDICES20_WIDTH = int(round(PATCH_SIZE_M / INDICES20_RES_M))
INDICES20_HEIGHT = INDICES20_WIDTH

WORKSPACE_BASE = _path(SETTINGS.get("workspace_base", "workspace"))
YEAR_ROOT = WORKSPACE_BASE / str(YEAR)
ACQUISITION_ROOT = YEAR_ROOT / "acquisition"
RAW_PROCESS_ROOT = ACQUISITION_ROOT / "process_eopatch"
PROCESSING_SUBDIR = str(SETTINGS.get('processing_subdir', 'processing'))
import re
if not re.fullmatch(r'[A-Za-z0-9_-]+', PROCESSING_SUBDIR) or PROCESSING_SUBDIR == 'acquisition':
    raise ValueError('processing_subdir must be a simple directory name distinct from acquisition.')
ROOT = YEAR_ROOT / PROCESSING_SUBDIR
PROCESS_ROOT = ROOT / "process_eopatch"

DATA_DIR = PROJECT_DIR / "data"
PATCH_GPKG = _path(SETTINGS.get("patch_gpkg", "data/EOPatch_PET.gpkg"))
PATCH_LAYER = SETTINGS.get("patch_layer") or None
PATCH_SELECTION_FILE = _path(SETTINGS.get("patch_selection_file", "data/patch_ids.txt"))

GRASS_ROOT = _path(SETTINGS.get("detection_grass_root", "data/grass_detection"))
GRASS_MANIFEST = _path(
    SETTINGS.get("detection_grass_manifest", str(GRASS_ROOT / "grass_masks_manifest.csv"))
)
STATS_GRASS_ROOT = _path(SETTINGS.get("stats_grass_root", "data/grass_stats"))
STATS_GRASS_MANIFEST = _path(
    SETTINGS.get("stats_grass_manifest", str(STATS_GRASS_ROOT / "grass_masks_manifest.csv"))
)

PATCH_CSV_ACQUISITION = ACQUISITION_ROOT / f"national_patches_{YEAR}_acquisition.csv"
PATCH_CSV = ROOT / f"national_patches_{YEAR}.csv"
OPENEO_LOG = ACQUISITION_ROOT / f"openeo_jobs_{YEAR}.csv"

AUTH_DIR = YEAR_ROOT / ".openeo_auth"
AUTH_TOKEN_FILE = AUTH_DIR / "refresh_tokens.json"
AUTH_ACCOUNT_FILE = AUTH_DIR / "account.json"
OPENEO_BACKEND_URL = str(
    SETTINGS.get("openeo_backend_url", "https://openeo.dataspace.copernicus.eu")
)
OPENEO_MAX_PARALLEL_JOBS = int(SETTINGS.get("openeo_parallel_jobs", 2))

LOAD_START = f"{YEAR}-{SETTINGS.get('load_start_mmdd', '03-15')}"
LOAD_END = f"{YEAR}-{SETTINGS.get('load_end_mmdd', '10-31')}"
MOWING_START = f"{YEAR}-{SETTINGS.get('mowing_start_mmdd', '04-01')}"
MOWING_END = f"{YEAR}-{SETTINGS.get('mowing_end_mmdd', '10-31')}"

PATCH_ID_FIELD = "PATCH_ID"
TYPE_FIELD = "TYPE"
NOTE_FIELD = "NOTE"
MAX_EVENTS = int(SETTINGS.get("max_events", 6))
if PATCH_SIZE_M != 5000 or NDVI_RES_M != 10 or INDICES20_RES_M != 20 or MAX_EVENTS != 6:
    raise ValueError('This implementation requires 5 km tiles, 10/20 m grids and max_events=6.')
if EXPECTED_PATCHES < 1 or not 2017 <= YEAR <= 2100:
    raise ValueError('Set a positive expected_patches and a Sentinel-2 L2A year (2017 onward).')
import re
if not re.fullmatch(r'[A-Za-z0-9_-]+', OUTPUT_PREFIX):
    raise ValueError('output_prefix must contain only letters, digits, underscores and hyphens.')

# Detector parameters frozen after the project pilot.
DROP_THR = float(SETTINGS.get("drop_thr", 0.093))
GROW_R = float(SETTINGS.get("grow_r", 1.3))
MIN_NDVI_VALID = SETTINGS.get("min_ndvi_valid", None)
MIN_EVENT_SPACING_DAYS = int(SETTINGS.get("min_event_spacing_days", 15))
RAPID_REBOUND_DAYS = int(SETTINGS.get("rapid_rebound_days", 5))
RAPID_REBOUND_THRESHOLD = float(SETTINGS.get("rapid_rebound_threshold", 0.15))
DATE_TOLERANCE = int(SETTINGS.get("date_tolerance_days", 5))

PRE_START_DAYS = int(SETTINGS.get("pre_start_days", 30))
PRE_END_DAYS = int(SETTINGS.get("pre_end_days", 5))
EVENT_BEFORE_DAYS = int(SETTINGS.get("event_before_days", 4))
EVENT_AFTER_DAYS = int(SETTINGS.get("event_after_days", 10))
CHILD_DATE_TOLERANCE = int(SETTINGS.get("child_date_tolerance_days", 5))

BASELINE_VERSION = str(SETTINGS.get("baseline_version", "v1_3"))
FINAL_VERSION = str(SETTINGS.get("final_version", "v0_3"))
if SETTINGS.get('detector_profile') != 'delivery_v0_3' or BASELINE_VERSION != 'v1_3' or FINAL_VERSION != 'v0_3':
    raise ValueError('This release requires detector_profile=delivery_v0_3, baseline_version=v1_3 and final_version=v0_3. Migrate an older project into a separate processing directory before recalculation.')
EVENT_CONFIRM_PRE_WINDOW_DAYS = int(SETTINGS.get('event_confirm_pre_window_days', 30))
EVENT_CONFIRM_PRE_MIN_OBS = int(SETTINGS.get('event_confirm_pre_min_obs', 2))
EVENT_CONFIRM_PRE_RANGE_MAX = float(SETTINGS.get('event_confirm_pre_range_max', 0.25))
EVENT_CONFIRM_START_MEDIAN_TOL = float(SETTINGS.get('event_confirm_start_median_tol', 0.06))
LOCAL_OBS_WINDOW_DAYS = int(SETTINGS.get('local_obs_window_days', 20))

PROCESS_LOG = Path(os.environ.get("MOWING_PROCESS_LOG", str(ROOT / f"national_processing_{YEAR}.csv")))
FINALIZE_LOG = Path(os.environ.get("MOWING_FINALIZE_LOG", str(ROOT / f"national_finalize_{YEAR}.csv")))
NATIONAL_QA = Path(os.environ.get("MOWING_NATIONAL_QA", str(ROOT / f"national_qa_{YEAR}.csv")))
OVERFLOW_ALL = Path(os.environ.get("MOWING_OVERFLOW_ALL", str(ROOT / f"national_overflow_gt6_detail_{YEAR}.csv")))
PILOT_SOURCE_ROOT = ROOT / "pilot_unused"

# Compatibility names from the older acquisition-only package.
OUTPUT_ROOT = ACQUISITION_ROOT
EXPECTED_USER_FILE = AUTH_ACCOUNT_FILE


def patch_dir(patch_id: str) -> Path:
    return PROCESS_ROOT / str(patch_id)


def raw_patch_dir(patch_id: str) -> Path:
    return RAW_PROCESS_ROOT / str(patch_id)


def grass_path(patch_id: str) -> Path:
    return GRASS_ROOT / f"{patch_id}_grass.tif"


def stats_grass_path(patch_id: str) -> Path:
    return STATS_GRASS_ROOT / f"{patch_id}_grass.tif"


def acquisition_patch_paths(patch_id: str) -> dict[str, Path]:
    base = raw_patch_dir(patch_id)
    return {
        "base": base,
        "ndvi": base / f"{patch_id}_NDVI_{YEAR}.nc",
        "indices20": base / f"{patch_id}_NDII_SWIRMean_explicit20m_{YEAR}.nc",
    }


def patch_paths(patch_id: str) -> dict[str, Path]:
    raw_base = raw_patch_dir(patch_id)
    out_base = patch_dir(patch_id)
    return {
        "base": out_base,
        "grass": grass_path(patch_id),
        "stats_grass": stats_grass_path(patch_id),
        "ndvi": raw_base / f"{patch_id}_NDVI_{YEAR}.nc",
        "indices20": raw_base / f"{patch_id}_NDII_SWIRMean_explicit20m_{YEAR}.nc",
        "baseline": out_base / f"{patch_id}_NDVI_baseline_{BASELINE_VERSION}_{YEAR}.tif",
        "dedup20": out_base / f"{patch_id}_NDVI_NDII_SWIR_event_responses_20m_v2_{YEAR}.csv",
        "overflow": out_base / f"{patch_id}_overflow_gt6_{YEAR}.csv",
        "final": out_base / f"{patch_id}_MOWING_FINAL_{FINAL_VERSION}_{YEAR}.tif",
        "diagnostics": out_base / f"{patch_id}_MOWING_DIAGNOSTICS_{FINAL_VERSION}_{YEAR}.tif",
        "events": out_base / f"{patch_id}_MOWING_FINAL_events_{FINAL_VERSION}_{YEAR}.csv",
    }
