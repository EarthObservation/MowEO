from __future__ import annotations

from pathlib import Path
import hashlib
import warnings

import numpy as np
import pandas as pd
import rasterio
import xarray as xr
from scipy.signal import find_peaks
from tqdm.auto import tqdm

from .config import (
    YEAR, TARGET_EPSG, MAX_EVENTS,
    DROP_THR, GROW_R, MIN_NDVI_VALID,
    MIN_EVENT_SPACING_DAYS, RAPID_REBOUND_DAYS, RAPID_REBOUND_THRESHOLD,
    DATE_TOLERANCE, MOWING_START, MOWING_END,
    PRE_START_DAYS, PRE_END_DAYS, EVENT_BEFORE_DAYS, EVENT_AFTER_DAYS,
    CHILD_DATE_TOLERANCE,
    EVENT_CONFIRM_PRE_WINDOW_DAYS, EVENT_CONFIRM_PRE_MIN_OBS,
    EVENT_CONFIRM_PRE_RANGE_MAX, EVENT_CONFIRM_START_MEDIAN_TOL,
    LOCAL_OBS_WINDOW_DAYS,
)

NODATA_OUT = -9999.0
NODATA_U16 = np.uint16(65535)


def _bool_series(series: pd.Series, field: str) -> pd.Series:
    """Parse booleans robustly after CSV round-trips."""
    if pd.api.types.is_bool_dtype(series):
        return series.astype(bool)
    mapping = {
        "true": True, "1": True, "yes": True, "y": True,
        "false": False, "0": False, "no": False, "n": False,
    }
    text = series.astype(str).str.strip().str.lower()
    bad = ~text.isin(mapping)
    if bad.any():
        raise ValueError(
            f"{field} contains invalid boolean values: "
            f"{sorted(text[bad].unique().tolist())}"
        )
    return text.map(mapping).astype(bool)


def read_patch_table(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path, dtype={"PATCH_ID": str})
    if "PATCH_ID" not in df.columns:
        raise KeyError(f"PATCH_ID missing from {csv_path}")
    df["PATCH_ID"] = df["PATCH_ID"].astype(str).str.strip()
    if "HAS_GRASS" not in df.columns:
        raise KeyError(
            f"HAS_GRASS missing from {csv_path}. "
            "Re-run moweo validate."
        )
    for field in [
        "HAS_GRASS", "HAS_DETECTION_GRASS", "HAS_STATS_GRASS",
        "DETECTION_MASK_READY", "STATS_MASK_READY", "STATS_CONFIGURED",
    ]:
        if field in df.columns:
            df[field] = _bool_series(df[field], field)
    for field in [
        "GRASS_PIXELS_MANIFEST", "DETECTION_GRASS_PIXELS_MANIFEST",
        "STATS_GRASS_PIXELS_MANIFEST",
    ]:
        if field in df.columns:
            df[field] = pd.to_numeric(df[field], errors="coerce")
    return df


def load_grass(grass_file: Path):
    if not grass_file.exists():
        raise FileNotFoundError(grass_file)

    with rasterio.open(grass_file) as src:
        grass = src.read(1)
        profile = src.profile.copy()
        transform = src.transform
        crs = src.crs
        bounds = src.bounds
        width = src.width
        height = src.height

    if crs is None or crs.to_epsg() != TARGET_EPSG:
        raise RuntimeError(f"Expected EPSG:{TARGET_EPSG}, got {crs}")
    if not np.isin(grass, [0, 1]).all():
        raise RuntimeError('Grass mask must contain only 0 and 1.')
    if transform.a != 10 or transform.e != -10 or transform.b != 0 or transform.d != 0:
        raise RuntimeError('Grass mask must be north-up with 10 m pixels and no rotation.')
    if (height, width) != (500, 500):
        raise RuntimeError(f"Expected 500x500 grass grid, got {height}x{width}")
    if not np.isclose(abs(transform.a), 10.0) or not np.isclose(abs(transform.e), 10.0):
        raise RuntimeError(f"Expected 10 m grid, got transform {transform}")
    if not np.isclose(bounds.right - bounds.left, 5000.0) or not np.isclose(bounds.top - bounds.bottom, 5000.0):
        raise RuntimeError(f"Expected 5 km x 5 km extent, got {bounds}")

    return grass, profile, transform, crs, bounds, width, height


def grass_mask_signature(grass_or_path) -> str:
    """Stable SHA256 of the binary 10 m detection mask."""
    if isinstance(grass_or_path, (str, Path)):
        grass, _, _, _, _, _, _ = load_grass(Path(grass_or_path))
    else:
        grass = np.asarray(grass_or_path)
    binary = (grass == 1).astype(np.uint8, copy=False)
    return hashlib.sha256(binary.tobytes(order="C")).hexdigest()


def load_ndvi_daily(ndvi_nc: Path, grass_transform, width: int, height: int):
    if not ndvi_nc.exists():
        raise FileNotFoundError(ndvi_nc)

    with xr.open_dataset(ndvi_nc, engine="h5netcdf") as ds:
        da = ds["NDVI"].transpose("t", "y", "x")
        da = da.sortby("x", ascending=True)
        da = da.sortby("y", ascending=False)

        expected_x = grass_transform.c + grass_transform.a / 2 + np.arange(width) * grass_transform.a
        expected_y = grass_transform.f + grass_transform.e / 2 + np.arange(height) * grass_transform.e

        if da.sizes["x"] != width or da.sizes["y"] != height:
            raise RuntimeError(f"NDVI grid size mismatch: {da.sizes['y']}x{da.sizes['x']}")
        if not np.allclose(da.x.values, expected_x, atol=0.01, rtol=0):
            raise RuntimeError("NDVI x-grid not aligned with grass grid.")
        if not np.allclose(da.y.values, expected_y, atol=0.01, rtol=0):
            raise RuntimeError("NDVI y-grid not aligned with grass grid.")

        dates_raw = pd.DatetimeIndex(pd.to_datetime(da.t.values)).normalize()
        ndvi_raw = da.values.astype(np.float32)

    dates = pd.DatetimeIndex(sorted(pd.unique(dates_raw)))
    ndvi_daily = np.full((len(dates), height, width), np.nan, dtype=np.float32)

    for i, day in enumerate(dates):
        idx = np.where(dates_raw == day)[0]
        if len(idx) == 1:
            ndvi_daily[i] = ndvi_raw[idx[0]]
        else:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                ndvi_daily[i] = np.nanmedian(ndvi_raw[idx], axis=0)

    return dates, ndvi_daily


def postprocess_mowing_events(events):
    if not events:
        return events

    cleaned = []
    for event in events:
        rebound_days = (event["event_end"] - event["bottom_timestamp"]).days
        rebound = event["value_end"] - event["value_bottom"]
        if rebound_days <= RAPID_REBOUND_DAYS and rebound > RAPID_REBOUND_THRESHOLD:
            continue
        cleaned.append(event)

    if len(cleaned) <= 1:
        return cleaned

    cleaned = sorted(cleaned, key=lambda e: e["bottom_timestamp"])
    filtered = [cleaned[0]]

    for event in cleaned[1:]:
        previous = filtered[-1]
        delta_days = (event["bottom_timestamp"] - previous["bottom_timestamp"]).days
        if delta_days > MIN_EVENT_SPACING_DAYS:
            filtered.append(event)
        elif event["value_drop"] > previous["value_drop"]:
            filtered[-1] = event

    return filtered


def _event_context_metrics(t: pd.DatetimeIndex, f: np.ndarray, pk: int, vl: int, rec: int) -> dict:
    """Event-local sampling and shape metrics used by v0.2 confirmation/confidence."""
    bottom = t[vl]
    gap_before = int((t[vl] - t[vl - 1]).days) if vl > 0 else np.nan
    gap_after = int((t[vl + 1] - t[vl]).days) if vl + 1 < len(t) else np.nan

    local_mask = (t >= bottom - pd.Timedelta(days=LOCAL_OBS_WINDOW_DAYS)) & (t <= bottom + pd.Timedelta(days=LOCAL_OBS_WINDOW_DAYS))
    local_obs_count = int(local_mask.sum())

    pre_mask = (t >= bottom - pd.Timedelta(days=EVENT_CONFIRM_PRE_WINDOW_DAYS)) & (t < bottom)
    pre_vals = f[pre_mask]
    pre_obs_count = int(len(pre_vals))
    if pre_obs_count:
        pre_range = float(np.nanmax(pre_vals) - np.nanmin(pre_vals))
        pre_median = float(np.nanmedian(pre_vals))
    else:
        pre_range = np.nan
        pre_median = np.nan

    drop = float(f[pk] - f[vl])
    recovery = float(f[rec] - f[vl])

    # Event confirmation evidence:
    # 1) the original detector already requires a recovery > DROP_THR*GROW_R;
    # 2) v0.3 records pre-event stability as evidence but does not reject
    #    otherwise valid drop + recovery candidates solely because the
    #    pre-event segment is sparse or moderately variable.
    pre_supported = pre_obs_count >= EVENT_CONFIRM_PRE_MIN_OBS
    pre_stable = (
        pre_supported
        and np.isfinite(pre_range)
        and pre_range <= EVENT_CONFIRM_PRE_RANGE_MAX
        and np.isfinite(pre_median)
        and float(f[pk]) >= pre_median - EVENT_CONFIRM_START_MEDIAN_TOL
    )

    return {
        "gap_before": gap_before,
        "gap_after": gap_after,
        "local_observation_count": local_obs_count,
        "pre_observation_count": pre_obs_count,
        "pre_range": pre_range,
        "pre_median": pre_median,
        "pre_stable": bool(pre_stable),
        "recovery_strength": recovery,
        "drop_strength": drop,
    }


def observation_density_score(gap_before, gap_after, local_observation_count) -> float:
    """0-100 event-local observation-density score using the agreed high/medium/low rules."""
    if not np.isfinite(gap_before) or not np.isfinite(gap_after):
        return 0.0
    if gap_before <= 10 and gap_after <= 10 and local_observation_count >= 4:
        return 100.0
    if gap_before <= 20 and gap_after <= 20 and local_observation_count >= 3:
        return 65.0
    return 30.0


def signal_quality_score(drop_strength, recovery_strength, pre_range, pre_stable) -> float:
    """0-100 shape score; diagnostic evidence, not a calibrated probability."""
    drop_scale = max(0.20 - DROP_THR, 1e-6)
    drop_score = 100.0 * np.clip((float(drop_strength) - DROP_THR) / drop_scale, 0.0, 1.0)
    rec_thr = DROP_THR * GROW_R
    rec_score = 100.0 * np.clip(float(recovery_strength) / max(rec_thr * 1.5, 1e-6), 0.0, 1.0)
    if np.isfinite(pre_range):
        stability_score = 100.0 * np.clip(1.0 - float(pre_range) / max(EVENT_CONFIRM_PRE_RANGE_MAX, 1e-6), 0.0, 1.0)
    else:
        stability_score = 0.0
    if not pre_stable:
        stability_score = min(stability_score, 25.0)
    return float(0.40 * drop_score + 0.35 * rec_score + 0.25 * stability_score)


def detect_mowing_sgs(
    dates,
    values,
    drop_thr=DROP_THR,
    grow_r=GROW_R,
    min_ndvi_valid=MIN_NDVI_VALID,
    mowing_start=None,
    mowing_end=None,
):
    values = np.asarray(values, dtype=np.float32)
    dates = pd.DatetimeIndex(dates)
    mowing_start = pd.Timestamp(mowing_start) if mowing_start is not None else None
    mowing_end = pd.Timestamp(mowing_end) if mowing_end is not None else None

    valid = np.isfinite(values)
    if min_ndvi_valid is not None:
        valid &= values > min_ndvi_valid

    f = values[valid]
    t = dates[valid]
    valid_obs = len(f)

    if valid_obs < 3:
        return [], valid_obs, 0

    gaps = np.diff(t.values).astype("timedelta64[D]").astype(int)
    max_gap_days = int(gaps.max()) if len(gaps) else 0

    thr1 = drop_thr
    thr2 = drop_thr * grow_r

    pk_ids = find_peaks(f)[0]
    vl_ids = find_peaks(-f)[0]

    if len(pk_ids) == 0 or len(vl_ids) == 0:
        return [], valid_obs, max_gap_days

    if vl_ids[0] < pk_ids[0]:
        pk_ids = np.insert(pk_ids, 0, 0)
    if vl_ids[-1] < pk_ids[-1]:
        vl_ids = np.append(vl_ids, len(f) - 1)

    pairs = []
    j = 0
    for pk in pk_ids:
        while j < len(vl_ids) and vl_ids[j] <= pk:
            j += 1
        if j >= len(vl_ids):
            break
        pairs.append([pk, vl_ids[j]])
        j += 1

    if not pairs:
        return [], valid_obs, max_gap_days

    pairs = np.asarray(pairs, dtype=int)
    new_pairs = [pairs[0].copy()]

    for pk2, vl2 in pairs[1:]:
        pk1, vl1 = new_pairs[-1]
        y11, y12 = f[pk1], f[vl1]
        y21, y22 = f[pk2], f[vl2]
        if (y21 - y12 < thr2) and (y22 < y12) and (y21 < y11):
            new_pairs[-1][1] = vl2
        else:
            new_pairs.append(np.array([pk2, vl2]))

    pairs = np.asarray(new_pairs, dtype=int)
    keep = (f[pairs[:, 0]] - f[pairs[:, 1]]) > thr1
    pairs = pairs[keep]

    if len(pairs) == 0:
        return [], valid_obs, max_gap_days

    markers = []
    for p_id, (pk0, vl0) in enumerate(pairs):
        pk = int(pk0)
        vl = int(vl0)

        for idx in range(vl, pk - 1, -1):
            if f[idx] - f[vl] > thr1:
                pk = idx
                break

        next_pk = pairs[p_id + 1][0] + 1 if p_id + 1 < len(pairs) else len(f)
        eligible = False
        rec = None

        for idx in range(vl, next_pk):
            if f[idx] - f[vl] > thr2:
                rec = idx
                eligible = True
                break
            if f[idx] < f[vl]:
                vl = idx

        if eligible:
            markers.append((pk, vl, rec))

    events = []
    for pk, vl, rec in markers:
        event_start = t[pk]
        bottom = t[vl]
        event_end = t[rec]
        duration = int((event_end - event_start).days)
        n_obs = int(rec - pk + 1)

        if n_obs == 3 and duration < 9:
            continue
        if mowing_start is not None and bottom < mowing_start:
            continue
        if mowing_end is not None and bottom > mowing_end:
            continue

        metrics = _event_context_metrics(t, f, pk, vl, rec)
        density_score = observation_density_score(
            metrics["gap_before"], metrics["gap_after"], metrics["local_observation_count"]
        )
        sig_score = signal_quality_score(
            metrics["drop_strength"], metrics["recovery_strength"], metrics["pre_range"], metrics["pre_stable"]
        )

        events.append({
            "event_start": event_start,
            "bottom_timestamp": bottom,
            "event_end": event_end,
            "value_start": float(f[pk]),
            "value_bottom": float(f[vl]),
            "value_end": float(f[rec]),
            "value_drop": float(f[pk] - f[vl]),
            "duration": duration,
            "num_observations": n_obs,
            "doy": int(bottom.dayofyear),
            **metrics,
            "observation_density_score": density_score,
            "signal_quality_score": sig_score,
        })

    return postprocess_mowing_events(events), valid_obs, max_gap_days

def write_baseline(
    baseline_tif: Path,
    grass,
    grass_profile,
    patch_id: str,
    mowing_count,
    mowing_doys,
    ndvi_drops,
    mowing_support,
    mean_event_support,
    low_support_events,
    valid_obs_arr,
    max_gap_arr,
    observation_support,
    provisional_confidence,
    event_gap_before,
    event_gap_after,
    event_local_obs_count,
    event_pre_range,
    event_pre_stable,
    event_recovery_strength,
    event_density_score,
    event_signal_quality,
):
    bands, names = [], []

    def add_band(arr, name, replace_65535=False):
        a = np.asarray(arr, dtype=np.float32).copy()
        if replace_65535:
            a[a == 65535] = NODATA_OUT
        a[~np.isfinite(a)] = NODATA_OUT
        a[grass != 1] = NODATA_OUT
        bands.append(a)
        names.append(name)

    add_band(mowing_count, "MOWING_COUNT", True)
    for i in range(MAX_EVENTS):
        add_band(mowing_doys[i], f"MOWING_DOY_{i+1}", True)
    for i in range(MAX_EVENTS):
        add_band(ndvi_drops[i], f"NDVI_DROP_{i+1}")
    for i in range(MAX_EVENTS):
        add_band(mowing_support[i], f"MOWING_SUPPORT_{i+1}", True)
    add_band(mean_event_support, "MEAN_EVENT_SUPPORT")
    add_band(low_support_events, "LOW_SUPPORT_EVENTS")
    add_band(valid_obs_arr, "VALID_OBS", True)
    add_band(max_gap_arr, "MAX_GAP_DAYS", True)
    add_band(observation_support, "OBSERVATION_SUPPORT")
    add_band(provisional_confidence, "PROVISIONAL_CONFIDENCE")
    for i in range(MAX_EVENTS): add_band(event_gap_before[i], f"EVENT_GAP_BEFORE_{i+1}")
    for i in range(MAX_EVENTS): add_band(event_gap_after[i], f"EVENT_GAP_AFTER_{i+1}")
    for i in range(MAX_EVENTS): add_band(event_local_obs_count[i], f"EVENT_LOCAL_OBS_{i+1}")
    for i in range(MAX_EVENTS): add_band(event_pre_range[i], f"EVENT_PRE_RANGE_{i+1}")
    for i in range(MAX_EVENTS): add_band(event_pre_stable[i], f"EVENT_PRE_STABLE_{i+1}")
    for i in range(MAX_EVENTS): add_band(event_recovery_strength[i], f"EVENT_RECOVERY_{i+1}")
    for i in range(MAX_EVENTS): add_band(event_density_score[i], f"EVENT_DENSITY_SCORE_{i+1}")
    for i in range(MAX_EVENTS): add_band(event_signal_quality[i], f"EVENT_SIGNAL_SCORE_{i+1}")

    profile = grass_profile.copy()
    profile.update(
        driver="GTiff", dtype="float32", count=len(bands), nodata=NODATA_OUT,
        compress="DEFLATE", predictor=3, tiled=True,
        blockxsize=256, blockysize=256, BIGTIFF="IF_SAFER",
    )

    baseline_tif.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(baseline_tif, "w", **profile) as dst:
        for i, (arr, name) in enumerate(zip(bands, names), start=1):
            dst.write(arr, i)
            dst.set_band_description(i, name)
        dst.update_tags(
            PRODUCT="NDVI mowing detection baseline",
            VERSION="v1.3",
            PATCH_ID=str(patch_id),
            YEAR=str(YEAR),
            METHOD="NDVI_SGS_baseline",
            NDVI_BANDS="B04,B08",
            DROP_THR=str(DROP_THR),
            GROW_R=str(GROW_R),
            RECOVERY_THR=str(round(DROP_THR * GROW_R, 4)),
            MOWING_START=str(pd.Timestamp(MOWING_START).date()),
            MOWING_END=str(pd.Timestamp(MOWING_END).date()),
            MIN_EVENT_SPACING_DAYS=str(MIN_EVENT_SPACING_DAYS),
            RAPID_REBOUND_DAYS=str(RAPID_REBOUND_DAYS),
            RAPID_REBOUND_THRESHOLD=str(RAPID_REBOUND_THRESHOLD),
            RAPID_REBOUND_NOTE="Experimental NDVI adaptation; not calibrated probability rule",
            NEIGHBOURHOOD="3x3 pixels",
            NEIGHBOUR_DATE_TOLERANCE_DAYS=str(DATE_TOLERANCE),
            LOW_SUPPORT_THRESHOLD_PERCENT="25",
            PROVISIONAL_CONFIDENCE_FORMULA=(
                "0.6*MEAN_EVENT_SUPPORT + 0.4*OBSERVATION_SUPPORT - 10*LOW_SUPPORT_EVENTS"
            ),
            PROVISIONAL_CONFIDENCE_NOTE="Diagnostic indicator only; not calibrated probability",
            EVENT_CONFIRMATION=(
                f"Recovery>{DROP_THR * GROW_R:.4f}; pre-event stability retained as confidence evidence: "
                f">= {EVENT_CONFIRM_PRE_MIN_OBS} observations "
                f"within {EVENT_CONFIRM_PRE_WINDOW_DAYS} d; pre-range <= {EVENT_CONFIRM_PRE_RANGE_MAX}; "
                f"event-start >= pre-median-{EVENT_CONFIRM_START_MEDIAN_TOL}; not a hard rejection filter"
            ),
            EVENT_DENSITY_HIGH="gap_before<=10 & gap_after<=10 & local_obs(+/-20d)>=4",
            EVENT_DENSITY_MEDIUM="gap_before<=20 & gap_after<=20 & local_obs(+/-20d)>=3",
            MAX_EVENTS_STORED=str(MAX_EVENTS),
            COUNT_POLICY="MOWING_COUNT capped to MAX_EVENTS_STORED for output consistency",
            GRASS_MASK_ROLE="detection +5m",
            GRASS_MASK_SHA256=grass_mask_signature(grass),
        )


def load_baseline_arrays(path: Path):
    with rasterio.open(path) as src:
        data = src.read().astype(np.float32)
        names = list(src.descriptions)
        nodata = src.nodata
        meta = {
            "crs": src.crs,
            "transform": src.transform,
            "bounds": src.bounds,
            "profile": src.profile.copy(),
        }

    required = (
        ["MOWING_COUNT"]
        + [f"MOWING_DOY_{i}" for i in range(1, MAX_EVENTS + 1)]
        + [f"NDVI_DROP_{i}" for i in range(1, MAX_EVENTS + 1)]
        + [f"MOWING_SUPPORT_{i}" for i in range(1, MAX_EVENTS + 1)]
        + ["MEAN_EVENT_SUPPORT", "LOW_SUPPORT_EVENTS", "VALID_OBS", "MAX_GAP_DAYS", "OBSERVATION_SUPPORT", "PROVISIONAL_CONFIDENCE"]
        + [f"EVENT_GAP_BEFORE_{i}" for i in range(1, MAX_EVENTS + 1)]
        + [f"EVENT_GAP_AFTER_{i}" for i in range(1, MAX_EVENTS + 1)]
        + [f"EVENT_LOCAL_OBS_{i}" for i in range(1, MAX_EVENTS + 1)]
        + [f"EVENT_PRE_RANGE_{i}" for i in range(1, MAX_EVENTS + 1)]
        + [f"EVENT_PRE_STABLE_{i}" for i in range(1, MAX_EVENTS + 1)]
        + [f"EVENT_RECOVERY_{i}" for i in range(1, MAX_EVENTS + 1)]
        + [f"EVENT_DENSITY_SCORE_{i}" for i in range(1, MAX_EVENTS + 1)]
        + [f"EVENT_SIGNAL_SCORE_{i}" for i in range(1, MAX_EVENTS + 1)]
    )
    missing = [x for x in required if x not in names]
    if missing:
        raise RuntimeError(f"Baseline missing bands: {missing}")

    by_name = {name: data[i] for i, name in enumerate(names)}

    def clean(name):
        a = by_name[name].copy()
        a[a == nodata] = np.nan
        return a

    arrays = {
        "mowing_count": clean("MOWING_COUNT"),
        "mowing_doys": np.stack([clean(f"MOWING_DOY_{i}") for i in range(1, MAX_EVENTS + 1)]),
        "ndvi_drops": np.stack([clean(f"NDVI_DROP_{i}") for i in range(1, MAX_EVENTS + 1)]),
        "mowing_support": np.stack([clean(f"MOWING_SUPPORT_{i}") for i in range(1, MAX_EVENTS + 1)]),
        "mean_event_support": clean("MEAN_EVENT_SUPPORT"),
        "low_support_events": clean("LOW_SUPPORT_EVENTS"),
        "valid_obs": clean("VALID_OBS"),
        "max_gap_days": clean("MAX_GAP_DAYS"),
        "observation_support": clean("OBSERVATION_SUPPORT"),
        "baseline_confidence": clean("PROVISIONAL_CONFIDENCE"),
        "event_gap_before": np.stack([clean(f"EVENT_GAP_BEFORE_{i}") for i in range(1, MAX_EVENTS + 1)]),
        "event_gap_after": np.stack([clean(f"EVENT_GAP_AFTER_{i}") for i in range(1, MAX_EVENTS + 1)]),
        "event_local_obs_count": np.stack([clean(f"EVENT_LOCAL_OBS_{i}") for i in range(1, MAX_EVENTS + 1)]),
        "event_pre_range": np.stack([clean(f"EVENT_PRE_RANGE_{i}") for i in range(1, MAX_EVENTS + 1)]),
        "event_pre_stable": np.stack([clean(f"EVENT_PRE_STABLE_{i}") for i in range(1, MAX_EVENTS + 1)]),
        "event_recovery_strength": np.stack([clean(f"EVENT_RECOVERY_{i}") for i in range(1, MAX_EVENTS + 1)]),
        "event_density_score": np.stack([clean(f"EVENT_DENSITY_SCORE_{i}") for i in range(1, MAX_EVENTS + 1)]),
        "event_signal_quality": np.stack([clean(f"EVENT_SIGNAL_SCORE_{i}") for i in range(1, MAX_EVENTS + 1)]),
    }
    return arrays, meta


def baseline_is_consistent(path: Path, grass_file: Path | None = None) -> bool:
    try:
        arrays, _ = load_baseline_arrays(path)
        vals = arrays["mowing_count"]
        if np.any(np.isfinite(vals) & (vals > MAX_EVENTS)):
            return False
        if grass_file is not None:
            with rasterio.open(path) as src:
                stored = src.tags().get("GRASS_MASK_SHA256", "")
            if not stored or stored != grass_mask_signature(grass_file):
                return False
        return True
    except Exception:
        return False


def build_baseline(patch_id: str, grass_file: Path, ndvi_nc: Path, baseline_tif: Path, force=False):
    grass, grass_profile, grass_transform, grass_crs, grass_bounds, width, height = load_grass(grass_file)

    if baseline_tif.exists() and not force and baseline_is_consistent(baseline_tif, grass_file):
        arrays, _ = load_baseline_arrays(baseline_tif)
        # Raw >MAX_EVENTS information cannot be reconstructed from the capped
        # baseline alone. New national runs store it separately in overflow CSV.
        return {
            "status": "existing",
            "overflow_pixels": np.nan,
            "overflow_records": [],
            "grass_pixels": int((grass == 1).sum()),
            **arrays,
        }

    dates, ndvi_daily = load_ndvi_daily(ndvi_nc, grass_transform, width, height)

    mowing_count = np.full((height, width), NODATA_U16, dtype=np.uint16)
    mowing_doys = np.full((MAX_EVENTS, height, width), NODATA_U16, dtype=np.uint16)
    valid_obs_arr = np.full((height, width), NODATA_U16, dtype=np.uint16)
    max_gap_arr = np.full((height, width), NODATA_U16, dtype=np.uint16)
    ndvi_drops = np.zeros((MAX_EVENTS, height, width), dtype=np.float32)
    event_gap_before = np.full((MAX_EVENTS, height, width), np.nan, dtype=np.float32)
    event_gap_after = np.full((MAX_EVENTS, height, width), np.nan, dtype=np.float32)
    event_local_obs_count = np.full((MAX_EVENTS, height, width), np.nan, dtype=np.float32)
    event_pre_range = np.full((MAX_EVENTS, height, width), np.nan, dtype=np.float32)
    event_pre_stable = np.full((MAX_EVENTS, height, width), np.nan, dtype=np.float32)
    event_recovery_strength = np.full((MAX_EVENTS, height, width), np.nan, dtype=np.float32)
    event_density_score = np.full((MAX_EVENTS, height, width), np.nan, dtype=np.float32)
    event_signal_quality = np.full((MAX_EVENTS, height, width), np.nan, dtype=np.float32)

    grass_pixels = grass == 1
    mowing_count[grass_pixels] = 0
    mowing_doys[:, grass_pixels] = 0
    valid_obs_arr[grass_pixels] = 0
    max_gap_arr[grass_pixels] = 0

    rows, cols = np.where(grass_pixels)
    overflow_pixels = 0
    overflow_records = []

    for r, c in tqdm(zip(rows, cols), total=len(rows), desc=f"{patch_id} NDVI detection"):
        events, valid_obs, max_gap = detect_mowing_sgs(
            dates=dates,
            values=ndvi_daily[:, r, c],
            mowing_start=MOWING_START,
            mowing_end=MOWING_END,
        )
        if len(events) > MAX_EVENTS:
            overflow_pixels += 1
            overflow_records.append({
                "PATCH_ID": patch_id,
                "row": int(r),
                "col": int(c),
                "RAW_EVENT_COUNT": int(len(events)),
                "RAW_DOYS": ";".join(str(int(e["doy"])) for e in events),
                "RAW_NDVI_DROPS": ";".join(
                    f"{float(e['value_drop']):.6f}" for e in events
                ),
            })

        kept = events[:MAX_EVENTS]
        mowing_count[r, c] = len(kept)
        valid_obs_arr[r, c] = valid_obs
        max_gap_arr[r, c] = max_gap

        for i, event in enumerate(kept):
            mowing_doys[i, r, c] = event["doy"]
            ndvi_drops[i, r, c] = event["value_drop"]
            event_gap_before[i, r, c] = event.get("gap_before", np.nan)
            event_gap_after[i, r, c] = event.get("gap_after", np.nan)
            event_local_obs_count[i, r, c] = event.get("local_observation_count", np.nan)
            event_pre_range[i, r, c] = event.get("pre_range", np.nan)
            event_pre_stable[i, r, c] = 1.0 if event.get("pre_stable", False) else 0.0
            event_recovery_strength[i, r, c] = event.get("recovery_strength", np.nan)
            event_density_score[i, r, c] = event.get("observation_density_score", np.nan)
            event_signal_quality[i, r, c] = event.get("signal_quality_score", np.nan)

    def get_pixel_doys(row, col):
        vals = mowing_doys[:, row, col]
        vals = vals[(vals > 0) & (vals < 65535)]
        return vals.astype(int)

    def event_neighbour_support(row, col, doy, tolerance=DATE_TOLERANCE):
        r0, r1 = max(0, row - 1), min(height, row + 2)
        c0, c1 = max(0, col - 1), min(width, col + 2)
        n_neighbours = 0
        n_matches = 0
        for rr in range(r0, r1):
            for cc in range(c0, c1):
                if rr == row and cc == col:
                    continue
                if grass[rr, cc] != 1:
                    continue
                n_neighbours += 1
                ndoys = get_pixel_doys(rr, cc)
                if len(ndoys) and np.any(np.abs(ndoys - doy) <= tolerance):
                    n_matches += 1
        return np.nan if n_neighbours == 0 else 100.0 * n_matches / n_neighbours

    mowing_support = np.full((MAX_EVENTS, height, width), 65535, dtype=np.uint16)
    mowing_support[:, grass == 1] = 0

    for r, c in tqdm(zip(rows, cols), total=len(rows), desc=f"{patch_id} neighbour support"):
        for i, doy in enumerate(get_pixel_doys(r, c)):
            support = event_neighbour_support(r, c, int(doy))
            if np.isfinite(support):
                mowing_support[i, r, c] = int(round(support))

    mean_event_support = np.full((height, width), np.nan, dtype=np.float32)
    low_support_events = np.full((height, width), np.nan, dtype=np.float32)

    for r, c in zip(rows, cols):
        n = int(mowing_count[r, c])
        if n == 0:
            low_support_events[r, c] = 0
            continue
        vals = mowing_support[:n, r, c].astype(float)
        vals = vals[vals != 65535]
        if len(vals):
            mean_event_support[r, c] = np.mean(vals)
            low_support_events[r, c] = np.sum(vals < 25)

    observation_support = np.full((height, width), np.nan, dtype=np.float32)
    for r, c in zip(rows, cols):
        obs_score = np.clip((float(valid_obs_arr[r, c]) - 20) / 20, 0, 1)
        gap_score = np.clip((40 - float(max_gap_arr[r, c])) / 25, 0, 1)
        observation_support[r, c] = 100 * (0.5 * obs_score + 0.5 * gap_score)

    provisional_confidence = np.full((height, width), np.nan, dtype=np.float32)
    for r, c in zip(rows, cols):
        if mowing_count[r, c] == 0:
            continue
        if np.isfinite(mean_event_support[r, c]) and np.isfinite(observation_support[r, c]) and np.isfinite(low_support_events[r, c]):
            score = 0.6 * mean_event_support[r, c] + 0.4 * observation_support[r, c] - 10 * low_support_events[r, c]
            provisional_confidence[r, c] = np.clip(score, 0, 100)

    write_baseline(
        baseline_tif, grass, grass_profile, patch_id,
        mowing_count, mowing_doys, ndvi_drops, mowing_support,
        mean_event_support, low_support_events, valid_obs_arr, max_gap_arr,
        observation_support, provisional_confidence,
        event_gap_before, event_gap_after, event_local_obs_count, event_pre_range, event_pre_stable,
        event_recovery_strength, event_density_score, event_signal_quality,
    )

    return {
        "status": "recomputed",
        "overflow_pixels": int(overflow_pixels),
        "overflow_records": overflow_records,
        "grass_pixels": int(grass_pixels.sum()),
        "mowing_count": mowing_count.astype(np.float32),
        "mowing_doys": mowing_doys.astype(np.float32),
        "ndvi_drops": ndvi_drops,
        "mowing_support": mowing_support.astype(np.float32),
        "mean_event_support": mean_event_support,
        "low_support_events": low_support_events,
        "valid_obs": valid_obs_arr.astype(np.float32),
        "max_gap_days": max_gap_arr.astype(np.float32),
        "observation_support": observation_support,
        "baseline_confidence": provisional_confidence,
        "event_gap_before": event_gap_before,
        "event_gap_after": event_gap_after,
        "event_local_obs_count": event_local_obs_count,
        "event_pre_range": event_pre_range,
        "event_pre_stable": event_pre_stable,
        "event_recovery_strength": event_recovery_strength,
        "event_density_score": event_density_score,
        "event_signal_quality": event_signal_quality,
    }


def load_indices20(indices20_nc: Path, grass_bounds):
    if not indices20_nc.exists():
        raise FileNotFoundError(indices20_nc)

    ds20 = xr.open_dataset(indices20_nc, engine="h5netcdf")
    if "NDII" not in ds20 or "SWIRMean" not in ds20:
        ds20.close()
        raise RuntimeError("20 m cube must contain NDII and SWIRMean.")
    if ds20.sizes["x"] != 250 or ds20.sizes["y"] != 250:
        ds20.close()
        raise RuntimeError(f"Expected 250x250 explicit 20 m cube, got {ds20.sizes['y']}x{ds20.sizes['x']}")

    expected_x20 = grass_bounds.left + 10 + np.arange(250) * 20
    expected_y20 = grass_bounds.top - 10 - np.arange(250) * 20
    if not np.allclose(ds20.x.values, expected_x20, atol=0.01, rtol=0):
        ds20.close()
        raise RuntimeError("20 m x-grid is not aligned.")
    if not np.allclose(ds20.y.values, expected_y20, atol=0.01, rtol=0):
        ds20.close()
        raise RuntimeError("20 m y-grid is not aligned.")

    ndii = ds20["NDII"].sortby("y", ascending=False)
    swir = ds20["SWIRMean"].sortby("y", ascending=False)
    valid_date = np.isfinite(ndii.values).any(axis=(1, 2)) & np.isfinite(swir.values).any(axis=(1, 2))
    ndii_values = ndii.isel(t=valid_date).values.astype(np.float32)
    swir_values = swir.isel(t=valid_date).values.astype(np.float32)
    dates20 = pd.DatetimeIndex(pd.to_datetime(ndii.isel(t=valid_date).t.values)).normalize()
    ds20.close()
    return dates20, ndii_values, swir_values


def spectral_response_from_arrays(row10, col10, event_date, dates20, ndii20, swir20):
    row20, col20 = row10 // 2, col10 // 2
    ndii_ts = ndii20[:, row20, col20]
    swir_ts = swir20[:, row20, col20] / 10000.0
    event_date = pd.Timestamp(event_date)

    pre_mask = (
        (dates20 >= event_date - pd.Timedelta(days=PRE_START_DAYS))
        & (dates20 <= event_date - pd.Timedelta(days=PRE_END_DAYS))
        & np.isfinite(ndii_ts) & np.isfinite(swir_ts)
    )
    pre_idx = np.where(pre_mask)[0]
    if len(pre_idx):
        pre_i = pre_idx[-1]
        pre_date = dates20[pre_i]
        ndii_pre = float(ndii_ts[pre_i])
        swir_pre = float(swir_ts[pre_i])
    else:
        pre_date, ndii_pre, swir_pre = pd.NaT, np.nan, np.nan

    event_mask = (
        (dates20 >= event_date - pd.Timedelta(days=EVENT_BEFORE_DAYS))
        & (dates20 <= event_date + pd.Timedelta(days=EVENT_AFTER_DAYS))
        & np.isfinite(ndii_ts) & np.isfinite(swir_ts)
    )
    event_idx = np.where(event_mask)[0]
    if len(event_idx):
        offsets = np.abs((dates20[event_idx] - event_date).days)
        event_i = event_idx[np.argmin(offsets)]
        spectral_date = dates20[event_i]
        ndii_event = float(ndii_ts[event_i])
        swir_event = float(swir_ts[event_i])
        spectral_offset = int((spectral_date - event_date).days)
    else:
        spectral_date, ndii_event, swir_event, spectral_offset = pd.NaT, np.nan, np.nan, np.nan

    ndii_drop = ndii_pre - ndii_event if np.isfinite(ndii_pre) and np.isfinite(ndii_event) else np.nan
    swir_rise = swir_event - swir_pre if np.isfinite(swir_pre) and np.isfinite(swir_event) else np.nan

    return {
        "spectral_pre_date": pre_date,
        "spectral_event_date": spectral_date,
        "spectral_offset_days": spectral_offset,
        "ndii_pre": ndii_pre,
        "ndii_event": ndii_event,
        "ndii_drop": ndii_drop,
        "swir_pre": swir_pre,
        "swir_event": swir_event,
        "swir_rise": swir_rise,
    }


def doy_to_date(doy):
    if not np.isfinite(doy) or doy <= 0:
        return pd.NaT
    return pd.Timestamp(f"{YEAR}-01-01") + pd.Timedelta(days=int(doy) - 1)


def build_dedup20_from_baseline(patch_id: str, grass_file: Path, baseline_tif: Path, indices20_nc: Path, out_csv: Path):
    grass, _, _, _, grass_bounds, _, _ = load_grass(grass_file)
    arrays, _ = load_baseline_arrays(baseline_tif)
    dates20, ndii20, swir20 = load_indices20(indices20_nc, grass_bounds)

    mowing_count = arrays["mowing_count"]
    mowing_doys = arrays["mowing_doys"]
    ndvi_drops = arrays["ndvi_drops"]
    mowing_support = arrays["mowing_support"]
    baseline_conf = arrays["baseline_confidence"]

    grass_fraction_20 = grass.reshape(250, 2, 250, 2).mean(axis=(1, 3)).astype(np.float32)
    pure20_on_10m = np.repeat(np.repeat(grass_fraction_20 == 1.0, 2, axis=0), 2, axis=1)

    rows, cols = np.where((grass == 1) & pure20_on_10m & np.isfinite(mowing_count) & (mowing_count > 0))
    records = []

    for r, c in tqdm(zip(rows, cols), total=len(rows), desc=f"{patch_id} NDII/SWIR responses"):
        for i in range(MAX_EVENTS):
            doy = mowing_doys[i, r, c]
            if not np.isfinite(doy) or doy <= 0:
                continue
            event_date = doy_to_date(doy)
            resp = spectral_response_from_arrays(r, c, event_date, dates20, ndii20, swir20)
            records.append({
                "row10": int(r), "col10": int(c),
                "row20": int(r // 2), "col20": int(c // 2),
                "event": i + 1,
                "mowing_date": event_date,
                "doy": int(doy),
                "ndvi_drop": float(ndvi_drops[i, r, c]) if np.isfinite(ndvi_drops[i, r, c]) else np.nan,
                "neighbour_support": float(mowing_support[i, r, c]) if np.isfinite(mowing_support[i, r, c]) else np.nan,
                "provisional_confidence": float(baseline_conf[r, c]) if np.isfinite(baseline_conf[r, c]) else np.nan,
                **resp,
            })

    responses = pd.DataFrame(records)
    if responses.empty:
        dedup20 = pd.DataFrame(columns=[
            "row20", "col20", "spectral_event_date", "n_10m_events", "n_10m_pixels",
            "mowing_date_min", "mowing_date_max", "ndvi_drop_median", "ndvi_drop_max",
            "neighbour_support_median", "provisional_confidence_median", "ndii_drop", "swir_rise",
            "spectral_offset_days", "child_pixel_support_percent", "spectral_response",
        ])
    else:
        responses["mowing_date"] = pd.to_datetime(responses["mowing_date"]).dt.normalize()
        responses["spectral_event_date"] = pd.to_datetime(responses["spectral_event_date"]).dt.normalize()
        responses["pixel10_id"] = responses["row10"].astype(str) + "_" + responses["col10"].astype(str)
        valid20 = responses.dropna(subset=["spectral_event_date", "ndii_drop", "swir_rise"]).copy()

        dedup20 = (
            valid20.groupby(["row20", "col20", "spectral_event_date"], as_index=False)
            .agg(
                n_10m_events=("event", "size"),
                n_10m_pixels=("pixel10_id", "nunique"),
                mowing_date_min=("mowing_date", "min"),
                mowing_date_max=("mowing_date", "max"),
                ndvi_drop_median=("ndvi_drop", "median"),
                ndvi_drop_max=("ndvi_drop", "max"),
                neighbour_support_median=("neighbour_support", "median"),
                provisional_confidence_median=("provisional_confidence", "median"),
                ndii_drop=("ndii_drop", "median"),
                swir_rise=("swir_rise", "median"),
                spectral_offset_days=("spectral_offset_days", "median"),
            )
        )
        dedup20["child_pixel_support_percent"] = 100.0 * dedup20["n_10m_pixels"] / 4.0
        conditions = [
            (dedup20["ndii_drop"] > 0) & (dedup20["swir_rise"] > 0),
            (dedup20["ndii_drop"] > 0) & (dedup20["swir_rise"] <= 0),
            (dedup20["ndii_drop"] <= 0) & (dedup20["swir_rise"] > 0),
        ]
        dedup20["spectral_response"] = np.select(
            conditions,
            ["NDII_down + SWIR_up", "NDII_down_only", "SWIR_up_only"],
            default="no_confirmation",
        )

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    dedup20.to_csv(out_csv, index=False)
    return dedup20


def get_child_support(grass, mowing_doys, row10, col10, event_doy, tolerance=CHILD_DATE_TOLERANCE):
    row20, col20 = row10 // 2, col10 // 2
    r0, c0 = row20 * 2, col20 * 2
    supporters = 0
    for rr in range(r0, r0 + 2):
        for cc in range(c0, c0 + 2):
            if grass[rr, cc] != 1:
                continue
            doys = mowing_doys[:, rr, cc]
            doys = doys[np.isfinite(doys) & (doys > 0)]
            if len(doys) and np.any(np.abs(doys - event_doy) <= tolerance):
                supporters += 1
    return 100.0 * supporters / 4.0


def spectral_score_from_code(code):
    if code == 3:
        return 100.0
    if code in (1, 2):
        return 50.0
    return 0.0


def confidence_class(score, spectral_code):
    if not np.isfinite(score):
        return 0
    if score >= 80 and spectral_code == 3:
        return 4
    if score >= 60:
        return 3
    if score >= 40:
        return 2
    return 1


def finalize_patch(patch_id: str, grass_file: Path, baseline_tif: Path, indices20_nc: Path, final_tif: Path, diagnostics_tif: Path, events_csv: Path):
    grass, grass_profile, grass_transform, grass_crs, grass_bounds, width, height = load_grass(grass_file)
    arrays, baseline_meta = load_baseline_arrays(baseline_tif)

    if np.any(np.isfinite(arrays["mowing_count"]) & (arrays["mowing_count"] > MAX_EVENTS)):
        raise RuntimeError("Baseline has MOWING_COUNT > MAX_EVENTS. Re-run moweo detect --force.")

    mowing_count = arrays["mowing_count"]
    mowing_doys = arrays["mowing_doys"]
    ndvi_drops = arrays["ndvi_drops"]
    neighbour_support = arrays["mowing_support"]
    valid_obs = arrays["valid_obs"]
    max_gap_days = arrays["max_gap_days"]
    observation_support = arrays["observation_support"]
    baseline_confidence = arrays["baseline_confidence"]
    event_gap_before = arrays["event_gap_before"]
    event_gap_after = arrays["event_gap_after"]
    event_local_obs_count = arrays["event_local_obs_count"]
    event_pre_range = arrays["event_pre_range"]
    event_pre_stable = arrays["event_pre_stable"]
    event_recovery_strength = arrays["event_recovery_strength"]
    event_density_score = arrays["event_density_score"]
    event_signal_quality = arrays["event_signal_quality"]

    dates20, ndii20, swir20 = load_indices20(indices20_nc, grass_bounds)

    shape3 = (MAX_EVENTS, height, width)
    ndii_drop_out = np.full(shape3, np.nan, dtype=np.float32)
    swir_rise_out = np.full(shape3, np.nan, dtype=np.float32)
    child_support_out = np.full(shape3, np.nan, dtype=np.float32)
    spectral_code_out = np.zeros(shape3, dtype=np.float32)
    miletich_flag_out = np.zeros(shape3, dtype=np.float32)
    event_conf_out = np.full(shape3, np.nan, dtype=np.float32)
    event_class_out = np.zeros(shape3, dtype=np.float32)
    spectral_offset_out = np.full(shape3, np.nan, dtype=np.float32)

    event_rows = []
    pixels_with_events = np.argwhere((grass == 1) & np.isfinite(mowing_count) & (mowing_count > 0))

    for r, c in tqdm(pixels_with_events, total=len(pixels_with_events), desc=f"{patch_id} final evidence"):
        r, c = int(r), int(c)
        base_conf = baseline_confidence[r, c]
        for i in range(MAX_EVENTS):
            doy = mowing_doys[i, r, c]
            if not np.isfinite(doy) or doy <= 0:
                continue

            event_date = doy_to_date(doy)
            child_support = get_child_support(grass, mowing_doys, r, c, doy)
            response = spectral_response_from_arrays(r, c, event_date, dates20, ndii20, swir20)

            ndii_drop = response["ndii_drop"]
            swir_rise = response["swir_rise"]
            spectral_offset = response["spectral_offset_days"]
            ndii_confirm = np.isfinite(ndii_drop) and ndii_drop > 0
            swir_confirm = np.isfinite(swir_rise) and swir_rise > 0

            if ndii_confirm and swir_confirm:
                spectral_code = 3
            elif ndii_confirm:
                spectral_code = 1
            elif swir_confirm:
                spectral_code = 2
            else:
                spectral_code = 0

            spectral_score = spectral_score_from_code(spectral_code)
            density_score = event_density_score[i, r, c]
            signal_score = event_signal_quality[i, r, c]
            # v0.3: event-specific confidence. This is an evidence score (0-100), not a calibrated probability.
            components = [
                (signal_score, 0.30),
                (density_score, 0.25),
                (child_support, 0.20),
                (spectral_score, 0.25),
            ]
            valid_components = [(float(v), w) for v, w in components if np.isfinite(v)]
            if valid_components:
                wsum = sum(w for _, w in valid_components)
                event_conf = sum(v * w for v, w in valid_components) / wsum
            else:
                event_conf = np.nan
            event_class = confidence_class(event_conf, spectral_code)
            ndvi_drop = ndvi_drops[i, r, c]
            miletich_flag = int(np.isfinite(ndvi_drop) and ndvi_drop > 0.20 and swir_confirm)

            ndii_drop_out[i, r, c] = ndii_drop
            swir_rise_out[i, r, c] = swir_rise
            child_support_out[i, r, c] = child_support
            spectral_code_out[i, r, c] = spectral_code
            miletich_flag_out[i, r, c] = miletich_flag
            event_conf_out[i, r, c] = event_conf
            event_class_out[i, r, c] = event_class
            spectral_offset_out[i, r, c] = spectral_offset

            event_rows.append({
                "PATCH_ID": patch_id,
                "row": r, "col": c, "event": i + 1,
                "doy": int(doy), "event_date": event_date,
                "ndvi_drop": float(ndvi_drop) if np.isfinite(ndvi_drop) else np.nan,
                "neighbour_support": float(neighbour_support[i, r, c]) if np.isfinite(neighbour_support[i, r, c]) else np.nan,
                "baseline_confidence": float(base_conf) if np.isfinite(base_conf) else np.nan,
                "child_support": float(child_support),
                "ndii_drop": ndii_drop, "swir_rise": swir_rise,
                "spectral_code": spectral_code,
                "spectral_offset_days": spectral_offset,
                "miletich_flag": miletich_flag,
                "gap_before": float(event_gap_before[i, r, c]) if np.isfinite(event_gap_before[i, r, c]) else np.nan,
                "gap_after": float(event_gap_after[i, r, c]) if np.isfinite(event_gap_after[i, r, c]) else np.nan,
                "local_observation_count": float(event_local_obs_count[i, r, c]) if np.isfinite(event_local_obs_count[i, r, c]) else np.nan,
                "pre_range": float(event_pre_range[i, r, c]) if np.isfinite(event_pre_range[i, r, c]) else np.nan,
                "pre_stable": float(event_pre_stable[i, r, c]) if np.isfinite(event_pre_stable[i, r, c]) else np.nan,
                "recovery_strength": float(event_recovery_strength[i, r, c]) if np.isfinite(event_recovery_strength[i, r, c]) else np.nan,
                "observation_density_score": float(density_score) if np.isfinite(density_score) else np.nan,
                "signal_quality_score": float(signal_score) if np.isfinite(signal_score) else np.nan,
                "event_confidence": event_conf,
                "event_class": event_class,
            })

    event_columns = [
        "PATCH_ID", "row", "col", "event", "doy", "event_date",
        "ndvi_drop", "neighbour_support", "baseline_confidence", "child_support",
        "ndii_drop", "swir_rise", "spectral_code", "spectral_offset_days",
        "miletich_flag", "gap_before", "gap_after", "local_observation_count",
        "pre_range", "pre_stable", "recovery_strength", "observation_density_score",
        "signal_quality_score", "event_confidence", "event_class",
    ]
    events_final = pd.DataFrame(event_rows, columns=event_columns)

    pixel_event_conf_mean = np.full((height, width), np.nan, dtype=np.float32)
    pixel_conf_class = np.zeros((height, width), dtype=np.float32)
    for r, c in pixels_with_events:
        r, c = int(r), int(c)
        vals = event_conf_out[:, r, c]
        vals = vals[np.isfinite(vals)]
        classes = event_class_out[:, r, c]
        classes = classes[classes > 0]
        if len(vals):
            pixel_event_conf_mean[r, c] = np.mean(vals)
        if len(classes):
            pixel_conf_class[r, c] = np.min(classes)

    core_bands, core_names = [], []
    def add_core(arr, name):
        a = np.asarray(arr, dtype=np.float32).copy()
        a[~np.isfinite(a)] = NODATA_OUT
        a[grass != 1] = NODATA_OUT
        core_bands.append(a)
        core_names.append(name)

    add_core(mowing_count, "MOWING_COUNT")
    for i in range(MAX_EVENTS):
        add_core(mowing_doys[i], f"MOWING_DOY_{i+1}")
    for i in range(MAX_EVENTS):
        add_core(event_conf_out[i], f"EVENT_CONF_{i+1}")
    add_core(valid_obs, "VALID_OBS")
    add_core(max_gap_days, "MAX_GAP_DAYS")
    add_core(pixel_event_conf_mean, "PIXEL_EVENT_CONF_MEAN")
    class_arr = pixel_conf_class.copy()
    class_arr[grass != 1] = NODATA_OUT
    core_bands.append(class_arr.astype(np.float32))
    core_names.append("PIXEL_CONF_CLASS")

    profile = grass_profile.copy()
    profile.update(driver="GTiff", dtype="float32", count=len(core_bands), nodata=NODATA_OUT,
                   compress="DEFLATE", predictor=3, tiled=True, blockxsize=256, blockysize=256, BIGTIFF="IF_SAFER")
    final_tif.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(final_tif, "w", **profile) as dst:
        for i, (arr, name) in enumerate(zip(core_bands, core_names), start=1):
            dst.write(arr.astype(np.float32), i)
            dst.set_band_description(i, name)
        dst.update_tags(
            PRODUCT="Grassland mowing final product", VERSION="v0.3",
            PATCH_ID=str(patch_id), YEAR=str(YEAR), STATUS="Operational evidence score; not calibrated probability",
            EVENT_CONFIDENCE_FORMULA="0.30*SIGNAL_QUALITY + 0.25*OBS_DENSITY + 0.20*CHILD_SUPPORT + 0.25*SPECTRAL_SCORE (renormalized if missing)",
            SPECTRAL_SCORE="0=none;50=one of NDII_down/SWIR_up;100=both",
            CONF_CLASS_1="low: score <40", CONF_CLASS_2="moderate: 40<=score<60",
            CONF_CLASS_3="high: 60<=score<80 or >=80 without both spectral directions",
            CONF_CLASS_4="very high: score>=80 and NDII_down+SWIR_up",
            CHILD_SUPPORT="percent of 4 child 10m pixels with event within +/-5 days",
            NDII_DEFINITION="(B8A-B11)/(B8A+B11), explicit 20m cube",
            SWIRMEAN_DEFINITION="(B11+B12)/2, explicit 20m cube",
            MILETICH_NOTE="NDVI drop >0.20 + SWIR rise retained as separate diagnostic flag",
            CALIBRATION_NOTE="Confidence formula/classes must be recalibrated with PlanetScope/reference data",
            GRASS_MASK_ROLE="detection +5m",
            GRASS_MASK_SHA256=grass_mask_signature(grass),
        )

    diag_bands, diag_names = [], []
    def add_diag(arr, name):
        a = np.asarray(arr, dtype=np.float32).copy()
        a[~np.isfinite(a)] = NODATA_OUT
        a[grass != 1] = NODATA_OUT
        diag_bands.append(a)
        diag_names.append(name)

    for i in range(MAX_EVENTS): add_diag(ndvi_drops[i], f"NDVI_DROP_{i+1}")
    for i in range(MAX_EVENTS): add_diag(neighbour_support[i], f"NEIGHBOUR_SUPPORT_{i+1}")
    for i in range(MAX_EVENTS): add_diag(ndii_drop_out[i], f"NDII_DROP_{i+1}")
    for i in range(MAX_EVENTS): add_diag(swir_rise_out[i], f"SWIR_RISE_{i+1}")
    for i in range(MAX_EVENTS): add_diag(child_support_out[i], f"CHILD_SUPPORT_{i+1}")
    for i in range(MAX_EVENTS): add_diag(spectral_code_out[i], f"SPECTRAL_CODE_{i+1}")
    for i in range(MAX_EVENTS): add_diag(miletich_flag_out[i], f"MILETICH_FLAG_{i+1}")
    add_diag(observation_support, "OBSERVATION_SUPPORT")
    add_diag(baseline_confidence, "BASELINE_CONFIDENCE")

    diag_profile = grass_profile.copy()
    diag_profile.update(driver="GTiff", dtype="float32", count=len(diag_bands), nodata=NODATA_OUT,
                        compress="DEFLATE", predictor=3, tiled=True, blockxsize=256, blockysize=256, BIGTIFF="IF_SAFER")
    with rasterio.open(diagnostics_tif, "w", **diag_profile) as dst:
        for i, (arr, name) in enumerate(zip(diag_bands, diag_names), start=1):
            dst.write(arr.astype(np.float32), i)
            dst.set_band_description(i, name)
        dst.update_tags(
            PRODUCT="Grassland mowing diagnostic companion", VERSION="v0.3", PATCH_ID=str(patch_id), YEAR=str(YEAR),
            SPECTRAL_CODE="0=none;1=NDII_down;2=SWIR_up;3=both",
            MILETICH_FLAG="1 when NDVI_DROP>0.20 and SWIR_RISE>0",
            EVENT_RESPONSE_PRE_WINDOW=f"{PRE_START_DAYS} to {PRE_END_DAYS} days before event",
            EVENT_RESPONSE_EVENT_WINDOW=f"-{EVENT_BEFORE_DAYS} to +{EVENT_AFTER_DAYS} days",
            STATUS="Diagnostic evidence; thresholds not probability-calibrated",
            GRASS_MASK_ROLE="detection +5m",
            GRASS_MASK_SHA256=grass_mask_signature(grass),
        )

    events_final.to_csv(events_csv, index=False)

    with rasterio.open(final_tif) as src:
        assert src.crs == grass_crs and src.transform == grass_transform and src.bounds == grass_bounds and src.count == 17
    with rasterio.open(diagnostics_tif) as src:
        assert src.crs == grass_crs and src.transform == grass_transform and src.bounds == grass_bounds and src.count == 44

    return events_final


def patch_qa(
    patch_id: str,
    patch_type,
    detection_grass_file: Path,
    stats_grass_file: Path | None,
    baseline_tif: Path,
    dedup20_csv: Path,
    events_csv: Path,
    *,
    has_stats_grass: bool,
):
    """QA/statistics using the strict -7.5 m mask.

    Detection/final rasters remain based on the +5 m processing mask. The
    statistics mask is independent and must be generated from the ORIGINAL
    grass polygons with an inner -7.5 m buffer. If the strict mask is empty,
    statistics are reported as empty/NaN rather than falling back to the
    detection mask.
    """
    detection, _, det_transform, det_crs, det_bounds, det_width, det_height = load_grass(
        detection_grass_file
    )

    if has_stats_grass:
        if stats_grass_file is None:
            raise RuntimeError("has_stats_grass=True but no statistics mask path was provided")
        stats, _, stat_transform, stat_crs, stat_bounds, stat_width, stat_height = load_grass(
            stats_grass_file
        )
        if (
            stat_transform != det_transform
            or stat_crs != det_crs
            or stat_bounds != det_bounds
            or stat_width != det_width
            or stat_height != det_height
        ):
            raise RuntimeError("Statistics mask grid does not match detection mask grid.")
        if np.any((stats == 1) & (detection != 1)):
            raise RuntimeError("Statistics mask is not a subset of the +5 m detection mask.")
        stats_mask = stats == 1
    else:
        stats_mask = np.zeros_like(detection, dtype=bool)

    arrays, _ = load_baseline_arrays(baseline_tif)
    mc = arrays["mowing_count"]
    vals = mc[stats_mask]
    vals = vals[np.isfinite(vals)]

    detection_pixels = int((detection == 1).sum())
    stats_pixels = int(stats_mask.sum())
    qa = {
        "PATCH_ID": patch_id,
        "TYPE": patch_type if pd.notna(patch_type) else "",
        # Backwards-compatible field now intentionally means STRICT statistics pixels.
        "GRASS_PIXELS": stats_pixels,
        "DETECTION_GRASS_PIXELS": detection_pixels,
        "STATS_GRASS_PIXELS": stats_pixels,
        "STATS_MASK_BUFFER_M": -7.5,
    }

    if stats_pixels:
        qa["VALID_OBS_MEDIAN"] = float(np.nanmedian(arrays["valid_obs"][stats_mask]))
        qa["MAX_GAP_MEDIAN"] = float(np.nanmedian(arrays["max_gap_days"][stats_mask]))
        denom = len(vals) if len(vals) else 1
        for k in range(5):
            qa[f"MOWING_{k}_PCT"] = 100.0 * np.sum(vals == k) / denom
        qa["MOWING_5PLUS_PCT"] = 100.0 * np.sum(vals >= 5) / denom
    else:
        qa["VALID_OBS_MEDIAN"] = np.nan
        qa["MAX_GAP_MEDIAN"] = np.nan
        for k in range(5):
            qa[f"MOWING_{k}_PCT"] = np.nan
        qa["MOWING_5PLUS_PCT"] = np.nan

    # 20 m QA: only cells whose four 10 m child pixels are all inside the
    # strict statistics mask. This avoids re-introducing boundary pixels.
    if dedup20_csv.exists():
        d_all = pd.read_csv(dedup20_csv)
        qa["EVENTS_20M_DETECTION"] = len(d_all)
        if stats_pixels and len(d_all):
            strict20 = stats_mask.reshape(250, 2, 250, 2).all(axis=(1, 3))
            rr = pd.to_numeric(d_all["row20"], errors="coerce").fillna(-1).astype(int)
            cc = pd.to_numeric(d_all["col20"], errors="coerce").fillna(-1).astype(int)
            valid_idx = (rr >= 0) & (rr < 250) & (cc >= 0) & (cc < 250)
            keep = np.zeros(len(d_all), dtype=bool)
            valid_positions = np.flatnonzero(valid_idx.to_numpy())
            keep[valid_positions] = strict20[
                rr.to_numpy()[valid_positions], cc.to_numpy()[valid_positions]
            ]
            d = d_all.loc[keep].copy()
        else:
            d = d_all.iloc[0:0].copy()

        qa["EVENTS_20M"] = len(d)
        if len(d):
            qa["NDII_DECREASE_PCT"] = 100.0 * (d["ndii_drop"] > 0).mean()
            qa["SWIR_INCREASE_PCT"] = 100.0 * (d["swir_rise"] > 0).mean()
            qa["BOTH_PCT"] = 100.0 * ((d["ndii_drop"] > 0) & (d["swir_rise"] > 0)).mean()
            strong = d[
                (d["neighbour_support_median"] >= 75)
                & (d["provisional_confidence_median"] >= 80)
            ]
            qa["STRONG_EVENTS_20M"] = len(strong)
            qa["STRONG_BOTH_PCT"] = (
                100.0 * ((strong["ndii_drop"] > 0) & (strong["swir_rise"] > 0)).mean()
                if len(strong) else np.nan
            )
        else:
            qa.update({
                "NDII_DECREASE_PCT": np.nan,
                "SWIR_INCREASE_PCT": np.nan,
                "BOTH_PCT": np.nan,
                "STRONG_EVENTS_20M": 0,
                "STRONG_BOTH_PCT": np.nan,
            })

    # Final event statistics: filter the 10 m event table with the strict mask.
    if events_csv.exists():
        e_all = pd.read_csv(events_csv)
        qa["FINAL_EVENTS_DETECTION"] = len(e_all)
        if stats_pixels and len(e_all):
            rr = pd.to_numeric(e_all["row"], errors="coerce").fillna(-1).astype(int)
            cc = pd.to_numeric(e_all["col"], errors="coerce").fillna(-1).astype(int)
            valid_idx = (
                (rr >= 0) & (rr < stats_mask.shape[0])
                & (cc >= 0) & (cc < stats_mask.shape[1])
            )
            keep = np.zeros(len(e_all), dtype=bool)
            valid_positions = np.flatnonzero(valid_idx.to_numpy())
            keep[valid_positions] = stats_mask[
                rr.to_numpy()[valid_positions], cc.to_numpy()[valid_positions]
            ]
            e = e_all.loc[keep].copy()
        else:
            e = e_all.iloc[0:0].copy()

        qa["FINAL_EVENTS"] = len(e)
        qa["EVENT_CONF_MEAN"] = float(e["event_confidence"].mean()) if len(e) else np.nan
        qa["EVENT_CONF_MEDIAN"] = float(e["event_confidence"].median()) if len(e) else np.nan
        qa["MILETICH_FLAG_PCT"] = 100.0 * e["miletich_flag"].mean() if len(e) else np.nan

    return qa
