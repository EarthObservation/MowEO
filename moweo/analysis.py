from __future__ import annotations

"""Mowing analysis for the configured study area and year.

Inputs
------
- national_qa_<YEAR>.csv
- per-EOPatch final event CSV files from process_eopatch/
- strict -7.5 m statistics masks
- optional EOPatch_PET.gpkg for patch-level maps

The event scan is streaming and resumable. Per-patch aggregates are cached under
national_analysis_<YEAR>/cache/events so the 42+ million strict-mask events never
need to be concatenated in memory.

Run from the national project folder:
    python 12_analyze.py

Fast QA-only run (no event CSV scan):
    python 12_analyze.py --skip-events

Force rebuilding event caches:
    python 12_analyze.py --rebuild-event-cache
"""

import argparse
import calendar
import math
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import rasterio

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover
    def tqdm(x, **kwargs):
        return x

try:
    import geopandas as gpd
except Exception:
    gpd = None

try:
    from . import config as cfg
except Exception:
    cfg = None


MOWING_COLS = [
    "MOWING_0_PCT", "MOWING_1_PCT", "MOWING_2_PCT",
    "MOWING_3_PCT", "MOWING_4_PCT", "MOWING_5PLUS_PCT",
]
MOWING_LABELS = ["0", "1", "2", "3", "4", "5+"]
SPECTRAL_COLS = ["NDII_DECREASE_PCT", "SWIR_INCREASE_PCT", "BOTH_PCT"]
SPECTRAL_LABELS = ["NDII decrease", "SWIR increase", "both"]


def atomic_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".part")
    df.to_csv(tmp, index=False, encoding="utf-8")
    tmp.replace(path)


def savefig(path: Path, dpi: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close()


def as_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def finite_number(value, default=np.nan):
    try:
        v = float(value)
        return v if np.isfinite(v) else default
    except Exception:
        return default


def weighted_percent(frame: pd.DataFrame, pct_col: str, weight_col: str) -> float:
    if pct_col not in frame.columns or weight_col not in frame.columns:
        return np.nan
    p = pd.to_numeric(frame[pct_col], errors="coerce")
    w = pd.to_numeric(frame[weight_col], errors="coerce")
    valid = p.notna() & w.notna() & (w > 0)
    if not valid.any():
        return np.nan
    return float(np.sum(p[valid] * w[valid]) / np.sum(w[valid]))


def weighted_mean(frame: pd.DataFrame, value_col: str, weight_col: str) -> float:
    if value_col not in frame.columns or weight_col not in frame.columns:
        return np.nan
    v = pd.to_numeric(frame[value_col], errors="coerce")
    w = pd.to_numeric(frame[weight_col], errors="coerce")
    valid = v.notna() & w.notna() & (w > 0)
    if not valid.any():
        return np.nan
    return float(np.sum(v[valid] * w[valid]) / np.sum(w[valid]))


def group_summary(frame: pd.DataFrame, group_name: str, group_value) -> dict:
    det_pixels = float(pd.to_numeric(frame["DETECTION_GRASS_PIXELS"], errors="coerce").fillna(0).sum())
    stats_pixels = float(pd.to_numeric(frame["STATS_GRASS_PIXELS"], errors="coerce").fillna(0).sum())
    final_det = float(pd.to_numeric(frame["FINAL_EVENTS_DETECTION"], errors="coerce").fillna(0).sum())
    final_stats = float(pd.to_numeric(frame["FINAL_EVENTS"], errors="coerce").fillna(0).sum())
    events20 = float(pd.to_numeric(frame["EVENTS_20M"], errors="coerce").fillna(0).sum())
    strong20 = float(pd.to_numeric(frame["STRONG_EVENTS_20M"], errors="coerce").fillna(0).sum())
    row = {
        group_name: group_value,
        "PATCHES": int(len(frame)),
        "DETECTION_GRASS_PIXELS": int(det_pixels),
        "STATS_GRASS_PIXELS": int(stats_pixels),
        "STATS_MASK_RETENTION_PCT": 100.0 * stats_pixels / det_pixels if det_pixels else np.nan,
        "FINAL_EVENTS_DETECTION": int(final_det),
        "FINAL_EVENTS": int(final_stats),
        "EVENT_RETENTION_PCT": 100.0 * final_stats / final_det if final_det else np.nan,
        "EVENTS_PER_STATS_PIXEL": final_stats / stats_pixels if stats_pixels else np.nan,
        "EVENTS_20M": int(events20),
        "STRONG_EVENTS_20M": int(strong20),
        "EVENT_CONF_MEAN": weighted_mean(frame, "EVENT_CONF_MEAN", "FINAL_EVENTS"),
        "MILETICH_FLAG_PCT": weighted_percent(frame, "MILETICH_FLAG_PCT", "FINAL_EVENTS"),
        "NDII_DECREASE_PCT": weighted_percent(frame, "NDII_DECREASE_PCT", "EVENTS_20M"),
        "SWIR_INCREASE_PCT": weighted_percent(frame, "SWIR_INCREASE_PCT", "EVENTS_20M"),
        "BOTH_PCT": weighted_percent(frame, "BOTH_PCT", "EVENTS_20M"),
        "STRONG_BOTH_PCT": weighted_percent(frame, "STRONG_BOTH_PCT", "STRONG_EVENTS_20M"),
        "OVERFLOW_PIXELS_GT6": int(pd.to_numeric(frame["OVERFLOW_PIXELS_GT6"], errors="coerce").fillna(0).sum()),
        "VALID_OBS_PATCH_MEDIAN": float(pd.to_numeric(frame["VALID_OBS_MEDIAN"], errors="coerce").median()),
        "MAX_GAP_PATCH_MEDIAN": float(pd.to_numeric(frame["MAX_GAP_MEDIAN"], errors="coerce").median()),
        "TOTAL_COST": float(pd.to_numeric(frame["TOTAL_COST"], errors="coerce").sum(min_count=1)),
    }
    for col, label in zip(MOWING_COLS, MOWING_LABELS):
        key = "MOWING_5PLUS_PCT" if label == "5+" else f"MOWING_{label}_PCT"
        row[key] = weighted_percent(frame, col, "STATS_GRASS_PIXELS")
    return row


def grouped_summary(ok: pd.DataFrame, col: str) -> pd.DataFrame:
    rows = []
    for value, group in ok.dropna(subset=[col]).groupby(col, dropna=True):
        rows.append(group_summary(group, col, value))
    return pd.DataFrame(rows)


def month_lookup(year: int) -> np.ndarray:
    lut = np.zeros(367, dtype=np.int16)
    start = pd.Timestamp(year=year, month=1, day=1)
    for doy in range(1, 367):
        try:
            lut[doy] = (start + pd.Timedelta(days=doy - 1)).month
        except Exception:
            lut[doy] = 0
    return lut


def cache_is_current(path: Path, event_file: Path, mask_file: Path) -> bool:
    if not path.exists():
        return False
    try:
        z = np.load(path, allow_pickle=False)
        return (
            int(z["event_size"]) == event_file.stat().st_size
            and int(z["event_mtime_ns"]) == event_file.stat().st_mtime_ns
            and int(z["mask_size"]) == mask_file.stat().st_size
            and int(z["mask_mtime_ns"]) == mask_file.stat().st_mtime_ns
        )
    except Exception:
        return False


def event_file_for(root: Path, patch_id: str, year: int) -> Path:
    return root / "process_eopatch" / patch_id / f"{patch_id}_MOWING_FINAL_events_{cfg.FINAL_VERSION}_{year}.csv"


def stats_mask_for(stats_root: Path, patch_id: str) -> Path:
    return stats_root / f"{patch_id}_grass.tif"


def aggregate_patch_events(
    patch_id: str,
    event_file: Path,
    mask_file: Path,
    cache_file: Path,
    year: int,
    rebuild: bool,
) -> dict:
    if not rebuild and cache_is_current(cache_file, event_file, mask_file):
        z = np.load(cache_file, allow_pickle=False)
        return {k: z[k] for k in z.files}

    with rasterio.open(mask_file) as src:
        stats_mask = src.read(1) == 1
    height, width = stats_mask.shape

    daily = np.zeros(367, dtype=np.int64)
    daily_order = np.zeros((7, 367), dtype=np.int64)
    monthly = np.zeros(13, dtype=np.int64)
    conf_hist = np.zeros(101, dtype=np.int64)
    class_counts = np.zeros(5, dtype=np.int64)
    spectral_counts = np.zeros(4, dtype=np.int64)
    order_counts = np.zeros(7, dtype=np.int64)
    miletich_counts = np.zeros(2, dtype=np.int64)
    month_lut = month_lookup(year)
    total_strict = 0
    invalid_rows = 0

    wanted = [
        "row", "col", "event", "doy", "event_confidence", "event_class",
        "spectral_code", "miletich_flag",
    ]
    header = pd.read_csv(event_file, nrows=0)
    missing = set(wanted) - set(header.columns)
    if missing:
        raise RuntimeError(f"{event_file}: missing event columns {sorted(missing)}")

    for chunk in pd.read_csv(event_file, usecols=wanted, chunksize=250_000):
        rr = pd.to_numeric(chunk["row"], errors="coerce").fillna(-1).astype(np.int32).to_numpy()
        cc = pd.to_numeric(chunk["col"], errors="coerce").fillna(-1).astype(np.int32).to_numpy()
        valid = (rr >= 0) & (rr < height) & (cc >= 0) & (cc < width)
        keep = np.zeros(len(chunk), dtype=bool)
        pos = np.flatnonzero(valid)
        if len(pos):
            keep[pos] = stats_mask[rr[pos], cc[pos]]
        invalid_rows += int((~valid).sum())
        if not keep.any():
            continue
        s = chunk.loc[keep]
        total_strict += len(s)

        doy = pd.to_numeric(s["doy"], errors="coerce").fillna(0).astype(np.int16).to_numpy()
        ev = pd.to_numeric(s["event"], errors="coerce").fillna(0).astype(np.int16).to_numpy()
        conf = pd.to_numeric(s["event_confidence"], errors="coerce").to_numpy(dtype=float)
        cls = pd.to_numeric(s["event_class"], errors="coerce").fillna(0).astype(np.int16).to_numpy()
        spec = pd.to_numeric(s["spectral_code"], errors="coerce").fillna(0).astype(np.int16).to_numpy()
        mil = pd.to_numeric(s["miletich_flag"], errors="coerce").fillna(0).astype(np.int16).to_numpy()

        good_doy = (doy >= 1) & (doy <= 366)
        np.add.at(daily, doy[good_doy], 1)
        months = month_lut[doy[good_doy]]
        np.add.at(monthly, months, 1)

        good_order = (ev >= 1) & (ev <= 6)
        np.add.at(order_counts, ev[good_order], 1)
        both = good_doy & good_order
        np.add.at(daily_order, (ev[both], doy[both]), 1)

        good_conf = np.isfinite(conf)
        conf_i = np.clip(np.floor(conf[good_conf]).astype(np.int16), 0, 100)
        np.add.at(conf_hist, conf_i, 1)

        good_cls = (cls >= 0) & (cls <= 4)
        np.add.at(class_counts, cls[good_cls], 1)
        good_spec = (spec >= 0) & (spec <= 3)
        np.add.at(spectral_counts, spec[good_spec], 1)
        mil_i = np.clip(mil, 0, 1)
        np.add.at(miletich_counts, mil_i, 1)

    data = {
        "event_size": np.int64(event_file.stat().st_size),
        "event_mtime_ns": np.int64(event_file.stat().st_mtime_ns),
        "mask_size": np.int64(mask_file.stat().st_size),
        "mask_mtime_ns": np.int64(mask_file.stat().st_mtime_ns),
        "total_strict": np.int64(total_strict),
        "invalid_rows": np.int64(invalid_rows),
        "daily": daily,
        "daily_order": daily_order,
        "monthly": monthly,
        "conf_hist": conf_hist,
        "class_counts": class_counts,
        "spectral_counts": spectral_counts,
        "order_counts": order_counts,
        "miletich_counts": miletich_counts,
    }
    cache_file.parent.mkdir(parents=True, exist_ok=True)
    tmp = cache_file.with_suffix(".part.npz")
    np.savez_compressed(tmp, **data)
    tmp.replace(cache_file)
    return data


def combine_event_caches(
    ok: pd.DataFrame,
    root: Path,
    stats_root: Path,
    cache_root: Path,
    year: int,
    rebuild: bool,
    max_patches: int | None,
) -> tuple[dict, pd.DataFrame]:
    totals = {
        "daily": np.zeros(367, dtype=np.int64),
        "daily_order": np.zeros((7, 367), dtype=np.int64),
        "monthly": np.zeros(13, dtype=np.int64),
        "conf_hist": np.zeros(101, dtype=np.int64),
        "class_counts": np.zeros(5, dtype=np.int64),
        "spectral_counts": np.zeros(4, dtype=np.int64),
        "order_counts": np.zeros(7, dtype=np.int64),
        "miletich_counts": np.zeros(2, dtype=np.int64),
    }
    checks = []
    rows = ok.sort_values("PATCH_ID")
    if max_patches is not None:
        rows = rows.head(max_patches)

    for row in tqdm(rows.itertuples(index=False), total=len(rows), desc="Strict-mask event scan"):
        patch_id = str(row.PATCH_ID)
        event_file = event_file_for(root, patch_id, year)
        mask_file = stats_mask_for(stats_root, patch_id)
        expected = int(finite_number(getattr(row, "FINAL_EVENTS", 0), 0))
        if not event_file.exists():
            checks.append({"PATCH_ID": patch_id, "STATUS": "MISSING_EVENT_CSV", "EXPECTED": expected, "SCANNED": np.nan})
            continue
        if not mask_file.exists():
            checks.append({"PATCH_ID": patch_id, "STATUS": "MISSING_STATS_MASK", "EXPECTED": expected, "SCANNED": np.nan})
            continue
        cache_file = cache_root / f"{patch_id}.npz"
        data = aggregate_patch_events(patch_id, event_file, mask_file, cache_file, year, rebuild)
        scanned = int(data["total_strict"])
        checks.append({
            "PATCH_ID": patch_id,
            "STATUS": "OK" if scanned == expected else "COUNT_MISMATCH",
            "EXPECTED": expected,
            "SCANNED": scanned,
            "INVALID_ROWCOL": int(data["invalid_rows"]),
        })
        for key in totals:
            totals[key] += data[key]
    return totals, pd.DataFrame(checks)


def write_event_tables(totals: dict, tables: Path, year: int) -> None:
    total_events = int(totals["daily"].sum())
    dates = pd.to_datetime(f"{year}-01-01") + pd.to_timedelta(np.arange(366), unit="D")
    daily = pd.DataFrame({
        "DOY": np.arange(1, 367),
        "DATE": dates.strftime("%Y-%m-%d"),
        "EVENTS": totals["daily"][1:367],
    })
    for ev in range(1, 7):
        daily[f"EVENT_{ev}"] = totals["daily_order"][ev, 1:367]
    atomic_csv(daily, tables / f"event_doy_daily_{year}.csv")

    month_names = [calendar.month_name[m] for m in range(1, 13)]
    monthly = pd.DataFrame({
        "MONTH": np.arange(1, 13),
        "MONTH_NAME": month_names,
        "EVENTS": totals["monthly"][1:13],
    })
    monthly["PCT"] = 100 * monthly["EVENTS"] / max(monthly["EVENTS"].sum(), 1)
    atomic_csv(monthly, tables / f"event_monthly_{year}.csv")

    conf = pd.DataFrame({"CONFIDENCE_FLOOR": np.arange(101), "EVENTS": totals["conf_hist"]})
    conf["PCT"] = 100 * conf["EVENTS"] / max(conf["EVENTS"].sum(), 1)
    atomic_csv(conf, tables / f"event_confidence_hist_{year}.csv")

    classes = pd.DataFrame({
        "EVENT_CLASS": np.arange(5),
        "LABEL": ["0/unavailable", "1 low", "2 moderate", "3 high", "4 very high"],
        "EVENTS": totals["class_counts"],
    })
    classes["PCT"] = 100 * classes["EVENTS"] / max(classes["EVENTS"].sum(), 1)
    atomic_csv(classes, tables / f"event_class_distribution_{year}.csv")

    spectral = pd.DataFrame({
        "SPECTRAL_CODE": np.arange(4),
        "LABEL": ["0 none", "1 NDII only", "2 SWIR only", "3 both"],
        "EVENTS": totals["spectral_counts"],
    })
    spectral["PCT"] = 100 * spectral["EVENTS"] / max(spectral["EVENTS"].sum(), 1)
    atomic_csv(spectral, tables / f"spectral_code_distribution_{year}.csv")

    order = pd.DataFrame({"EVENT_ORDER": np.arange(1, 7), "EVENTS": totals["order_counts"][1:7]})
    order["PCT"] = 100 * order["EVENTS"] / max(order["EVENTS"].sum(), 1)
    atomic_csv(order, tables / f"event_order_distribution_{year}.csv")

    mil = pd.DataFrame({"MILETICH_FLAG": [0, 1], "EVENTS": totals["miletich_counts"]})
    mil["PCT"] = 100 * mil["EVENTS"] / max(mil["EVENTS"].sum(), 1)
    atomic_csv(mil, tables / f"miletich_flag_distribution_{year}.csv")


def make_qa_figures(ok: pd.DataFrame, out: Path, dpi: int, year: int) -> None:
    # 01 national mowing distribution
    vals = [weighted_percent(ok, c, "STATS_GRASS_PIXELS") for c in MOWING_COLS]
    plt.figure(figsize=(8, 5))
    plt.bar(MOWING_LABELS, vals)
    plt.xlabel("Number of detected mowings")
    plt.ylabel("Share of inner grassland pixels (%)")
    plt.title(f"Mowing-count distribution, {year}")
    savefig(out / "01_mowing_distribution_national.png", dpi)

    # 02 by biogeographic region
    bio = grouped_summary(ok, "BIOGEO_R")
    if len(bio):
        x = np.arange(len(MOWING_LABELS))
        width = 0.8 / max(len(bio), 1)
        plt.figure(figsize=(10, 5.5))
        for i, row in bio.reset_index(drop=True).iterrows():
            y = [row["MOWING_5PLUS_PCT" if lab == "5+" else f"MOWING_{lab}_PCT"] for lab in MOWING_LABELS]
            plt.bar(x + (i - (len(bio)-1)/2) * width, y, width=width, label=str(row["BIOGEO_R"]))
        plt.xticks(x, MOWING_LABELS)
        plt.xlabel("Number of detected mowings")
        plt.ylabel("Share of inner grassland pixels (%)")
        plt.title(f"Mowing counts by region, {year}")
        plt.legend()
        savefig(out / "02_mowing_distribution_biogeo.png", dpi)

        x = np.arange(len(SPECTRAL_LABELS))
        width = 0.8 / max(len(bio), 1)
        plt.figure(figsize=(9, 5.5))
        for i, row in bio.reset_index(drop=True).iterrows():
            y = [row[c] for c in SPECTRAL_COLS]
            plt.bar(x + (i - (len(bio)-1)/2) * width, y, width=width, label=str(row["BIOGEO_R"]))
        plt.xticks(x, SPECTRAL_LABELS)
        plt.ylabel("Share of 20 m events (%)")
        plt.title(f"Spectral confirmation of events by region, {year}")
        plt.legend()
        savefig(out / "03_spectral_confirmation_biogeo.png", dpi)

    # 04 PET spectral confirmation
    pet = grouped_summary(ok, "PET")
    if len(pet):
        pet = pet.sort_values("PET")
        x = np.arange(len(pet))
        width = 0.25
        plt.figure(figsize=(12, 5.5))
        for i, (col, label) in enumerate(zip(SPECTRAL_COLS, SPECTRAL_LABELS)):
            plt.bar(x + (i - 1) * width, pet[col], width=width, label=label)
        plt.xticks(x, pet['PET'].astype(str).tolist())
        plt.xlabel("PET")
        plt.ylabel("Share of 20 m events (%)")
        plt.title(f"Spectral confirmation by PET group, {year}")
        plt.legend()
        savefig(out / "04_spectral_confirmation_pet.png", dpi)

    # 05 confidence by biogeo
    if len(bio):
        plt.figure(figsize=(7, 5))
        plt.bar(bio["BIOGEO_R"].astype(str), bio["EVENT_CONF_MEAN"])
        plt.ylim(0, 100)
        plt.ylabel("Mean evidence score")
        plt.title(f"Mean evidence score by region, {year}")
        savefig(out / "05_event_confidence_biogeo.png", dpi)

    # 06 observation quality
    plt.figure(figsize=(8, 6))
    for name, group in ok.groupby("BIOGEO_R", dropna=False):
        plt.scatter(group["VALID_OBS_MEDIAN"], group["MAX_GAP_MEDIAN"], s=14, alpha=0.55, label=str(name))
    plt.xlabel("Median valid observations per tile")
    plt.ylabel("Median maximum observation gap (days)")
    plt.title(f"Time-series quality by tile, {year}")
    plt.legend()
    savefig(out / "06_observation_quality.png", dpi)

    # 07 detection vs strict mask
    det = pd.to_numeric(ok["DETECTION_GRASS_PIXELS"], errors="coerce")
    sta = pd.to_numeric(ok["STATS_GRASS_PIXELS"], errors="coerce")
    plt.figure(figsize=(7, 6))
    plt.scatter(det, sta, s=12, alpha=0.5)
    mx = float(np.nanmax(det)) if len(det) else 1
    plt.plot([0, mx], [0, mx], linestyle="--", linewidth=1)
    plt.xlabel("Pixels in the +5 m detection mask")
    plt.ylabel("Pixels in the -7.5 m statistics mask")
    plt.title("Effect of the inner statistics mask")
    savefig(out / "07_detection_vs_stats_mask.png", dpi)

    # 08 mask retention distribution
    retention = 100 * sta / det.replace(0, np.nan)
    plt.figure(figsize=(8, 5))
    plt.hist(retention.dropna(), bins=30)
    plt.xlabel("Pixels retained after the -7.5 m inner buffer (%)")
    plt.ylabel("Number of tiles")
    plt.title("Grassland pixel retention in the statistics mask")
    savefig(out / "08_stats_mask_retention.png", dpi)

    # 09 costs
    costs = pd.to_numeric(ok["TOTAL_COST"], errors="coerce").dropna()
    plt.figure(figsize=(8, 5))
    plt.hist(costs, bins=np.arange(math.floor(costs.min()) if len(costs) else 0,
                                   math.ceil(costs.max()) + 2 if len(costs) else 2) - 0.5)
    plt.xlabel("openEO credits per tile")
    plt.ylabel("Number of tiles")
    plt.title(f"Acquisition cost distribution, {year}")
    savefig(out / "09_acquisition_costs.png", dpi)


def make_event_figures(totals: dict, out: Path, dpi: int, year: int) -> None:
    # 10 monthly
    months = np.arange(1, 13)
    counts = totals["monthly"][1:13]
    plt.figure(figsize=(10, 5))
    plt.bar(months, counts)
    plt.xticks(months, [calendar.month_abbr[m] for m in months])
    plt.ylabel("Number of events")
    plt.title(f"Monthly detected mowing events, {year} (-7.5 m mask)")
    savefig(out / "10_event_dates_monthly.png", dpi)

    # 11 daily
    doys = np.arange(1, 367)
    plt.figure(figsize=(12, 5))
    plt.plot(doys, totals["daily"][1:367], linewidth=1)
    plt.xlim(1, 366)
    plt.xlabel("Day of year (DOY)")
    plt.ylabel("Number of events")
    plt.title(f"Detected mowing dates, {year}")
    savefig(out / "11_event_dates_doy.png", dpi)

    # 12 daily by event order
    plt.figure(figsize=(12, 6))
    for ev in range(1, 7):
        if totals["order_counts"][ev] > 0:
            plt.plot(doys, totals["daily_order"][ev, 1:367], linewidth=1, label=f"event {ev}")
    plt.xlim(1, 366)
    plt.xlabel("Day of year (DOY)")
    plt.ylabel("Number of events")
    plt.title(f"Mowing date by event order, {year}")
    plt.legend(ncol=3)
    savefig(out / "12_event_dates_by_order.png", dpi)

    # 13 confidence histogram
    x = np.arange(101)
    plt.figure(figsize=(10, 5))
    plt.bar(x, totals["conf_hist"], width=1.0)
    plt.xlim(0, 100)
    plt.xlabel("Evidence score")
    plt.ylabel("Number of events")
    plt.title(f"Event evidence score distribution, {year}")
    savefig(out / "13_event_confidence_distribution.png", dpi)

    # 14 event classes
    classes = np.arange(1, 5)
    vals = totals["class_counts"][1:5]
    plt.figure(figsize=(8, 5))
    plt.bar(["Low", "Moderate", "High", "Very high"], vals)
    plt.ylabel("Number of events")
    plt.title(f"Evidence score classes, {year}")
    savefig(out / "14_event_class_distribution.png", dpi)

    # 15 spectral codes
    plt.figure(figsize=(8, 5))
    plt.bar(["none", "NDII", "SWIR", "both"], totals["spectral_counts"])
    plt.ylabel("Number of events")
    plt.title(f"Spectral confirmation of final events, {year}")
    savefig(out / "15_spectral_code_distribution.png", dpi)

    # 16 event order
    plt.figure(figsize=(8, 5))
    plt.bar(np.arange(1, 7), totals["order_counts"][1:7])
    plt.xlabel("Event order within a pixel")
    plt.ylabel("Number of events")
    plt.title(f"Events by order, {year}")
    savefig(out / "16_event_order_distribution.png", dpi)


def make_maps(ok: pd.DataFrame, gpkg: Path, out: Path, dpi: int, year: int) -> None:
    if gpd is None or gpkg is None or not gpkg.exists():
        return
    try:
        grid = gpd.read_file(gpkg)
    except Exception as exc:
        print(f"WARNING: map skipped; cannot read {gpkg}: {exc}")
        return
    if "PATCH_ID" not in grid.columns:
        print("WARNING: map skipped; GPKG has no PATCH_ID column")
        return

    cols = ["PATCH_ID", "MOWING_0_PCT", "MOWING_1_PCT", "MOWING_2_PCT", "MOWING_3_PCT",
            "MOWING_4_PCT", "MOWING_5PLUS_PCT", "EVENT_CONF_MEAN", "STATS_GRASS_PIXELS"]
    a = ok[cols].copy()
    a["MOWING_2PLUS_PCT"] = a[["MOWING_2_PCT", "MOWING_3_PCT", "MOWING_4_PCT", "MOWING_5PLUS_PCT"]].sum(axis=1, min_count=1)
    m = grid.merge(a[["PATCH_ID", "MOWING_2PLUS_PCT", "EVENT_CONF_MEAN", "STATS_GRASS_PIXELS"]], on="PATCH_ID", how="left")

    fig, ax = plt.subplots(figsize=(10, 8))
    m.plot(column="MOWING_2PLUS_PCT", ax=ax, legend=True, missing_kwds={"color": "lightgrey"})
    ax.set_axis_off()
    ax.set_title(f"Share of inner grassland pixels with at least two mowings, {year}")
    fig.tight_layout()
    fig.savefig(out / "17_map_mowing_2plus_pct.png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 8))
    m.plot(column="EVENT_CONF_MEAN", ax=ax, legend=True, missing_kwds={"color": "lightgrey"})
    ax.set_axis_off()
    ax.set_title(f"Mean evidence score by tile, {year}")
    fig.tight_layout()
    fig.savefig(out / "18_map_event_confidence.png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def write_report_summary(path: Path, national: dict, event_total: int | None, year: int) -> None:
    lines = [
        f"NATIONAL MOWING ANALYSIS {year}",
        "=" * 72,
        f"OK EOPatches: {national['PATCHES']}",
        f"Detection grass pixels (+5 m): {national['DETECTION_GRASS_PIXELS']:,}",
        f"Statistics grass pixels (-7.5 m): {national['STATS_GRASS_PIXELS']:,}",
        f"Statistics-mask retention: {national['STATS_MASK_RETENTION_PCT']:.2f}%",
        f"Final detection events: {national['FINAL_EVENTS_DETECTION']:,}",
        f"Final strict-mask events: {national['FINAL_EVENTS']:,}",
        f"Event retention: {national['EVENT_RETENTION_PCT']:.2f}%",
        f"Events per strict grass pixel: {national['EVENTS_PER_STATS_PIXEL']:.3f}",
        f"NDII decrease: {national['NDII_DECREASE_PCT']:.2f}%",
        f"SWIR increase: {national['SWIR_INCREASE_PCT']:.2f}%",
        f"Both spectral directions: {national['BOTH_PCT']:.2f}%",
        f"Strong-event BOTH: {national['STRONG_BOTH_PCT']:.2f}%",
        f"Mean event evidence score: {national['EVENT_CONF_MEAN']:.2f}",
        f"Miletich flag: {national['MILETICH_FLAG_PCT']:.2f}%",
        f"Overflow pixels >6 events: {national['OVERFLOW_PIXELS_GT6']:,}",
        f"Total acquisition cost: {national['TOTAL_COST']:.0f} credits",
    ]
    if event_total is not None:
        lines.append(f"Event rows independently rescanned inside strict mask: {event_total:,}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    if cfg is None:
        raise RuntimeError("Could not import config.py. Run 00_setup_project.py first.")
    default_root = Path(cfg.ROOT)
    default_year = int(cfg.YEAR)
    default_stats = Path(cfg.STATS_GRASS_ROOT)
    default_gpkg = Path(cfg.PATCH_GPKG)

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=default_root)
    parser.add_argument("--year", type=int, default=default_year)
    parser.add_argument("--qa", type=Path)
    parser.add_argument("--stats-root", type=Path, default=default_stats)
    parser.add_argument("--gpkg", type=Path, default=default_gpkg)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--skip-events", action="store_true", help="Create QA summaries/figures without scanning event CSVs.")
    parser.add_argument("--rebuild-event-cache", action="store_true")
    parser.add_argument("--max-patches", type=int, help="Testing only: scan first N OK patches.")
    parser.add_argument("--dpi", type=int, default=220)
    parser.add_argument("--no-maps", action="store_true")
    args = parser.parse_args()

    root = args.root
    year = args.year
    qa_path = args.qa or root / f"national_qa_{year}.csv"
    out = args.output or root / f"national_analysis_{year}"
    tables = out / "tables"
    figures = out / "figures"
    cache_root = out / "cache" / "events"
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    if not qa_path.exists():
        raise FileNotFoundError(f"QA file not found: {qa_path}")
    qa = pd.read_csv(qa_path)
    required = {
        "PATCH_ID", "STATUS", "DETECTION_GRASS_PIXELS", "STATS_GRASS_PIXELS",
        "FINAL_EVENTS_DETECTION", "FINAL_EVENTS",
    }
    missing = required - set(qa.columns)
    if missing:
        raise RuntimeError(f"QA file missing columns: {sorted(missing)}")

    qa["HAS_GRASS"] = as_bool(qa["HAS_GRASS"]) if "HAS_GRASS" in qa.columns else qa["STATUS"].eq("OK")
    qa["HAS_STATS_GRASS"] = as_bool(qa["HAS_STATS_GRASS"]) if "HAS_STATS_GRASS" in qa.columns else qa["STATS_GRASS_PIXELS"].fillna(0).gt(0)
    ok = qa[qa["STATUS"].astype(str).eq("OK")].copy()
    no_grass = qa[qa["STATUS"].astype(str).eq("NO_GRASS")].copy()
    other = qa[~qa["STATUS"].astype(str).isin(["OK", "NO_GRASS"])].copy()

    expected_patches = int(getattr(cfg, "EXPECTED_PATCHES", 1))
    if len(qa) != expected_patches:
        print(f"WARNING: QA has {len(qa)} rows, expected {expected_patches}.")
    if len(other):
        print(f"WARNING: {len(other)} rows have statuses other than OK/NO_GRASS.")

    national = group_summary(ok, "SCOPE", str(getattr(cfg, "STUDY_AREA_NAME", "Study area")).upper())
    national["NO_GRASS_PATCHES"] = int(len(no_grass))
    national["ALL_PATCHES"] = int(len(qa))
    national_df = pd.DataFrame([national])
    atomic_csv(national_df, tables / f"national_summary_{year}.csv")

    bio = grouped_summary(ok, "BIOGEO_R")
    pet = grouped_summary(ok, "PET")
    atomic_csv(bio, tables / f"summary_by_biogeo_{year}.csv")
    atomic_csv(pet, tables / f"summary_by_pet_{year}.csv")

    # Patch-level table + transparent derived fields.
    patch = qa.copy()
    patch["STATS_MASK_RETENTION_PCT"] = 100 * pd.to_numeric(patch["STATS_GRASS_PIXELS"], errors="coerce") / pd.to_numeric(patch["DETECTION_GRASS_PIXELS"], errors="coerce").replace(0, np.nan)
    patch["EVENT_RETENTION_PCT"] = 100 * pd.to_numeric(patch["FINAL_EVENTS"], errors="coerce") / pd.to_numeric(patch["FINAL_EVENTS_DETECTION"], errors="coerce").replace(0, np.nan)
    patch["MOWING_2PLUS_PCT"] = patch[["MOWING_2_PCT", "MOWING_3_PCT", "MOWING_4_PCT", "MOWING_5PLUS_PCT"]].sum(axis=1, min_count=1)
    atomic_csv(patch, tables / f"national_patch_table_{year}.csv")

    # QA flags are review flags, not automatic errors.
    flags = ok[
        ok["PET"].isna()
        | ok["BIOGEO_R"].isna()
        | (pd.to_numeric(ok["AREA_PCT"], errors="coerce") < 50)
        | (pd.to_numeric(ok["STATS_GRASS_PIXELS"], errors="coerce") < 100)
        | (pd.to_numeric(ok["OVERFLOW_PIXELS_GT6"], errors="coerce") > 0)
    ].copy()
    reasons = []
    for r in flags.itertuples(index=False):
        rr = []
        if pd.isna(getattr(r, "PET", np.nan)) or pd.isna(getattr(r, "BIOGEO_R", np.nan)):
            rr.append("missing_PET_or_BIOGEO")
        if finite_number(getattr(r, "AREA_PCT", np.nan)) < 50:
            rr.append("AREA_PCT_lt_50")
        if finite_number(getattr(r, "STATS_GRASS_PIXELS", np.nan), 0) < 100:
            rr.append("stats_pixels_lt_100")
        if finite_number(getattr(r, "OVERFLOW_PIXELS_GT6", np.nan), 0) > 0:
            rr.append("overflow_gt6")
        reasons.append(";".join(rr))
    flags["REVIEW_REASON"] = reasons
    atomic_csv(flags, tables / f"qa_review_flags_{year}.csv")

    make_qa_figures(ok, figures, args.dpi, year)
    if not args.no_maps:
        make_maps(ok, args.gpkg, figures, args.dpi, year)

    event_total = None
    if not args.skip_events:
        totals, checks = combine_event_caches(
            ok, root, args.stats_root, cache_root, year,
            args.rebuild_event_cache, args.max_patches,
        )
        atomic_csv(checks, tables / f"event_scan_checks_{year}.csv")
        write_event_tables(totals, tables, year)
        make_event_figures(totals, figures, args.dpi, year)
        event_total = int(totals["daily"].sum())
        mismatches = checks[checks["STATUS"] != "OK"]
        if len(mismatches):
            print(f"WARNING: {len(mismatches)} event-scan checks are not OK. See event_scan_checks_{year}.csv")
            raise RuntimeError('Event table scan is incomplete or inconsistent with QA; inspect event_scan_checks before publishing reports.')

    write_report_summary(out / f"national_summary_{year}.txt", national, event_total, year)

    print("=" * 78)
    print(f"NATIONAL ANALYSIS {year} COMPLETE")
    print("=" * 78)
    print(f"QA rows: {len(qa)} | OK: {len(ok)} | NO_GRASS: {len(no_grass)} | other: {len(other)}")
    print(f"Strict grass pixels: {national['STATS_GRASS_PIXELS']:,}")
    print(f"Strict-mask final events from QA: {national['FINAL_EVENTS']:,}")
    if event_total is not None:
        print(f"Strict-mask final events rescanned from event CSVs: {event_total:,}")
    print(f"Output: {out}")
    print(f"Tables: {tables}")
    print(f"Figures: {figures}")
    print("Event scan is resumable via per-patch caches in:", cache_root)


if __name__ == "__main__":
    main()
