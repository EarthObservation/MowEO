from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
import traceback

import numpy as np
import pandas as pd
import rasterio

from .config import (
    FINALIZE_LOG,
    NATIONAL_QA,
    OPENEO_LOG,
    PATCH_CSV,
    PATCH_ID_FIELD,
    PROCESS_LOG,
    STATS_GRASS_ROOT,
    STATS_GRASS_MANIFEST,
    TYPE_FIELD,
    patch_paths,
)
from .core import (
    baseline_is_consistent, finalize_patch, grass_mask_signature, load_grass, patch_qa, read_patch_table
)


def atomic_csv(df: pd.DataFrame, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.stem + ".part" + path.suffix)
    df.to_csv(tmp, index=False)
    tmp.replace(path)


def append_log(path: Path, row: dict):
    new = pd.DataFrame([row])
    if path.exists():
        try:
            old = pd.read_csv(path)
            new = pd.concat([old, new], ignore_index=True, sort=False)
        except Exception:
            pass
    atomic_csv(new, path)


def costs_by_patch():
    if not OPENEO_LOG.exists():
        return pd.DataFrame(columns=["PATCH_ID", "NDVI_COST", "INDICES20_COST", "TOTAL_COST"])
    j = pd.read_csv(OPENEO_LOG)
    if "costs" not in j.columns:
        return pd.DataFrame(columns=["PATCH_ID", "NDVI_COST", "INDICES20_COST", "TOTAL_COST"])
    j["costs"] = pd.to_numeric(j["costs"], errors="coerce")
    if "timestamp" in j.columns:
        j = j.sort_values("timestamp")
    j = j[j["status"].astype(str).str.startswith("finished")]
    j = j.drop_duplicates(subset=["PATCH_ID", "product"], keep="last")
    piv = j.pivot_table(index="PATCH_ID", columns="product", values="costs", aggfunc="last").reset_index()
    piv["NDVI_COST"] = piv["NDVI10"] if "NDVI10" in piv.columns else np.nan
    piv["INDICES20_COST"] = piv["NDII_SWIR20"] if "NDII_SWIR20" in piv.columns else np.nan
    piv["TOTAL_COST"] = piv[["NDVI_COST", "INDICES20_COST"]].sum(axis=1, min_count=1)
    return piv[["PATCH_ID", "NDVI_COST", "INDICES20_COST", "TOTAL_COST"]]


def latest_overflow_count(patch_id: str):
    if not PROCESS_LOG.exists():
        return np.nan
    try:
        d = pd.read_csv(PROCESS_LOG)
        d = d[d["PATCH_ID"].astype(str) == str(patch_id)]
        vals = pd.to_numeric(d.get("overflow_pixels_gt6"), errors="coerce").dropna()
        return float(vals.iloc[-1]) if len(vals) else np.nan
    except Exception:
        return np.nan


def row_value(row, field, default=""):
    if field not in row.index or pd.isna(row[field]):
        return default
    return row[field]


def add_region_metadata(qa: dict, row):
    qa["PET"] = row_value(row, "PET", np.nan)
    qa["BIOGEO_R"] = row_value(row, "BIOGEO_R", "")
    qa["AREA_PCT"] = row_value(row, "AREA_PCT", np.nan)
    return qa


def validate_final_outputs(final_tif: Path, diagnostics_tif: Path, events_csv: Path, grass_file: Path):
    _, _, transform, crs, bounds, width, height = load_grass(grass_file)
    expected_mask_hash = grass_mask_signature(grass_file)
    with rasterio.open(final_tif) as src:
        if src.count != 17 or src.width != width or src.height != height:
            raise RuntimeError("Final TIFF dimensions/band count invalid.")
        if src.crs != crs or src.transform != transform or src.bounds != bounds:
            raise RuntimeError("Final TIFF grid invalid.")
        if src.tags().get("GRASS_MASK_SHA256", "") != expected_mask_hash:
            raise RuntimeError("Final TIFF was created with a different detection grass mask.")
    with rasterio.open(diagnostics_tif) as src:
        if src.count != 44 or src.width != width or src.height != height:
            raise RuntimeError("Diagnostics TIFF dimensions/band count invalid.")
        if src.crs != crs or src.transform != transform or src.bounds != bounds:
            raise RuntimeError("Diagnostics TIFF grid invalid.")
        if src.tags().get("GRASS_MASK_SHA256", "") != expected_mask_hash:
            raise RuntimeError("Diagnostics TIFF was created with a different detection grass mask.")
    e = pd.read_csv(events_csv)
    required = {"PATCH_ID", "row", "col", "event", "event_confidence", "miletich_flag"}
    missing = required - set(e.columns)
    if missing:
        raise RuntimeError(f"Final events CSV missing columns: {sorted(missing)}")
    return e


def update_national_qa(row: dict):
    new = pd.DataFrame([row])
    if NATIONAL_QA.exists():
        try:
            old = pd.read_csv(NATIONAL_QA)
            old = old[old["PATCH_ID"].astype(str) != str(row["PATCH_ID"])]
            new = pd.concat([old, new], ignore_index=True, sort=False)
        except Exception:
            pass

    preferred = [
        "PATCH_ID", "PET", "BIOGEO_R", "AREA_PCT", "TYPE",
        "HAS_GRASS", "HAS_STATS_GRASS", "DETECTION_GRASS_PIXELS", "STATS_GRASS_PIXELS",
    ]
    front = [c for c in preferred if c in new.columns]
    rest = [c for c in new.columns if c not in front]
    new = new[front + rest].sort_values("PATCH_ID")
    atomic_csv(new, NATIONAL_QA)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "NATIONAL: create final/diagnostic mowing products using the +5 m detection mask, "
            "then calculate QA/statistics using the independent -7.5 m inner mask."
        )
    )
    parser.add_argument("--force", action="store_true", help="Re-create final outputs.")
    parser.add_argument("--patch", action="append", help="Process only this PATCH_ID (repeatable).")
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()

    if not PATCH_CSV.exists():
        raise FileNotFoundError(f"Missing {PATCH_CSV}. Run 09_prepare_processing.py first.")

    df = read_patch_table(PATCH_CSV)
    if args.patch:
        df = df[df[PATCH_ID_FIELD].isin(args.patch)]
    if args.limit:
        df = df.head(args.limit)

    grass_rows = df[df["HAS_GRASS"]].copy()
    if len(grass_rows):
        required_cols = {"STATS_CONFIGURED", "STATS_MASK_READY", "HAS_STATS_GRASS"}
        missing = required_cols - set(df.columns)
        if missing or not grass_rows["STATS_CONFIGURED"].all() or not grass_rows["STATS_MASK_READY"].all():
            raise RuntimeError(
                "Strict -7.5 m statistics masks are not ready in national_patches CSV.\n"
                f"Expected folder: {STATS_GRASS_ROOT}\n"
                f"Expected manifest: {STATS_GRASS_MANIFEST}\n"
                "Create the statistics masks from ORIGINAL polygons, then re-run:\n"
                "  python 09_prepare_processing.py --require-stats\n"
                "01 acquisition and 02 detection do not need to be repeated."
            )

    cost_table = costs_by_patch()
    costs = cost_table.set_index("PATCH_ID") if len(cost_table) else pd.DataFrame()
    print("National patches selected:", len(df))
    print("Final rasters: +5 m detection mask")
    print("QA/statistics: -7.5 m strict inner mask")

    n_ok = 0
    n_no_grass = 0
    n_err = 0

    for _, row in df.iterrows():
        patch_id = str(row[PATCH_ID_FIELD])
        patch_type = row_value(row, TYPE_FIELD, "")
        print("\n" + "=" * 78)
        print("FINALIZE:", patch_id)

        if not bool(row["HAS_GRASS"]):
            print("NO_GRASS -> no final rasters/events required")
            qa = {
                "PATCH_ID": patch_id,
                "TYPE": patch_type,
                "HAS_GRASS": False,
                "HAS_STATS_GRASS": False,
                "GRASS_PIXELS": 0,
                "DETECTION_GRASS_PIXELS": 0,
                "STATS_GRASS_PIXELS": 0,
                "STATS_MASK_BUFFER_M": -7.5,
                "VALID_OBS_MEDIAN": np.nan,
                "MAX_GAP_MEDIAN": np.nan,
                "MOWING_0_PCT": np.nan,
                "MOWING_1_PCT": np.nan,
                "MOWING_2_PCT": np.nan,
                "MOWING_3_PCT": np.nan,
                "MOWING_4_PCT": np.nan,
                "MOWING_5PLUS_PCT": np.nan,
                "EVENTS_20M_DETECTION": 0,
                "EVENTS_20M": 0,
                "NDII_DECREASE_PCT": np.nan,
                "SWIR_INCREASE_PCT": np.nan,
                "BOTH_PCT": np.nan,
                "STRONG_EVENTS_20M": 0,
                "STRONG_BOTH_PCT": np.nan,
                "FINAL_EVENTS_DETECTION": 0,
                "FINAL_EVENTS": 0,
                "EVENT_CONF_MEAN": np.nan,
                "EVENT_CONF_MEDIAN": np.nan,
                "MILETICH_FLAG_PCT": np.nan,
                "OVERFLOW_PIXELS_GT6": 0,
                "NDVI_COST": 0,
                "INDICES20_COST": 0,
                "TOTAL_COST": 0,
                "STATUS": "NO_GRASS",
            }
            qa = add_region_metadata(qa, row)
            update_national_qa(qa)
            append_log(FINALIZE_LOG, {
                "timestamp": datetime.now().isoformat(timespec="seconds"),
                "PATCH_ID": patch_id,
                "status": "NO_GRASS",
                "final_events": 0,
            })
            n_no_grass += 1
            continue

        p = patch_paths(patch_id)
        append_log(FINALIZE_LOG, {
            "timestamp": datetime.now().isoformat(timespec="seconds"),
            "PATCH_ID": patch_id,
            "status": "STARTED",
        })

        try:
            for key in ("grass", "baseline", "indices20", "dedup20"):
                if not p[key].exists():
                    raise FileNotFoundError(p[key])
            if not baseline_is_consistent(p["baseline"], p["grass"]):
                raise RuntimeError(
                    "Baseline is missing the current +5 m mask signature or is otherwise inconsistent. "
                    "Run 10_process_mowing.py first."
                )
            has_stats_grass = bool(row["HAS_STATS_GRASS"])
            if has_stats_grass and not p["stats_grass"].exists():
                raise FileNotFoundError(p["stats_grass"])

            valid_existing = False
            if p["final"].exists() and p["diagnostics"].exists() and p["events"].exists() and not args.force:
                try:
                    events = validate_final_outputs(
                        p["final"], p["diagnostics"], p["events"], p["grass"]
                    )
                    valid_existing = True
                    print("Final outputs exist + validate -> skipped creation")
                except Exception as exc:
                    print("Existing final outputs invalid -> recompute:", exc)

            if not valid_existing:
                temp = {}
                for key in ("final", "diagnostics", "events"):
                    path = p[key]
                    temp[key] = path.with_name(path.stem + ".part" + path.suffix)
                    if temp[key].exists():
                        temp[key].unlink()

                events = finalize_patch(
                    patch_id,
                    p["grass"],
                    p["baseline"],
                    p["indices20"],
                    temp["final"],
                    temp["diagnostics"],
                    temp["events"],
                )
                validate_final_outputs(
                    temp["final"], temp["diagnostics"], temp["events"], p["grass"]
                )
                temp["final"].replace(p["final"])
                temp["diagnostics"].replace(p["diagnostics"])
                temp["events"].replace(p["events"])
                print("Final detection events:", len(events))
                print("Saved:", p["final"])
                print("Saved:", p["diagnostics"])
                print("Saved:", p["events"])

            qa = patch_qa(
                patch_id,
                patch_type,
                p["grass"],
                p["stats_grass"] if has_stats_grass else None,
                p["baseline"],
                p["dedup20"],
                p["events"],
                has_stats_grass=has_stats_grass,
            )
            qa = add_region_metadata(qa, row)
            qa["HAS_GRASS"] = True
            qa["HAS_STATS_GRASS"] = has_stats_grass
            qa["OVERFLOW_PIXELS_GT6"] = latest_overflow_count(patch_id)
            qa["STATUS"] = "OK" if has_stats_grass else "OK_NO_STRICT_STATS_PIXELS"

            qa.update(NDVI_COST=np.nan, INDICES20_COST=np.nan, TOTAL_COST=np.nan)
            if isinstance(costs, pd.DataFrame) and patch_id in costs.index:
                c = costs.loc[patch_id]
                qa["NDVI_COST"] = c.get("NDVI_COST", np.nan)
                qa["INDICES20_COST"] = c.get("INDICES20_COST", np.nan)
                qa["TOTAL_COST"] = c.get("TOTAL_COST", np.nan)

            update_national_qa(qa)
            append_log(FINALIZE_LOG, {
                "timestamp": datetime.now().isoformat(timespec="seconds"),
                "PATCH_ID": patch_id,
                "status": qa["STATUS"],
                "final_events_detection": qa.get("FINAL_EVENTS_DETECTION", np.nan),
                "final_events_stats": qa.get("FINAL_EVENTS", np.nan),
            })
            n_ok += 1

        except Exception as exc:
            print("ERROR:", exc)
            traceback.print_exc()
            qa = {
                "PATCH_ID": patch_id,
                "TYPE": patch_type,
                "STATUS": f"ERROR: {exc}",
            }
            qa = add_region_metadata(qa, row)
            qa["OVERFLOW_PIXELS_GT6"] = latest_overflow_count(patch_id)
            update_national_qa(qa)
            append_log(FINALIZE_LOG, {
                "timestamp": datetime.now().isoformat(timespec="seconds"),
                "PATCH_ID": patch_id,
                "status": f"ERROR: {exc}",
            })
            n_err += 1

    print("\n" + "=" * 78)
    print("National finalize finished.")
    print("OK:", n_ok, "NO_GRASS:", n_no_grass, "ERROR:", n_err)
    print("National QA:", NATIONAL_QA)
    print("Finalize log:", FINALIZE_LOG)
    if n_err:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
