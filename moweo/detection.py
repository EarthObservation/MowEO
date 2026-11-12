from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
import traceback

import numpy as np
import pandas as pd

from .config import (
    OVERFLOW_ALL,
    PATCH_CSV,
    PATCH_ID_FIELD,
    PROCESS_LOG,
    PROCESS_ROOT,
    YEAR,
    patch_paths,
)
from .core import (
    baseline_is_consistent,
    build_baseline,
    build_dedup20_from_baseline,
    load_grass,
    read_patch_table,
)

OVERFLOW_COLUMNS = [
    "PATCH_ID", "row", "col", "RAW_EVENT_COUNT", "RAW_DOYS", "RAW_NDVI_DROPS"
]


def atomic_csv(df: pd.DataFrame, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.stem + ".part" + path.suffix)
    df.to_csv(tmp, index=False)
    tmp.replace(path)


def append_process_log(row: dict):
    new = pd.DataFrame([row])
    if PROCESS_LOG.exists():
        try:
            old = pd.read_csv(PROCESS_LOG)
            new = pd.concat([old, new], ignore_index=True, sort=False)
        except Exception:
            pass
    atomic_csv(new, PROCESS_LOG)


def previous_overflow_count(patch_id: str):
    if not PROCESS_LOG.exists():
        return np.nan
    try:
        d = pd.read_csv(PROCESS_LOG)
        if "PATCH_ID" not in d or "overflow_pixels_gt6" not in d:
            return np.nan
        d = d[d["PATCH_ID"].astype(str) == str(patch_id)]
        vals = pd.to_numeric(d["overflow_pixels_gt6"], errors="coerce").dropna()
        return float(vals.iloc[-1]) if len(vals) else np.nan
    except Exception:
        return np.nan


def validate_dedup20(path: Path):
    required = {
        "row20", "col20", "spectral_event_date", "ndii_drop", "swir_rise",
        "neighbour_support_median", "provisional_confidence_median",
    }
    d = pd.read_csv(path)
    missing = required - set(d.columns)
    if missing:
        raise RuntimeError(f"dedup20 missing columns: {sorted(missing)}")
    return d


def write_overflow_detail(path: Path, records: list[dict]):
    df = pd.DataFrame(records, columns=OVERFLOW_COLUMNS)
    atomic_csv(df, path)


def rebuild_national_overflow_detail():
    frames = []
    for path in PROCESS_ROOT.glob(f"*/**/*_overflow_gt6_{YEAR}.csv"):
        try:
            d = pd.read_csv(path)
            if len(d):
                frames.append(d)
        except Exception:
            continue
    out = pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame(columns=OVERFLOW_COLUMNS)
    atomic_csv(out, OVERFLOW_ALL)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "NATIONAL: local NDVI mowing detection, baseline and independent 20 m event table. "
            "Checkpointed per EOPatch."
        )
    )
    parser.add_argument("--force", action="store_true", help="Recompute baseline and dedup20.")
    parser.add_argument("--patch", action="append", help="Process only this PATCH_ID (repeatable).")
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument(
        "--backfill-overflow",
        action="store_true",
        help=(
            "If an imported baseline exists but no detailed >6-event CSV exists, "
            "recompute the baseline locally to recover row/col/raw event details."
        ),
    )
    args = parser.parse_args()

    if not PATCH_CSV.exists():
        raise FileNotFoundError(f"Missing {PATCH_CSV}. Run 09_prepare_processing.py first.")

    df = read_patch_table(PATCH_CSV)
    if args.patch:
        df = df[df[PATCH_ID_FIELD].isin(args.patch)]
    if args.limit:
        df = df.head(args.limit)

    print("National patches selected:", len(df))
    print("Local processing is intentionally sequential; openEO acquisition is the two-job step.")

    n_ok = 0
    n_no_grass = 0
    n_err = 0

    for _, row in df.iterrows():
        patch_id = str(row[PATCH_ID_FIELD])
        print("\n" + "=" * 78)
        print("PROCESS:", patch_id)

        if not bool(row["HAS_GRASS"]):
            print("NO_GRASS -> skipped local processing")
            append_process_log({
                "timestamp": datetime.now().isoformat(timespec="seconds"),
                "PATCH_ID": patch_id,
                "status": "NO_GRASS",
                "baseline_status": "not_applicable",
                "grass_pixels": 0,
                "overflow_pixels_gt6": 0,
                "events20": 0,
            })
            n_no_grass += 1
            continue

        p = patch_paths(patch_id)
        p["base"].mkdir(parents=True, exist_ok=True)

        append_process_log({
            "timestamp": datetime.now().isoformat(timespec="seconds"),
            "PATCH_ID": patch_id,
            "status": "STARTED",
        })

        try:
            for key in ("grass", "ndvi", "indices20"):
                if not p[key].exists():
                    raise FileNotFoundError(p[key])

            grass, _, _, _, _, _, _ = load_grass(p["grass"])
            grass_pixels = int((grass == 1).sum())

            baseline_valid = p["baseline"].exists() and baseline_is_consistent(p["baseline"], p["grass"])
            need_overflow_backfill = args.backfill_overflow and not p["overflow"].exists()
            recompute_baseline = args.force or not baseline_valid or need_overflow_backfill

            if recompute_baseline:
                tmp_baseline = p["baseline"].with_name(
                    p["baseline"].stem + ".part" + p["baseline"].suffix
                )
                if tmp_baseline.exists():
                    tmp_baseline.unlink()
                result = build_baseline(
                    patch_id, p["grass"], p["ndvi"], tmp_baseline, force=True
                )
                if not baseline_is_consistent(tmp_baseline, p["grass"]):
                    raise RuntimeError("New baseline failed consistency validation.")
                tmp_baseline.replace(p["baseline"])
                write_overflow_detail(p["overflow"], result.get("overflow_records", []))
                overflow_count = int(result.get("overflow_pixels", 0))
                baseline_status = "recomputed"
                print("Baseline recomputed:", p["baseline"])
                print("Pixels with >6 raw events:", overflow_count)
            else:
                baseline_status = "existing"
                if p["overflow"].exists():
                    try:
                        overflow_count = len(pd.read_csv(p["overflow"]))
                    except Exception:
                        overflow_count = previous_overflow_count(patch_id)
                else:
                    overflow_count = previous_overflow_count(patch_id)
                print("Baseline exists + consistent -> skipped")

            recompute_dedup = args.force or recompute_baseline or not p["dedup20"].exists()
            if not recompute_dedup:
                try:
                    d = validate_dedup20(p["dedup20"])
                    print("dedup20 exists + valid -> skipped | events:", len(d))
                except Exception as exc:
                    print("dedup20 invalid -> recompute:", exc)
                    recompute_dedup = True

            if recompute_dedup:
                tmp_csv = p["dedup20"].with_name(
                    p["dedup20"].stem + ".part" + p["dedup20"].suffix
                )
                if tmp_csv.exists():
                    tmp_csv.unlink()
                d = build_dedup20_from_baseline(
                    patch_id, p["grass"], p["baseline"], p["indices20"], tmp_csv
                )
                validate_dedup20(tmp_csv)
                tmp_csv.replace(p["dedup20"])
                print("Independent 20 m events:", len(d))
                print("Saved:", p["dedup20"])

            append_process_log({
                "timestamp": datetime.now().isoformat(timespec="seconds"),
                "PATCH_ID": patch_id,
                "status": "OK",
                "baseline_status": baseline_status,
                "grass_pixels": grass_pixels,
                "overflow_pixels_gt6": overflow_count,
                "events20": len(d),
            })
            n_ok += 1

        except Exception as exc:
            print("ERROR:", exc)
            traceback.print_exc()
            append_process_log({
                "timestamp": datetime.now().isoformat(timespec="seconds"),
                "PATCH_ID": patch_id,
                "status": f"ERROR: {exc}",
            })
            n_err += 1

    rebuild_national_overflow_detail()

    print("\n" + "=" * 78)
    print("National local processing finished.")
    print("OK:", n_ok, "NO_GRASS:", n_no_grass, "ERROR:", n_err)
    print("Processing log:", PROCESS_LOG)
    print("Detailed >6-event table:", OVERFLOW_ALL)
    if n_err:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
