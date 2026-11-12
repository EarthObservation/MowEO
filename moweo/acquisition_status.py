from __future__ import annotations

import pandas as pd

from .config import OPENEO_LOG, PATCH_CSV_ACQUISITION, YEAR, acquisition_patch_paths


def main():
    print("=" * 78)
    print(f"ACQUISITION STATUS {YEAR}")
    print("=" * 78)

    if not PATCH_CSV_ACQUISITION.exists():
        print("Prepared patch table missing:", PATCH_CSV_ACQUISITION)
        return

    patches = pd.read_csv(PATCH_CSV_ACQUISITION)
    ids = patches["PATCH_ID"].astype(str).tolist()

    ndvi_ok = sum(acquisition_patch_paths(pid)["ndvi"].exists() for pid in ids)
    idx_ok = sum(acquisition_patch_paths(pid)["indices20"].exists() for pid in ids)
    both_ok = sum(
        acquisition_patch_paths(pid)["ndvi"].exists()
        and acquisition_patch_paths(pid)["indices20"].exists()
        for pid in ids
    )

    pct = 100.0 * both_ok / len(ids) if ids else 0.0
    print("Selected EOPatches:", len(ids))
    print("NDVI10 files       :", ndvi_ok)
    print("NDII/SWIR20 files  :", idx_ok)
    print(f"Complete EOPatches : {both_ok} / {len(ids)} ({pct:.1f}%)")
    print("Remaining complete :", len(ids) - both_ok)

    if not OPENEO_LOG.exists():
        print("\nNo job log yet:", OPENEO_LOG)
        return

    log = pd.read_csv(OPENEO_LOG)
    if log.empty:
        return

    if "timestamp" in log.columns:
        log = log.sort_values("timestamp")
    latest = log.drop_duplicates(["PATCH_ID", "product"], keep="last")
    print("\nLatest task statuses:")
    print(latest["status"].value_counts(dropna=False).to_string())

    if "costs" in log.columns:
        cost_rows = log[pd.to_numeric(log["costs"], errors="coerce").notna()].copy()
        if not cost_rows.empty:
            cost_rows["costs"] = pd.to_numeric(cost_rows["costs"], errors="coerce")
            charged = cost_rows.drop_duplicates("job_id", keep="last")
            print("\nRecorded openEO costs:", float(charged["costs"].sum()))

    if "timestamp" in log.columns:
        print("Last log timestamp:", log["timestamp"].iloc[-1])


if __name__ == "__main__":
    main()
