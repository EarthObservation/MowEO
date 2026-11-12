from __future__ import annotations

import argparse
from collections import deque
from datetime import datetime
from pathlib import Path
import time
import traceback

import numpy as np
import pandas as pd
import xarray as xr
from tqdm.auto import tqdm

from .auth import connect_project_account

from .config import (
    INDICES20_HEIGHT,
    INDICES20_WIDTH,
    LOAD_END,
    LOAD_START,
    NDVI_HEIGHT,
    NDVI_WIDTH,
    OPENEO_LOG,
    PATCH_CSV_ACQUISITION,
    TARGET_EPSG,
    OPENEO_MAX_PARALLEL_JOBS,
    YEAR,
    acquisition_patch_paths,
)

DEFAULT_PARALLEL_JOBS = min(2, OPENEO_MAX_PARALLEL_JOBS)
DEFAULT_POLL_SECONDS = 30
START_GAP_SECONDS = 6


def validate_ndvi(path: Path, west: float, north: float):
    with xr.open_dataset(path, engine="h5netcdf") as ds:
        if "NDVI" not in ds:
            raise RuntimeError("NDVI variable missing")
        if ds.sizes["x"] != NDVI_WIDTH or ds.sizes["y"] != NDVI_HEIGHT:
            raise RuntimeError(
                f"NDVI expected {NDVI_HEIGHT}x{NDVI_WIDTH}, "
                f"got {ds.sizes['y']}x{ds.sizes['x']}"
            )

        expected_x = west + 5.0 + np.arange(NDVI_WIDTH) * 10.0
        expected_y = north - 5.0 - np.arange(NDVI_HEIGHT) * 10.0
        if not np.allclose(ds.x.values, expected_x, atol=0.01, rtol=0):
            raise RuntimeError("NDVI x-grid misaligned")
        if not np.allclose(ds.y.values, expected_y, atol=0.01, rtol=0):
            raise RuntimeError("NDVI y-grid misaligned")


def validate_indices20(path: Path, west: float, north: float):
    with xr.open_dataset(path, engine="h5netcdf") as ds:
        if "NDII" not in ds or "SWIRMean" not in ds:
            raise RuntimeError("NDII/SWIRMean variable missing")
        if ds.sizes["x"] != INDICES20_WIDTH or ds.sizes["y"] != INDICES20_HEIGHT:
            raise RuntimeError(
                f"20 m expected {INDICES20_HEIGHT}x{INDICES20_WIDTH}, "
                f"got {ds.sizes['y']}x{ds.sizes['x']}"
            )
        expected_x = west + 10.0 + np.arange(INDICES20_WIDTH) * 20.0
        expected_y = north - 10.0 - np.arange(INDICES20_HEIGHT) * 20.0
        if not np.allclose(ds.x.values, expected_x, atol=0.01, rtol=0):
            raise RuntimeError("20 m x-grid misaligned")
        if not np.allclose(ds.y.values, expected_y, atol=0.01, rtol=0):
            raise RuntimeError("20 m y-grid misaligned")


def validate_task_output(task: dict, path: Path):
    if task["product"] == "NDVI10":
        validate_ndvi(path, task["west"], task["north"])
    elif task["product"] == "NDII_SWIR20":
        validate_indices20(path, task["west"], task["north"])
    else:
        raise ValueError(f"Unknown product: {task['product']}")


def append_log(row: dict):
    new = pd.DataFrame([row])
    if OPENEO_LOG.exists():
        try:
            old = pd.read_csv(OPENEO_LOG)
            if not old.empty:
                new = pd.concat([old, new], ignore_index=True)
        except Exception:
            pass
    OPENEO_LOG.parent.mkdir(parents=True, exist_ok=True)
    new.to_csv(OPENEO_LOG, index=False)


def log_status(task: dict, status: str, job_id: str = "", costs=None, message: str = ""):
    append_log(
        {
            "timestamp": datetime.now().isoformat(timespec="seconds"),
            "PATCH_ID": task["PATCH_ID"],
            "product": task["product"],
            "job_id": job_id,
            "costs": costs,
            "status": status,
            "output": str(task.get("output", "")),
            "message": message,
        }
    )


def latest_job_records() -> dict[tuple[str, str], dict]:
    if not OPENEO_LOG.exists():
        return {}
    try:
        df = pd.read_csv(OPENEO_LOG)
    except Exception:
        return {}
    required = {"PATCH_ID", "product", "job_id", "status"}
    if not required.issubset(df.columns):
        return {}
    if "timestamp" in df.columns:
        df = df.sort_values("timestamp")
    df = df.drop_duplicates(subset=["PATCH_ID", "product"], keep="last")
    return {
        (str(row["PATCH_ID"]), str(row["product"])): row.to_dict()
        for _, row in df.iterrows()
    }



def connect_openeo(*, reuse_auth: bool = False):
    return connect_project_account(force_interactive=not reuse_auth, confirm=True)


def build_ndvi_cube(connection, spatial_extent: dict):
    scl = connection.load_collection(
        "SENTINEL2_L2A",
        spatial_extent=spatial_extent,
        temporal_extent=[LOAD_START, LOAD_END],
        bands=["SCL"],
        max_cloud_cover=95,
    )
    cloud_mask = scl.process(
        "to_scl_dilation_mask",
        data=scl,
        kernel1_size=17,
        kernel2_size=77,
        mask1_values=[2, 4, 5, 6, 7],
        mask2_values=[3, 8, 9, 10, 11],
        erosion_kernel_size=3,
    )
    s2 = connection.load_collection(
        "SENTINEL2_L2A",
        spatial_extent=spatial_extent,
        temporal_extent=[LOAD_START, LOAD_END],
        bands=["B04", "B08"],
        max_cloud_cover=95,
    )
    return (
        s2.mask(cloud_mask).resample_spatial(
            resolution=10, projection=TARGET_EPSG, method='near', align='upper-left'
        )
        .ndvi(red="B04", nir="B08", target_band="NDVI")
        .filter_bands(["NDVI"])
    )


def build_indices20_cube(connection, spatial_extent: dict):
    featureflags = {"no_resample_on_read": True}

    s2_20 = connection.load_collection(
        "SENTINEL2_L2A",
        spatial_extent=spatial_extent,
        temporal_extent=[LOAD_START, LOAD_END],
        bands=["B8A", "B11", "B12"],
        max_cloud_cover=95,
    )
    s2_20.result_node().update_arguments(featureflags=featureflags)
    s2_20 = s2_20.resample_spatial(
        resolution=20,
        projection=TARGET_EPSG,
        method="near",
        align="upper-left",
    )

    scl_20 = connection.load_collection(
        "SENTINEL2_L2A",
        spatial_extent=spatial_extent,
        temporal_extent=[LOAD_START, LOAD_END],
        bands=["SCL"],
        max_cloud_cover=95,
    )
    scl_20.result_node().update_arguments(featureflags=featureflags)
    scl_20 = scl_20.resample_spatial(
        resolution=20,
        projection=TARGET_EPSG,
        method="near",
        align="upper-left",
    )

    cloud_mask_20 = scl_20.process(
        "to_scl_dilation_mask",
        data=scl_20,
        kernel1_size=17,
        kernel2_size=77,
        mask1_values=[2, 4, 5, 6, 7],
        mask2_values=[3, 8, 9, 10, 11],
        erosion_kernel_size=3,
    )

    clean = s2_20.mask(cloud_mask_20)
    b8a = clean.band("B8A")
    b11 = clean.band("B11")
    b12 = clean.band("B12")

    ndii = ((b8a - b11) / (b8a + b11)).add_dimension(
        name="bands", label="NDII", type="bands"
    )
    swir = ((b11 + b12) / 2).add_dimension(
        name="bands", label="SWIRMean", type="bands"
    )
    return ndii.merge_cubes(swir)


def create_job_for_task(connection, task: dict):
    if task["product"] == "NDVI10":
        cube = build_ndvi_cube(connection, task["spatial_extent"])
        title = f"{task['PATCH_ID']} NDVI {YEAR}"
        description = (
            f"Sentinel-2 L2A NDVI 10 m, {task['PATCH_ID']}, "
            f"{LOAD_START} to {LOAD_END}, SCL dilation mask."
        )
    elif task["product"] == "NDII_SWIR20":
        cube = build_indices20_cube(connection, task["spatial_extent"])
        title = f"{task['PATCH_ID']} NDII + SWIRMean explicit20m {YEAR}"
        description = (
            f"Sentinel-2 L2A B8A/B11/B12 explicit 20 m cube, "
            f"{task['PATCH_ID']}, {LOAD_START} to {LOAD_END}."
        )
    else:
        raise ValueError(task["product"])

    # Persist the actual process graph before submitting any remote job.
    import json
    from .config import YEAR_ROOT
    graph_dir = YEAR_ROOT / 'provenance' / 'process_graphs'
    graph_dir.mkdir(parents=True, exist_ok=True)
    (graph_dir / f"{task['PATCH_ID']}_{task['product']}.json").write_text(
        json.dumps(cube.flat_graph(), indent=2) + '\n', encoding='utf-8')
    return cube.create_job(out_format="NetCDF", title=title, description=description)


def make_task(row: pd.Series, product: str):
    patch_id = str(row["PATCH_ID"])
    west = float(row["WEST"])
    south = float(row["SOUTH"])
    east = float(row["EAST"])
    north = float(row["NORTH"])
    paths = acquisition_patch_paths(patch_id)
    output = paths["ndvi"] if product == "NDVI10" else paths["indices20"]
    return {
        "PATCH_ID": patch_id,
        "product": product,
        "output": output,
        "west": west,
        "south": south,
        "east": east,
        "north": north,
        "spatial_extent": {
            "west": west,
            "south": south,
            "east": east,
            "north": north,
            "crs": f"EPSG:{TARGET_EPSG}",
        },
        "retry_count": 0,
    }


def download_and_validate(job, task: dict):
    output = Path(task["output"])
    output.parent.mkdir(parents=True, exist_ok=True)
    tmp = output.with_name(output.stem + ".part" + output.suffix)
    if tmp.exists():
        tmp.unlink()
    print(f"Downloading {task['PATCH_ID']} {task['product']} -> {output.name}")
    job.get_results().download_file(tmp)
    validate_task_output(task, tmp)
    tmp.replace(output)


def job_cost(job):
    try:
        return job.describe().get("costs")
    except Exception:
        return None


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Acquire Sentinel-2 NDVI 10 m and NDII/SWIRMean 20 m cubes "
            "for EOPatches listed in the prepared CSV."
        )
    )
    parser.add_argument("--run", action="store_true", help="Submit/resume openEO jobs. Default is dry run.")
    parser.add_argument("--force", action="store_true", help="Re-create cubes even if valid outputs exist.")
    parser.add_argument("--patch", action="append", help="Only this PATCH_ID; can be repeated.")
    parser.add_argument("--limit", type=int, default=None, help="Only first N selected patches.")
    parser.add_argument("--parallel", type=int, default=DEFAULT_PARALLEL_JOBS, help="Concurrent openEO jobs (1 or 2).")
    parser.add_argument("--poll-seconds", type=int, default=DEFAULT_POLL_SECONDS, help="Status polling interval, >=10 s.")
    parser.add_argument("--max-retries", type=int, default=1, help="Retries after backend/download error.")
    parser.add_argument("--no-resume", action="store_true", help="Do not resume jobs from the persistent log.")
    parser.add_argument("--reuse-auth", action="store_true", help="Reuse this year/project local refresh token instead of forcing a fresh login. Default: fresh interactive login.")
    args = parser.parse_args()

    if args.parallel < 1 or args.parallel > OPENEO_MAX_PARALLEL_JOBS:
        raise ValueError(f"--parallel must be between 1 and {OPENEO_MAX_PARALLEL_JOBS}.")
    if args.poll_seconds < 10:
        raise ValueError("Use --poll-seconds >= 10.")
    if not PATCH_CSV_ACQUISITION.exists():
        raise FileNotFoundError(
            f"Missing prepared table: {PATCH_CSV_ACQUISITION}\nRun 01_prepare_acquisition.py first."
        )

    df = pd.read_csv(PATCH_CSV_ACQUISITION)
    required = {"PATCH_ID", "WEST", "SOUTH", "EAST", "NORTH"}
    missing = required - set(df.columns)
    if missing:
        raise RuntimeError(f"Prepared table missing columns: {sorted(missing)}")

    if args.patch:
        requested = set(args.patch)
        unknown = sorted(requested - set(df["PATCH_ID"].astype(str)))
        if unknown:
            raise RuntimeError(f"Unknown PATCH_ID: {unknown}")
        df = df[df["PATCH_ID"].astype(str).isin(requested)]
    if args.limit is not None:
        if args.limit < 1:
            raise ValueError("--limit must be >= 1")
        df = df.head(args.limit)

    print("=" * 78)
    print(f"{YEAR} SENTINEL-2 EOPATCH ACQUISITION")
    print("=" * 78)
    print("Selected EOPatches:", len(df))
    print("Temporal extent:", LOAD_START, "to", LOAD_END)
    print("Mode:", f"RUN, parallel={args.parallel}" if args.run else "DRY RUN")
    print("No grass masks are used. Acquisition covers full 5 x 5 km EOPatches.")

    tasks = []
    complete_patch_ids = set()
    skipped = 0
    invalid_existing = 0

    for _, row in df.iterrows():
        patch_id = str(row["PATCH_ID"])
        paths = acquisition_patch_paths(patch_id)
        paths["base"].mkdir(parents=True, exist_ok=True)

        print("\n" + "-" * 78)
        print(patch_id)
        print(
            "Extent:",
            float(row["WEST"]), float(row["SOUTH"]),
            float(row["EAST"]), float(row["NORTH"]),
        )

        ndvi_needed = args.force or not paths["ndvi"].exists()
        if not ndvi_needed:
            try:
                validate_ndvi(paths["ndvi"], float(row["WEST"]), float(row["NORTH"]))
                print("NDVI exists + grid OK -> skipped")
                skipped += 1
            except Exception as exc:
                print("NDVI exists but invalid -> re-create:", exc)
                ndvi_needed = True
                invalid_existing += 1
        if ndvi_needed:
            tasks.append(make_task(row, "NDVI10"))
            print("NDVI10 ->", "QUEUED" if args.run else "WOULD SUBMIT")

        idx_needed = args.force or not paths["indices20"].exists()
        if not idx_needed:
            try:
                validate_indices20(paths["indices20"], float(row["WEST"]), float(row["NORTH"]))
                print("NDII/SWIR 20 m exists + grid OK -> skipped")
                skipped += 1
            except Exception as exc:
                print("20 m cube exists but invalid -> re-create:", exc)
                idx_needed = True
                invalid_existing += 1
        if idx_needed:
            tasks.append(make_task(row, "NDII_SWIR20"))
            print("NDII_SWIR20 ->", "QUEUED" if args.run else "WOULD SUBMIT")

        if not ndvi_needed and not idx_needed:
            complete_patch_ids.add(patch_id)

    print("\nProduct tasks to run:", len(tasks))
    print("Existing valid products skipped:", skipped)
    if invalid_existing:
        print("Invalid existing products to replace:", invalid_existing)

    if not args.run:
        print("\nDry run finished. No openEO jobs submitted.")
        return
    if not tasks:
        print("\nNothing to run. All selected outputs already exist and validate.")
        return

    connection = connect_openeo(reuse_auth=args.reuse_auth)
    latest = latest_job_records() if not args.no_resume and not args.force else {}

    pending = deque()
    active = {}
    completed = 0
    failed = 0
    progress = tqdm(
        total=len(df),
        initial=len(complete_patch_ids),
        desc=f"{YEAR} complete EOPatches",
        unit="patch",
        dynamic_ncols=True,
    )

    def mark_patch_complete(patch_id: str):
        if patch_id in complete_patch_ids:
            return
        paths = acquisition_patch_paths(patch_id)
        if paths["ndvi"].exists() and paths["indices20"].exists():
            complete_patch_ids.add(patch_id)
            progress.update(1)

    def refresh_progress(active_count: int, pending_count: int, failed_count: int):
        progress.set_postfix(
            active=active_count, pending=pending_count, failed=failed_count,
            products_done=completed, refresh=False
        )

    for task in tasks:
        key = (task["PATCH_ID"], task["product"])
        rec = latest.get(key)
        if rec and str(rec.get("job_id", "")).strip() not in {"", "nan", "None"}:
            job_id = str(rec["job_id"])
            try:
                job = connection.job(job_id)
                status = job.status()
                if status == "finished":
                    print(f"RESUME finished {job_id}: {task['PATCH_ID']} {task['product']}")
                    download_and_validate(job, task)
                    costs = job_cost(job)
                    log_status(task, "finished", job_id=job_id, costs=costs)
                    completed += 1
                    mark_patch_complete(task["PATCH_ID"])
                    refresh_progress(len(active), len(pending), failed)
                    continue
                if status in {"created", "queued", "running"}:
                    print(f"RESUME {status} {job_id}: {task['PATCH_ID']} {task['product']}")
                    active[job_id] = {"job": job, "task": task, "last_status": status}
                    continue
                print(f"Previous job {job_id} status={status}; will submit a new job.")
            except Exception as exc:
                print(f"Could not resume job {job_id}: {exc}")
        pending.append(task)

    print(f"\nScheduler start: pending={len(pending)}, resumed_active={len(active)}, parallel={args.parallel}")
    refresh_progress(len(active), len(pending), failed)

    try:
        while pending or active:
            while pending and len(active) < args.parallel:
                task = pending.popleft()
                try:
                    job = create_job_for_task(connection, task)
                    job.start()
                    active[job.job_id] = {"job": job, "task": task, "last_status": "created"}
                    log_status(task, "created", job_id=job.job_id)
                    print(
                        f"STARTED {job.job_id}: {task['PATCH_ID']} {task['product']} "
                        f"| active {len(active)}/{args.parallel}"
                    )
                    refresh_progress(len(active), len(pending), failed)
                    if pending and len(active) < args.parallel:
                        time.sleep(START_GAP_SECONDS)
                except Exception as exc:
                    task["retry_count"] += 1
                    log_status(task, "start_failed", message=f"start failed: {exc}")
                    print(f"START ERROR {task['PATCH_ID']} {task['product']}: {exc}")
                    if task["retry_count"] <= args.max_retries:
                        pending.append(task)
                    else:
                        failed += 1
                    refresh_progress(len(active), len(pending), failed)
                    time.sleep(args.poll_seconds)

            if not active:
                continue

            print(
                f"\nPolling {len(active)} active job(s) | pending={len(pending)} "
                f"| completed={completed} | failed={failed}"
            )
            finished_ids = []

            for job_id, info in list(active.items()):
                job = info["job"]
                task = info["task"]
                try:
                    status = job.status()
                except Exception as exc:
                    print(f"STATUS ERROR {job_id}: {exc}")
                    continue

                print(f"  {status:8s} {job_id} {task['PATCH_ID']} {task['product']}")
                if status != info["last_status"]:
                    log_status(task, status, job_id=job_id)
                    info["last_status"] = status

                if status == "finished":
                    try:
                        download_and_validate(job, task)
                        costs = job_cost(job)
                        log_status(task, "finished", job_id=job_id, costs=costs)
                        print(
                            f"FINISHED {task['PATCH_ID']} {task['product']} "
                            f"job={job_id} costs={costs}"
                        )
                        completed += 1
                        mark_patch_complete(task["PATCH_ID"])
                        refresh_progress(len(active), len(pending), failed)
                    except Exception as exc:
                        task["retry_count"] += 1
                        log_status(
                            task,
                            "download_error",
                            job_id=job_id,
                            message=f"download/validation failed: {exc}",
                        )
                        print(f"DOWNLOAD/VALIDATION ERROR {task['PATCH_ID']} {task['product']}: {exc}")
                        traceback.print_exc()
                        if task["retry_count"] <= args.max_retries:
                            pending.append(task)
                        else:
                            failed += 1
                    finished_ids.append(job_id)

                elif status in {"error", "canceled"}:
                    costs = job_cost(job)
                    log_status(task, status, job_id=job_id, costs=costs)
                    task["retry_count"] += 1
                    if task["retry_count"] <= args.max_retries:
                        print(
                            f"JOB {status.upper()} -> retry "
                            f"({task['retry_count']}/{args.max_retries}): "
                            f"{task['PATCH_ID']} {task['product']}"
                        )
                        pending.append(task)
                    else:
                        print(f"JOB {status.upper()} -> giving up: {task['PATCH_ID']} {task['product']}")
                        failed += 1
                    finished_ids.append(job_id)

            for job_id in finished_ids:
                active.pop(job_id, None)
            refresh_progress(len(active), len(pending), failed)

            if pending or active:
                time.sleep(args.poll_seconds)

    except KeyboardInterrupt:
        progress.close()
        print("\nInterrupted locally with Ctrl+C.")
        print("Queued/running openEO jobs were NOT canceled.")
        print(f"Run the same command again to resume from {OPENEO_LOG.name}.")
        raise SystemExit(130)

    progress.close()
    print("\n" + "=" * 78)
    print("ACQUISITION SCHEDULER FINISHED")
    print("=" * 78)
    print("Completed product tasks:", completed)
    print("Failed product tasks   :", failed)
    print("Job log:", OPENEO_LOG)
    if failed:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
