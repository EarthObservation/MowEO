from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pandas as pd

from . import config as cfg

HERE = cfg.PROJECT_DIR
from .commands import command


def atomic_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.stem + ".part" + path.suffix)
    df.to_csv(tmp, index=False)
    tmp.replace(path)


def split_round_robin(items: list[str], n: int) -> list[list[str]]:
    groups = [[] for _ in range(n)]
    for i, item in enumerate(items):
        groups[i % n].append(item)
    return [g for g in groups if g]


def read_patch_ids() -> list[str]:
    if not cfg.PATCH_CSV.exists():
        raise FileNotFoundError(f"Missing patch table: {cfg.PATCH_CSV}")
    df = pd.read_csv(cfg.PATCH_CSV)
    if cfg.PATCH_ID_FIELD not in df.columns:
        raise RuntimeError(f"{cfg.PATCH_CSV} has no {cfg.PATCH_ID_FIELD} column")
    return df[cfg.PATCH_ID_FIELD].astype(str).tolist()


def merge_history(canonical: Path, worker_files: list[Path]) -> None:
    frames = []
    if canonical.exists():
        try:
            frames.append(pd.read_csv(canonical))
        except Exception:
            pass
    for path in worker_files:
        if path.exists():
            try:
                frames.append(pd.read_csv(path))
            except Exception as exc:
                print(f"WARNING: could not merge {path}: {exc}")
    if frames:
        atomic_csv(pd.concat(frames, ignore_index=True, sort=False), canonical)


def merge_qa(canonical: Path, worker_files: list[Path]) -> None:
    frames = []
    if canonical.exists():
        try:
            frames.append(pd.read_csv(canonical))
        except Exception:
            pass
    for path in worker_files:
        if path.exists():
            try:
                frames.append(pd.read_csv(path))
            except Exception as exc:
                print(f"WARNING: could not merge {path}: {exc}")
    if not frames:
        return
    out = pd.concat(frames, ignore_index=True, sort=False)
    if "PATCH_ID" in out.columns:
        out["PATCH_ID"] = out["PATCH_ID"].astype(str)
        out = out.drop_duplicates(subset=["PATCH_ID"], keep="last").sort_values("PATCH_ID")
    atomic_csv(out, canonical)


def rebuild_overflow(canonical: Path) -> None:
    frames = []
    for path in cfg.PROCESS_ROOT.glob(f"*/**/*_overflow_gt6_{cfg.YEAR}.csv"):
        try:
            d = pd.read_csv(path)
            if len(d):
                frames.append(d)
        except Exception:
            continue
    if frames:
        out = pd.concat(frames, ignore_index=True, sort=False)
    else:
        out = pd.DataFrame(columns=["PATCH_ID", "row", "col", "RAW_EVENT_COUNT", "RAW_DOYS", "RAW_NDVI_DROPS"])
    atomic_csv(out, canonical)


def run_parallel(step: str, workers: int) -> None:
    patch_ids = read_patch_ids()
    workers = max(1, min(workers, len(patch_ids)))

    if workers == 1:
        script = "10_process_mowing.py" if step == "process" else "11_finalize.py"
        subprocess.run(command([script]), cwd=HERE, check=True)
        return

    script = "10_process_mowing.py" if step == "process" else "11_finalize.py"
    groups = split_round_robin(patch_ids, workers)

    temp_root = cfg.ROOT / "_parallel_logs" / step
    if temp_root.exists():
        shutil.rmtree(temp_root)
    temp_root.mkdir(parents=True, exist_ok=True)

    print("=" * 78)
    print(f"PARALLEL {step.upper()} | year={cfg.YEAR} | workers={len(groups)} | patches={len(patch_ids)}")
    print("Each worker processes its own EOPatch subset; raw NetCDF inputs are read-only.")
    print("=" * 78)

    procs: list[tuple[int, subprocess.Popen, Path, Path | None]] = []
    for i, group in enumerate(groups, start=1):
        env = os.environ.copy()
        if step == "process":
            worker_log = temp_root / f"process_worker_{i}.csv"
            worker_aux = temp_root / f"overflow_worker_{i}.csv"
            env["MOWING_PROCESS_LOG"] = str(worker_log)
            env["MOWING_OVERFLOW_ALL"] = str(worker_aux)
        else:
            worker_log = temp_root / f"finalize_worker_{i}.csv"
            worker_aux = temp_root / f"qa_worker_{i}.csv"
            env["MOWING_FINALIZE_LOG"] = str(worker_log)
            env["MOWING_NATIONAL_QA"] = str(worker_aux)

        cmd = command([script])
        for patch_id in group:
            cmd.extend(["--patch", patch_id])

        print(f"Starting worker {i}/{len(groups)}: {len(group)} patches")
        proc = subprocess.Popen(cmd, cwd=HERE, env=env)
        procs.append((i, proc, worker_log, worker_aux))

    failures = []
    for i, proc, _, _ in procs:
        rc = proc.wait()
        if rc != 0:
            failures.append((i, rc))
        print(f"Worker {i} finished with exit code {rc}")

    if step == "process":
        merge_history(cfg.ROOT / f"national_processing_{cfg.YEAR}.csv", [x[2] for x in procs])
        rebuild_overflow(cfg.ROOT / f"national_overflow_gt6_detail_{cfg.YEAR}.csv")
    else:
        merge_history(cfg.ROOT / f"national_finalize_{cfg.YEAR}.csv", [x[2] for x in procs])
        merge_qa(cfg.ROOT / f"national_qa_{cfg.YEAR}.csv", [x[3] for x in procs if x[3] is not None])

    if failures:
        raise SystemExit(f"Parallel {step} failed in worker(s): {failures}")

    print(f"Parallel {step} complete.")


def main() -> None:
    p = argparse.ArgumentParser(description="Run patch-local processing/finalization in parallel without shared-log write conflicts.")
    p.add_argument("step", choices=["process", "finalize"])
    p.add_argument("--workers", type=int, default=4)
    args = p.parse_args()
    if args.workers < 1:
        raise SystemExit("--workers must be >= 1")
    run_parallel(args.step, args.workers)


if __name__ == "__main__":
    main()
