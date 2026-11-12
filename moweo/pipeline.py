from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys

from . import config as cfg

HERE = cfg.PROJECT_DIR
from .commands import command


def run(label: str, args: list[str], dry_run: bool = False) -> None:
    cmd = command(args)
    print("\n" + "=" * 78)
    print(label)
    print("COMMAND:", " ".join(cmd))
    print("=" * 78)
    if not dry_run:
        subprocess.run(cmd, cwd=HERE, check=True)


def main() -> None:
    p = argparse.ArgumentParser(description="Run the complete local mowing-processing pipeline for the configured year.")
    p.add_argument("--workers", type=int, default=4, help="Parallel EOPatch workers for detection/finalization.")
    p.add_argument("--skip-analysis", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()
    if args.workers < 1:
        raise SystemExit("--workers must be >= 1")

    print("Configured year:", cfg.YEAR)
    print("Raw NetCDF input:", cfg.RAW_PROCESS_ROOT)
    print("Processing output:", cfg.ROOT)
    print("Workers:", args.workers)

    steps = [
        ("Prepare/validate processing table and masks", ["09_prepare_processing.py", "--require-stats"]),
        ("Detect mowing events", ["parallel_step.py", "process", "--workers", str(args.workers)]),
        ("Finalize patch products and QA", ["parallel_step.py", "finalize", "--workers", str(args.workers)]),
    ]
    if not args.skip_analysis:
        steps.append(("National statistics and figures", ["12_analyze.py"]))
    steps.append(("Final national mowing-count mosaic", ["13_mosaic.py"]))
    steps.append(('GeoPackage of tile QA summaries', ['vector_export.py']))

    for label, cmd in steps:
        run(label, cmd, args.dry_run)

    print("\nProcessing pipeline complete.")
    print("Final mosaic:", cfg.PROCESS_ROOT.parent / "national_maps" / f"{cfg.OUTPUT_PREFIX}_MOWING_COUNT_{cfg.YEAR}_detection_plus5m.tif")


if __name__ == "__main__":
    main()
