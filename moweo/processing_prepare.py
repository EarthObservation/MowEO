from __future__ import annotations

import argparse
import re

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio

from .config import (
    EXPECTED_PATCHES,
    GRASS_MANIFEST,
    GRASS_ROOT,
    PATCH_CSV,
    PATCH_GPKG,
    PATCH_LAYER,
    PROCESS_ROOT,
    ROOT,
    STATS_GRASS_MANIFEST,
    STATS_GRASS_ROOT,
    YEAR,
    patch_paths,
)
from .core import load_grass
from importlib import import_module


def _clean_int_series(series: pd.Series, field_name: str) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    if numeric.isna().any():
        bad = series[numeric.isna()].tolist()
        raise ValueError(f"Field '{field_name}' contains invalid values: {bad}")
    return numeric.astype(int).astype(str)


def _bool_series(series: pd.Series, field: str) -> pd.Series:
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


def _patch_id(frame: pd.DataFrame) -> pd.Series:
    if 'PATCH_ID' in frame.columns:
        ids = frame['PATCH_ID'].astype(str).str.strip()
        if not ids.map(lambda value: bool(re.fullmatch(r'[A-Za-z0-9_-]+', value))).all():
            raise ValueError('PATCH_ID must contain only letters, digits, underscores and hyphens.')
        return ids
    return (
        frame["name"].astype(str).str.strip()
        + "_" + _clean_int_series(frame["split_x"], "split_x")
        + "_" + _clean_int_series(frame["split_y"], "split_y")
    )


def _read_manifest(path, prefix: str, expected_patches: int) -> pd.DataFrame:
    manifest = pd.read_csv(path, keep_default_na=False)
    required = {"filename", "has_grass", "grass_pixels"}
    missing = required - set(manifest.columns)
    if missing:
        raise KeyError(f"{path} missing columns: {sorted(missing)}")
    if len(manifest) != expected_patches:
        raise RuntimeError(
            f"Expected {expected_patches} rows in {path.name}, got {len(manifest)}."
        )
    manifest["PATCH_ID"] = _patch_id(manifest)
    if manifest["PATCH_ID"].duplicated().any():
        dup = manifest.loc[manifest["PATCH_ID"].duplicated(keep=False), "PATCH_ID"]
        raise RuntimeError(f"Duplicate PATCH_ID in {path}: {sorted(dup.unique())}")
    manifest[f"HAS_{prefix}_GRASS"] = _bool_series(manifest["has_grass"], "has_grass")
    manifest[f"{prefix}_GRASS_PIXELS_MANIFEST"] = pd.to_numeric(
        manifest["grass_pixels"], errors="raise"
    ).astype(int)
    manifest[f"{prefix}_GRASS_FILENAME"] = manifest["filename"].astype(str).str.strip()
    return manifest


def _validate_manifest_files(table: pd.DataFrame, prefix: str, root, path_key: str) -> list[str]:
    problems = []
    has_col = f"HAS_{prefix}_GRASS"
    pixels_col = f"{prefix}_GRASS_PIXELS_MANIFEST"
    filename_col = f"{prefix}_GRASS_FILENAME"

    for _, row in table.iterrows():
        patch_id = str(row.PATCH_ID)
        path = patch_paths(patch_id)[path_key]
        expected_name = f"{patch_id}_grass.tif"
        has_grass = bool(row[has_col])
        pixels = int(row[pixels_col])
        filename = str(row[filename_col]).strip()

        if has_grass:
            if pixels <= 0:
                problems.append(f"{patch_id}: {has_col}=True but pixels={pixels}")
            if filename != expected_name:
                problems.append(
                    f"{patch_id}: manifest filename={filename!r}, expected={expected_name!r}"
                )
            if not path.exists():
                problems.append(f"{patch_id}: required mask missing: {path}")
            else:
                try:
                    arr, _, _, _, bounds, _, _ = load_grass(path)
                    if int((arr == 1).sum()) != pixels:
                        problems.append(f'{patch_id}: manifest grass_pixels differs from raster')
                    if not np.allclose([bounds.left, bounds.bottom, bounds.right, bounds.top],
                                       [row.WEST, row.SOUTH, row.EAST, row.NORTH], rtol=0, atol=0.01):
                        problems.append(f'{patch_id}: mask bounds differ from acquisition tile')
                except Exception as exc:
                    problems.append(f'{patch_id}: invalid mask: {exc}')
        else:
            if pixels != 0:
                problems.append(f"{patch_id}: {has_col}=False but pixels={pixels}")
            if filename:
                problems.append(
                    f"{patch_id}: {has_col}=False but manifest filename={filename!r}"
                )
            if path.exists():
                problems.append(f"{patch_id}: empty-tile mask should not exist: {path}")
    return problems


def _validate_stats_subset(table: pd.DataFrame) -> list[str]:
    """Strict -7.5 m statistics mask must be a same-grid subset of +5 m detection mask."""
    problems = []
    for _, row in table.iterrows():
        if not bool(row["HAS_STATS_GRASS"]):
            continue
        patch_id = str(row.PATCH_ID)
        if not bool(row["HAS_DETECTION_GRASS"]):
            problems.append(f"{patch_id}: statistics mask has grass but detection mask does not")
            continue
        det_path = patch_paths(patch_id)["grass"]
        stat_path = patch_paths(patch_id)["stats_grass"]
        try:
            with rasterio.open(det_path) as det, rasterio.open(stat_path) as stat:
                if det.crs != stat.crs or det.transform != stat.transform or det.shape != stat.shape:
                    problems.append(f"{patch_id}: statistics mask grid differs from detection mask")
                    continue
                det_arr = det.read(1)
                stat_arr = stat.read(1)
                if not np.isin(stat_arr, [0, 1]).all():
                    problems.append(f"{patch_id}: statistics mask is not binary 0/1")
                if np.any((stat_arr == 1) & (det_arr != 1)):
                    problems.append(
                        f"{patch_id}: -7.5 m statistics mask is not a subset of +5 m detection mask"
                    )
        except Exception as exc:
            problems.append(f"{patch_id}: statistics-mask validation failed: {exc}")
    return problems


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Prepare the national EOPatch table. Detection uses +5 m masks. "
            "Statistics use separate -7.5 m masks when available."
        )
    )
    parser.add_argument(
        "--require-stats",
        action="store_true",
        help="Fail if the -7.5 m statistics-mask folder/manifest is not ready.",
    )
    args = parser.parse_args()

    print("=" * 78)
    print("PREPARE NATIONAL EOPATCH TABLE")
    print("=" * 78)
    print("YEAR:", YEAR)
    print("ROOT:", ROOT)
    print("PROCESS_ROOT:", PROCESS_ROOT)
    print("DETECTION GRASS ROOT (+5 m):", GRASS_ROOT)
    print("DETECTION MANIFEST:", GRASS_MANIFEST)
    print("STATS GRASS ROOT (-7.5 m):", STATS_GRASS_ROOT)
    print("STATS MANIFEST:", STATS_GRASS_MANIFEST)
    print("PATCH_GPKG:", PATCH_GPKG)
    print("PATCH_LAYER:", PATCH_LAYER or "<first layer>")

    ROOT.mkdir(parents=True, exist_ok=True)
    PROCESS_ROOT.mkdir(parents=True, exist_ok=True)

    if not PATCH_GPKG.exists():
        raise FileNotFoundError(f"National GeoPackage not found:\n{PATCH_GPKG}")
    if not GRASS_ROOT.exists():
        raise FileNotFoundError(f"Detection mask directory not found:\n{GRASS_ROOT}")
    if not GRASS_MANIFEST.exists():
        raise FileNotFoundError(
            f"Detection mask manifest not found:\n{GRASS_MANIFEST}\n"
            "Re-create the +5 m masks with grass_masks_manifest.csv."
        )

    gdf = gpd.read_file(PATCH_GPKG, layer=PATCH_LAYER) if PATCH_LAYER else gpd.read_file(PATCH_GPKG)
    if len(gdf) != EXPECTED_PATCHES:
        raise RuntimeError(f"Expected {EXPECTED_PATCHES} EOPatches, got {len(gdf)}.")

    if "PATCH_ID" in gdf.columns:
        gdf["PATCH_ID"] = _patch_id(gdf)
    else:
        required = ["name", "split_x", "split_y"]
        missing = [f for f in required if f not in gdf.columns]
        if missing:
            raise KeyError(f"Cannot create PATCH_ID; missing fields: {missing}")
        gdf["PATCH_ID"] = _patch_id(gdf)

    if gdf["PATCH_ID"].duplicated().any() or gdf["PATCH_ID"].eq("").any():
        raise RuntimeError("PATCH_ID must be unique and non-empty.")

    for field, default in [('PET', 'all'), ('BIOGEO_R', 'all'), ('AREA_PCT', 100.0)]:
        if field not in gdf.columns:
            gdf[field] = default
    if "TYPE" not in gdf.columns:
        gdf["TYPE"] = ""
    if "NOTE" not in gdf.columns:
        gdf["NOTE"] = ""

    det = _read_manifest(GRASS_MANIFEST, "DETECTION", EXPECTED_PATCHES)
    if set(det.PATCH_ID) != set(gdf.PATCH_ID):
        raise RuntimeError("PATCH_ID mismatch between EOPatch_PET.gpkg and detection mask manifest.")

    keep_det = [
        "PATCH_ID", "HAS_DETECTION_GRASS", "DETECTION_GRASS_PIXELS_MANIFEST",
        "DETECTION_GRASS_FILENAME",
    ]
    keep_gdf = ["PATCH_ID", "PET", "BIOGEO_R", "AREA_PCT", "TYPE", "NOTE"]
    for col in ["name", "split_x", "split_y"]:
        if col in gdf.columns:
            keep_gdf.append(col)

    table = gdf[keep_gdf].merge(det[keep_det], on="PATCH_ID", validate="one_to_one")
    # Reuse the same tile reconstruction as acquisition, even for projected input grids.
    from .config import TARGET_EPSG
    prepare = import_module('moweo.acquisition_prepare')
    projected = gdf.to_crs(epsg=TARGET_EPSG)
    bounds = {str(row.PATCH_ID): prepare._snap_tile_from_centroid(row.geometry, str(row.PATCH_ID))[:4]
              for _, row in projected.iterrows()}
    for index, name in enumerate(['WEST', 'SOUTH', 'EAST', 'NORTH']):
        table[name] = table.PATCH_ID.map(lambda pid: bounds[pid][index])
    det_problems = _validate_manifest_files(table, "DETECTION", GRASS_ROOT, "grass")
    if det_problems:
        print("\nDetection-mask problems:")
        for problem in det_problems:
            print(" -", problem)
        raise RuntimeError("Detection-mask validation failed. Acquisition/processing is blocked.")

    # Backwards-compatible name used by 01/02/03.
    table["HAS_GRASS"] = table["HAS_DETECTION_GRASS"].astype(bool)
    table["GRASS_PIXELS_MANIFEST"] = table["DETECTION_GRASS_PIXELS_MANIFEST"]
    table["GRASS_FILE"] = [
        str(patch_paths(pid)["grass"]) if has else ""
        for pid, has in zip(table["PATCH_ID"], table["HAS_DETECTION_GRASS"])
    ]
    table["DETECTION_MASK_READY"] = True

    # Statistics masks are intentionally optional for 01/02. 03 will refuse to
    # calculate national statistics until these are ready and 00 is re-run.
    stats_configured = STATS_GRASS_ROOT.exists() and STATS_GRASS_MANIFEST.exists()
    if stats_configured:
        stats = _read_manifest(STATS_GRASS_MANIFEST, "STATS", EXPECTED_PATCHES)
        if set(stats.PATCH_ID) != set(gdf.PATCH_ID):
            raise RuntimeError("PATCH_ID mismatch between EOPatch_PET.gpkg and statistics mask manifest.")
        keep_stats = [
            "PATCH_ID", "HAS_STATS_GRASS", "STATS_GRASS_PIXELS_MANIFEST",
            "STATS_GRASS_FILENAME",
        ]
        table = table.merge(stats[keep_stats], on="PATCH_ID", validate="one_to_one")
        stats_problems = _validate_manifest_files(table, "STATS", STATS_GRASS_ROOT, "stats_grass")
        stats_problems += _validate_stats_subset(table)
        if stats_problems:
            print("\nStatistics-mask problems:")
            for problem in stats_problems:
                print(" -", problem)
            raise RuntimeError("Statistics-mask validation failed.")
        table["STATS_MASK_READY"] = True
        table["STATS_GRASS_FILE"] = [
            str(patch_paths(pid)["stats_grass"]) if has else ""
            for pid, has in zip(table.PATCH_ID, table.HAS_STATS_GRASS)
        ]
    else:
        if args.require_stats:
            raise FileNotFoundError(
                "Statistics masks are required but not ready. Expected:\n"
                f"{STATS_GRASS_ROOT}\n{STATS_GRASS_MANIFEST}"
            )
        table["HAS_STATS_GRASS"] = False
        table["STATS_GRASS_PIXELS_MANIFEST"] = np.nan
        table["STATS_GRASS_FILENAME"] = ""
        table["STATS_GRASS_FILE"] = ""
        table["STATS_MASK_READY"] = False

    table["STATS_CONFIGURED"] = bool(stats_configured)
    table["PATCH_DIR"] = table["PATCH_ID"].map(lambda x: str(patch_paths(x)["base"]))
    table = table.sort_values("PATCH_ID").reset_index(drop=True)

    n_total = len(table)
    n_detection = int(table.HAS_DETECTION_GRASS.sum())
    n_no_detection = n_total - n_detection
    print("\nNational patches:", n_total)
    print("Detection patches WITH grass (+5 m):", n_detection)
    print("Detection patches WITHOUT grass:", n_no_detection)
    print("Required +5 m masks found:", n_detection, "/", n_detection)

    if stats_configured:
        n_stats = int(table.HAS_STATS_GRASS.sum())
        print("Statistics masks configured: YES")
        print("Statistics patches WITH robust inner pixels (-7.5 m):", n_stats)
        print("Statistics patches WITHOUT robust inner pixels:", n_total - n_stats)
        print("Statistics mask grid/subset validation: OK")
    else:
        print("Statistics masks configured: NO")
        print("NOTE: 01 acquisition and 02 detection can run now.")
        print("      Before 11_finalize.py, create -7.5 m masks and re-run 09_prepare_processing.py.")

    tmp = PATCH_CSV.with_name(PATCH_CSV.stem + ".part" + PATCH_CSV.suffix)
    table.to_csv(tmp, index=False, encoding="utf-8-sig")
    tmp.replace(PATCH_CSV)
    print("\nSaved:", PATCH_CSV)
    print("READY for acquisition/processing with +5 m detection masks.")


if __name__ == "__main__":
    main()
