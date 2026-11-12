from __future__ import annotations

from pathlib import Path
import re

import geopandas as gpd
import numpy as np
import pandas as pd
from pyproj import CRS
from shapely.geometry import box

from .config import (
    EXPECTED_PATCHES,
    PATCH_CSV_ACQUISITION,
    PATCH_GPKG,
    PATCH_LAYER,
    PATCH_SELECTION_FILE,
    PATCH_SIZE_M,
    RAW_PROCESS_ROOT,
    OUTPUT_ROOT,
    TARGET_EPSG,
    YEAR,
)


def _clean_int_series(series: pd.Series, field_name: str) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    if numeric.isna().any():
        bad = series[numeric.isna()].tolist()
        raise ValueError(f"Field '{field_name}' contains invalid values: {bad[:10]}")
    return numeric.astype(int).astype(str)


def _make_patch_id(frame: pd.DataFrame) -> pd.Series:
    if "PATCH_ID" in frame.columns:
        out = frame["PATCH_ID"].astype(str).str.strip()
        if out.ne("").all():
            return out

    required = ["name", "split_x", "split_y"]
    missing = [f for f in required if f not in frame.columns]
    if missing:
        raise KeyError(
            "Cannot create PATCH_ID. GeoPackage must contain PATCH_ID or "
            f"the fields {required}. Missing: {missing}"
        )

    return (
        frame["name"].astype(str).str.strip()
        + "_" + _clean_int_series(frame["split_x"], "split_x")
        + "_" + _clean_int_series(frame["split_y"], "split_y")
    )


def _snap_tile_from_centroid(geom, patch_id: str) -> tuple[float, float, float, float, float]:
    """Recover the exact 5 km UTM tile from the EOPatch geometry centroid.

    The national EOPatch grid is aligned to 5 km UTM coordinates. Snapping the
    centroid is more robust than trusting floating bounds after a CRS round-trip.
    The overlap check prevents silently assigning a geometry to the wrong tile.
    """
    if geom is None or geom.is_empty:
        raise ValueError(f"{patch_id}: empty geometry")

    c = geom.centroid
    size = PATCH_SIZE_M
    half = size / 2.0

    west = round((c.x - half) / size) * size
    south = round((c.y - half) / size) * size
    east = west + size
    north = south + size

    tile = box(west, south, east, north)
    if not tile.contains(c) and not tile.touches(c):
        raise RuntimeError(f"{patch_id}: centroid is outside reconstructed 5 km tile")

    geom_area = float(geom.area)
    if geom_area <= 0:
        raise ValueError(f"{patch_id}: invalid geometry area")

    overlap = float(geom.intersection(tile).area)
    overlap_ratio = overlap / min(geom_area, size * size)

    # EOPatch polygons should almost exactly match the reconstructed tile.
    # A generous threshold also tolerates small projection/rounding differences.
    if overlap_ratio < 0.90:
        raise RuntimeError(
            f"{patch_id}: reconstructed 5 km tile overlaps geometry by only "
            f"{overlap_ratio:.3f}. Check the GeoPackage CRS/grid."
        )

    return west, south, east, north, overlap_ratio


def _read_selection(path: Path) -> list[str] | None:
    if not path.exists():
        return None
    ids = []
    for line in path.read_text(encoding="utf-8").splitlines():
        value = line.strip()
        if value and not value.startswith("#"):
            ids.append(value)
    if not ids:
        raise RuntimeError(f"{path} exists but contains no PATCH_ID values.")
    if len(ids) != len(set(ids)):
        raise RuntimeError(f"Duplicate PATCH_ID values in {path}.")
    return ids


def main():
    crs = CRS.from_epsg(TARGET_EPSG)
    if not crs.is_projected or any(axis.unit_name != 'metre' for axis in crs.axis_info):
        raise ValueError('target_epsg must be a projected CRS with metre units.')
    print("=" * 78)
    print("PREPARE ACQUISITION TABLE FROM EOPatch GeoPackage")
    print("=" * 78)
    print("YEAR:", YEAR)
    print("GPKG:", PATCH_GPKG)
    print("TARGET CRS: EPSG:", TARGET_EPSG, sep="")
    print("Acquisition does not require grass masks; local detection does.")
    print()

    if not PATCH_GPKG.exists():
        raise FileNotFoundError(
            f"GeoPackage not found:\n{PATCH_GPKG}\n\n"
            "Copy EOPatch_PET.gpkg into the data folder."
        )

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    RAW_PROCESS_ROOT.mkdir(parents=True, exist_ok=True)

    gdf = (
        gpd.read_file(PATCH_GPKG, layer=PATCH_LAYER)
        if PATCH_LAYER
        else gpd.read_file(PATCH_GPKG)
    )

    print("Input rows:", len(gdf))
    print("Input CRS:", gdf.crs)

    if len(gdf) != EXPECTED_PATCHES:
        raise RuntimeError(
            f"Expected {EXPECTED_PATCHES} national EOPatches, got {len(gdf)}."
        )
    if gdf.crs is None:
        raise RuntimeError("GeoPackage layer has no CRS.")

    gdf = gdf.copy()
    gdf["PATCH_ID"] = _make_patch_id(gdf)
    if not gdf['PATCH_ID'].map(lambda value: bool(re.fullmatch(r'[A-Za-z0-9_-]+', value))).all():
        raise ValueError('PATCH_ID must contain only letters, digits, underscores and hyphens.')
    if gdf["PATCH_ID"].duplicated().any() or gdf["PATCH_ID"].eq("").any():
        raise RuntimeError("PATCH_ID values must be unique and non-empty.")

    gdf_utm = gdf.to_crs(epsg=TARGET_EPSG)

    records = []
    for (_, src_row), (_, utm_row) in zip(gdf.iterrows(), gdf_utm.iterrows()):
        patch_id = str(src_row["PATCH_ID"])
        west, south, east, north, overlap_ratio = _snap_tile_from_centroid(
            utm_row.geometry, patch_id
        )

        rec = {
            "PATCH_ID": patch_id,
            "WEST": west,
            "SOUTH": south,
            "EAST": east,
            "NORTH": north,
            "EPSG": TARGET_EPSG,
            "WIDTH10": int(round(PATCH_SIZE_M / 10.0)),
            "HEIGHT10": int(round(PATCH_SIZE_M / 10.0)),
            "WIDTH20": int(round(PATCH_SIZE_M / 20.0)),
            "HEIGHT20": int(round(PATCH_SIZE_M / 20.0)),
            "GEOM_TILE_OVERLAP": overlap_ratio,
        }
        for col in ["PET", "BIOGEO_R", "AREA_PCT", "TYPE", "NOTE", "name", "split_x", "split_y"]:
            if col in src_row.index:
                rec[col] = src_row[col]
        records.append(rec)

    table = pd.DataFrame(records)

    # A duplicated 5 km extent would mean the grid reconstruction failed.
    extent_cols = ["WEST", "SOUTH", "EAST", "NORTH"]
    if table.duplicated(extent_cols).any():
        bad = table.loc[table.duplicated(extent_cols, keep=False), ["PATCH_ID"] + extent_cols]
        raise RuntimeError(f"Duplicate reconstructed tile extents:\n{bad.head(20)}")

    selection = _read_selection(PATCH_SELECTION_FILE)
    if selection is not None:
        missing = sorted(set(selection) - set(table["PATCH_ID"]))
        if missing:
            raise RuntimeError(
                f"patch_ids.txt contains {len(missing)} unknown IDs. First: {missing[:10]}"
            )
        table = table[table["PATCH_ID"].isin(selection)].copy()
        order = {pid: i for i, pid in enumerate(selection)}
        table["_order"] = table["PATCH_ID"].map(order)
        table = table.sort_values("_order").drop(columns="_order")
        print(f"Selection file found: {PATCH_SELECTION_FILE}")
        print("Selected patches:", len(table))
    else:
        table = table.sort_values("PATCH_ID").reset_index(drop=True)
        print("No patch_ids.txt -> all national EOPatches will be acquired.")
        print("Selected patches:", len(table))

    table.to_csv(PATCH_CSV_ACQUISITION, index=False)

    print("Minimum geometry/tile overlap:", f"{table['GEOM_TILE_OVERLAP'].min():.4f}")
    print("Saved:", PATCH_CSV_ACQUISITION)
    print()
    print("READY FOR DRY RUN / ACQUISITION.")
    print("Each selected EOPatch will produce two NetCDF cubes:")
    print("  1) NDVI 10 m")
    print("  2) NDII + SWIRMean explicit 20 m")


if __name__ == "__main__":
    main()
