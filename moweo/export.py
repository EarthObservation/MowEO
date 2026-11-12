from pathlib import Path
import math
import re

import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window

from .config import (
    YEAR,
    PATCH_CSV,
    PROCESS_ROOT,
    GRASS_ROOT,
    FINAL_VERSION,
    OUTPUT_PREFIX,
)

# =============================================================================
# NATIONAL MOSAIC: NUMBER OF DETECTED MOWINGS ON THE +5 m DETECTION MASK
# -----------------------------------------------------------------------------
# Purpose:
#   Build a national mosaic of MOWING_COUNT, but display it on the +5 m
#   detection mask instead of the stricter -7.5 m statistics mask.
#
# Why this script exists:
#   This is a diagnostic / visualization product to check whether white seams
#   seen between neighboring EOPatches are caused by the strict statistics mask
#   or by the mosaicking itself.
#
# Output semantics:
#   - 0   = grass pixel in the +5 m detection mask, but no detected mowing
#   - 1-6 = number of detected mowings
#   - 255 = outside the +5 m detection mask / no data
# =============================================================================

OUT_DIR = PROCESS_ROOT.parent / "national_maps"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_TIF = OUT_DIR / f"{OUTPUT_PREFIX}_MOWING_COUNT_{YEAR}_detection_plus5m.tif"

NODATA = 255
OUT_DTYPE = "uint8"

TARGET_BAND_NAMES = {
    "MOWING_COUNT",
    "MOWINGCOUNT",
}


def natural_key(value):
    return [
        int(x) if x.isdigit() else x.lower()
        for x in re.split(r"(\d+)", str(value))
    ]


def final_tif_path(patch_id: str) -> Path:
    return (
        PROCESS_ROOT
        / str(patch_id)
        / f"{patch_id}_MOWING_FINAL_{FINAL_VERSION}_{YEAR}.tif"
    )


def det_mask_path(patch_id: str) -> Path:
    return GRASS_ROOT / f"{patch_id}_grass.tif"


def find_mowing_count_band(src: rasterio.io.DatasetReader) -> int:
    """
    Find MOWING_COUNT from band descriptions.
    Falls back to band 1, which is the expected band in the final product.
    """
    descriptions = src.descriptions or ()

    for i, desc in enumerate(descriptions, start=1):
        if not desc:
            continue
        key = desc.strip().upper().replace(" ", "_")
        if key in TARGET_BAND_NAMES:
            return i

    return 1


def same_grid(a: rasterio.io.DatasetReader, b: rasterio.io.DatasetReader) -> bool:
    return (
        a.crs == b.crs
        and a.width == b.width
        and a.height == b.height
        and np.allclose(tuple(a.transform), tuple(b.transform), atol=1e-9, rtol=0)
    )


def window_from_aligned_transform(src_transform, src_width, src_height, out_transform):
    """
    Compute an integer-aligned write window from affine transforms.

    This avoids occasional 1-pixel seams that can appear when using
    floating-point bounds -> window conversion.
    """
    xres = out_transform.a
    yres = abs(out_transform.e)

    col_off = int(round((src_transform.c - out_transform.c) / xres))
    row_off = int(round((out_transform.f - src_transform.f) / yres))

    return Window(col_off=col_off, row_off=row_off, width=src_width, height=src_height)


def main():
    print("=" * 78)
    print("NATIONAL MOWING_COUNT MOSAIC (+5 m DETECTION MASK)")
    print("=" * 78)
    print(f"Year: {YEAR}")
    print(f"Process root: {PROCESS_ROOT}")
    print(f"Detection masks (+5 m): {GRASS_ROOT}")
    print(f"Output: {OUT_TIF}")
    print()

    if not PATCH_CSV.exists():
        raise FileNotFoundError(f"PATCH_CSV not found: {PATCH_CSV}")

    from .core import read_patch_table
    df = read_patch_table(PATCH_CSV)
    df = df[df['HAS_GRASS']].copy()
    if "PATCH_ID" not in df.columns:
        raise ValueError("PATCH_CSV must contain column PATCH_ID")

    patch_ids = sorted(df["PATCH_ID"].astype(str).tolist(), key=natural_key)

    # -------------------------------------------------------------------------
    # 1. Collect existing final rasters and determine the national grid extent.
    # -------------------------------------------------------------------------
    final_files = []
    missing_final = []

    for patch_id in patch_ids:
        p = final_tif_path(patch_id)
        if p.exists():
            final_files.append((patch_id, p))
        else:
            missing_final.append(patch_id)

    if not final_files:
        raise RuntimeError("No final mowing rasters found.")
    if missing_final:
        raise RuntimeError(f'Missing final rasters for grassland tiles: {missing_final[:20]}')

    with rasterio.open(final_files[0][1]) as src0:
        crs = src0.crs
        xres = src0.transform.a
        yres = abs(src0.transform.e)
        band_idx = find_mowing_count_band(src0)

        print(f"Final rasters found: {len(final_files)}")
        print(f"Missing final rasters: {len(missing_final)}")
        print(f"CRS: {crs}")
        print(f"Resolution: {xres} x {yres} m")
        print(f"MOWING_COUNT band: {band_idx}")
        print(f"Band descriptions: {src0.descriptions}")
        print()

    left = math.inf
    bottom = math.inf
    right = -math.inf
    top = -math.inf

    for i, (_, path) in enumerate(final_files, start=1):
        with rasterio.open(path) as src:
            if src.crs != crs:
                raise RuntimeError(f"CRS mismatch: {path}")
            if not np.isclose(src.transform.a, xres) or not np.isclose(
                abs(src.transform.e), yres
            ):
                raise RuntimeError(f"Resolution mismatch: {path}")

            b = src.bounds
            left = min(left, b.left)
            bottom = min(bottom, b.bottom)
            right = max(right, b.right)
            top = max(top, b.top)

        if i % 100 == 0:
            print(f"Checked extents: {i}/{len(final_files)}")

    left = math.floor(left / xres) * xres
    right = math.ceil(right / xres) * xres
    bottom = math.floor(bottom / yres) * yres
    top = math.ceil(top / yres) * yres

    width = int(round((right - left) / xres))
    height = int(round((top - bottom) / yres))
    out_transform = rasterio.transform.from_origin(left, top, xres, yres)

    print()
    print(f"National raster size: {width} x {height}")
    print(f"Bounds: {left}, {bottom}, {right}, {top}")
    print()

    profile = {
        "driver": "GTiff",
        "height": height,
        "width": width,
        "count": 1,
        "dtype": OUT_DTYPE,
        "crs": crs,
        "transform": out_transform,
        "nodata": NODATA,
        "compress": "DEFLATE",
        "predictor": 2,
        "tiled": True,
        "blockxsize": 512,
        "blockysize": 512,
        "BIGTIFF": "YES",
    }

    written = 0
    skipped_no_mask = 0

    # -------------------------------------------------------------------------
    # 2. Write EOPatch by EOPatch on the +5 m detection mask.
    # -------------------------------------------------------------------------
    with rasterio.open(OUT_TIF, "w", **profile) as dst:
        # Explicitly initialize gaps between non-adjacent tiles to NoData.
        for _, window in dst.block_windows(1):
            dst.write(np.full((int(window.height), int(window.width)), NODATA, dtype=np.uint8),
                      1, window=window)
        dst.set_band_description(1, "MOWING_COUNT")
        dst.update_tags(
            YEAR=str(YEAR),
            VARIABLE="Number of detected mowings",
            MASK="Detection grass mask buffered +5 m",
            CLASS_0="No detected mowing",
            CLASS_1="1 mowing",
            CLASS_2="2 mowings",
            CLASS_3="3 mowings",
            CLASS_4="4 mowings",
            CLASS_5="5 mowings",
            CLASS_6="6 mowings",
            NODATA=str(NODATA),
        )

        for i, (patch_id, final_path) in enumerate(final_files, start=1):
            mask_path = det_mask_path(patch_id)

            if not mask_path.exists():
                raise FileNotFoundError(mask_path)

            with rasterio.open(final_path) as src, rasterio.open(mask_path) as msrc:
                if not same_grid(src, msrc):
                    raise RuntimeError(
                        "Grid mismatch between final raster and +5 m detection mask:\n"
                        f"  final: {final_path}\n"
                        f"  mask : {mask_path}"
                    )

                band_idx = find_mowing_count_band(src)
                data = src.read(band_idx)
                mask = msrc.read(1) > 0

                out = np.full(data.shape, NODATA, dtype=np.uint8)
                valid = mask & np.isfinite(data) & (data >= 0) & (data <= 6)
                out[valid] = data[valid].astype(np.uint8)

                win = window_from_aligned_transform(
                    src.transform, src.width, src.height, out_transform
                )

                if (
                    int(win.col_off) < 0
                    or int(win.row_off) < 0
                    or int(win.col_off + win.width) > width
                    or int(win.row_off + win.height) > height
                ):
                    raise RuntimeError(
                        f"Computed window falls outside output raster for {final_path}: {win}"
                    )

                dst.write(out, 1, window=win)
                written += 1

            if i % 25 == 0 or i == len(final_files):
                print(f"Written: {i}/{len(final_files)}")

    print()
    print("DONE")
    print(f"Output written: {OUT_TIF}")
    print(f"Patches written: {written}")
    print(f"Patches skipped (no +5 m mask): {skipped_no_mask}")
    print()
    print("Interpretation:")
    print(f"  0   = inside +5 m grass mask, no detected mowing")
    print(f"  1-6 = number of detected mowings")
    print(f"  {NODATA} = outside +5 m grass mask / no data")


if __name__ == "__main__":
    main()
