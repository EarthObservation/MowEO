from __future__ import annotations

from pathlib import Path
import math
import re

import numpy as np
import pandas as pd
import rasterio
from rasterio.features import sieve
from rasterio.windows import Window
from scipy import ndimage
import xarray as xr

from .config import YEAR, PATCH_CSV, PROCESS_ROOT, GRASS_ROOT, FINAL_VERSION, OUTPUT_PREFIX, patch_paths
from .core import read_patch_table, load_grass

OUT_DIR = PROCESS_ROOT.parent / "national_maps"
OUT_DIR.mkdir(parents=True, exist_ok=True)

HARMONIZED_CONFIDENCE_THRESHOLD = 60.0
MAJORITY_MIN_SUPPORT_3X3 = 5
DOY_SEIVE_MIN_SIZE_PIXELS = 25


def natural_key(value):
    return [int(x) if x.isdigit() else x.lower() for x in re.split(r"(\d+)", str(value))]


def final_tif_path(patch_id: str) -> Path:
    return PROCESS_ROOT / str(patch_id) / f"{patch_id}_MOWING_FINAL_{FINAL_VERSION}_{YEAR}.tif"


def det_mask_path(patch_id: str) -> Path:
    return GRASS_ROOT / f"{patch_id}_grass.tif"


def same_grid(a, b) -> bool:
    return (
        a.crs == b.crs
        and a.width == b.width
        and a.height == b.height
        and np.allclose(tuple(a.transform), tuple(b.transform), atol=1e-9, rtol=0)
    )


def window_from_aligned_transform(src_transform, src_width, src_height, out_transform):
    xres = out_transform.a
    yres = abs(out_transform.e)
    col_off = int(round((src_transform.c - out_transform.c) / xres))
    row_off = int(round((out_transform.f - src_transform.f) / yres))
    return Window(col_off=col_off, row_off=row_off, width=src_width, height=src_height)


def band_index(src, name: str) -> int:
    for i, desc in enumerate(src.descriptions or (), start=1):
        if desc and desc.strip().upper() == name.upper():
            return i
    raise RuntimeError(f"Band {name!r} not found in {src.name}. Available: {src.descriptions}")


def unique_ndvi_dates(ndvi_nc: Path) -> int:
    if not ndvi_nc.exists():
        return 0
    with xr.open_dataset(ndvi_nc, engine="h5netcdf") as ds:
        if "t" not in ds.coords:
            return 0
        dates = pd.DatetimeIndex(pd.to_datetime(ds.t.values)).normalize()
    return int(len(pd.unique(dates)))


def output_specs():
    specs = [
        ("MOWING_COUNT", f"{OUTPUT_PREFIX}_MOWING_COUNT_{YEAR}.tif", "uint8", 255),
        ("MOWING_COUNT_SPATIAL_HARMONIZED", f"{OUTPUT_PREFIX}_MOWING_COUNT_SPATIAL_HARMONIZED_{YEAR}.tif", "uint8", 255),
        ("MOWING_COUNT_SEIVE", f"{OUTPUT_PREFIX}_MOWING_COUNT_SEIVE_{YEAR}.tif", "uint8", 255),
    ]
    for i in range(1, 7):
        specs.append((f"MOWING_DOY_{i}", f"{OUTPUT_PREFIX}_MOWING_DOY_{i}_{YEAR}.tif", "uint16", 65535))
    specs.extend([
        ("PIXEL_EVENT_CONF_MEAN", f"{OUTPUT_PREFIX}_MOWING_CONFIDENCE_{YEAR}.tif", "uint8", 255),
        ("VALID_OBS", f"{OUTPUT_PREFIX}_VALID_OBSERVATIONS_{YEAR}.tif", "uint16", 65535),
        ("SCL_MASKED_PERCENT", f"{OUTPUT_PREFIX}_SCL_MASKED_PERCENT_{YEAR}.tif", "float32", -9999.0),
    ])
    return specs


def majority_value(values: np.ndarray) -> tuple[int, int]:
    values = values[(values >= 0) & (values <= 6)]
    if values.size == 0:
        return 255, 0
    counts = np.bincount(values.astype(np.uint8), minlength=7)
    value = int(np.argmax(counts))
    return value, int(counts[value])


def harmonize_mowing_count(count: np.ndarray, confidence: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Create confidence-aware spatially harmonized mowing count for one patch.

    The raw count remains unchanged elsewhere. Harmonization only changes
    low-confidence pixels (<60) inside the grass mask.
    """
    valid = mask & np.isfinite(count) & (count >= 0) & (count <= 6)
    low_conf = valid & (~np.isfinite(confidence) | (confidence < HARMONIZED_CONFIDENCE_THRESHOLD))
    out = np.full(count.shape, 255, dtype=np.uint8)
    out[valid] = count[valid].astype(np.uint8)

    # Step 1: local 3x3 majority within the grass mask. A low-confidence pixel
    # is changed only when at least five valid pixels in its 3x3 neighbourhood
    # agree on the same count.
    padded_count = np.pad(out, 1, mode="constant", constant_values=255)
    padded_valid = np.pad(valid, 1, mode="constant", constant_values=False)
    local = out.copy()
    rows, cols = np.where(low_conf)
    for r, c in zip(rows, cols):
        window_counts = padded_count[r : r + 3, c : c + 3]
        window_valid = padded_valid[r : r + 3, c : c + 3]
        mode, support = majority_value(window_counts[window_valid])
        if mode != 255 and support >= MAJORITY_MIN_SUPPORT_3X3 and mode != out[r, c]:
            local[r, c] = mode

    # Step 2: object-level consistency. In the absence of parcel IDs in this
    # raster workflow, objects are connected grass components inside the patch.
    # Low-confidence pixels that differ from the component-dominant count are
    # set to that dominant count.
    labels, n_labels = ndimage.label(mask, structure=np.ones((3, 3), dtype=np.uint8))
    harmonized = local.copy()
    for label_id in range(1, n_labels + 1):
        obj = labels == label_id
        obj_valid = obj & valid
        if not obj_valid.any():
            continue
        dominant, _ = majority_value(local[obj_valid])
        if dominant == 255:
            continue
        replace = obj_valid & low_conf & (local != dominant)
        harmonized[replace] = dominant

    harmonized[~valid] = 255
    return harmonized


def seive_doy_count(src, mask: np.ndarray) -> np.ndarray:
    """Sieve each DOY event layer and count remaining non-zero DOY events."""
    count = np.zeros(mask.shape, dtype=np.uint8)
    for i in range(1, 7):
        doy = src.read(band_index(src, f"MOWING_DOY_{i}")).astype(np.int32)
        valid = mask & (doy >= 1) & (doy <= 366)
        layer = np.zeros(doy.shape, dtype=np.int32)
        layer[valid] = doy[valid]

        # GRAME applies an MMU of 0.25 ha, equivalent to 25 10-m pixels.
        # Rasterio/GDAL sieve removes connected regions below this size from
        # the multi-valued date layer. The count is then derived from the
        # remaining non-zero event-date layers.
        sieved = sieve(
            layer,
            size=DOY_SEIVE_MIN_SIZE_PIXELS,
            connectivity=8,
            mask=mask,
        )
        count += ((sieved >= 1) & (sieved <= 366)).astype(np.uint8)

    out = np.full(mask.shape, 255, dtype=np.uint8)
    out[mask] = count[mask]
    return out


def main():
    print("=" * 78)
    print("NATIONAL v0.3 PRODUCTS")
    print("=" * 78)
    print("YEAR:", YEAR)
    print("PATCH TABLE:", PATCH_CSV)
    print("OUTPUT:", OUT_DIR)

    if not PATCH_CSV.exists():
        raise FileNotFoundError(PATCH_CSV)
    df = read_patch_table(PATCH_CSV)
    df = df[df['HAS_GRASS']].copy()
    patch_ids = sorted(df["PATCH_ID"].astype(str).tolist(), key=natural_key)

    final_files = [(pid, final_tif_path(pid)) for pid in patch_ids if final_tif_path(pid).exists()]
    missing = [pid for pid in patch_ids if not final_tif_path(pid).exists()]
    if missing:
        raise RuntimeError(f'Missing final rasters for grassland tiles: {missing[:20]}')
    for pid, path in final_files:
        load_grass(det_mask_path(pid))
        with rasterio.open(path) as src:
            if src.tags().get('VERSION') != 'v0.3':
                raise RuntimeError(f'{pid}: recalculate with the delivery v0.3 detector before export.')
            for name in ['MOWING_COUNT', 'PIXEL_EVENT_CONF_MEAN', 'VALID_OBS'] + [f'MOWING_DOY_{i}' for i in range(1, 7)]:
                band_index(src, name)
        if unique_ndvi_dates(patch_paths(pid)['ndvi']) < 1:
            raise RuntimeError(f'{pid}: raw NDVI dates are required for the unavailable-observation layer.')
    if not final_files:
        raise RuntimeError("No final v0.3 patch rasters found.")

    with rasterio.open(final_files[0][1]) as src0:
        crs = src0.crs
        xres = src0.transform.a
        yres = abs(src0.transform.e)

    left, bottom, right, top = math.inf, math.inf, -math.inf, -math.inf
    for _, path in final_files:
        with rasterio.open(path) as src:
            if src.crs != crs:
                raise RuntimeError(f"CRS mismatch: {path}")
            b = src.bounds
            left, bottom = min(left, b.left), min(bottom, b.bottom)
            right, top = max(right, b.right), max(top, b.top)

    left = math.floor(left / xres) * xres
    right = math.ceil(right / xres) * xres
    bottom = math.floor(bottom / yres) * yres
    top = math.ceil(top / yres) * yres
    width = int(round((right - left) / xres))
    height = int(round((top - bottom) / yres))
    out_transform = rasterio.transform.from_origin(left, top, xres, yres)

    specs = output_specs()
    writers = {}
    try:
        for variable, filename, dtype, nodata in specs:
            profile = {
                "driver": "GTiff", "height": height, "width": width, "count": 1,
                "dtype": dtype, "crs": crs, "transform": out_transform, "nodata": nodata,
                "compress": "DEFLATE", "tiled": True, "blockxsize": 512, "blockysize": 512,
                "BIGTIFF": "YES",
            }
            if dtype == "float32":
                profile["predictor"] = 3
            else:
                profile["predictor"] = 2
            path = OUT_DIR / (filename + ".part.tif")
            dst = rasterio.open(path, "w", **profile)
            dst.set_band_description(1, variable)
            for _, window in dst.block_windows(1):
                dst.write(np.full((int(window.height), int(window.width)), nodata, dtype=dtype), 1, window=window)
            tags = {
                "YEAR": str(YEAR),
                "MASK": "Detection grass mask buffered +5 m",
                "VERSION": "v0.3",
                "NODATA": str(nodata),
                "DETECTOR_PROFILE": "delivery_v0_3",
                "PROCESSING_SCOPE": "Per-tile filtering before mosaicking; connectivity does not cross tile boundaries",
            }
            if variable == "PIXEL_EVENT_CONF_MEAN":
                tags.update({
                    "VARIABLE": "Mean confidence score of detected mowing events (0-100)",
                    "NOTE": "Evidence score, not a calibrated probability; UInt8 rounded 0-100, NoData=255",
                })
            elif variable == "MOWING_COUNT_SPATIAL_HARMONIZED":
                tags.update({
                    "VARIABLE": "Spatially harmonized mowing count",
                    "SOURCE": "MOWING_COUNT + PIXEL_EVENT_CONF_MEAN + detection grass mask",
                    "RULE_1": (
                        "3x3 majority filter inside grass mask; low-confidence pixels "
                        f"(<{HARMONIZED_CONFIDENCE_THRESHOLD:g}) change only if majority support >= {MAJORITY_MIN_SUPPORT_3X3}"
                    ),
                    "RULE_2": (
                        "Connected grass object dominant count replaces differing low-confidence pixels"
                    ),
                    "NOTE": "Raw SI_MOWING_COUNT is also written unchanged.",
                })
            elif variable == "MOWING_COUNT_SEIVE":
                tags.update({
                    "VARIABLE": "Mowing count derived from sieved DOY event layers",
                    "SOURCE": "MOWING_DOY_1..6 + detection grass mask",
                    "RULE": (
                        f"Each DOY layer is sieved with size={DOY_SEIVE_MIN_SIZE_PIXELS} "
                        "10-m pixels and 8-connectivity; count is number of remaining non-zero DOY layers"
                    ),
                    "NOTE": "Name intentionally uses SEIVE for easy visual distinction from other test outputs.",
                })
            elif variable == "SCL_MASKED_PERCENT":
                tags.update({
                    "VARIABLE": "Percent of NDVI acquisition dates unavailable after preprocessing",
                    "FORMULA": "100*(1-VALID_OBS/N_UNIQUE_NDVI_DATES)",
                    "IMPORTANT": "Not pure meteorological cloud cover; may include non-coverage and other masked/invalid observations",
                })
            dst.update_tags(**tags)
            writers[variable] = dst

        for i, (patch_id, final_path) in enumerate(final_files, start=1):
            mask_path = det_mask_path(patch_id)
            if not mask_path.exists():
                raise FileNotFoundError(mask_path)
            with rasterio.open(final_path) as src, rasterio.open(mask_path) as msrc:
                if not same_grid(src, msrc):
                    raise RuntimeError(f"Grid mismatch for {patch_id}")
                mask = msrc.read(1) == 1
                win = window_from_aligned_transform(src.transform, src.width, src.height, out_transform)

                # Direct bands from the v0.3 final product.
                for variable, _, dtype, nodata in specs:
                    if variable == "SCL_MASKED_PERCENT":
                        continue
                    if variable == "MOWING_COUNT_SPATIAL_HARMONIZED":
                        count = src.read(band_index(src, "MOWING_COUNT")).astype(np.float32)
                        conf = src.read(band_index(src, "PIXEL_EVENT_CONF_MEAN")).astype(np.float32)
                        conf[conf == src.nodata] = np.nan
                        out = harmonize_mowing_count(count, conf, mask)
                        writers[variable].write(out, 1, window=win)
                        continue
                    if variable == "MOWING_COUNT_SEIVE":
                        out = seive_doy_count(src, mask)
                        writers[variable].write(out, 1, window=win)
                        continue
                    idx = band_index(src, variable)
                    data = src.read(idx).astype(np.float32)
                    out = np.full(data.shape, nodata, dtype=np.dtype(dtype))
                    valid = mask & np.isfinite(data) & (data != src.nodata)
                    if variable == "MOWING_COUNT":
                        valid &= (data >= 0) & (data <= 6)
                    elif variable.startswith("MOWING_DOY_"):
                        valid &= (data >= 0) & (data <= 366)
                    elif variable == "PIXEL_EVENT_CONF_MEAN":
                        valid &= (data >= 0) & (data <= 100)
                    elif variable == "VALID_OBS":
                        valid &= data >= 0
                    if variable == 'PIXEL_EVENT_CONF_MEAN':
                        out[valid] = np.rint(np.clip(data[valid], 0, 100)).astype(dtype)
                    else:
                        out[valid] = data[valid].astype(dtype)
                    writers[variable].write(out, 1, window=win)

                # Post-SCL unavailable fraction derived from the same raw NDVI cube.
                valid_obs = src.read(band_index(src, "VALID_OBS")).astype(np.float32)
                n_dates = unique_ndvi_dates(patch_paths(patch_id)["ndvi"])
                scl_out = np.full(valid_obs.shape, -9999.0, dtype=np.float32)
                if n_dates > 0:
                    valid = mask & np.isfinite(valid_obs) & (valid_obs != src.nodata)
                    scl_out[valid] = np.clip(100.0 * (1.0 - valid_obs[valid] / float(n_dates)), 0.0, 100.0)
                writers["SCL_MASKED_PERCENT"].write(scl_out, 1, window=win)

            if i % 25 == 0 or i == len(final_files):
                print(f"Mosaicked: {i}/{len(final_files)}")
    finally:
        for dst in writers.values():
            dst.close()

    for _, filename, _, _ in specs:
        (OUT_DIR / (filename + '.part.tif')).replace(OUT_DIR / filename)
    import json
    manifest = [{'file': name, 'variable': variable, 'dtype': dtype, 'nodata': nodata} for variable, name, dtype, nodata in specs]
    (OUT_DIR / 'products_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
    print("\nDONE. National products:")
    for _, filename, _, _ in specs:
        print(" -", OUT_DIR / filename)


if __name__ == "__main__":
    main()
