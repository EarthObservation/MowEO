from __future__ import annotations

import json
from pathlib import Path

HERE = Path.cwd().resolve()
OUT = HERE / "project_config.json"
EXAMPLE = Path(__file__).with_name("project_config.example.json")


def ask(prompt: str, default: str) -> str:
    value = input(f"{prompt} [{default}]: ").strip()
    return value or default


def main() -> None:
    defaults = json.loads(EXAMPLE.read_text(encoding="utf-8"))

    print("=" * 78)
    print("ZRSVN MOWING WORKFLOW - PROJECT SETUP")
    print("=" * 78)
    print("This writes local project_config.json. It is ignored by Git.")
    print()

    defaults["year"] = int(ask("Processing year", str(defaults["year"])))
    defaults['target_epsg'] = int(ask('Local metric projected EPSG code', str(defaults['target_epsg'])))
    defaults['study_area_name'] = ask('Study area name', defaults['study_area_name'])
    defaults['output_prefix'] = ask('Output filename prefix', defaults['output_prefix'])
    defaults["expected_patches"] = int(ask("Expected number of EOPatches", str(defaults["expected_patches"])))
    defaults["patch_gpkg"] = ask("EOPatch GeoPackage", defaults["patch_gpkg"])
    defaults["detection_grass_root"] = ask(
        "Detection grass-mask folder (+5 m)", defaults["detection_grass_root"]
    )
    defaults["detection_grass_manifest"] = ask(
        "Detection-mask manifest", defaults["detection_grass_manifest"]
    )
    defaults["stats_grass_root"] = ask(
        "Statistics grass-mask folder (-7.5 m)", defaults["stats_grass_root"]
    )
    defaults["stats_grass_manifest"] = ask(
        "Statistics-mask manifest", defaults["stats_grass_manifest"]
    )
    defaults["workspace_base"] = ask("Workspace/output root", defaults["workspace_base"])

    print()
    print("Year:", defaults["year"])
    print("GeoPackage:", defaults["patch_gpkg"])
    print("Workspace:", defaults["workspace_base"])
    answer = input("Write this configuration? Type YES: ").strip()
    if answer != "YES":
        raise SystemExit("Not saved.")

    OUT.write_text(json.dumps(defaults, indent=2) + "\n", encoding="utf-8")
    print("Saved:", OUT)
    print("Next: python 01_prepare_acquisition.py")


if __name__ == "__main__":
    main()
