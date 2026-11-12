"""Compatibility launcher; prefer the moweo command."""
import runpy
if __name__ == "__main__":
    runpy.run_module("moweo.acquisition_status", run_name="__main__")
