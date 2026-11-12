"""Compatibility launcher; prefer the moweo command."""
import runpy
if __name__ == "__main__":
    runpy.run_module("moweo.export", run_name="__main__")
