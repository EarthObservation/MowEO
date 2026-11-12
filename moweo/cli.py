"""Single entry point; importing the CLI does not require a configured project."""
import argparse
import os
from pathlib import Path
import runpy
import sys

STAGES = {
    'init-aoi': 'area', 'masks': 'masks',
    'init': 'inputs', 'setup': 'setup_project', 'check': 'environment_check',
    'prepare': 'acquisition_prepare', 'login': 'authenticate',
    'acquire': 'acquisition', 'status': 'acquisition_status',
    'process': 'pipeline', 'detect': 'detection', 'finalize': 'finalization',
    'analyze': 'analysis', 'export': 'export',
    'validate': 'processing_prepare',
    'geoexport': 'vector_export',
}


def main():
    p = argparse.ArgumentParser(description='MowEO - Automated Grassland Mowing Detection Using Sentinel-2 and openEO')
    from . import __version__
    p.add_argument('--version', action='version', version='MowEO ' + __version__)
    p.add_argument('--project', type=Path, default=Path.cwd(), help='Project directory with local configuration and inputs.')
    p.add_argument('stage', choices=[*STAGES, 'run'])
    args, rest = p.parse_known_args()
    args.project = args.project.resolve()
    args.project.mkdir(parents=True, exist_ok=True)
    os.chdir(args.project)
    os.environ['MOWEO_PROJECT_DIR'] = str(args.project.resolve())
    if args.stage == 'run':
        r = argparse.ArgumentParser(description='Download annual cubes, then run the complete local pipeline.')
        r.add_argument('--run', action='store_true', help='Required to submit cloud jobs; otherwise show instructions only.')
        r.add_argument('--parallel', type=int, default=2)
        r.add_argument('--workers', type=int, default=4)
        options = r.parse_args(rest)
        if not options.run:
            print('Configure with moweo init or moweo setup first. Run moweo run --run to acquire and process.')
            return
        stages = [('acquisition_prepare', []), ('acquisition', ['--run', '--parallel', str(options.parallel)]),
                  ('pipeline', ['--workers', str(options.workers)])]
    else:
        stages = [(STAGES[args.stage], rest)]
    from .provenance import record
    for module, options in stages:
        if module not in ('inputs', 'setup_project', 'area'):
            record(module)
        sys.argv = ['moweo ' + module, *options]
        runpy.run_module('moweo.' + module, run_name='__main__')


if __name__ == '__main__':
    main()
