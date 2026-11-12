import sys
MODULES = {'00_setup_project.py': 'setup_project', '00_check_environment.py': 'environment_check', '01_prepare_acquisition.py': 'acquisition_prepare', '02_authenticate.py': 'authenticate', '03_acquire.py': 'acquisition', '04_acquisition_status.py': 'acquisition_status', '09_prepare_processing.py': 'processing_prepare', '10_process_mowing.py': 'detection', '11_finalize.py': 'finalization', '12_analyze.py': 'analysis', '13_mosaic.py': 'export', '14_run_processing.py': 'pipeline', 'prepare_inputs.py': 'inputs', 'parallel_step.py': 'parallel', 'mowing_core.py': 'core', 'project_auth.py': 'auth', 'config.py': 'config'}
MODULES['vector_export.py'] = 'vector_export'
def command(args):
    return [sys.executable, "-m", "moweo." + MODULES[args[0]], *args[1:]]
