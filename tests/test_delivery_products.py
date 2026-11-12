"""Behavioral checks for delivery postprocessing and event-local quality."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


class DeliveryProducts(unittest.TestCase):
    def test_sieve_harmonization_and_quality_rules(self):
        source = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as tmp:
            project = Path(tmp)
            if os.environ.get('MOWEO_TEST_INSTALLED') != '1':
                shutil.copytree(source / 'moweo', project / 'moweo', ignore=shutil.ignore_patterns('__pycache__'))
            shutil.copy2(source / 'project_config.example.json', project / 'project_config.json')
            code = '''
import numpy as np
from moweo.export import seive_doy_count, harmonize_mowing_count
from moweo.core import observation_density_score, signal_quality_score
mask=np.ones((50,50), dtype=bool)
mask[0,:]=False
layers=np.zeros((6,50,50), dtype=np.int32)
layers[0,15:35,15:35]=100
layers[0,3:5,3:5]=150
class DateRaster:
    descriptions=tuple(f'MOWING_DOY_{i}' for i in range(1,7))
    name='synthetic'
    def read(self, band):
        return layers[band-1]
result=seive_doy_count(DateRaster(),mask)
assert (result[15:35,15:35]==1).all()
assert (result[3:5,3:5]==0).all()
assert (result[~mask]==255).all()
count=np.zeros((50,50),dtype=np.float32)
count[5,5]=3
count[10,10]=3
confidence=np.full((50,50),80,dtype=np.float32)
confidence[5,5]=30
harm=harmonize_mowing_count(count,confidence,mask)
assert harm[5,5]==0
assert harm[10,10]==3
assert (harm[~mask]==255).all()
assert observation_density_score(10,10,4)==100
assert observation_density_score(20,20,3)==65
assert observation_density_score(21,20,3)==30
assert observation_density_score(np.nan,10,5)==0
assert np.isclose(signal_quality_score(1,1,0,True),100)
assert np.isclose(signal_quality_score(1,1,0,False),81.25)
print('Delivery quality and postprocessing rules passed.')
'''
            result = subprocess.run([sys.executable, '-c', code], cwd=project, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == '__main__':
    unittest.main()
