#!/usr/bin/env python3
"""Run TIDE prediction on glioma expression data."""
import os
PROJ = os.getcwd()
import pandas as pd
try:
    import tidepy
except ImportError:
    import subprocess, sys
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'tidepy'])
    import tidepy

expr = pd.read_csv(os.path.join(PROJ, 'results', 'tide_input.txt'), sep='\t', index_col=0)
result = tidepy.predict_tide(expr, cancer_type='GBM')
result.to_csv(os.path.join(PROJ, "results", "tide_result.csv"))
print('TIDE results saved successfully')
print(result.head())
