#!/usr/bin/env python
import sys, os
import itertools
import subprocess
import numpy as np

runScriptName = 'run_chi2.sbatch'
nToys = 10000

n_vals = [1000, 2000, 5000, 10000]
s_vals = np.concatenate([
#            np.logspace(-4, 0, 100),
           np.arange(1, 10)/10,
           np.arange(1, 10)/100,
           np.arange(1, 10)/1000,
         ])
m_vals = np.concatenate([
#            np.logspace(-2, 1, 100),
           np.arange(1, 10)/100,
           np.arange(1, 10)/10,
           np.arange(1, 10)/1,
         ])

if not os.path.exists('results'):
  os.makedirs('results')

s_str = ','.join([f'{x:g}' for x in s_vals])
s_str2 = s_str.replace(',', '_')
for i, (m, n) in enumerate(itertools.product(m_vals, n_vals)):
  suffix = f'iter_{i}__dm41_{m}__nSignal_{n}'.replace('.', 'p')
  jobName = f'Chi2_{suffix}'
  fName = f'results/result__{suffix}.root'
  cmd = ['sbatch', f'--job-name={jobName}',
         f'--export=V_SIN14="{s_str}",V_DM41={m},V_NSIG={n},V_FNAME={fName},V_NTOYS={nToys}',
         runScriptName]

  if os.path.exists(fName):
    print(f"{fName} already exists. Skip this point...")
  else:
    print(' '.join(cmd))
    res = subprocess.run(cmd, capture_output=True, text=True, check=True)
    
