#!/usr/bin/env python
import sys, os
import itertools
import subprocess
import numpy as np

runScriptName = 'run_hypoTest.sbatch'

n_vals = [1000, 2000, 5000,
          10000, 20000, 50000,
          100000]
s_vals = np.concatenate([
           np.arange(1, 10)/10,
           np.arange(1, 10)/100,
           np.arange(1, 10)/1000,
         ])
m_vals = np.concatenate([
           np.arange(1, 10)/100,
           np.arange(1, 10)/10,
           np.arange(1, 10)/1,
         ])

if not os.path.exists('results'):
  os.makedirs('results')

for s, m, n in itertools.product(s_vals, m_vals, n_vals):
  suffix = f's_{s}__m_{m}__n_{n}'.replace('.', 'p')
  jobName = f'HypoTest_{suffix}'
  fName = f'results/result__{suffix}.root'
  cmd = ['sbatch', f'--job-name={jobName}',
         f'--export=V_SIN14={s},V_DM41={m},V_NSIG={s},V_FNAME={fName}',
         runScriptName]

  if os.path.exists(fName):
    print(f"{fName} already exists. Skip this point...")
  else:
    print(' '.join(cmd))
    res = subprocess.run(cmd, capture_output=True, text=True, check=True)
    
