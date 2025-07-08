#!/usr/bin/env python
import sys, os
import itertools
import subprocess
import numpy as np

runScriptName = 'run_hypoTest.sbatch'

n_vals = [1000, 10000, 100000]
#n_vals = [1000, 2000, 5000,
#          10000, 20000, 50000,
#          100000]
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

#(ds4hep) jhgoh@lugia:test$ python -i ./run_hypoTest.py -s 1e-4,1e-3,1e-2,1e-1,0.5 -m 2 -n 10000 -o res.root^C

s_str = ','.join([f'{x:g}' for x in s_vals])
s_str2 = s_str.replace(',', '_')
for m, n in itertools.product(m_vals, n_vals):
  suffix = f's_{s_str2}__m_{m}__n_{n}'.replace('.', 'p')
  jobName = f'HypoTest_{suffix}'
  fName = f'results/result__{suffix}.root'
  cmd = ['sbatch', f'--job-name={jobName}',
         f'--export=V_SIN14={s_str},V_DM41={m},V_NSIG={n},V_FNAME={fName}',
         runScriptName]

  if os.path.exists(fName):
    print(f"{fName} already exists. Skip this point...")
  else:
    print(' '.join(cmd))
    res = subprocess.run(cmd, capture_output=True, text=True, check=True)
    
