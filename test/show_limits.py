#!/usr/bin/env python
import sys
from glob import glob
from tqdm import tqdm

import numpy as np
import pandas as pd

import ROOT

pars = []
for fName0 in tqdm(glob('results/result_*.root')):
    dName, fName = fName0.rsplit('/', 1)

    par = fName[len('result_'):-len('.root')].strip('_')
    par = [x.rsplit('_',1) for x in par.split('__')]
    par = [[x[0], float(x[1].replace('p','.'))] for x in par]
    par = dict(par)

    f = ROOT.TFile(fName0)
    r = f.Get("result")
    cls = r.CLs()
    sig = r.Significance()

    par['fName'] = fName0
    par['CLs'] = cls
    par['Signif'] = sig
    pars.append(par)

df = pd.DataFrame(pars)

print(df)
