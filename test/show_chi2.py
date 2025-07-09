#!/usr/bin/env python
import sys
from glob import glob
from tqdm import tqdm

import numpy as np
import pandas as pd

import ROOT

pars = []
for fName0 in tqdm(glob('results/result_*__n_1000.root')):
    dName, fName = fName0.rsplit('/', 1)

    par = fName[len('result_'):-len('.root')].strip('_')
    par = [x.split('_') for x in par.split('__')]
    par = [[x[0], float(x[1].replace('p','.'))] for x in par if len(x) == 2]
    par = dict(par)

    f = ROOT.TFile(fName0)
    for key in f.GetListOfKeys():
      obj = f.Get(key.GetName())
      if not obj.IsA().InheritsFrom("TGraph"): continue

      nPoint = obj.GetN()
      xvals = np.array([obj.GetX()[i] for i in range(nPoint)])
      yvals = np.array([obj.GetY()[i] for i in range(nPoint)])

      idxs = np.argsort(xvals)
      xvals = xvals[idxs]
      yvals = yvals[idxs]

      par['fName'] = fName0
      for xval, yval in zip(xvals, yvals):
        par['sin14'] = xval
        par['nll'] = yval
      
        pars.append(dict(par))

df = pd.DataFrame(pars)

nll_min = df['nll'].min()
df['dnll'] = df['nll'] - nll_min

g2 = ROOT.TGraph2D()
for i, row in df.iterrows():
    g2.SetPoint(i, row['sin14'], row['dm41'], row['dnll'])
g2.SetNpx(500)
g2.SetNpy(500)
ROOT.gStyle.SetNumberContours(256)
#h2 = g2.GetHistogram()
#h2.SetTitle("Delta NLL scan;sin^{2}(2#theta_{14});#Delta m^{2}_{41} (eV);#Delta NLL")
#h2.SetContour(3)
#h2.SetContourLevel(0, 4.61)  # 90% CL
#h2.SetContourLevel(1, 11.83) # 3sigma
#h2.SetContourLevel(2, 28.74) # 5sigma
#h2.SetMaximum(30)

c = ROOT.TCanvas("c", "cNLL", 500, 500)
c.SetLogx()
c.SetLogy()
g2.Draw("CONT4 Z")
#h2.Draw("CONT Z LIST")
#h2.Draw("CONT5 LIST")
c.Update()
