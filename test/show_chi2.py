#!/usr/bin/env python
import sys
from glob import glob
from tqdm import tqdm

import numpy as np
import pandas as pd

import ROOT

import argparse

parser = argparse.ArgumentParser(description='Run the chi2 tests for sin14 and dm41')
parser.add_argument('-o', '--output', type=str, required=True, help='output ROOT file name to store results')
args = parser.parse_args()

chain = ROOT.TChain("limit")
chain.Add(args.output)
#chain.Add("results/result_*_nSignal_1000.root")
#chain.Draw("nll:sin14", "", "COLZ")

branches = [bn.GetName() for bn in chain.GetListOfBranches()]
pars = {name: [] for name in branches}
n_entries = chain.GetEntries()
for i in range(n_entries):
    chain.GetEntry(i)
    for name in branches:
        pars[name].append(getattr(chain, name))
for name in pars:
    pars[name] = np.array(pars[name])

df = pd.DataFrame(pars)
print(df)

nll_min = df['nll'].min()
df['dnll'] = df['nll'] - nll_min
df = df[df['dnll'] < 30]

g2 = ROOT.TGraph2D()
for i, row in df.iterrows():
    g2.SetPoint(i, row['sin14'], row['dm41'], row['dnll'])
g2.SetNpx(512)
g2.SetNpy(512)
#ROOT.gStyle.SetNumberContours(256)
h2 = g2.GetHistogram()
h2.SetTitle("Delta NLL scan;sin^{2}(2#theta_{14});#Delta m^{2}_{41} (eV);#Delta NLL")
h2.SetContour(3)
h2.SetContourLevel(0, 4.61)  # 90% CL
h2.SetContourLevel(1, 11.83) # 3sigma
h2.SetContourLevel(2, 28.74) # 5sigma
h2.SetMaximum(30)

c = ROOT.TCanvas("c", "cNLL", 500, 500)
c.SetLogx()
c.SetLogy()
#g2.Draw("CONT4 Z")
#h2.Draw("CONT Z LIST")
h2.Draw("CONT5 Z LIST")
c.Update()
