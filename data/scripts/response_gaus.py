#!/usr/bin/env python
## Build a response matrix, E_neutrino vs E_positron
## This is not conventional, but I think it is worth to split
## kinematics of IBD from the detector effect
import numpy as np
import ROOT
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

binW = 0.1
xmin, xmax = 0, 10
ymin, ymax = 0, 10

## Constants
#a, b, c = 0.042, 0.021, 0.010 ## NEOS ~5%/sqrt(E)
a, b, c = 0.06, 0.01, 0.0005 ## RENO? ~5-6%/sqrt(E)
#a, b, c = 0.08, 0.01, 0.0005 ## Daya Bay? ~6-8%/sqrt(E)
#a, b, c = 0.025, 0.01, 0.005 ## JUNO? 3% at 1MeV

def gaus(x, mu, sigma):
  y = np.sqrt(2*np.pi)*np.exp(-0.5*((x-mu)/sigma)**2)
  y /= np.sum(y)
  return y

xvals = np.arange(xmin, xmax, binW)+binW/2
yvals = np.arange(ymin, ymax, binW)+binW/2
sumW = []
for eTrue in xvals:
  sigma = np.sqrt((a**2)*eTrue + (b**2)*(eTrue**2) + (c**2))
  _sumW = gaus(yvals, eTrue, sigma)
  sumW.append(_sumW)

## Copy to the TH2D
h = ROOT.TH2D("hResp_ETrue_vs_EReco", "Detector response matrix;True energy (MeV);Measured energy (MeV)",
              len(xvals), xmin, xmax, len(yvals), ymin, ymax)
for i in range(len(xvals)):
  for j in range(len(yvals)):
    h.SetBinContent(i+1, j+1, sumW[i][j])

h.Draw("COLZ")

f = ROOT.TFile("response_gaus.root", "recreate")
h.Write()
f.Close()

