#!/usr/bin/env python
## Build a response matrix, E_neutrino vs E_positron
## This is not conventional, but I think it is worth to split
## kinematics of IBD from the detector effect
import argparse
import os

parser = argparse.ArgumentParser(description="Generate a gaussian detector response matrix")
parser.add_argument('--a', type=float, default=0.06,
                    help='Parameter "a" for energy resolution: sigma=sqrt(a^2*E + b^2*E^2 + c^2')
parser.add_argument('--b', type=float, default=0.01,
                    help='Parameter "b" for energy resolution: sigma=sqrt(a^2*E + b^2*E^2 + c^2')
parser.add_argument('--c', type=float, default=0.0005,
                    help='Parameter "c" for energy resolution: sigma=sqrt(a^2*E + b^2*E^2 + c^2')
parser.add_argument('--xmin', type=float, default=0.0,
                    help='Minimum energy for both true and reconstructed energy (MeV)')
parser.add_argument('--xmax', type=float, default=10.0,
                    help='Maximum energy for both true and reconstructed energy (MeV)')
parser.add_argument('--dx', type=float, default=0.1,
                    help='Bin width for both true and reconstructed energy (MeV)')
parser.add_argument('-o', '--output', type=str, default='response_gaus.root',
                    help='Output ROOT file name for response matrix')
parser.add_argument('-g', '--gui', action='store_true',
                    help='Display the response matrix using ROOT TCanvas')
args = parser.parse_args()

import numpy as np
from scipy.stats import norm
import ROOT

## Set parameters
a, b, c = args.a, args.b, args.c

xmin, xmax = args.xmin, args.xmax
dx = args.dx
nbins = int(np.ceil((xmax-xmin)/dx))
if not np.isclose(xmax, xmin + nbins*dx):
  print(f'@@@ Warning: Actual maximum energy differs from estimated range from the input dx.')
  print(f'@@@          The x-axis range is adjusted: ({xmax}) -> ({xmin+nbins*dx})')
  xmax = xmin + nbins*dx
ymin, ymax = xmin, xmax

binEdges = np.linspace(xmin, xmax, nbins+1)
binCenters = (binEdges[:-1]+binEdges[1:])/2

def get_sigma(eTrue, a, b, c):
  s2 = (a**2 * eTrue) + (c**2 * eTrue**2) + (b**2)
  return np.sqrt(s2)

## Build the response matrix
respMatrix = np.zeros((nbins, nbins), dtype=float)
for j in range(nbins):
  eTrue = binCenters[j]
  sigma = get_sigma(eTrue, a, b, c)

  probs = norm.cdf(binEdges[1:], loc=eTrue, scale=sigma) - \
          norm.cdf(binEdges[:-1], loc=eTrue, scale=sigma)
  respMatrix[:, j] = probs

print(respMatrix)

## Normalize the matrix
sums = np.sum(respMatrix, axis=0, keepdims=True)
respMatrix /= sums

## Copy to the TH2D
h = ROOT.TH2D("hResp_ETrue_vs_EReco", "Detector response matrix;True energy (MeV);Measured energy (MeV)",
              nbins, xmin, xmax, nbins, ymin, ymax)
for j in range(nbins):
  for i in range(nbins):
    h.SetBinContent(j+1, i+1, respMatrix[i, j])

if args.gui:
  ROOT.gStyle.SetOptStat(0)
  ROOT.gStyle.SetOptTitle(0)
  c = ROOT.TCanvas("c", "c", 500, 500)
  h.Draw("COLZ")
  c.Update()
  input('Press return to exit')

f = ROOT.TFile(args.output, "recreate")
h.Write()
f.Close()


