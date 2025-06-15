#!/usr/bin/env python
## Script to build convoluted response matrix
import numpy as np
import ROOT
import argparse
import sys

parser = argparse.ArgumentParser(description='Convolute two ROOT TH2D response matrices')
parser.add_argument('input_eNu_eTrue', type=str,
                    help='Path to the ROOT file containing P(E_true_e+ | E_nu). '
                         'Format: "filePath" or "filePath:histName"')
parser.add_argument('input_eTrue_eReco', type=str,
                    help='Path to the ROOT file containing P(E_reco | E_true_e+). '
                         'Format: "filePath" or "filePath:histName"')
parser.add_argument('-o', '--output', type=str,
                    help='Path to the output. Format: "filePath" or "filePath:histName". default histName="hResp_Ereco_vs_ENu"')
parser.add_argument('-g', '--gui', action='store_true', default=False,
                    help='Display the combined response matrix')
args = parser.parse_args()

################################################################################
## Get histograms from the input files
################################################################################
def getHistogram(path):
  ## Extract input file and histogram paths
  path = path.split(':', 1)
  fPath = path[0]
  hPath = path[1] if len(path) > 1 else None

  ## Open files
  f = ROOT.TFile(fPath)
  if not f or not f.IsOpen():
    print(f'!!! Cannot open ROOT file: {fPath}')
    sys.exit(1)

  ## Load histograms
  h = None
  if hPath != None:
    h = f.Get(hPath)
    if not h or not isinstance(h, ROOT.TH2):
      print(f'!!! Requested histogram {hPath} does not exist or not a TH2D type')
      sys.exit(2)
  else:
    for key in f.GetListOfKeys():
      obj = key.ReadObj()
      if isinstance(obj, ROOT.TH2):
        h = obj
        print(f'@@@ Using the first-appearing TH2 object: {key}')
        break
    if not h:
      print(f'!!! Cannot find TH2 object in the ROOT file')

  return f, h

f1, h1 = getHistogram(args.input_eNu_eTrue)
f2, h2 = getHistogram(args.input_eTrue_eReco)
################################################################################

################################################################################
## Convert TH2 histogram to numpy 2d array
################################################################################
def convertTH2toNP(h):
  ## Get binning
  nbinsX = h.GetNbinsX()
  nbinsY = h.GetNbinsY()
  xmin, xmax = h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax()
  ymin, ymax = h.GetYaxis().GetXmin(), h.GetYaxis().GetXmax()
  xtitle = h.GetXaxis().GetTitle()
  ytitle = h.GetYaxis().GetTitle()

  hNP = np.zeros((nbinsY, nbinsX), dtype=float)
  for i in range(nbinsX):
    for j in range(nbinsY):
      hNP[j,i] = h.GetBinContent(i+1,j+1)

  return (hNP, (nbinsX, xmin, xmax, xtitle), (nbinsY, ymin, ymax, ytitle))

hNP1, xaxis1, yaxis1 = convertTH2toNP(h1)
hNP2, xaxis2, yaxis2 = convertTH2toNP(h2)
################################################################################

################################################################################
## Check consistency of matrices
################################################################################
## Before proceed, check consistency of two matrices
## x-axis of h2 and y-axis of h1 should be consistent
## True energy of e+ corresponds to x-axis of h2 and y axis of h1
xbinW2 = (xaxis2[2]-xaxis2[1])/xaxis2[0]
ybinW1 = (yaxis1[2]-yaxis1[1])/yaxis1[0]
if not np.isclose(xbinW2, ybinW1):
  print(f'!!! Mismatch in intermediate axis (True E_e+) bin widths')
  sys.exit(3)

## Adjust bin ranges, find common ranges
cmin = max(xaxis2[1], yaxis1[1])
cmax = min(xaxis2[2], yaxis1[2])
if cmin >= cmax:
  print(f'!!! No overapping range found for intermediate axis (True E_e+)')
  sys.exit(3)
## find bin index
startBin1 = int(np.floor((cmin-yaxis1[1])/ybinW1 + 1e-9))
startBin2 = int(np.floor((cmin-xaxis2[1])/xbinW2 + 1e-9))
endBin1 = int(np.ceil((cmax-yaxis1[1])/ybinW1 - 1e-9))
endBin2 = int(np.ceil((cmax-xaxis2[1])/ybinW1 - 1e-9))

hhNP1 = hNP1[startBin1:endBin1, :]
hhNP2 = hNP2[:, startBin2:endBin2]
################################################################################

################################################################################
## Do the convolution
################################################################################
hNP3 = hhNP2 @ hhNP1
xaxis3 = xaxis1
yaxis3 = yaxis2
print(f"Matrix1: {hNP1.shape} -> {hhNP1.shape}")
print(f"Matrix2: {hNP2.shape} -> {hhNP2.shape}")
print(f"Convolved matrix: {hNP3.shape}")
print(f"Axis info: x:{xaxis3}, y:{yaxis3}")
################################################################################

################################################################################
## Prepare output
################################################################################
outPath = args.output.split(':', 1)
outFilePath = outPath[0]
outHistName = outPath[1] if len(outPath) > 1 else "hResp_Ereco_vs_ENu"

h3 = ROOT.TH2D(outHistName, "Response matrix;"+xaxis3[-1]+';'+yaxis3[-1],
               xaxis3[0], xaxis3[1], xaxis3[2],
               yaxis3[0], yaxis3[1], yaxis3[2])
for i in range(xaxis3[0]):
  for j in range(yaxis3[0]):
    h3.SetBinContent(i+1, j+1, hNP3[j,i])
################################################################################

################################################################################
## Draw histogram if requested
################################################################################
if args.gui:
  ROOT.gStyle.SetOptStat(0)
  ROOT.gStyle.SetOptTitle(0)

  c = ROOT.TCanvas("c", "c", 500, 500)
  h3.Draw("COLZ")
  c.Update()

fout = ROOT.TFile(outFilePath, "recreate")
h3.Write()
fout.Close()

