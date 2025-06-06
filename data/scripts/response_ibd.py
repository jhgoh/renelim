#!/usr/bin/env python
## Build a response matrix, E_neutrino vs E_positron
## This is not conventional, but I think it is worth to split
## kinematics of IBD from the detector effect
import numpy as np
import ROOT
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

nbins = 80 
xmin, xmax = 0, 8
n = 1000000
m = 1000

## Constants
qIBD = 1.806 ## Q-value of IBD
mE = 0.511   ## Electron mass
mN = 939.6   ## Neutron mass
mP = 938.3   ## Proton mass
lam = -1.27  ## g_A/g_V

def dSdE(eAe, eNu): ## neutrino energy, positron energy
  ## Mask to select kinematically allowed regions
  mask = (eNu >= qIBD) & (eAe >= mE) & (eAe <= (eNu - qIBD + mE))

  ## skip xsec calculation if none of samples meet criteria
  xsec = np.zeros(eNu.shape)
  if not np.any(mask): return xsec

  ## Momentum of positrons
  pAe = np.sqrt(eAe**2 - mE**2)
  kAe = eAe - mE

  ## Constants in diff xsec.
  fSq = 1.0**2
  gSq = lam**2
  coefDelta = (mN**2 - mP**2 - mE**2)/2
  coefEpsil = (mN**2 - mP**2 + mE**2)/2

  ## Compute terms
  term1 = (fSq + gSq)*kAe*eAe
  term2 = (fSq - gSq)*coefDelta/2*mE**2 / eNu
  term3 = (fSq + gSq)*(coefDelta/eNu*kAe + coefEpsil/eNu*eAe)

  ## Fill into the results, but skip for the forbidden regions
  xsec[mask] = (term1 + term2 + term3)[mask]
  xsec[~mask] = 0

  return xsec

## Do the numerical integration, for each bins.
sumW = np.zeros([nbins, nbins])
import tqdm
for _ in tqdm.tqdm(range(m)):
  xs = np.random.uniform(xmin, xmax, n)
  ys = np.random.uniform(xmin, xmax, n)
  xsec = dSdE(xs, ys)/n

  _sumW, _, _ = np.histogram2d(xs, ys, weights=xsec,
                              bins=[nbins, nbins], range=[[xmin, xmax], [xmin, xmax]])
  sumW += _sumW
## Normalize histogram
sums = np.sum(sumW, axis=0, keepdims=True)
sumW = np.where(sums > 0, sumW / sums, 0)

## Copy to the TH2D
h = ROOT.TH2D("hResp_ENu_vs_EAe", "Response matrix of Kinematics;Neutrino energy (MeV);Positron energy (MeV)",
              nbins, xmin, xmax, nbins, xmin, xmax)
for i in range(nbins):
  for j in range(nbins):
    h.SetBinContent(i+1, j+1, sumW[j,i])
#h.Scale(1/m)

del(xs)
del(ys)
del(xsec)

h.Draw("COLZ")

f = ROOT.TFile("response_ibd.root", "recreate")
h.Write()
f.Close()

