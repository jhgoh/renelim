#!/usr/bin/env python
import ROOT
import numpy as np
import sys

sys.path.append("python")
from RENElib import SterileNuOscChi2

################################################################################
## Load the RENO and NEOS data
fRENO = ROOT.TFile("data/reno.root")
hRENO = fRENO.Get("hSpectrum")
covRENO = fRENO.Get("hCovariance")

fNEOS = ROOT.TFile("data/neos.root")
hNEOS = fNEOS.Get("hSpectrum")
covNEOS = fNEOS.Get("hCovariance")
################################################################################

computeChi2 = SterileNuOscChi2(419, 23.7, hRENO, hNEOS, covRENO, covNEOS)

sin14vals = np.power(10.0, np.arange(-4, 0, 0.01))
dm14vals = np.power(10, np.arange(-1, 2, 0.01))
chi2vals = np.zeros([len(sin14vals), len(dm14vals)])
for i, sin14 in enumerate(sin14vals):
    for j, dm14 in enumerate(dm14vals):
        chi2 = computeChi2.chi2(2.2E-2, 2.45E-3, sin14, dm14)
        chi2vals[i, j] = chi2

