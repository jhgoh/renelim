#!/usr/bin/env python
import numpy as np
#from ROOT import TGraph
from ROOT import TH1
from ROOT import TCanvas

class SterileNuOscChi2:
    def __init__(self, baseline1, baseline2, spect1, spect2, cov1, cov2):
        ## Initialize member variables
        self.isValid = False
        self.xbins = None
        self.baseline1 = baseline1
        self.baseline2 = baseline2
        self.spect1 = None
        self.spect2 = None
        self.cov1 = None
        self.cov2 = None

        ## Read input spectrums
        ## They can have different x-ranges. Therefore we have to find common range
        ## We assume spect1 and spect2 are ROOT.TH1 type
        nx1, nx2 = spect1.GetNbinsX(), spect2.GetNbinsX()
        xbins1 = np.array([spect1.GetXaxis().GetBinLowEdge(i+1) for i in range(nx1)]
                         +[spect1.GetXaxis().GetBinUpEdge(nx1)])
        xbins2 = np.array([spect2.GetXaxis().GetBinLowEdge(i+1) for i in range(nx2)]
                         +[spect2.GetXaxis().GetBinUpEdge(nx2)])
        self.xbins = np.intersect1d(xbins1, xbins2)
        nx = len(self.xbins)-1

        self.spect1 = np.zeros(nx)
        self.spect2 = np.zeros(nx)

        for i in range(nx):
            x = (self.xbins[i+1] + self.xbins[i])/2

            bin1 = spect1.GetXaxis().FindBin(x)
            bin2 = spect2.GetXaxis().FindBin(x)

            self.spect1[i] = spect1.GetBinContent(bin1)
            self.spect2[i] = spect2.GetBinContent(bin2)
        ## Normalize the rates
        area1 = (self.spect1 * (self.xbins[1:]-self.xbins[:-1])).sum()
        area2 = (self.spect2 * (self.xbins[1:]-self.xbins[:-1])).sum()
        self.spect1 /= area1
        self.spect2 /= area2

        ## Read the covariance matrices, build covariance matrices with common binning
        ## By definition, covariance matrix has to be symmetric, therefore rectangular
        ## We assume cov1 and cov2 are ROOT.TH2 type
        ## FIXME: How can we treat overflow bins?
        self.cov1 = np.zeros([nx, nx])
        self.cov2 = np.zeros([nx, nx])

        for i in range(nx):
            x = (self.xbins[i+1]+self.xbins[i])/2
            ibin1 = cov1.GetXaxis().FindBin(x)
            ibin2 = cov2.GetXaxis().FindBin(x)

            for j in range(nx):
                y = (self.xbins[j+1]+self.xbins[j])/2
                jbin1 = cov1.GetXaxis().FindBin(y)
                jbin2 = cov2.GetXaxis().FindBin(y)

                self.cov1[i,j] = cov1.GetBinContent(ibin1, jbin1)
                self.cov2[i,j] = cov2.GetBinContent(ibin2, jbin2)
        ## Normalize the covariance matrices
        self.cov1 /= area1
        self.cov2 /= area2

        ## Validate the object
        self.isValid = True

    def chi2(self, sin13, dm13, sin14, dm14, alpha=1):
        xx = (self.xbins[1:] + self.xbins[:-1])/2

        prob1 = self.survivalProb(self.baseline1, xx, sin13, dm13, sin14, dm14)
        prob2 = self.survivalProb(self.baseline2, xx, sin13, dm13, sin14, dm14)

        relProb = prob1/prob2
        vec = self.spect1 - alpha * relProb * self.spect2
        vec = np.matrix(vec)

        scaleMat = np.multiply.outer(relProb, relProb)
        cov = self.cov1 + alpha * alpha * scaleMat * self.cov2
        cov = np.matrix(cov)

        return np.array(vec * cov * vec.T).flatten()[0]

    def survivalProb(self, baseline, energy, sin13, dm13, sin14, dm14):
        ## [baseline] = [meter], [energy] = [MeV]
        ## or [baseline] = [km], [energy] = [GeV]
        coef = 1.267 ## 1/hbar/c/4 with hbarc = 197.3 MeV*fm

        ## Note: energy can be numpy array
        ##       baseline must be a scalar
        ##       or should be in the same shape of energy array

        sinD13 = np.sin(coef*dm13*baseline/energy)
        sinD14 = np.sin(coef*dm14*baseline/energy)

        prob13 = 1 - sin13 * sinD13 * sinD13
        prob14 = prob13 - sin14 * sinD14 * sinD14

        return prob14
