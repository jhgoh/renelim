#!/usr/bin/env python
import ROOT
import numpy as np
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
#ROOT.gROOT.ProcessLineSync(".x src/PiecewiseLinearPdf.cxx+")
ROOT.gROOT.ProcessLineSync(".x src/NuOscIBDPdf.cxx+")

####################################################################################################
## Constants, common variables
elemNames = ["U235", "Pu239", "Pu241"]
## The Observable is v_NuE
v_NuE = ROOT.RooRealVar("v_NuE", "Neutrino Energy", 0, 10, unit="MeV")
v_NuE.setBins(200)

v_sin13 = ROOT.RooRealVar("v_sin13", "sin^{2}(2#theta_{13})", 0, 1)
v_sin14 = ROOT.RooRealVar("v_sin14", "sin^{2}(2#theta_{14})", 0, 1)
v_dm13 = ROOT.RooRealVar("v_dm13", "#Delta^{2}_{13}", 0, 5, unit="eV")
v_dm14 = ROOT.RooRealVar("v_dm14", "#Delta^{2}_{14}", 0, 5, unit="eV")
v_L = ROOT.RooRealVar("v_L", "L", 0, 10000, unit="meter") ## up to 10km for now.
####################################################################################################

####################################################################################################
## Set default values
v_sin13.setVal(0.093) ## 0.093 +- 0.008
v_dm13.setVal(24.4E-4) ## 24.4 +- 0.6 * 10^-4 eV^2 in normal hierarchy
v_L.setVal(20) ## 20 meter
#v_L.setVal(294) ## 294 meter
####################################################################################################
## Our parameter of interests
v_sin14.setVal(0.1) ## To be measured
v_dm14.setVal(2.0) ## To be measured, in eV^2
####################################################################################################

####################################################################################################
## (Maybe later) Consider the burn-up effect, 
## the fuel composition is a subject to be changed in time.
v_elemFracs = [1., 0., 0.] ## FIXME: (1,0,0) for the initial version
for i in range(len(v_elemFracs)):
    elemName = elemNames[i]
    elemFrac = v_elemFracs[i]
    v_elemFracs[i] = ROOT.RooRealVar("v_"+elemName, "v_"+elemName, elemFrac, elemFrac)
####################################################################################################

####################################################################################################
## Load the Neutrino flux model, such as Huber-Mueller, incorporating the fuel compositions
fHM = ROOT.TFile("data/hubermueller.root")
grps_HM = []
for elemName in elemNames:
    grp = fHM.Get("g_"+elemName)
    grps_HM.append(grp)
####################################################################################################

####################################################################################################
## Load the IBD cross section table
fXsec = ROOT.TFile("data/ibdxsec.root")
grp_Xsec = fXsec.Get("g_LowE")

pdf_NuE = ROOT.NuOscIBDPdf("pdf_NuE", "pdf_NuE", v_NuE, v_L,
                           v_sin13, v_dm13, v_sin14, v_dm14,
                           v_elemFracs[0], v_elemFracs[1], v_elemFracs[2],
                           grps_HM[0], grps_HM[1], grps_HM[2],
                           grp_Xsec)

####################################################################################################

####################################################################################################
## Remake them as binned-pdfs
v_NuE.setBins(25)
pdf_BNuE = ROOT.RooBinSamplingPdf("pdf_BNuE", "pdf_NuE", v_NuE, pdf_NuE)
####################################################################################################

### Test drawing
cNuE = ROOT.TCanvas("cNuE", "cNuE", 500, 500)
frameNuE = v_NuE.frame()
pdf_NuE.plotOn(frameNuE, ROOT.RooFit.LineColor(ROOT.kBlack), ROOT.RooFit.LineWidth(2))
pdf_BNuE.plotOn(frameNuE, ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.LineWidth(1))
frameNuE.Draw()

