#!/usr/bin/env python
import ROOT
import numpy as np
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
#ROOT.gROOT.ProcessLineSync(".x src/PiecewiseLinearPdf.cxx+")
ROOT.gROOT.ProcessLineSync(".x src/NuOscIBDPdf.cxx+")

####################################################################################################
## Constants, common variables
elemNames = ["U235", "Pu239", "Pu241", "U238"]
## The Observable is v_NuE
v_NuE = ROOT.RooRealVar("v_NuE", "Neutrino Energy", 2, 8, unit="MeV")
v_NuE.setBins(500)

v_sin13 = ROOT.RooRealVar("v_sin13", "sin^{2}(2#theta_{13})", 0.093-3*0.008, 0.093+3*0.008)
v_sin14 = ROOT.RooRealVar("v_sin14", "sin^{2}(2#theta_{14})", 0, 1)
v_dm13 = ROOT.RooRealVar("v_dm13", "#Delta^{2}_{13}", 24.4E-4-3*0.6E-4, 24.4E-4+3*0.6E-4, unit="eV")
v_dm14 = ROOT.RooRealVar("v_dm14", "#Delta^{2}_{14}", 0, 5, unit="eV")
v_L = ROOT.RooRealVar("v_L", "L", 0, 10000, unit="meter") ## up to 10km for now.
####################################################################################################

####################################################################################################
## Set default values
v_sin13.setVal(0.093) ## 0.093 +- 0.008
v_dm13.setVal(24.4E-4) ## 24.4 +- 0.6 * 10^-4 eV^2 in normal hierarchy
v_L.setVal(20) ## 20 meter
#v_L.setVal(294) ## 294 meter

v_sin14.setVal(0.1) ## To be measured
v_dm14.setVal(2.0) ## To be measured, in eV^2
####################################################################################################

####################################################################################################
## (Maybe later) Consider the burn-up effect, 
## the fuel composition is a subject to be changed in time.
v_elemFracs = [0.8, 0.1, 0.1] ## U238 is to be set automatically to satisfy sum=1
formula = "1"
for i in range(len(elemNames)-1):
    elemName = elemNames[i]
    elemFrac = v_elemFracs[i]
    v_elemFracs[i] = ROOT.RooRealVar("v_"+elemName, "v_"+elemName, 0, 1)
    formula += ("-v_" + elemName)
print(formula)
v_elemFracs.append(ROOT.RooFormulaVar("v_"+elemNames[-1], "v_"+elemNames[-1], formula, ROOT.RooArgList(*v_elemFracs)))
del(formula)
####################################################################################################

####################################################################################################
## Load the Neutrino flux model, such as Huber-Mueller, incorporating the fuel compositions
fHuber = ROOT.TFile("data/huber.root")
fMueller = ROOT.TFile("data/mueller.root")
grps_HM = []
for elemName in elemNames:
    grp = fHuber.Get("g_"+elemName)
    if grp is None or grp == None:
        grp = fMueller.Get("g_"+elemName)
    grps_HM.append(grp)
####################################################################################################

####################################################################################################
## Load the IBD cross section table
fXsec = ROOT.TFile("data/ibdxsec.root")
grp_Xsec = fXsec.Get("g_LowE")
####################################################################################################

####################################################################################################
## Build the PDF
pdf_NuE = ROOT.NuOscIBDPdf("pdf_NuE", "pdf_NuE", v_NuE, v_L,
                           v_sin13, v_dm13, v_sin14, v_dm14,
                           v_elemFracs[0], v_elemFracs[1], v_elemFracs[2], v_elemFracs[3],
                           grps_HM[0], grps_HM[1], grps_HM[2], grps_HM[3],
                           grp_Xsec)

## Remake binned-pdfs
#v_NuE.setBins(25)
#pdf_BNuE = ROOT.RooBinSamplingPdf("pdf_BNuE", "pdf_NuE", v_NuE, pdf_NuE)

## Clean up temporary objects
del(grps_HM)
del(grp_Xsec)
del(fHuber)
del(fMueller)
del(fXsec)
####################################################################################################

####################################################################################################
## Load the RENO and NEOS data
fRENO = ROOT.TFile("data/reno.root")
hRENO = fRENO.Get("hSpectrum")
dh_RENO = ROOT.RooDataHist("dh_RENO", "dh_RENO", ROOT.RooArgSet(v_NuE), hRENO)

fNEOS = ROOT.TFile("data/neos.root")
hNEOS = fNEOS.Get("hSpectrum")
dh_NEOS = ROOT.RooDataHist("dh_NEOS", "dh_NEOS", ROOT.RooArgSet(v_NuE), hNEOS)
####################################################################################################

v_sin13.setConstant(False)
v_dm13.setConstant(False)
v_sin14.setVal(0)
v_dm14.setVal(0)
v_sin14.setConstant(True)
v_dm14.setConstant(True)
v_L.setVal(294) ## 294 meter
v_L.setConstant(True)
#for v_elemFrac in v_elemFracs:
#    v_elemFrac.setConstant(False)
fitResult_RENO = pdf_NuE.fitTo(dh_RENO, ROOT.RooFit.Save())
#fitResult_NEOS = pdf_NuE.fitTo(dh_NEOS, ROOT.RooFit.Save())

fitResult_RENO.Print()
#fitResult_NEOS.Print()

####################################################################################################
## Test drawing
cNuE = ROOT.TCanvas("cNuE", "cNuE", 500, 500)
frameNuE = v_NuE.frame()
#dh_NEOS.plotOn(frameNuE, ROOT.RooFit.LineColor(ROOT.kBlack), ROOT.RooFit.LineWidth(1))
dh_RENO.plotOn(frameNuE, ROOT.RooFit.LineColor(ROOT.kBlack), ROOT.RooFit.LineWidth(1))

#pdf_NuE.plotOn(frameNuE, ROOT.RooFit.ProjWData(v_NuE, dh_NEOS),
#               ROOT.RooFit.LineColor(ROOT.kBlue+1), ROOT.RooFit.LineWidth(2))
pdf_NuE.plotOn(frameNuE, ROOT.RooFit.ProjWData(v_NuE, dh_RENO),
               ROOT.RooFit.LineColor(ROOT.kBlue+1), ROOT.RooFit.LineWidth(2))
#pdf_BNuE.plotOn(frameNuE, ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.LineWidth(1))

frameNuE.Draw()
####################################################################################################
