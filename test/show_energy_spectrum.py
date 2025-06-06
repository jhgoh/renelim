#!/usr/bin/env python
import ROOT
import sys
import numpy as np
from array import array

###############################################################################
## Common constants
minNuE = 2.0 ## The Huber Mueller table is available only from 2MeV
maxNuE = 8.0 ## The Huber Mueller table is available only up to 8MeV
###############################################################################

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
#ROOT.gROOT.ProcessLineSync(".x src/PiecewiseLinearPdf.cxx+")
ROOT.gROOT.ProcessLineSync(".x src/NuOscIBDPdf.cxx+")

sys.path.append('python')
from Config import ConfigRENE
config = ConfigRENE('config.yaml')

#print(config.get('reactors'))
#print(config.getReactors('name'))
reaPoses = config.getReactors('position')
reaPoses = np.array(reaPoses)
#print(config.getReactors('elements'))
#print(config.getReactors('power'))

#print(config.get('detectors'))
detPoses = config.getDetectors('position')
detPoses = np.array(detPoses)
#print(config.getDetectors('efficiency'))
baselines = config.getBaselines('RENE0')

###############################################################################
## Draw the reactor & detector layout
###############################################################################
xmin, xmax = np.min(reaPoses), np.max(reaPoses)
gReaLayout = ROOT.TGraph(len(reaPoses), array('d', reaPoses[:,0]), array('d', reaPoses[:,1]))
gDetLayout = ROOT.TGraph(len(detPoses), array('d', detPoses[:,0]), array('d', detPoses[:,1]))

cLayout = ROOT.TCanvas("cLayout", "cLayout", 500, 500)
hFrameLayout = ROOT.TH2F("hFrameLayout", ";Relative position (m);Relative position (m)",
                         100, xmin-100, xmax+100, 100, xmin-100, xmax+100)
hFrameLayout.Draw()

gReaLayout.SetMarkerColor(ROOT.kRed)
gReaLayout.SetMarkerStyle(ROOT.kCircle)
gReaLayout.SetMarkerSize(2)

gDetLayout.SetMarkerColor(ROOT.kBlue)
gDetLayout.SetMarkerStyle(ROOT.kFullCircle)
gDetLayout.SetMarkerSize(1)

gReaLayout.SetEditable(False)
gDetLayout.SetEditable(False)
gReaLayout.Draw("P")
gDetLayout.Draw("P")
del(xmin)
del(xmax)

cLayout.Update()
###############################################################################

###############################################################################
## Draw Neutrino energy spectrums
###############################################################################
## (Maybe later) Consider the burn-up effect, 
## the fuel composition is a subject to be changed in time.
## We choose the 6th core only (1st one in the configuration)
elemNames = config.getReactors('elements')[0]['name']
v_elemFracs = config.getReactors('elements')[0]['fraction']
formula = "1"
for i in range(len(elemNames)-1):
    elemName = elemNames[i]
    elemFrac = v_elemFracs[i]
    v_elemFracs[i] = ROOT.RooRealVar("v_"+elemName, "v_"+elemName, 0, 1)
    formula += ("-v_" + elemName)
#print(formula)
v_elemFracs[-1] = ROOT.RooFormulaVar("v_"+elemNames[-1], "v_"+elemNames[-1], formula, ROOT.RooArgList(*v_elemFracs[:-1]))
del(formula)

## Load the Neutrino flux model, such as Huber-Mueller, incorporating the fuel compositions
grps_HM = []
for elemName in elemNames:
    fPath, gPath = config.get(f'physics.isotope_flux.{elemName}').split(':')
    grp = ROOT.TFile(fPath).Get(gPath)
    grps_HM.append(grp)

cHM = ROOT.TCanvas("cHM", "Huber Mueller Energy spectrum", 500, 500)
hFrameHM = ROOT.TH1D("hFrameHM", ";Neutrino energy (MeV);Arbitrary Unit", 100, minNuE, maxNuE)
legHM = ROOT.TLegend(0.60, 0.65, 0.85, 0.85)
legHM.SetBorderSize(0)
legHM.SetFillStyle(0)
maxY = 0
colors = [ROOT.kRed+1, ROOT.kBlue+2, ROOT.kGreen+3, ROOT.kViolet+5, ROOT.kOrange+7]
for grp, elemName, color in zip(grps_HM, elemNames, colors):
  grp.SetLineColor(color)
  maxY = max(maxY, max(grp.GetY()))
  grp.SetEditable(False)
  legHM.AddEntry(grp, elemName)

hFrameHM.SetMaximum(1.1*maxY)
hFrameHM.Draw()
for grp in grps_HM:
  grp.Draw("L")
legHM.Draw()
del(maxY)

cHM.Update()
################################################################################

###############################################################################
## Draw IBD cross section
###############################################################################
fPath, gPath = config.get('physics.ibd_xsec').split(':')
grp_Xsec = ROOT.TFile(fPath).Get(gPath)
idx = np.searchsorted(grp_Xsec.GetX(), maxNuE)
maxY = grp_Xsec.GetY()[idx]

cIBDXsec = ROOT.TCanvas("cIBDXsec", "IBD cross section", 500, 500)
hFrameIBDXsec = ROOT.TH1D("hFrameIBDXsec", ";Neutrino energy (MeV);Arbitrary Unit", 100, minNuE, maxNuE)
hFrameIBDXsec.SetMaximum(1.1*maxY)
hFrameIBDXsec.Draw()
grp_Xsec.Draw("Lsame")

cIBDXsec.Update()
################################################################################

################################################################################
## build the pdf
v_NuE = ROOT.RooRealVar("v_NuE", "Neutrino Energy", minNuE, maxNuE, unit="MeV")
v_NuE.setBins(500)

v_sin13 = ROOT.RooRealVar("v_sin13", "sin^{2}(#theta_{13})", 0, 1)
v_sin14 = ROOT.RooRealVar("v_sin14", "sin^{2}(#theta_{14})", 0, 1)
v_dm31 = ROOT.RooRealVar("v_dm31", "#Delta^{2}_{31}", 0, 1, unit="eV^{2}")
v_dm41 = ROOT.RooRealVar("v_dm41", "#Delta^{2}_{41}", 0, 5, unit="eV^{2}")

v_sin13.setVal(config.get("physics.oscillation.sin13"))
v_dm31.setVal(config.get("physics.oscillation.dm31"))

v_sin14.setVal(0.0) ## No oscillation
v_dm41.setVal(1.0) ## To be measured, in eV^2

v_L = ROOT.RooRealVar("v_L", "L", 0, 10000, unit="meter") ## up to 10km for now.
v_L.setVal(baselines[0])

v_elemFracs[0].setVal(1)
v_elemFracs[1].setVal(0)
v_elemFracs[2].setVal(0)
pdf_NuE = ROOT.NuOscIBDPdf("pdf_NuE", "pdf_NuE", v_NuE, v_L,
                           v_sin13, v_dm31, v_sin14, v_dm41,
                           v_elemFracs[0], v_elemFracs[1], v_elemFracs[2], v_elemFracs[3],
                           grps_HM[0], grps_HM[1], grps_HM[2], grps_HM[3],
                           grp_Xsec)

cNuE = ROOT.TCanvas("cNuE", "cNuE", 500, 500)
legNuE = ROOT.TLegend(0.55, 0.6, 0.85, 0.85)
legNuE.SetBorderSize(0)
legNuE.SetFillStyle(0)
frameNuE = v_NuE.frame()

v_sin13.setVal(0.0)
v_sin14.setVal(0.0)
pdf_NuE.plotOn(frameNuE, ROOT.RooFit.Name("pdf_NoOsc"),
               ROOT.RooFit.LineColor(ROOT.kBlue+2), ROOT.RooFit.LineWidth(2))
legNuE.AddEntry(frameNuE.findObject("pdf_NoOsc"), "No Osc.")

v_sin13.setVal(config.get("physics.oscillation.sin13"))
pdf_NuE.plotOn(frameNuE, ROOT.RooFit.Name("pdf_noNu4"),
               ROOT.RooFit.LineColor(ROOT.kRed+1), ROOT.RooFit.LineWidth(2))
legNuE.AddEntry(frameNuE.findObject("pdf_noNu4"), "with #nu_{3}")

v_sin14.setVal(0.2)
v_dm41.setVal(1.0)
pdf_NuE.plotOn(frameNuE, ROOT.RooFit.Name("pdf_A"), 
               ROOT.RooFit.LineColor(ROOT.kGreen+3), ROOT.RooFit.LineWidth(2))
legNuE.AddEntry(frameNuE.findObject("pdf_A"),
                   "sin_{14}=%.2f, #Delta m_{41}=%.2f" % (v_sin14.getVal(), v_dm41.getVal()))

v_sin14.setVal(0.2)
v_dm41.setVal(2.0)
pdf_NuE.plotOn(frameNuE, ROOT.RooFit.Name("pdf_B"),
               ROOT.RooFit.LineColor(ROOT.kOrange+5), ROOT.RooFit.LineWidth(2))
legNuE.AddEntry(frameNuE.findObject("pdf_B"),
                   "sin_{14}=%.2f, #Delta m_{41}=%.2f" % (v_sin14.getVal(), v_dm41.getVal()))

v_sin14.setVal(0.2)
v_dm41.setVal(5.0)
pdf_NuE.plotOn(frameNuE, ROOT.RooFit.Name("pdf_C"),
               ROOT.RooFit.LineColor(ROOT.kViolet+1), ROOT.RooFit.LineWidth(2))
legNuE.AddEntry(frameNuE.findObject("pdf_C"),
                   "sin_{14}=%.2f, #Delta m_{41}=%.2f" % (v_sin14.getVal(), v_dm41.getVal()))

frameNuE.Draw()
legNuE.Draw()

cNuE.Update()

################################################################################

input("Press Enter to exit...")
