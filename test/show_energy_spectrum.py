#!/usr/bin/env python
import ROOT
import sys
import numpy as np
from array import array

###############################################################################
## Common constants
minENu = 2.0 ## The Huber Mueller table is available only from 2MeV
maxENu = 8.0 ## The Huber Mueller table is available only up to 8MeV
###############################################################################

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
#ROOT.gROOT.ProcessLineSync(".x src/PiecewiseLinearPdf.cxx+")
ROOT.gROOT.ProcessLineSync(".x src/NuOscIBDPdf.cxx+")

ws = ROOT.RooWorkspace("ws", "ws")

sys.path.append('python')
from Config import *
config = ConfigRENE('config.yaml')

#print(config.get('reactors'))
#print(config.getReactors('name'))
reaPoses = config.getReactors('position')
reaPoses = np.array(reaPoses)
#print(config.getReactors('elements'))
#print(config.getReactors('power'))

#print(config.get('detectors'))
detNames = config.getDetectors('name')
detPoses = np.array(config.getDetectors('position'))
detEffs = config.getDetectors('efficiency')
responses = config.getDetectors('response')
allBaselines = []
for detName in detNames:
  allBaselines.append(config.getBaselines(detName))

###############################################################################
## Draw the reactor & detector layout
## We take only one detector for this version.
## We also guess range of neutrino energy distribution from the response matrix
###############################################################################
detName = detNames[0] ## Maybe not used
defEff = detEffs[0] ## Maybe not used
baselines = allBaselines[0]

xmin, xmax = np.min(reaPoses), np.max(reaPoses)
ROOT.gROOT.cd()
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
## Detector informations
## We take only one detector for this version.
## We also guess range of neutrino energy distribution from the response matrix
###############################################################################
detName = detNames[0] ## Maybe not used
defEff = detEffs[0] ## Maybe not used
baselines = allBaselines[0]

_fResp, _hResp = getFileAndObj(responses[0])
ROOT.gROOT.cd()

nbinsENu = _hResp.GetNbinsX()
minENu = _hResp.GetXaxis().GetXmin()
maxENu = _hResp.GetXaxis().GetXmax()

nbinsEReco = _hResp.GetNbinsY()
minEReco = _hResp.GetYaxis().GetXmin()
maxEReco = _hResp.GetYaxis().GetXmax()

reactorIdx = np.argmin(baselines)
baseline = baselines[reactorIdx]
###############################################################################

###############################################################################
## Draw Neutrino energy spectrums
###############################################################################
## (Maybe later) Consider the burn-up effect, 
## the fuel composition is a subject to be changed in time.
## We choose the 6th core only (1st one in the configuration)
elemNames = config.getReactors('elements')[reactorIdx]['name']
v_elemFracs = config.getReactors('elements')[reactorIdx]['fraction']
formula = "1"
formulaVars = []
for i in range(len(elemNames)-1):
    elemName = elemNames[i]
    elemFrac = v_elemFracs[i]
    ws.factory(f'v_{elemName}[{elemFrac}, 0, 1]')
    v_elemFracs[i] = ws.var(f'v_{elemName}')
    v_elemFracs[i].setConstant(True)
    formula += f"-@{i}"
    formulaVars.append(f'v_{elemName}')
formulaVars = ','.join(formulaVars)
ws.factory(f'EXPR::v_{elemNames[-1]}("{formula}", {{ {formulaVars} }})')
del(formula)
del(formulaVars)
for i in range(len(elemNames)-1):
  v_elemFracs[i] = ws.var(f'v_{elemNames[i]}')
v_elemFracs[-1] = ws.function(f'v_{elemNames[-1]}')

## Load the Neutrino flux model, such as Huber-Mueller, incorporating the fuel compositions
grps_HM = []
for elemName in elemNames:
    _, grp = getFileAndObj(config.get(f'physics.isotope_flux.{elemName}'))
    ROOT.gROOT.cd()
    grps_HM.append(grp)

ROOT.gROOT.cd()
cHM = ROOT.TCanvas("cHM", "Huber Mueller Energy spectrum", 500, 500)
hFrameHM = ROOT.TH1D("hFrameHM", ";Neutrino energy (MeV);Arbitrary Unit", 100, minENu, maxENu)
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
_, grp_Xsec = getFileAndObj(config.get('physics.ibd_xsec'))
ROOT.gROOT.cd()
idx = np.searchsorted(grp_Xsec.GetX(), maxENu)
maxY = grp_Xsec.GetY()[idx]

ROOT.gROOT.cd()
cIBDXsec = ROOT.TCanvas("cIBDXsec", "IBD cross section", 500, 500)
hFrameIBDXsec = ROOT.TH1D("hFrameIBDXsec", ";Neutrino energy (MeV);Arbitrary Unit", 100, minENu, maxENu)
hFrameIBDXsec.SetMaximum(1.1*maxY)
hFrameIBDXsec.Draw()
grp_Xsec.Draw("Lsame")

cIBDXsec.Update()
################################################################################

################################################################################
## build the pdf
v_ENu = ROOT.RooRealVar("v_ENu", "Neutrino Energy", minENu, maxENu, unit="MeV")
v_ENu.setBins(1000)
getattr(ws, 'import')(v_ENu)
v_ENu = ws.var('v_ENu')
v_ENu.setConstant(True)

v_sin13 = ROOT.RooRealVar("v_sin13", "sin^{2}(#theta_{13})", 0, 1)
v_sin14 = ROOT.RooRealVar("v_sin14", "sin^{2}(#theta_{14})", 0, 1)
v_dm31 = ROOT.RooRealVar("v_dm31", "#Delta^{2}_{31}", 0, 1, unit="eV^{2}")
v_dm41 = ROOT.RooRealVar("v_dm41", "#Delta^{2}_{41}", 0, 5, unit="eV^{2}")
getattr(ws, 'import')(v_sin13)
getattr(ws, 'import')(v_sin14)
getattr(ws, 'import')(v_dm31)
getattr(ws, 'import')(v_dm41)
v_sin13 = ws.var('v_sin13')
v_sin14 = ws.var('v_sin14')
v_dm31 = ws.var('v_dm31')
v_dm41 = ws.var('v_dm41')

v_sin13.setVal(config.get("physics.oscillation.sin13")[0])
v_dm31.setVal(config.get("physics.oscillation.dm31")[0])
v_sin13.setConstant(True)
v_dm31.setConstant(True)

v_sin14.setVal(0.0) ## No oscillation
v_dm41.setVal(1.0) ## To be measured, in eV^2

v_L = ROOT.RooRealVar("v_L", "L", 0, 10000, unit="meter") ## up to 10km for now.
v_L.setVal(baseline)
getattr(ws, 'import')(v_L)
v_L = ws.var('v_L')
v_L.setConstant(True)

#v_elemFracs[0].setVal(1)
#v_elemFracs[1].setVal(0)
#v_elemFracs[2].setVal(0)
pdf_ENu = ROOT.NuOscIBDPdf("pdf_ENu", "pdf_ENu", v_ENu, v_L,
                           v_sin13, v_dm31, v_sin14, v_dm41,
                           v_elemFracs[0], v_elemFracs[1], v_elemFracs[2], v_elemFracs[3],
                           grps_HM[0], grps_HM[1], grps_HM[2], grps_HM[3],
                           grp_Xsec)
getattr(ws, 'import')(pdf_ENu)
pdf_ENu = ws.pdf('pdf_ENu')

cENu = ROOT.TCanvas("cENu", "cENu", 500, 500)
legENu = ROOT.TLegend(0.55, 0.6, 0.85, 0.85)
legENu.SetBorderSize(0)
legENu.SetFillStyle(0)
frameENu = v_ENu.frame()

v_sin13.setVal(0.0)
v_sin14.setVal(0.0)
pdf_ENu.plotOn(frameENu, ROOT.RooFit.Name("pdf_ENu_NoOsc"),
               ROOT.RooFit.LineColor(ROOT.kBlue+2), ROOT.RooFit.LineWidth(2))
legENu.AddEntry(frameENu.findObject("pdf_ENu_NoOsc"), "No Osc.")

v_sin13.setVal(config.get("physics.oscillation.sin13")[0])
pdf_ENu.plotOn(frameENu, ROOT.RooFit.Name("pdf_ENu_noNu4"),
               ROOT.RooFit.LineColor(ROOT.kRed+1), ROOT.RooFit.LineWidth(2))
legENu.AddEntry(frameENu.findObject("pdf_ENu_noNu4"), "with #nu_{3}")

v_sin14.setVal(0.5)
v_dm41.setVal(1.0)
pdf_ENu.plotOn(frameENu, ROOT.RooFit.Name("pdf_ENu_A"),
               ROOT.RooFit.LineColor(ROOT.kGreen+3), ROOT.RooFit.LineWidth(2))
legENu.AddEntry(frameENu.findObject("pdf_ENu_A"),
                "sin_{14}=%.2f, #Delta m_{41}=%.2f" % (v_sin14.getVal(), v_dm41.getVal()))

v_sin14.setVal(0.5)
v_dm41.setVal(2.0)
pdf_ENu.plotOn(frameENu, ROOT.RooFit.Name("pdf_ENu_B"),
               ROOT.RooFit.LineColor(ROOT.kOrange+5), ROOT.RooFit.LineWidth(2))
legENu.AddEntry(frameENu.findObject("pdf_ENu_B"),
                "sin_{14}=%.2f, #Delta m_{41}=%.2f" % (v_sin14.getVal(), v_dm41.getVal()))

v_sin14.setVal(0.5)
v_dm41.setVal(5.0)
pdf_ENu.plotOn(frameENu, ROOT.RooFit.Name("pdf_ENu_C"),
               ROOT.RooFit.LineColor(ROOT.kViolet+1), ROOT.RooFit.LineWidth(2))
legENu.AddEntry(frameENu.findObject("pdf_ENu_C"),
                "sin_{14}=%.2f, #Delta m_{41}=%.2f" % (v_sin14.getVal(), v_dm41.getVal()))

frameENu.Draw()
legENu.Draw()

cENu.Update()

################################################################################

################################################################################
## Build the convoluted PDF
## to do the convolution with conditional variable, the response matrix 
## has to be transposed.
################################################################################
v_EReco = ROOT.RooRealVar("v_EReco", "Reconstructed Energy", minEReco, maxEReco, unit="MeV")
v_EReco.setBins(nbinsEReco)
getattr(ws, 'import')(v_EReco)
v_EReco = ws.var('v_EReco')

_hRespT = ROOT.TH2D("hRespT", _hResp.GetTitle(),
                    _hResp.GetNbinsY(), _hResp.GetYaxis().GetXmin(), _hResp.GetYaxis().GetXmax(),
                    _hResp.GetNbinsX(), _hResp.GetXaxis().GetXmin(), _hResp.GetXaxis().GetXmax())
for i in range(_hResp.GetNbinsX()):
  for j in range(_hResp.GetNbinsY()):
    _hRespT.SetBinContent(j+1, i+1, _hResp.GetBinContent(i+1, j+1))
_dhRespT = ROOT.RooDataHist("dhRespT", _hRespT.GetTitle(),
                            ROOT.RooArgList(v_EReco, v_ENu), _hRespT)
_pdf_RespT = ROOT.RooHistPdf("pdf_RespT", "pdf_RespT",
                             ROOT.RooArgList(v_EReco, v_ENu), _dhRespT)

_pdf_Joint = ROOT.RooProdPdf("pdf_Joint", "Joint pdf",
                             ROOT.RooArgList(pdf_ENu, _pdf_RespT))
pdf_EReco = _pdf_Joint.createProjection(ROOT.RooArgSet(v_ENu))
pdf_EReco.SetName("pdf_EReco")
pdf_EReco.SetTitle("PDF of reconstructed energy")
getattr(ws, 'import')(pdf_EReco)
pdf_EReco = ws.pdf('pdf_EReco')

cEReco = ROOT.TCanvas("cEReco", "cEReco", 500, 500)
legEReco = ROOT.TLegend(0.55, 0.6, 0.85, 0.85)
legEReco.SetBorderSize(0)
legEReco.SetFillStyle(0)
frameEReco = v_EReco.frame()

v_sin13.setVal(0.0)
v_sin14.setVal(0.0)
pdf_EReco.plotOn(frameEReco, ROOT.RooFit.Name("pdf_EReco_NoOsc"),
               ROOT.RooFit.LineColor(ROOT.kBlue+2), ROOT.RooFit.LineWidth(2))
legEReco.AddEntry(frameEReco.findObject("pdf_EReco_NoOsc"), "No Osc.")

v_sin13.setVal(config.get("physics.oscillation.sin13")[0])
pdf_EReco.plotOn(frameEReco, ROOT.RooFit.Name("pdf_EReco_noNu4"),
               ROOT.RooFit.LineColor(ROOT.kRed+1), ROOT.RooFit.LineWidth(2))
legEReco.AddEntry(frameEReco.findObject("pdf_EReco_noNu4"), "with #nu_{3}")

v_sin14.setVal(0.5)
v_dm41.setVal(1.0)
pdf_EReco.plotOn(frameEReco, ROOT.RooFit.Name("pdf_EReco_A"),
               ROOT.RooFit.LineColor(ROOT.kGreen+3), ROOT.RooFit.LineWidth(2))
legEReco.AddEntry(frameEReco.findObject("pdf_EReco_A"),
                "sin_{14}=%.2f, #Delta m_{41}=%.2f" % (v_sin14.getVal(), v_dm41.getVal()))

v_sin14.setVal(0.5)
v_dm41.setVal(2.0)
pdf_EReco.plotOn(frameEReco, ROOT.RooFit.Name("pdf_EReco_B"),
               ROOT.RooFit.LineColor(ROOT.kOrange+5), ROOT.RooFit.LineWidth(2))
legEReco.AddEntry(frameEReco.findObject("pdf_EReco_B"),
                "sin_{14}=%.2f, #Delta m_{41}=%.2f" % (v_sin14.getVal(), v_dm41.getVal()))

v_sin14.setVal(0.5)
v_dm41.setVal(5.0)
pdf_EReco.plotOn(frameEReco, ROOT.RooFit.Name("pdf_EReco_C"),
               ROOT.RooFit.LineColor(ROOT.kViolet+1), ROOT.RooFit.LineWidth(2))
legEReco.AddEntry(frameEReco.findObject("pdf_EReco_C"),
                "sin_{14}=%.2f, #Delta m_{41}=%.2f" % (v_sin14.getVal(), v_dm41.getVal()))

frameEReco.Draw()
legEReco.Draw()

cEReco.Update()


################################################################################

input("Press Enter to exit...")
