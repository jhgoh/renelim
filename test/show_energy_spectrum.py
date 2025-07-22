#!/usr/bin/env python
import ROOT
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
import sys
import numpy as np
from array import array

sys.path.append('python')
from ModelConfig import load_model

model = load_model('config.yaml')
ws = model['ws']
config = model['config']
v_sin13 = model['v_sin13']
v_sin14 = model['v_sin14']
v_dm31 = model['v_dm31']
v_dm41 = model['v_dm41']
v_L = model['v_L']
v_ENu = model['v_ENu']
v_EReco = model['v_EReco']
pdf_ENu = model['pdf_ENu']
pdf_EReco = model['pdf_EReco']
constrs = model['constrs']

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
reactorIdx = int(np.argmin(baselines))
baseline = baselines[reactorIdx]

xmin, xmax = np.min(reaPoses), np.max(reaPoses)
minENu = v_ENu.getMin()
maxENu = v_ENu.getMax()
minEReco = v_EReco.getMin()
maxEReco = v_EReco.getMax()
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

###############################################################################
## Neutrino energy spectrums
###############################################################################
## (Maybe later) Consider the burn-up effect, 
## the fuel composition is a subject to be changed in time.
## We choose the 6th core only (1st one in the configuration)
elemNames = config.getReactors('elements')[reactorIdx]['name']
_elemFracs = config.getReactors('elements')[reactorIdx]['fraction']
formula = "1"
formulaVars = []
for i in range(len(elemNames)-1):
    elemName = elemNames[i]
    ws.factory(f'v_{elemName}[{_elemFracs[i]}, 0, 1]')
    ws.var(f'v_{elemName}').setConstant(True)
    formula += f"-@{i}"
    formulaVars.append(f'v_{elemName}')
formulaVars = ','.join(formulaVars)
ws.factory(f'EXPR::v_{elemNames[-1]}("{formula}", {{ {formulaVars} }})')
del(formula)
del(formulaVars)

v_elemFracs = ROOT.RooArgList()
for i in range(len(elemNames)-1):
  v_elemFracs.add(ws.var(f'v_{elemNames[i]}'))
v_elemFracs.add(ws.function(f'v_{elemNames[-1]}'))

## Load the Neutrino flux model, such as Huber-Mueller, incorporating the fuel compositions
grps_HM = ROOT.std.vector('TGraph')()
for elemName in elemNames:
    _, grp = getFileAndObj(config.get(f'physics.isotope_flux.{elemName}'))
    ROOT.gROOT.cd()
    grps_HM.push_back(grp.Clone())
    del(grp)
################################################################################

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
## IBD cross section
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

## Reduce precision of integrator for speed up
#epsAbs = v_ENu.getIntegratorConfig().epsAbs()
#epsRel = v_ENu.getIntegratorConfig().epsRel()
#v_ENu.getIntegratorConfig().setEpsAbs(1e-3)
#v_ENu.getIntegratorConfig().setEpsRel(1e-3)
#v_ENu.getIntegratorConfig().setEpsAbs(5e-7)
#v_ENu.getIntegratorConfig().setEpsRel(5e-7)
#v_ENu.getIntegratorConfig().setEpsAbs(epsAbs)
#v_ENu.getIntegratorConfig().setEpsRel(epsRel)

################################################################################

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
