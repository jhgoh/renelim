#!/usr/bin/env python
#from tqdm import tqdm
import argparse

parser = argparse.ArgumentParser(description='Run the hypothesis tests for sin14 and dm41')
parser.add_argument('-s', '--sin14', type=float, required=True,
                    help='Oscillation parameter sin^2(2 theta_14)')
parser.add_argument('-m', '--dm41', type=float, required=True,
                    help='Oscillation parameter Delta m^2_14')
parser.add_argument('-n', '--n_signal', type=float, required=True,
                    help='Number of expected signals')
parser.add_argument('-o', '--output', type=str, required=True,
                    help='output ROOT file name to store HypoTestResult (result)')
args = parser.parse_args()

import sys
import ROOT
import numpy as np
import resource

print("@@@ Setting up...")
print("@@@ MAXRSS now = ", resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)

###############################################################################
## Set contants
###############################################################################
nSignal = args.n_signal
sin14 = args.sin14
dm41 = args.dm41
foutName = args.output
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
detEffs = config.getDetectors('efficiency')
responses = config.getDetectors('response')
allBaselines = []
for detName in detNames:
  allBaselines.append(config.getBaselines(detName))

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
## Neutrino energy spectrums
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
_grps_HM = []
for elemName in elemNames:
    _, grp = getFileAndObj(config.get(f'physics.isotope_flux.{elemName}'))
    ROOT.gROOT.cd()
    _grps_HM.append(grp.Clone())
################################################################################

###############################################################################
## IBD cross section
###############################################################################
_fXsec, _grp_Xsec = getFileAndObj(config.get('physics.ibd_xsec'))
ROOT.gROOT.cd()
################################################################################

################################################################################
## Build the Oscillated neutrino energy spectrum
################################################################################
v_ENu = ROOT.RooRealVar("v_ENu", "Neutrino Energy", minENu, maxENu, unit="MeV")
v_ENu.setBins(nbinsENu)
getattr(ws, 'import')(v_ENu)
v_ENu = ws.var('v_ENu')
v_ENu.setConstant(True)

v_EReco = ROOT.RooRealVar("v_EReco", "Reconstructed Energy", minEReco, maxEReco, unit="MeV")
v_EReco.setBins(nbinsEReco)
getattr(ws, 'import')(v_EReco)
v_EReco = ws.var('v_EReco')

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

v_sin13.setVal(config.get("physics.oscillation.sin13"))
v_dm31.setVal(config.get("physics.oscillation.dm31"))
v_sin13.setConstant(True)
v_dm31.setConstant(True)

v_sin14.setVal(0.0) ## No oscillation
v_dm41.setVal(1.0) ## To be measured, in eV^2

v_L = ROOT.RooRealVar("v_L", "L", 0, 10000, unit="meter") ## up to 10km for now.
v_L.setVal(baseline)
getattr(ws, 'import')(v_L)
v_L = ws.var('v_L')
v_L.setConstant(True)

pdf_ENu = ROOT.NuOscIBDPdf("pdf_ENu", "pdf_ENu", v_ENu, v_L,
                           v_sin13, v_dm31, v_sin14, v_dm41,
                           v_elemFracs[0], v_elemFracs[1], v_elemFracs[2], v_elemFracs[3],
                           _grps_HM[0], _grps_HM[1], _grps_HM[2], _grps_HM[3],
                           _grp_Xsec)
getattr(ws, 'import')(pdf_ENu)
pdf_ENu = ws.pdf('pdf_ENu')
################################################################################

################################################################################
## Build the convoluted PDF
## to do the convolution with conditional variable, the response matrix 
## has to be transposed.
################################################################################
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

## Reduce precision of integrator for speed up
#epsAbs = v_ENu.getIntegratorConfig().epsAbs()
#epsRel = v_ENu.getIntegratorConfig().epsRel()
#v_ENu.getIntegratorConfig().setEpsAbs(epsAbs)
#v_ENu.getIntegratorConfig().setEpsRel(epsRel)
#v_ENu.getIntegratorConfig().setEpsAbs(5e-7)
#v_ENu.getIntegratorConfig().setEpsRel(5e-7)

################################################################################
## Start setting up RooStats
################################################################################
ws.factory('v_nSignal[0, 1e9]')
v_nSignal = ws.var('v_nSignal')
_model = ROOT.RooExtendPdf("model", "Signal-only PDF", pdf_EReco, v_nSignal) 
getattr(ws, 'import')(_model)

## Set the null hypothesis
mcNull = ROOT.RooStats.ModelConfig("mcNull", ws)
mcNull.SetPdf(ws.pdf("model"))
mcNull.SetObservables(ROOT.RooArgSet(v_EReco))
mcNull.SetParametersOfInterest(ROOT.RooArgSet(v_sin14, v_dm41))

v_nSignal.setVal(nSignal)
v_nSignal.setConstant(True)
v_sin14.setVal(0.0)
v_dm41.setVal(0.0)
v_sin14.setConstant(True)
v_dm41.setConstant(True)

poi_nuis_params = ROOT.RooArgSet()
poi_nuis_params.add(mcNull.GetParametersOfInterest())
#poi_nuis_params.add(mcNull.GetNuisanceParameters())
mcNull.SetSnapshot(poi_nuis_params)

## Set the alternative hypothesis
mcAlt = ROOT.RooStats.ModelConfig("mcAlt", ws)
mcAlt.SetPdf(ws.pdf("model"))
mcAlt.SetObservables(ROOT.RooArgSet(v_EReco))
mcAlt.SetParametersOfInterest(ROOT.RooArgSet(v_sin14, v_dm41))

v_nSignal.setVal(nSignal)
v_nSignal.setConstant(False)
v_sin14.setVal(sin14)
v_dm41.setVal(dm41)
v_sin14.setConstant(False)
v_dm41.setConstant(False)

mcAlt.SetSnapshot(poi_nuis_params)

## Create Asimov dataset
#asimovData = mcNull.GetPdf().generateAsimovData(ROOT.RooArgSet(mcNull.GetObservables()))
asimovData = ROOT.RooStats.AsymptoticCalculator.MakeAsimovData(mcNull, mcNull.GetObservables(), ROOT.RooArgSet())

## Scan over the parameters of interests
calc = ROOT.RooStats.AsymptoticCalculator(asimovData, mcAlt, mcNull)

result = calc.GetHypoTest()
#signif = result.Significance()
#signif = calc.ExpectedSignificance(mcAlt.GetParametersOfInterest())

result.Print()

fout = ROOT.TFile(foutName, 'recreate')
result.SetName("result")
result.SetTitle("Test result sin^{2}(2#theta_{14})="+f"{sin14}"+"#Delta m^{2}_{41}="+f"{dm41}")
result.Write()
fout.Write()
fout.Close()

print("@@@ Finishing...")
print("@@@ MAXRSS now = ", resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
