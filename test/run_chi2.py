#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description='Run the chi2 tests for sin14 and dm41')
parser.add_argument('-s', '--sin14', type=str,
                    help='Oscillation parameter sin^2(2 theta_14)')
parser.add_argument('-m', '--dm41', type=float, required=True,
                    help='Oscillation parameter Delta m^2_14')
parser.add_argument('-n', '--n_signal', type=float, required=True,
                    help='Number of expected signals')
parser.add_argument('-o', '--output', type=str, required=True,
                    help='output ROOT file name to store results')
parser.add_argument('--tol', type=float, default=1e-4,
                    help='Tolerance factor for adaptive NLL scan')
parser.add_argument('--toys', type=int, default=1000,
                    help='Number of ToyMC samples')
parser.add_argument('-g', '--gui', action='store_true',
                    help='Display results using ROOT TCanvas')
args = parser.parse_args()

import sys
import ROOT
import numpy as np
from array import array
import resource
from tqdm import tqdm

print("@@@ Setting up...")
print("@@@ MAXRSS now = ", resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)

###############################################################################
## Set constants
###############################################################################
nSignal = args.n_signal
dm41 = args.dm41
foutName = args.output

sin14_toScan = [0.0]
if args.sin14:
  sin14_toScan = [float(x) for x in args.sin14.split(',')]
  sin14_toScan = np.array(sin14_toScan)
else:
  sin14_toScan = np.concatenate([
                   np.zeros(1),
                   np.arange(1e-4, 1e-3, 1e-4), ## We are still interested to see NLL landscape near the null hypothesis
                   #np.arange(1e-3, 1e-2, 5e-4),
                   #np.arange(1e-2, 1e-1, 5e-3),
                   np.arange(1e-1, 5e-1, 1e-1), ## Roughly scan over a wide range
                   np.ones(1)*0.5,
                 ])
  print("@@@ sin14 values are not given. Perform scanning with default points (may take 1 hour):")
  print(sin14_toScan)

sin14 = sin14_toScan[-1] ## Set a dummy value
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
## Important parameters of interests
## Define these variables in advance, to avoid possible buggy behaviours
###############################################################################
v_sin13 = ROOT.RooRealVar("v_sin13", "sin^{2}(#theta_{13})", 0, 0.5)
v_sin14 = ROOT.RooRealVar("v_sin14", "sin^{2}(#theta_{14})", 0, 0.5)
v_dm31 = ROOT.RooRealVar("v_dm31", "#Delta^{2}_{31}", 0, 1, unit="eV^{2}")
v_dm41 = ROOT.RooRealVar("v_dm41", "#Delta^{2}_{41}", 0, 5, unit="eV^{2}")
## Note: mind the ordering
ws.Import(v_sin14)
ws.Import(v_dm41)
ws.Import(v_sin13)
ws.Import(v_dm31)
v_sin13 = ws.var('v_sin13')
v_sin14 = ws.var('v_sin14')
v_dm31 = ws.var('v_dm31')
v_dm41 = ws.var('v_dm41')
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
    del(grp)
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
ws.Import(v_ENu)
v_ENu = ws.var('v_ENu')

v_EReco = ROOT.RooRealVar("v_EReco", "Reconstructed Energy", minEReco, maxEReco, unit="MeV")
v_EReco.setBins(nbinsEReco)
ws.Import(v_EReco)
v_EReco = ws.var('v_EReco')

_sin13, _sin13err = config.get("physics.oscillation.sin13")
_dm31, _dm31err = config.get("physics.oscillation.dm31")
v_sin13.setVal(_sin13)
v_dm31.setVal(_dm31)
#v_sin13.setConstant(True)
#v_dm31.setConstant(True)
_v_constr_sin13_m = ROOT.RooConstVar("v_constr_sin13_m", "sin13_mean", _sin13)
_v_constr_sin13_s = ROOT.RooConstVar("v_constr_sin13_s", "sin13_sigma", _sin13err)
_v_constr_dm31_m = ROOT.RooConstVar("v_constr_dm31_m", "dm31_mean", _dm31)
_v_constr_dm31_s = ROOT.RooConstVar("v_constr_dm31_s", "dm31_sigma", _dm31err)
v_constr_sin13 = ROOT.RooGaussian("v_constr_sin13", "Constraint on sin13", v_sin13, _v_constr_sin13_m, _v_constr_sin13_s)
v_constr_dm31 = ROOT.RooGaussian("v_constr_dm31", "Constraint on dm31", v_dm31, _v_constr_dm31_m, _v_constr_dm31_s)
ws.Import(v_constr_sin13)
ws.Import(v_constr_dm31)
v_constr_sin13 = ws.pdf('v_constr_sin13')
v_constr_dm31 = ws.pdf('v_constr_dm31')

constrs = ROOT.RooArgSet(v_constr_sin13, v_constr_dm31)

v_sin14.setVal(sin14)
v_dm41.setVal(dm41)

v_L = ROOT.RooRealVar("v_L", "L", 0, 10000, unit="meter") ## up to 10km for now.
v_L.setVal(baseline)
ws.Import(v_L)
v_L = ws.var('v_L')
v_L.setConstant(True)

pdf_ENu = ROOT.NuOscIBDPdf("pdf_ENu", "pdf_ENu", v_ENu, v_L,
                           v_sin13, v_dm31, v_sin14, v_dm41,
                           v_elemFracs[0], v_elemFracs[1], v_elemFracs[2], v_elemFracs[3],
                           _grps_HM[0], _grps_HM[1], _grps_HM[2], _grps_HM[3],
                           _grp_Xsec)
ws.Import(pdf_ENu)
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
                             ROOT.RooArgList(pdf_ENu, _pdf_RespT),
                             ROOT.RooFit.Conditional(ROOT.RooArgSet(_pdf_RespT), ROOT.RooArgSet(v_EReco)))
pdf_EReco = _pdf_Joint.createProjection(ROOT.RooArgSet(v_ENu))
pdf_EReco.SetName("pdf_EReco")
pdf_EReco.SetTitle("PDF of reconstructed energy")
ws.Import(pdf_EReco)
pdf_EReco = ws.pdf('pdf_EReco')

## Reduce precision of integrator for speed up
#epsAbs = v_ENu.getIntegratorConfig().epsAbs()
#epsRel = v_ENu.getIntegratorConfig().epsRel()
#v_ENu.getIntegratorConfig().setEpsAbs(1e-2)
#v_ENu.getIntegratorConfig().setEpsRel(1e-2)
#v_ENu.getIntegratorConfig().setEpsAbs(5e-7)
#v_ENu.getIntegratorConfig().setEpsRel(5e-7)
#v_ENu.getIntegratorConfig().setEpsAbs(epsAbs)
#v_ENu.getIntegratorConfig().setEpsRel(epsRel)

################################################################################
## Start setting up RooStats
################################################################################
v_ENu.setConstant(True)

ws.factory('v_nSignal[0, 1e9]')
v_nSignal = ws.var('v_nSignal')
_model = ROOT.RooExtendPdf("model", "Signal-only PDF", pdf_EReco, v_nSignal)
ws.Import(_model)
model = ws.pdf('model')

vs_obs = ROOT.RooArgSet(v_EReco)
vs_poi = ROOT.RooArgSet(v_sin14) #, v_dm41)
vs_nui = ROOT.RooArgSet()
vs_nui.add(v_nSignal)
vs_nui.add(v_sin13)
vs_nui.add(v_dm31)
vs_poi_nui = ROOT.RooArgSet()
vs_poi_nui.add(vs_poi)
vs_poi_nui.add(vs_nui)

## Set the null hypothesis and keep the snapshot
mcNull = ROOT.RooStats.ModelConfig("mcNull", ws)
mcNull.SetPdf(model)
mcNull.SetObservables(vs_obs)
mcNull.SetParametersOfInterest(vs_poi)
mcNull.SetNuisanceParameters(vs_nui)

v_nSignal.setVal(nSignal)
#v_nSignal.setConstant(True)
v_sin14.setVal(0.0)
v_dm41.setVal(0.0)
v_sin14.setConstant(True)
v_dm41.setConstant(True)

mcNull.SetSnapshot(vs_poi)

## Create Asimov dataset
#asimovData = mcNull.GetPdf().generateBinned(mcNull.GetObservables(), ROOT.RooFit.Extended(), ROOT.RooFit.Asimov())
asimovData = ROOT.RooStats.AsymptoticCalculator.MakeAsimovData(mcNull, mcNull.GetObservables(), ROOT.RooArgSet())

## Set the alternative hypothesis and keep the snapshot
mcAlt = ROOT.RooStats.ModelConfig("mcAlt", ws)
mcAlt.SetPdf(model)
mcAlt.SetObservables(vs_obs)
mcAlt.SetParametersOfInterest(vs_poi)
mcAlt.SetNuisanceParameters(vs_nui)

v_nSignal.setVal(nSignal)
#v_nSignal.setConstant(True)
v_sin14.setVal(sin14)
v_dm41.setVal(dm41)
v_sin14.setConstant(True)
v_dm41.setConstant(True)

mcAlt.SetSnapshot(vs_poi)

## Calculate the likelihood
nll = model.createNLL(asimovData, ROOT.RooFit.Constrain(constrs)) 

minimizer = ROOT.RooMinuit(nll)
minimizer.setStrategy(2)
#minimizer.setPrintLevel(3)
minimizer.setEps(1e-6)

sin14_scanned = []
nll_scanned = []
for _sin14 in tqdm(sin14_toScan):
    v_sin14.setVal(_sin14)
    v_sin14.setConstant(True)

    _status = minimizer.migrad()
    _nll = nll.getVal()

    if _status == -1:
      print(f"Fit failure ({_status}): sin14={_sin14}, nll={_nll}")
      continue

    sin14_scanned.append(_sin14)
    nll_scanned.append(_nll)

sin14_scanned = np.array(sin14_scanned)
nll_scanned = np.array(nll_scanned)

if len(sin14_scanned) == 0:
  print("No successful fits")
  exit()

_idxs = np.argsort(sin14_scanned)
sin14_scanned = sin14_scanned[_idxs]
nll_scanned = nll_scanned[_idxs]

## Perform minimization to find a global minimum of this search
_imin = np.argmin(nll_scanned)
v_sin14.setVal(sin14_scanned[_imin])
v_sin14.setConstant(False)
minimizer.setStrategy(2)
_status = minimizer.migrad()
_sin14 = v_sin14.getVal()
if _status != -1 and _sin14 not in sin14_scanned:
  _idx = np.searchsorted(sin14_scanned, _sin14, side='left')
  sin14_scanned = np.insert(sin14_scanned, _idx, _sin14)
  nll_scanned = np.insert(nll_scanned, _idx, nll.getVal())

## Proceed to the Feldman-Cousins tests
#pl = ROOT.RooStats.ProfileLikelihoodTestStat(mcNull.GetPdf())
#pl.SetOneSided(True)

t90_scanned = np.zeros(len(sin14_scanned))
pNull_scanned = np.zeros(len(sin14_scanned))
pAlt_scanned = np.zeros(len(sin14_scanned))
for i, _sin14 in tqdm(enumerate(sin14_scanned)):
  v_sin14.setVal(_sin14)
  v_sin14.setConstant(True)
  mcAlt.SetSnapshot(vs_poi)

  t90 = 0
  #limitFC = ROOT.RooStats.FeldmanCousins(asimovData, mcNull)
  #limitFC.SetTestSize(0.1)
  #limitFC.SetTestStatistic(pl)
  #limitFC.SetNuisanceParameters(vs_nui)
  #limitFC.SetNBins(20)
  #limitFC.UseAdaptiveSampling(True)
  #limitFC.FluctuateNumDataEntries(False)
  #interval = limitFC.GetInterval()
  #upper_limit = interval.UpperLimit(vs_poi)

  ## Create test statistic and calculators
  calcFreq = ROOT.RooStats.FrequentistCalculator(asimovData, mcNull, mcAlt)
  calcFreq.SetToys(args.toys, args.toys)
  resFreq = calcFreq.GetHypoTest()

  t90_scanned[i] = t90
  pNull_scanned[i] = resFreq.NullPValue()
  pAlt_scanned[i] = resFreq.AlternatePValue()

fout = ROOT.TFile(foutName, 'recreate')
tree = ROOT.TTree("limit", "limit tree")

ptr_dm41 = array('d', [0])
ptr_sin14 = array('d', [0])
ptr_nll = array('d', [0])
ptr_t90 = array('d', [0])
ptr_pNull = array('d', [0])
ptr_pAlt = array('d', [0])

tree.Branch("dm41", ptr_dm41, "dm41/D")
tree.Branch("sin14", ptr_sin14, "sin14/D")
tree.Branch("nll", ptr_nll, "nll/D")
tree.Branch("t90", ptr_t90, "t90/D")
tree.Branch("pNull", ptr_pNull, "pNull/D")
tree.Branch("pAlt", ptr_pAlt, "pAlt/D")

for i in range(len(sin14_scanned)):
  ptr_dm41[0] = dm41
  ptr_sin14[0] = sin14_scanned[i]
  ptr_nll[0] = nll_scanned[i]
  ptr_t90[0] = t90_scanned[i]
  ptr_pNull[0] = pNull_scanned[i]
  ptr_pAlt[0] = pAlt_scanned[i]

  tree.Fill()
tree.Write()

if args.gui:
  grpNLL = ROOT.TGraph(len(nll_scanned), array('d', sin14_scanned), array('d', nll_scanned))
  grpNLL.SetLineWidth(2)
  grpNLL.SetEditable(False)

  cNLL = ROOT.TCanvas("cNLL", "cNLL", 500, 500)
  grpNLL.Draw("ALP")
  cNLL.Update()

  cAS = ROOT.TCanvas("cAS", "cAS", 500, 500)
  frameAS = v_EReco.frame()
  asimovData.plotOn(frameAS, ROOT.RooFit.LineColor(ROOT.kBlack), ROOT.RooFit.MarkerSize(1), ROOT.RooFit.DataError(getattr(ROOT.RooAbsData, "None")))
  mcAlt.LoadSnapshot()
  mcAlt.GetWS().pdf("model").plotOn(frameAS, ROOT.RooFit.LineColor(ROOT.kRed))
  mcNull.LoadSnapshot()
  mcNull.GetWS().pdf("model").plotOn(frameAS, ROOT.RooFit.LineColor(ROOT.kBlue))
  frameAS.Draw()

  input("")

  grpNLL.SetName(f"grpNLL_dm41_{dm41:g}".replace('-', 'm').replace('.','p'))
  grpNLL.SetTitle("NLL scan for #Delta m^{2}_{41} = "+f"{dm41}"+";sin^{2}(2#theta_{14});-log(L)")
  grpNLL.Write()

fout.Write()
fout.Close()

print("@@@ Finishing...")
print(f"@@@ Scanned for dm41={dm41}, sin14={sin14_scanned}")
print("@@@ MAXRSS now = ", resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
