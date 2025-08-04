#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description='Run the chi2 tests for sin14 and dm41')
parser.add_argument('-s', '--sin14', type=str,
                    help='Oscillation parameter sin^2(2 theta_14)')
parser.add_argument('-m', '--dm41', type=float, required=True,
                    help='Oscillation parameter Delta m^2_14')
parser.add_argument('-n', '--nsignal', type=float, required=True,
                    help='Number of expected signals')
parser.add_argument('-o', '--output', type=str, required=True,
                    help='output ROOT file name to store results')
parser.add_argument('--toys', type=int, default=1000,
                    help='Number of ToyMC samples')
parser.add_argument('--seed', type=int, default=0,
                    help='Random number seed')
parser.add_argument('-g', '--gui', action='store_true',
                    help='Display results using ROOT TCanvas')
args = parser.parse_args()

import sys
import ROOT
if args.gui:
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
import numpy as np
from array import array
import resource
import time
from tqdm import tqdm

print("@@@ Setting up...")
print("@@@ MAXRSS now = ", resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)

###############################################################################
## Set constants
###############################################################################
nSignal = args.nsignal
dm41 = args.dm41
foutName = args.output
ROOT.RooRandom.randomGenerator().SetSeed(args.seed)

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

sys.path.append('python')
from ModelConfig import load_model

model = load_model('config.yaml')
ws = model['ws']
config = model['config']
v_sin13 = ws.var('v_sin13')
v_sin14 = model['v_sin14']
v_dm31 = ws.var('v_dm31')
v_dm41 = model['v_dm41']
v_ENu = model['v_ENu']
v_EReco = model['v_EReco']
pdf_ENu = model['pdf_ENu']
pdf_EReco = model['pdf_EReco']

v_sin14.setVal(sin14)
v_dm41.setVal(dm41)

# Build Gaussian constraints for sin13 and dm31
s13, s13err = config.get("physics.oscillation.sin13")
dm31, dm31err = config.get("physics.oscillation.dm31")
c_s13_m = ROOT.RooConstVar("v_constr_sin13_m", "sin13_mean", s13)
c_s13_s = ROOT.RooConstVar("v_constr_sin13_s", "sin13_sigma", s13err)
c_dm31_m = ROOT.RooConstVar("v_constr_dm31_m", "dm31_mean", dm31)
c_dm31_s = ROOT.RooConstVar("v_constr_dm31_s", "dm31_sigma", dm31err)
v_constr_sin13 = ROOT.RooGaussian("v_constr_sin13", "Constraint on sin13",
                                  v_sin13, c_s13_m, c_s13_s)
v_constr_dm31 = ROOT.RooGaussian("v_constr_dm31", "Constraint on dm31",
                                 v_dm31, c_dm31_m, c_dm31_s)
constrs = ROOT.RooArgSet(v_constr_sin13, v_constr_dm31)

# Save default integrator precision so we can restore it later if needed
#epsAbs_default = v_ENu.getIntegratorConfig().epsAbs()
#epsRel_default = v_ENu.getIntegratorConfig().epsRel()
#v_ENu.getIntegratorConfig().setEpsAbs(5e-7)
#v_ENu.getIntegratorConfig().setEpsRel(5e-7)
# To restore the original values:
# v_ENu.getIntegratorConfig().setEpsAbs(epsAbs_default)
# v_ENu.getIntegratorConfig().setEpsRel(epsRel_default)

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

status_scanned = []
covQual_scanned = []
edm_scanned = []
time_scanned = []

for _sin14 in tqdm(sin14_toScan):
    v_sin14.setVal(_sin14)
    v_sin14.setConstant(True)

    t0 = time.time()
    _status = minimizer.migrad()
    minimizer.hesse()
    _nll = nll.getVal()
    t1 = time.time()

    if _status == -1:
      print(f"Fit failure ({_status}): sin14={_sin14}, nll={_nll}")
      continue

    result = minimizer.save()

    sin14_scanned.append(_sin14)
    nll_scanned.append(_nll)

    status_scanned.append(_status)
    covQual_scanned.append(result.covQual())
    edm_scanned.append(result.edm())
    time_scanned.append(t1-t0)

sin14_scanned = np.array(sin14_scanned)
nll_scanned = np.array(nll_scanned)

status_scanned = np.array(status_scanned)
covQual_scanned = np.array(covQual_scanned)
edm_scanned = np.array(edm_scanned)
time_scanned = np.array(time_scanned)

if len(sin14_scanned) == 0:
  print("No successful fits")
  exit()

_idxs = np.argsort(sin14_scanned)

sin14_scanned = sin14_scanned[_idxs]
nll_scanned = nll_scanned[_idxs]

status_scanned = status_scanned[_idxs]
covQual_scanned = covQual_scanned[_idxs]
edm_scanned = edm_scanned[_idxs]
time_scanned = time_scanned[_idxs]

## Perform minimization to find a global minimum of this search
_imin = np.argmin(nll_scanned)
v_sin14.setVal(sin14_scanned[_imin])
v_sin14.setConstant(False)
minimizer.setStrategy(2)
t0 = time.time()
_status = minimizer.migrad()
minimizer.hesse()
t1 = time.time()
_sin14 = v_sin14.getVal()
result = minimizer.save()
if _status != -1 and _sin14 not in sin14_scanned:
  _idx = np.searchsorted(sin14_scanned, _sin14, side='left')

  sin14_scanned = np.insert(sin14_scanned, _idx, _sin14)
  nll_scanned = np.insert(nll_scanned, _idx, nll.getVal())

  status_scanned = np.insert(status_scanned, _idx, _status)
  covQual_scanned = np.insert(covQual_scanned, _idx, result.covQual())
  edm_scanned = np.insert(edm_scanned, _idx, result.edm())
  time_scanned = np.insert(time_scanned, _idx, t1-t0)

pNull_scanned = np.zeros(len(sin14_scanned))
pAlt_scanned = np.zeros(len(sin14_scanned))
for i, _sin14 in tqdm(enumerate(sin14_scanned)):
  v_sin14.setVal(_sin14)
  v_sin14.setConstant(True)
  mcAlt.SetSnapshot(vs_poi)

  ## Create test statistic and calculators
  calcFreq = ROOT.RooStats.FrequentistCalculator(asimovData, mcNull, mcAlt)
  calcFreq.SetToys(args.toys, args.toys)
  sampler = calcFreq.GetTestStatSampler()
  resFreq = calcFreq.GetHypoTest()

  pNull_scanned[i] = resFreq.NullPValue()
  pAlt_scanned[i] = resFreq.AlternatePValue()

fout = ROOT.TFile(foutName, 'recreate')
tree = ROOT.TTree("limit", "limit tree")

ptr_dm41 = array('d', [0])
ptr_sin14 = array('d', [0])
ptr_nll = array('d', [0])
ptr_pNull = array('d', [0])
ptr_pAlt = array('d', [0])

ptr_status = array('i', [0])
ptr_covQual = array('i', [0])
ptr_edm = array('d', [0])
ptr_time = array('d', [0])

tree.Branch("dm41", ptr_dm41, "dm41/D")
tree.Branch("sin14", ptr_sin14, "sin14/D")
tree.Branch("nll", ptr_nll, "nll/D")
tree.Branch("pNull", ptr_pNull, "pNull/D")
tree.Branch("pAlt", ptr_pAlt, "pAlt/D")

for i in range(len(sin14_scanned)):
  ptr_dm41[0] = dm41
  ptr_sin14[0] = sin14_scanned[i]

  ptr_nll[0] = nll_scanned[i]
  ptr_pNull[0] = pNull_scanned[i]
  ptr_pAlt[0] = pAlt_scanned[i]

  ptr_status[0] = status_scanned[i]
  ptr_covQual[0] = covQual_scanned[i]
  ptr_edm[0] = edm_scanned[i]
  ptr_time[0] = time_scanned[i]

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
