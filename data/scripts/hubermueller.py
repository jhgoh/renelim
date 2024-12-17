#!/usr/bin/env python
import ROOT
import numpy as np
ROOT.gStyle.SetOptStat(0)

def buildHist(fName):
    grp = ROOT.TGraphAsymmErrors()
    with open(fName) as f:
        for l in f.readlines():
            l = l.strip()
            if len(l) == 0 or l[0] == '#': continue

            energy, flux, _, _, _, _, _, _, _, _, _, _, errLo, errHi = l.split()
            energy, flux, errLo, errHi = float(energy), float(flux), float(errLo), float(errHi)

            n = grp.GetN()
            grp.SetPoint(n, energy, flux)
            grp.SetPointError(n, 0, 0, abs(errLo), errHi)
    return grp

def cpp2Arr(n, obj):
    a = np.zeros(n)
    for i in range(n):
        a[i] = obj[i]
    return a

## Build Huber-Mueller spectrum in TGraphAsymmErrors
fout = ROOT.TFile("hubermueller.root", "recreate")

elements = ["U235", "Pu239", "Pu241"]
grps = []
for element in elements:
    grp = buildHist("original/"+element+"-anti-neutrino-flux-250keV.dat")
    grps.append(grp)

    grp.SetName("g_"+element)
    grp.SetTitle(element)

    #grp.SetDirectory(fout)
    grp.Write()

## In addtion, build histograms
hists = []
for grp, element in zip(grps, elements):
    n = grp.GetN()

    binCenters = np.array([grp.GetX()[i] for i in range(n)])
    dxs = binCenters[1:]-binCenters[:-1]
    bins = np.zeros(n+1)
    bins[:-2] = binCenters[:-1]-dxs/2
    bins[-2:] = binCenters[-2:]+dxs[-2:]/2

    h = ROOT.TH1D("h_"+element, element+";Neutrino energy (MeV);Flux", len(bins)-1, bins)
    hists.append(h)
    for i in range(n):
        val = grp.GetY()[i]
        err = (grp.GetEYhigh()[i]+grp.GetEYlow()[i])/2
        h.SetBinContent(i+1, val)
        h.SetBinError(i+1, err)
    #h.Write()

fout.Write()

## Draw them
colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+1]
for color, grp, hist in zip(colors, grps, hists):
    grp.SetLineColor(color)
    hist.SetLineColor(color)

xmax = max([grp.GetX()[grp.GetN()-1] for grp in grps])
ymax = max([np.max(cpp2Arr(grp.GetN(), grp.GetY())+cpp2Arr(grp.GetN(), grp.GetEYhigh())) for grp in grps])
hFrame = ROOT.TH1D("hFrame", ";Neutrino energy (MeV)", 100, 0, xmax)
hFrame.SetMaximum(1.2*ymax)

hFrame.Draw()
leg = ROOT.TLegend(0.6, 0.7, 0.9, 0.9)
leg.SetBorderSize(0)
#for grp in grps:
#    leg.AddEntry(grp)
#    grp.Draw("same")
for hist in hists:
    leg.AddEntry(hist)
    hist.Draw("same")
leg.Draw("brNDC")

