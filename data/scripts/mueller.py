#!/usr/bin/env python
import ROOT
import numpy as np
ROOT.gStyle.SetOptStat(0)

def tab2data(s):
    vals = {}
    for l in s.split('\n'):
        l = l.strip()
        if l == '': continue

        xx = l.split()
        vals[float(xx[0])] = float(xx[1])
    es = sorted(vals.keys())
    xs = [vals[e] for e in es]
    return np.array(es), np.array(xs)

## KE(MeV) N_nuebar(/fission) TotalErr
u235 = """
2.00 	1.31E0 	2.1
2.25 	1.11E0 	2.1
2.50 	9.27E-1 	2.1
2.75 	7.75E-1 	2.1
3.00 	6.51E-1 	2.1
3.25 	5.47E-1 	2.1
3.50 	4.49E-1 	2.3
3.75 	3.63E-1 	2.3
4.00 	2.88E-1 	2.3
4.25 	2.27E-1 	2.6
4.50 	1.77E-1 	2.6
4.75 	1.37E-1 	2.6
5.00 	1.09E-1 	2.8
5.25 	8.54E-2 	2.8
5.50 	6.56E-2 	3.8
5.75 	4.99E-2 	4.1
6.00 	3.68E-2 	4.1
6.25 	2.74E-2 	4.1
6.50 	2.07E-2 	4.1
6.75 	1.56E-2 	4.1
7.00 	1.11E-2 	4.1
7.25 	6.91E-3 	4.1
7.50 	4.30E-3 	4.4
7.75 	2.78E-3 	4.4
8.00 	1.49E-3 	4.7
"""

## KE(MeV) N_nuebar(/fission) TotalErr
Pu239 = """
2.00 1.13 	2.3 
2.25 9.19E-1 	2.3 
2.50 7.28E-1 	2.4 
2.75 6.13E-1 	2.4 
3.00 5.04E-1 	2.4 
3.25 4.10E-1 	2.4 
3.50 3.21E-1 	2.6 
3.75 2.54E-1 	2.6 
4.00 2.00E-1 	2.7 
4.25 1.51E-1 	2.9 
4.50 1.10E-1 	3.0 
4.75 7.97E-2 	3.0 
5.00 6.15E-2 	3.3 
5.25 4.68E-2 	3.3 
5.50 3.50E-2 	4.4 
5.75 2.55E-2 	4.6 
6.00 1.82E-2 	4.9 
6.25 1.32E-2 	5.0 
6.50 9.82E-3 	5.2 
6.75 7.32E-3 	5.2 
7.00 5.13E-3 	7.1 
7.25 3.15E-3 	9.2 
7.50 1.83E-3 	11.1
7.75 1.03E-3 	15.7
8.00 4.91E-4 	20.6
"""

## KE(MeV) N_nuebar(/fission) TotalErr
Pu241 = """
2.00  1.27 	2.2
2.25  1.07 	2.2
2.50  9.06E-1 	2.2
2.75  7.63E-1 	2.2
3.00  6.39E-1 	2.2
3.25  5.31E-1 	2.2
3.50  4.33E-1 	2.4
3.75  3.51E-1 	2.4
4.00  2.82E-1 	2.5
4.25  2.18E-1 	2.7
4.50  1.65E-1 	2.8
4.75  1.22E-1 	2.8
5.00  9.59E-2 	3.1
5.25  7.36E-2 	3.1
5.50  5.52E-2 	4.3
5.75  4.01E-2 	4.5
6.00  2.81E-2 	4.7
6.25  2.04E-2 	4.7
6.50  1.50E-2 	4.9
6.75  1.07E-2 	4.9
7.00  7.20E-3 	5.3
7.25  4.47E-3 	5.3
7.50 	2.54E-3 	5.7
7.75 	1.65E-3 	5.7
8.00 	9.63E-4 	7.0
"""

u238 = """
2.00 	1.43
2.25 	1.26
2.50 	1.12
2.75 	9.80E-1
3.00 	8.70E-1
3.25 	7.57E-1
3.50 	6.40E-1
3.75 	5.39E-1
4.00 	4.50E-1
4.25 	3.67E-1
4.50 	2.94E-1
4.75 	2.32E-1
5.00 	1.83E-1
5.25 	1.43E-1
5.50 	1.10E-1
5.75 	8.35E-2
6.00 	6.21E-2
6.25 	4.70E-2
6.50 	3.58E-2
6.75 	2.71E-2
7.00 	1.95E-2
7.25 	1.32E-2
7.50 	8.65E-3
7.75 	6.01E-3
8.00 	3.84E-3
"""

u238_450d = """
2.00 	1.48
2.25 	1.30
2.50 	1.15
2.75 	1.00
3.00 	8.76E-1
3.25 	7.59E-1
3.50 	6.42E-1
3.75 	5.39E-1
4.00 	4.51E-1
4.25 	3.67E-1
4.50 	2.93E-1
4.75 	2.32E-1
5.00 	1.83E-1
5.25 	1.43E-1
5.50 	1.10E-1
5.75 	8.35E-2
6.00 	6.21E-2
6.25 	4.70E-2
6.50 	3.58E-2
6.75 	2.71E-2
7.00 	1.95E-2
7.25 	1.33E-2
7.50 	8.65E-3
7.75 	6.01E-3
8.00 	3.84E-3
"""

## Build graphs
u235 = tab2data(u235)
Pu239 = tab2data(Pu239)
Pu241 = tab2data(Pu241)
u238 = tab2data(u238)

elemNames = ["U235", "Pu239", "Pu241", "U238"]
tables = [u235, Pu239, Pu241, u238]

fout = ROOT.TFile("mueller.root", "recreate")
grps = []
for elemName, table in zip(elemNames, tables):
    grp = ROOT.TGraph()
    for i, (energy, val) in enumerate(zip(table[0], table[1])):
        grp.SetPoint(i, energy, val)
    grp.SetName("g_"+elemName)
    grp.Write()

    grps.append(grp)

fout.Write()

## Draw them
ymax = max([np.max(table[1]) for table in tables])
xmax = max([np.max(table[0]) for table in tables])
hFrame = ROOT.TH1D("hFrame", ";Neutrino energy (MeV);N per fission", 100, 0, xmax)
hFrame.SetMaximum(1.2*ymax)

hFrame.Draw()
colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+1, ROOT.kMagenta+1]
for grp, color in zip(grps, colors):
    grp.SetLineColor(color)
    grp.Draw("same")

