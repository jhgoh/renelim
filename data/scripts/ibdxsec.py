#!/usr/bin/env python
import ROOT
import numpy as np
ROOT.gStyle.SetOptStat(0)

def tab2data(s):
    vals = {}
    for l in s.split('\n'):
        l = l.strip()
        if l == '': continue

        e1, x1, e2, x2, e3, x3 = l.split()
        vals[float(e1)] = float(x1)
        vals[float(e2)] = float(x2)
        vals[float(e3)] = float(x3)
    es = sorted(vals.keys())
    xs = [vals[e] for e in es]
    return np.array(es), np.array(xs)

xsecHiE = """
2. 0.00331709 35. 8.42907 68. 25.9144
3. 0.0265181 36. 8.86955 69. 26.5082
4. 0.0680334 37. 9.31767 70. 27.1041
5. 0.127583 38. 9.77319 71. 27.7021
6. 0.204738 39. 10.2359 72. 28.302
7. 0.299076 40. 10.7055 73. 28.9038
8. 0.41018 41. 11.1819 74. 29.5074
9. 0.537644 42. 11.6648 75. 30.1126
10. 0.681069 43. 12.1541 76. 30.7195
11. 0.840063 44. 12.6494 77. 31.3278
12. 1.01424 45. 13.1507 78. 31.9375
13. 1.20324 46. 13.6577 79. 32.5485
14. 1.40667 47. 14.1702 80. 33.1607
15. 1.62418 48. 14.688 81. 33.7741
16. 1.85542 49. 15.211 82. 34.3886
17. 2.10003 50. 15.739 83. 35.004
18. 2.35767 51. 16.2718 84. 35.6204
19. 2.62801 52. 16.8093 85. 36.2376
20. 2.91071 53. 17.3512 86. 36.8556
21. 3.20546 54. 17.8974 87. 37.4743
22. 3.51193 55. 18.4478 88. 38.0937
23. 3.82982 56. 19.0022 89. 38.7136
24. 4.15881 57. 19.5604 90. 39.3341
25. 4.49861 58. 20.1223 91. 39.955
26. 4.84893 59. 20.6877 92. 40.5763
27. 5.20946 60. 21.2566 93. 41.1979
28. 5.57995 61. 21.8287 94. 41.8199
29. 5.96009 62. 22.404 95. 42.4421
30. 6.34963 63. 22.9823 96. 43.0645
31. 6.7483 64. 23.5635 97. 43.687
32. 7.15583 65. 24.1474 98. 44.3096
33. 7.57197 66. 24.7339 99. 44.9322
34. 7.99647 67. 25.323 100. 45.5549
"""

xsecLoE = """
1.9 0.00190183 5.3 0.1489 8.7 0.497711
2. 0.00331709 5.4 0.156356 8.8 0.510862
2.1 0.00484225 5.5 0.163987 8.9 0.524173
2.2 0.00652675 5.6 0.171791 9. 0.537644
2.3 0.00838533 5.7 0.179769 9.1 0.551275
2.4 0.0104239 5.8 0.187919 9.2 0.565064
2.5 0.0126452 5.9 0.196243 9.3 0.579013
2.6 0.0150505 6. 0.204738 9.4 0.59312
2.7 0.0176403 6.1 0.213406 9.5 0.607386
2.8 0.0204149 6.2 0.222245 9.6 0.621809
2.9 0.0233742 6.3 0.231255 9.7 0.636389
3. 0.0265181 6.4 0.240435 9.8 0.651126
3.1 0.0298462 6.5 0.249786 9.9 0.666019
3.2 0.0333584 6.6 0.259306 10. 0.681069
3.3 0.0370542 6.7 0.268996 10.1 0.696274
3.4 0.0409332 6.8 0.278854 10.2 0.711634
3.5 0.0449951 6.9 0.288881 10.3 0.72715
3.6 0.0492395 7. 0.299076 10.4 0.74282
3.7 0.0536659 7.1 0.309439 10.5 0.758644
3.8 0.058274 7.2 0.319969 10.6 0.774622
3.9 0.0630633 7.3 0.330665 10.7 0.790753
4. 0.0680334 7.4 0.341529 10.8 0.807037
4.1 0.0731839 7.5 0.352558 10.9 0.823474
4.2 0.0785142 7.6 0.363753 11. 0.840063
4.3 0.0840241 7.7 0.375113 11.1 0.856804
4.4 0.089713 7.8 0.386638 11.2 0.873697
4.5 0.0955806 7.9 0.398327 11.3 0.890741
4.6 0.101626 8. 0.41018 11.4 0.907935
4.7 0.10785 8.1 0.422197 11.5 0.92528
4.8 0.114251 8.2 0.434377 11.6 0.942774
4.9 0.120829 8.3 0.44672 11.7 0.960418
5. 0.127583 8.4 0.459225 11.8 0.978212
5.1 0.134513 8.5 0.471892 11.9 0.996154
5.2 0.141619 8.6 0.484721 12. 1.01424
"""

## Build graphs
energyHi, xsecHi = tab2data(xsecHiE)
energyLo, xsecLo = tab2data(xsecLoE)

## Find bin edges
binsHi = np.zeros(len(energyHi)+1)
dxsHi = energyHi[1:]-energyHi[:-1]
binsHi[:-2] = energyHi[:-1]-dxsHi/2
binsHi[-2:] = energyHi[-2:]+dxsHi[-2:]/2

binsLo = np.zeros(len(energyLo)+1)
dxsLo = energyLo[1:]-energyLo[:-1]
binsLo[:-2] = energyLo[:-1]-dxsLo/2
binsLo[-2:] = energyLo[-2:]+dxsLo[-2:]/2

fout = ROOT.TFile("ibdxsec.root", "recreate")
gHiE = ROOT.TGraph()
for i, (energy, xsec) in enumerate(zip(energyHi, xsecHi)):
    gHiE.SetPoint(i, energy, xsec)
gLoE = ROOT.TGraph()
for i, (energy, xsec) in enumerate(zip(energyLo, xsecLo)):
    gLoE.SetPoint(i, energy, xsec)
gHiE.SetName("g_HighE")
gLoE.SetName("g_LowE")
gHiE.Write()
gLoE.Write()

hHiE = ROOT.TH1D("h_HighE", "High E;Neutrino energy (MeV);Cross section (10^{-41}cm^{2})", len(binsHi)-1, binsHi)
hLoE = ROOT.TH1D("h_LowE", "Low E;Neutrino energy (MeV);Cross section (10^{-41}cm^{2})", len(binsLo)-1, binsLo)
for i in range(hHiE.GetNbinsX()):
    hHiE.SetBinContent(i+1, xsecHi[i])
for i in range(hLoE.GetNbinsX()):
    hLoE.SetBinContent(i+1, xsecLo[i])

fout.Write()

## Draw them
colors = [ROOT.kRed, ROOT.kBlue]
ymax = max(np.max(xsecHi), np.max(xsecLo))
hFrame = ROOT.TH1D("hFrame", ";Neutrino energy (MeV);Cross section (10^{-41}cm^{2})", 100, 0, np.max(energyHi))
hFrame.SetMaximum(1.2*ymax)

hFrame.Draw()
hHiE.Draw("same")
hLoE.Draw("same")

