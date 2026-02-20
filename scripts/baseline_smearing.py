#!/usr/bin/env python
## Build a baseline distribution due to the reactor core and detector volume
import argparse
import os

parser = argparse.ArgumentParser(description="Generate a baseline distribution due to core and detector")
parser.add_argument('--core-height', type=float, default=3.8,
                    help='Reactor core height')
parser.add_argument('--core-radius', type=float, default=1.75,
                    help='Reactor core radius')
parser.add_argument('--core-axial-shape', choices=['flat', 'cosine'], default='cosine',
                    help='Axial shape of power distribution')
parser.add_argument('--core-radial-shape', choices=['flat', 'quadratic'], default='quadratic',
                    help='Radial shape of power distribution')

parser.add_argument('--det-length', type=float, default=1.2,
                    help='Length of the detector')
parser.add_argument('--det-radius', type=float, default=0.534/2,
                    help='Radius of the detector')
parser.add_argument('--det-orientation', choices=['vertical', 'horizontal'], default='horizontal',
                    help='Orientation of the detector')

parser.add_argument('--core-z', type=float, default=-3.8/2,
                    help='z-position of the reactor core with respect to the ground')
parser.add_argument('--det-z', type=float, default=-11.5+1.5,
                    help='z-position of the detector center with respect to the ground')
parser.add_argument('--det-y', type=float, default=20+1.5,
                    help='horizontal distance between the detector and the core')

parser.add_argument('--hist-nbins', type=int, default=100,
                    help='Number of bins')
parser.add_argument('--hist-xmin', type=float, default=20,
                    help='Minimum of histogram x-axis')
parser.add_argument('--hist-xmax', type=float, default=30,
                    help='Maximum of histogram x-axis')

parser.add_argument('-n', type=int, default=100000,
                    help='Number of MC samples')
parser.add_argument('-o', '--output', type=str, default='baseline.root',
                    help='Output ROOT file name for baseline distribution')
parser.add_argument('-g', '--gui', action='store_true',
                    help='Display the baseline distribution using ROOT TCanvas')
args = parser.parse_args()

import numpy as np
import ROOT

## Set parameters
n = args.n

## Generate points
def generate_cylinder(radius, height, n):
  z = np.random.uniform(-height/2, height/2, n)
  theta = np.random.uniform(0, 2*np.pi, n)

  u = np.random.uniform(0, 1, n)
  r = radius*np.sqrt(u)

  x = r*np.cos(theta)
  y = r*np.sin(theta)

  return x, y, z

coreX, coreY, coreZ = generate_cylinder(args.core_radius, args.core_height, n)
coreR = np.hypot(coreX, coreY)
coreZ += args.core_z

detX, detY, detZ = generate_cylinder(args.det_radius, args.det_length, n)
if args.det_orientation == 'horizontal':
  detZ, detX = detX, detZ 
detY += args.det_y
detZ += args.det_z

## Estimate weight due to the fuel power
wZ = np.ones_like(coreZ)
wR = np.ones_like(coreR)
if args.core_axial_shape == 'cosine':
  wZ = np.cos(np.pi * coreZ / args.core_height )
  wZ = np.clip(wZ, 0, None)
if args.core_radial_shape == 'quadratic':
  wR = 1.0 - (coreR / args.core_radius)**2
  wR = np.clip(wR, 0, None)
w = wZ*wR

## Get the baseline
dx, dy, dz = coreX-detX, coreY-detY, coreZ-detZ
l = np.sqrt(dx*dx + dy*dy + dz*dz)

## Save them
nbins = args.hist_nbins
xmin, xmax = args.hist_xmin, args.hist_xmax
h = ROOT.TH1D("hBaseline", "Baseline;Baseline (m);Arbitrary", nbins, xmin, xmax)
for i in range(n):
  print(f'{i+1}/{n}', end='\r')
  h.Fill(l[i], w[i])
print()
if h.Integral() != 0:
  h.Scale(1./h.Integral())

if args.gui:
  #ROOT.gStyle.SetOptStat(0)
  ROOT.gStyle.SetOptTitle(0)
  c = ROOT.TCanvas("c", "c", 500, 500)
  h.Draw("COLZ")
  c.Update()
  input('Press return to exit')

f = ROOT.TFile(args.output, "recreate")
h.Write()
f.Close()

