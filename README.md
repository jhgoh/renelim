# RENELim: Limit calculation for RENE project
## Initial setup
One has to install pyROOT to run this tool.
For convenience, numpy is required as well.
```bash
conda create -n ds4hep -c conda-forge
conda activate ds4hep

conda install root -c conda-forge -y
conda install numpy -c conda-forge -y
#conda install matplotlib -c conda-forge -y
```

Prepare neutrino spectrum and cross section table.
```bash
cd data
python ./scripts/huber.py
python ./scripts/mueller.py
python ./scripts/ibdxsec.py
cd ..
```

Then you will have input root files
```
$ ls data
huber.root  ibdxsec.root  mueller.root  neos.root  reno.root  scripts
```

Export data of past experiments, RENO and NEOS data.
This output file is to be used later for cross-experiment results.
```bash
cd data
python ./scripts/reno_neos_to_root.py
cd ..
```

Fit functions has to be compiled using the (py)ROOT.
```bash
python -i fit.py
```

## Reference of original data files
Huber-Mueller flux
- Paper: https://journals.aps.org/prc/abstract/10.1103/PhysRevC.84.024617
- Title: Determination of antineutrino spectra from nuclear reactors
- Authors: Patrick Huber
- URL: https://www1.phys.vt.edu/~pahuber/reactorfluxes/

Mueller ab initio approach
- Paper: https://journals.aps.org/prc/abstract/10.1103/PhysRevC.83.054615
- Title: Improved predictions of reactor antineutrino spectra
- Authors: Th. A. Mueller, D. Lhuillier, M. Fallot, A. Letourneau, S. Cormon, M. Fechner, L. Giot, T. Lasserre, and J. Martino et al.

IBD cross section
- Paper: https://link.springer.com/article/10.1007/JHEP08(2022)212
- Title: An accurate evaluation of electron (anti-)neutrino scattering on nucleons
- Authors: Giulia Ricciardi, Natascia Vignaroli & Francesco Vissani 
