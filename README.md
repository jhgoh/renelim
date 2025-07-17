# RENELim: Limit calculation for RENE project
## Initial setup
One has to install pyROOT to run this tool.
For convenience, numpy is required as well.
```bash
conda create -n ds4hep -c conda-forge
conda activate ds4hep

conda install root -c conda-forge -y
conda install numpy pandas tqdm pyyaml -c conda-forge -y
#conda install matplotlib -c conda-forge -y
```

Prepare neutrino spectrum and cross section table.
```bash
cd data
mkdir original
for ELEM in U235 Pu239 Pu241; do
  curl https://www1.phys.vt.edu/~pahuber/reactorfluxes/$ELEM-anti-neutrino-flux-250keV.dat > original/$ELEM-anti-neutrino-flux-250keV.dat
done

python ./scripts/huber.py
python ./scripts/mueller.py
python ./scripts/ibdxsec.py
cd ..
```

Then you will have input root files
```
$ ls data
huber.root  ibdxsec.root  mueller.root  scripts
```

Export data of past experiments, RENO and NEOS data.
This output file is to be used later for cross-experiment results.
```bash
cd data
python ./scripts/reno_neos_to_root.py
cd ..
```

## Running the NLL scan
The script `test/run_chi2.py` automatically compiles the RooFit PDF and performs a negative log-likelihood (NLL) scan.  P-values are obtained using toy Monte Carlo samples.
```bash
python test/run_chi2.py -m 1.0 -n 1000 -o results/result_dm41_1.root --toys 1000
```
The resulting ROOT files can be visualised with:
```bash
python test/show_chi2.py
```
This script draws the expected limits together with 3σ and 5σ contours.

For a full grid scan on a batch system you can submit many jobs using
```bash
python test/submit_chi2.py
```
This script creates a `results` directory and dispatches `sbatch`
jobs defined in `test/run_chi2.sbatch`.

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
