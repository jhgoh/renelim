# RENELim: Limit calculation for RENE project

## Initial setup

Install [conda](https://docs.conda.io) and create the environment:

```bash
conda create -n ds4hep -c conda-forge
conda activate ds4hep
conda install root numpy pandas tqdm pyyaml -c conda-forge -y
```

### Response matrix

The detector response matrix must be generated before running the fit.
A Gaussian smearing response can be created with:

```bash
python scripts/response_gaus.py
```

This produces `data/response_gaus.root` used by the default `config.yaml`.

## Configuration

The behaviour of RENELim is controlled via `config.yaml`. Detector
positions, reactor properties and response matrices can be adjusted by
editing this file to match your setup.

Physics inputs (neutrino flux and IBD cross section) are stored as YAML
tables in `data/`:

| File | Contents |
|---|---|
| `data/huber.yaml` | Huber flux for U235, Pu239, Pu241 |
| `data/mueller.yaml` | Mueller flux for U235, U238, Pu239, Pu241 |
| `data/ibdxsec.yaml` | IBD cross section |

Select which model to use by editing the `physics.isotope_flux` and
`physics.ibd_xsec` keys in `config.yaml`.

## Running the NLL scan

The script `test/run_chi2.py` automatically compiles the RooFit PDF and
performs a negative log-likelihood (NLL) scan. P-values are obtained using
toy Monte Carlo samples.

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

## Visualising the energy spectrum

To inspect the detector smearing and flux models interactively run:

```bash
python test/show_energy_spectrum.py
```

The script compiles the necessary ROOT classes on the fly and opens
several canvases displaying the spectrum. Press `Enter` to close them.

## Python modules

| Module | Description |
|---|---|
| `python/Config.py` | YAML configuration loader (`Config`, `ConfigRENE`). Also provides `loadYamlData()` for reading flux/xsec tables. |
| `python/ModelConfig.py` | Builds the full RooFit workspace from a configuration file (`load_model()`). |

## Reference of original data files

Huber flux (U235, Pu239, Pu241)
- Paper: https://journals.aps.org/prc/abstract/10.1103/PhysRevC.84.024617
- Title: Determination of antineutrino spectra from nuclear reactors
- Authors: Patrick Huber

Mueller flux (U235, U238, Pu239, Pu241)
- Paper: https://journals.aps.org/prc/abstract/10.1103/PhysRevC.83.054615
- Title: Improved predictions of reactor antineutrino spectra
- Authors: Th. A. Mueller et al.

IBD cross section
- Paper: https://link.springer.com/article/10.1007/JHEP08(2022)212
- Title: An accurate evaluation of electron (anti-)neutrino scattering on nucleons
- Authors: Giulia Ricciardi, Natascia Vignaroli & Francesco Vissani
