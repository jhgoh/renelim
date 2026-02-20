# RENELim: Limit calculation for RENE project

RENELim computes exclusion limits on sterile neutrino oscillation parameters
($\sin^2 2\theta_{14}$, $\Delta m^2_{41}$) using reactor antineutrinos detected
via inverse beta decay (IBD). The analysis targets the Hanbit nuclear power plant
reactors and the RENE near-detector and is built on top of RooFit.

## Project structure

```
config.yaml          Main configuration (detector, reactors, physics)
data/                Physics input tables (flux models, IBD cross section)
python/
  Config.py          YAML configuration loader
  ModelConfig.py     RooFit workspace builder
scripts/
  response_gaus.py   Generate a Gaussian detector response matrix
  baseline_smearing.py  Generate a baseline distribution from core/detector geometry
src/
  NuOscIBDPdf        Unbinned neutrino oscillation IBD RooFit PDF
  BinnedNuOscIBDPdf  Binned version with detector response matrix
test/
  run_chi2.py        NLL scan and p-value calculation via toy MC
  show_chi2.py       Visualise limits (3σ/5σ contours)
  show_energy_spectrum.py  Inspect energy spectrum interactively
  submit_chi2.py     Dispatch grid scan jobs to a SLURM batch system
```

See [`src/README.md`](src/README.md) for the mathematical details and API of the
RooFit PDF classes.

## Initial setup

Install [mamba](https://mamba.readthedocs.io) (or [conda](https://docs.conda.io))
and create the environment:

```bash
mamba create -n hep2026.01 -c conda-forge root numpy pandas scipy tqdm pyyaml -y
mamba activate hep2026.01
```

### ROOT runtime notes

- The scripts in `test/` compile and load C++ classes from `src/` through PyROOT.
- Run commands from the repository root (`/workspace/renelim`) so relative paths like
  `python/`, `data/`, and `config.yaml` are resolved correctly.

## Typical analysis workflow

### 1. Generate the response matrix

The detector response matrix must be created before running the fit.
The default Gaussian smearing matrix is built with:

```bash
python scripts/response_gaus.py
```

This produces `data/response_gaus.root` used by the default `config.yaml`.

### 2. Generate the baseline distribution

To account for the finite sizes of the reactor core and the detector volume,
generate the baseline smearing histogram with:

```bash
python scripts/baseline_smearing.py
```

### 3. Configure

Edit `config.yaml` to match your setup. The file has three main sections:

- **`physics`** – flux model (`huber` or `mueller`) and oscillation parameters
- **`detectors`** – position, efficiency, number of target protons, response matrix
- **`reactors`** – Hanbit units 1–6 positions and thermal powers

### 4. Inspect the energy spectrum

```bash
python test/show_energy_spectrum.py
```

Compiles the required ROOT classes on the fly and opens several canvases
displaying the detector-smeared spectrum. Press `Enter` to close them.

### 5. Run the NLL scan

```bash
mkdir -p results
python test/run_chi2.py -m 1.0 -n 1000 -o results/result_dm41_1.root --toys 1000
```

Useful options:

- `-s/--sin14`: comma-separated custom scan points, e.g. `-s 0,0.001,0.01,0.1`
- `--seed`: random seed for reproducible toy studies
- `-g/--gui`: display ROOT canvases while running

### 6. Visualise results

```bash
python test/show_chi2.py
```

This draws the expected exclusion limits together with 3σ and 5σ contours
from the ROOT files matching `results/result_*_nSignal_1000.root`.

### 7. Batch grid scan

To submit a full ($\sin^2 2\theta_{14}$, $\Delta m^2_{41}$) grid to a SLURM
cluster:

```bash
python test/submit_chi2.py
```

This creates the `results/` directory and dispatches jobs defined in
`test/run_chi2.sbatch`.

## Configuration

The default `config.yaml` assumes:

- Hanbit units 1–6 as reactor sources
- one near detector with a Gaussian response matrix in `data/response_gaus.root`
- neutrino flux/cross-section tables in `data/*.yaml`

Before long production runs, verify file paths and detector/reactor coordinates in `config.yaml`.

### Physics inputs

| File | Contents |
|---|---|
| `data/huber.yaml` | Huber flux for U235, Pu239, Pu241 |
| `data/mueller.yaml` | Mueller flux for U235, U238, Pu239, Pu241 |
| `data/ibdxsec.yaml` | IBD cross section |

Select the model via `physics.isotope_flux` and `physics.ibd_xsec` in `config.yaml`.

### Python modules

| Module | Description |
|---|---|
| `python/Config.py` | YAML configuration loader (`Config`, `ConfigRENE`). Also provides `loadYamlData()` for reading flux/cross-section tables. |
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
