# scripts/: input-generation utilities

This directory contains helper scripts that generate ROOT histogram inputs used by
`BinnedNuOscIBDPdf` in the main analysis.

## Overview

The binned RooFit PDF depends on two external histograms:

1. **Detector response matrix** (`TH2`) mapping true neutrino energy to reconstructed energy.
2. **Baseline-smearing distribution** (`TH1`) describing the effective baseline spread from
   finite reactor-core and detector sizes.

These are produced by:

- `response_gaus.py`
- `baseline_smearing.py`

Both scripts are intended as *pre-processing* tools. Re-run them whenever detector
resolution assumptions or reactor/detector geometry settings are changed.

---

## `response_gaus.py`

Generate a Gaussian detector response matrix in ROOT format.

### What it does

- Defines true-energy and reconstructed-energy binning on a common range.
- Builds a response matrix where each true-energy bin is smeared by a Gaussian:

  \[
  \sigma(E) = \sqrt{a^2 E + b^2 E^2 + c^2}
  \]

- Normalizes each true-energy column to 1.
- Stores the matrix as `TH2D` named `hResp_ETrue_vs_EReco`.

### Default usage

```bash
python scripts/response_gaus.py -o data/response_gaus.root
```

### Important options

- `--a`, `--b`, `--c`: energy-resolution coefficients.
- `--xmin`, `--xmax`, `--dx`: binning range and step in MeV.
- `-o/--output`: output ROOT filename.
- `-g/--gui`: show matrix on a ROOT canvas.

### Notes

- If `dx` does not divide `(xmax - xmin)` exactly, the script adjusts the actual `xmax`
  and prints a warning.
- The output is consumed in `config.yaml` via detector response-matrix path settings.

---

## `baseline_smearing.py`

Generate a baseline distribution histogram from reactor-core and detector geometry.

### What it does

- Samples points in a cylindrical reactor core and cylindrical detector volume.
- Supports configurable detector orientation (`vertical`/`horizontal`).
- Applies optional axial/radial core power profiles.
- Computes event-by-event baseline distances and fills `TH1D` `hBaseline`.
- Normalizes the histogram to unit integral.

### Default usage

```bash
python scripts/baseline_smearing.py -o data/baseline.root
```

### Important options

- Core geometry/profile:
  - `--core-height`, `--core-radius`
  - `--core-axial-shape` (`flat`, `cosine`)
  - `--core-radial-shape` (`flat`, `quadratic`)
- Detector geometry:
  - `--det-length`, `--det-radius`
  - `--det-orientation` (`vertical`, `horizontal`)
- Relative placement:
  - `--core-z`, `--det-z`, `--det-y`
- Histogram control:
  - `--hist-nbins`, `--hist-xmin`, `--hist-xmax`
- Sampling/statistics:
  - `-n` number of Monte Carlo samples
- Output/visualization:
  - `-o/--output`, `-g/--gui`

### Notes

- Larger `-n` gives a smoother baseline histogram but increases runtime.
- Regenerate this file whenever geometry parameters change.
- In the binned PDF path, baseline dependence is driven by this histogram.

---

## Suggested workflow

1. Generate response matrix:
   - `python scripts/response_gaus.py -o data/response_gaus.root`
2. Generate baseline histogram:
   - `python scripts/baseline_smearing.py -o data/baseline.root`
3. Reference both outputs from `config.yaml`.
4. Run fitting and scan scripts from `test/`.
