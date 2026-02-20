# Neutrino Oscillation PDFs

This directory provides RooFit probability density functions for neutrino oscillation analyses.

## NuOscIBDPdf

`NuOscIBDPdf` models the inverse beta decay (IBD) energy spectrum without detector smearing. The prediction for a neutrino of energy $E$ at baseline $L$ is

$$
P(E) = \phi(E)\sigma(E)
\left[1 - \sin^2 2\theta_{13}\sin^2\left(\frac{1.27\Delta m^2_{31}L}{E}\right)
      - \sin^2 2\theta_{14}\sin^2\left(\frac{1.27\Delta m^2_{41}L}{E}\right)\right],
$$

where $\phi(E)$ is the neutrino flux (summed over fuel isotopes weighted by their fractions), $\sigma(E)$ the IBD cross section, $L$ is in metres, and the mass splittings are in $\text{eV}^2$.

### Constructors

Two constructor overloads are provided:

```cpp
// Primary: spectra and cross section as pre-built vectors
NuOscIBDPdf(name, title, x, l, sin13, dm31, sin14, dm41, elemFracs,
            vector<vector<double>> elemSpectsX, vector<vector<double>> elemSpectsY,
            vector<double> ibdXsecX, vector<double> ibdXsecY);

// Legacy: load directly from TGraph objects
NuOscIBDPdf(name, title, x, l, sin13, dm31, sin14, dm41, elemFracs,
            vector<const TGraph*> elemSpects, const TGraph* grpXsec);
```

The legacy constructor converts the TGraph data to vectors internally and delegates to the primary constructor.

## BinnedNuOscIBDPdf

`BinnedNuOscIBDPdf` extends the above by incorporating an energy-dependent detector response described by a response matrix $R(\tilde{E},E)$. The reconstructed energy distribution is

$$
P(\tilde{E}) = \int dE R(\tilde{E},E)\phi(E)\sigma(E)
\left[1 - \sin^2 2\theta_{13}\sin^2\left(\frac{K_{31}}{E}\right)
      - \sin^2 2\theta_{14}\sin^2\left(\frac{K_{41}}{E}\right)\right],
$$

with $K_{ij} = 1.27\Delta m^2_{ij}L$.

Because the response matrix $R$ is provided as a histogram and the flux and cross section are treated as piecewise linear functions, the integral over each true-energy bin can be performed analytically:

$$
P(\tilde{E}_{i}) = \sum_{j} R_{ij} \int_{E_j}^{E_{j+1}} dE \phi(E)\sigma(E)
\left[1 - \sin^2 2\theta_{13}\sin^2\left(\frac{K_{31}}{E}\right)
      - \sin^2 2\theta_{14}\sin^2\left(\frac{K_{41}}{E}\right)\right].
$$

Expanding $\phi$ and $\sigma$ as linear functions within each bin and separating the oscillating term gives integrals of the form $\int E^n \cos(2K/E) dE$, which are expressed in closed form using the sine integral $\text{Si}(x)$ and cosine integral $\text{Ci}(x)$ available in ROOT as `Math::sinint(x)` and `Math::cosint(x)`.

### Constructors

The same two-overload pattern as `NuOscIBDPdf` is available, with an additional `TH2` response matrix argument:

```cpp
// Primary: spectra and cross section as pre-built vectors
BinnedNuOscIBDPdf(name, title, x, xInt, l, sin13, dm31, sin14, dm41, elemFracs,
                  vector<vector<double>> elemSpectsX, vector<vector<double>> elemSpectsY,
                  vector<double> ibdXsecX, vector<double> ibdXsecY,
                  const TH2* hResp);

// Legacy: load directly from TGraph objects
BinnedNuOscIBDPdf(name, title, x, xInt, l, sin13, dm31, sin14, dm41, elemFracs,
                  vector<const TGraph*> elemSpects, const TGraph* grpXsec,
                  const TH2* hResp);
```
