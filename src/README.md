# Neutrino Oscillation PDFs

This directory provides RooFit probability density functions for neutrino oscillation analyses.

## NuOscIBDPdf

`NuOscIBDPdf` models the inverse beta decay (IBD) energy spectrum without detector smearing. The prediction for a neutrino of energy $E$ at baseline $L$ is

$$
P(E) = F(E) S(E)
\left[1 - \sin^2 2\theta_{13}\sin^2\left(\frac{1.27\,\Delta m^2_{31} L}{E}\right)
      - \sin^2 2\theta_{14}\sin^2\left(\frac{1.27\,\Delta m^2_{41} L}{E}\right)\right],
$$

where $F(E)$ is the neutrino flux, $S(E)$ the IBD cross section, $L$ is in metres and the mass splittings are in $\text{eV}^2$.

## SmearedNuOscIBDPdf

`SmearedNuOscIBDPdf` extends the above by incorporating an energy-dependent detector response described by a response matrix $R(E',E)$. The reconstructed energy distribution is

$$
P(E') = \int dE\, R(E',E) F(E) S(E)
\left[1 - \sin^2 2\theta_{13}\sin^2\left(\frac{K_{31}}{E}\right)
      - \sin^2 2\theta_{14}\sin^2\left(\frac{K_{41}}{E}\right)\right],
$$

with $K_{ij} = 1.27\,\Delta m^2_{ij} L$. The flux and cross section are treated as piecewise linear functions, allowing the integrals over each true-energy bin to be performed analytically using sine and cosine integral functions.
