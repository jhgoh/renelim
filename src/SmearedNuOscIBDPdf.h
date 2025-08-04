#ifndef SmearedNuOscIBDPdf_H
#define SmearedNuOscIBDPdf_H

// SmearedNuOscIBDPdf: Energy spectum of sterile neutrinos with detector resolution
// This class gives neutrino energy spectrum with the detector smearing.
//
// In many cases, smeared energy distribution can be obtained with convolution (FFT, etc)
// with a certain kernel K(E-E'), which is not applicable in this case because the
// resolution itself is energy dependent.
//
// Therefore, smearing has to be done by multiplying the response matrix then apply projection
// on the reco-energy axis.
//
/// The reco energy distribution can be written as:
//   P(E') = \int d(E) R(E',E) * F(E) * S(E) ( 1 - sin14*sin^2(K/E) )
// where E' = reconstructed (smeared) energy, E = true energy
//       R(E',E) = response 'matrix'
//       F(E) = neutrino energy flux (Huber-Mueller)
//       S(E) = IBD cross section
//
// Natually, this operation involves numerical integral, which takes high computing cost.
// The RooFit provides multiplication and projection with RooProdPdf() and createProjection().
// The numeric integral has to be performed for each evaluation unless the analytic integral
// is implemented.
//
// The idea of this class is to perform the integrals in advance, since the response matrix
// is given as a 2D histogram with finite number of bins, and the flux, cross sections are
// given as piecewise linear curves where their analytic integral is elementary.
// The coefficients can be computed and cached at the constructor, and only the disappearance
// part is computed with existing (semi) analytic formula,
//   P(E') = P_i = R_ij \sum_j \int_{E_i}^{E_{i+1}) F(E)*S(E) (1-sin14*sin^2(K/E))
// by splitting integration ranges according to the response matrix then extracting response matrix.
// and inside of the integral, the flux and cross section can be written as linear functions
//   F(E) = F(E_i) + (dF/dE) * (E-E_i) and S(E) = S(E_i) + (dS/dE) * (E-E_i)
// then the analytic integration can be obtained. It involves integral of E^n sin^2(K/E),
// where the solution can be found with sine-intetral (TMath::sinint) and cosine-integrals
// (TMath::cosint).

#include "RooAbsCategory.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooCategoryProxy.h"
#include "RooRealProxy.h"

#include "NuOscIBDPdf.h"

#include "TH2.h"
#include "TMatrixD.h"
#include <vector>

class SmearedNuOscIBDPdf : public NuOscIBDPdf {
public:
  SmearedNuOscIBDPdf() = default;
  SmearedNuOscIBDPdf(const char *name, const char *title, RooAbsReal &x, RooAbsReal &xInt,
                     RooAbsReal &l,                       // baseline L (meter)
                     RooAbsReal &sin13, RooAbsReal &dm31, // sin^2(2theta_13), Delta m_13^2
                     RooAbsReal &sim14, RooAbsReal &dm41, // sin^2(2theta_14), Delta m_14^2
                     const RooArgList &elemFracs,
                     const std::vector<const TGraph *> elemSpects, // Fuel composition and spectrum
                     const TGraph *grpXsec,                        // IBD Cross-section curve
                     const TH2 *hResp // Response matrix, (x,y) = (ETrue, EReco)
  );
  SmearedNuOscIBDPdf(const SmearedNuOscIBDPdf &other, const char *name = 0);
  virtual TObject *clone(const char *newname) const override {
    return new SmearedNuOscIBDPdf(*this, newname);
  }
  inline virtual ~SmearedNuOscIBDPdf() override = default;

protected:
  RooRealProxy xr_;
  TMatrixD respMat_;
  std::vector<double> binsT_, binsR_; // bin edges, trueE and recoE

  double evaluate() const override;
  int getAnalyticalIntegral(RooArgSet &allVars, RooArgSet &analVars,
                            const char *rangeName = 0) const override;
  double analyticalIntegral(int code, const char *rangeName = 0) const override;

private:
  ClassDef(SmearedNuOscIBDPdf, 1) // Your description goes here...
};

#endif
