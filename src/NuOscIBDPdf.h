#ifndef NuOscIBDPdf_H
#define NuOscIBDPdf_H

#include "RooAbsCategory.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooArgList.h"
#include "RooCategoryProxy.h"
#include "RooRealProxy.h"

class NuOscIBDPdf : public RooAbsPdf {
public:
  /**
   * Construct the unsmeared neutrino oscillation spectrum.
   *
   * @param x      Neutrino energy.
   * @param l      Baseline length in metres.
   * @param sin13  \f$\sin^2 2\theta_{13}\f$.
   * @param dm31   \f$\Delta m^2_{31}\f$.
   * @param sin14  \f$\sin^2 2\theta_{14}\f$.
   * @param dm41   \f$\Delta m^2_{41}\f$.
   * @param elemFracs  Fuel composition fractions.
   * @param elemSpects Fuel spectra for each element.
   * @param grpXsec    IBD cross-section curve.
   */
  NuOscIBDPdf(const char *name, const char *title, RooAbsReal &x, RooAbsReal &l, RooAbsReal &sin13,
              RooAbsReal &dm31, RooAbsReal &sin14, RooAbsReal &dm41, const RooArgList &elemFracs,
              const std::vector<const TGraph *> elemSpects, const TGraph *grpXsec);
  NuOscIBDPdf(const NuOscIBDPdf &other, const char *name = 0);
  virtual TObject *clone(const char *newname) const override {
    return new NuOscIBDPdf(*this, newname);
  }
  inline ~NuOscIBDPdf() override = default;

protected:
  RooRealProxy x_;                               //!< Energy variable
  RooRealProxy l_, sin13_, dm31_, sin14_, dm41_; //!< Oscillation parameters
  RooArgList elemFracs_;                         //!< Fuel fractions

  std::vector<std::vector<double>> elemSpectsX_; //!< Spectrum energies
  std::vector<std::vector<double>> elemSpectsY_; //!< Spectrum values
  std::vector<double> ibdXsecX_, ibdXsecY_;      //!< Cross-section curve
  std::vector<double> xEdges_;                   //!< Combined x bin edges

  double evaluate() const override;
  int getAnalyticalIntegral(RooArgSet &allVars, RooArgSet &analVars,
                            const char *rangeName = 0) const override;
  double analyticalIntegral(int code, const char *rangeName = 0) const override;

  void loadFromTGraph(const TGraph *grp, std::vector<double> &xx, std::vector<double> &yy);
  double interpolate(const double x, const std::vector<double> &xx,
                     const std::vector<double> &yy) const;
  double subIntegral(const double e0, const double e1, const double f0, const double f1,
                     const double s0, const double s1, const double s13, const double k31,
                     const double s14, const double k41) const;

private:
  ClassDef(NuOscIBDPdf, 1) // RooFit class definition
};

#endif
