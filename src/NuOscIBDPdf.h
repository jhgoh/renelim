#ifndef NuOscIBDPdf_H
#define NuOscIBDPdf_H

#include "RooAbsCategory.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooArgList.h"
#include "RooCategoryProxy.h"
#include "RooRealProxy.h"

#include <vector>

class NuOscIBDPdf : public RooAbsPdf {
public:
  NuOscIBDPdf() = default;
  NuOscIBDPdf(const char *name, const char *title, RooAbsReal &x,
              RooAbsReal &l,                       // baseline L (meter)
              RooAbsReal &sin13, RooAbsReal &dm31, // sin^2(2theta_13), Delta m_13^2
              RooAbsReal &sim14, RooAbsReal &dm41, // sin^2(2theta_14), Delta m_14^2
              const RooArgList &elemFracs,
              const std::vector<const TGraph *> elemSpects, // Fuel composition and spectrum
              const TGraph *grpXsec                         // IBD Cross-section curve
  );
  NuOscIBDPdf(const NuOscIBDPdf &other, const char *name = 0);
  virtual TObject *clone(const char *newname) const override {
    return new NuOscIBDPdf(*this, newname);
  }
  inline virtual ~NuOscIBDPdf() override = default;

protected:
  RooRealProxy x_;                               // Variable
  RooRealProxy l_, sin13_, dm31_, sin14_, dm41_; // Parameters
  RooArgList elemFracs_;

  std::vector<std::vector<double>> elemSpectsX_;
  std::vector<std::vector<double>> elemSpectsY_;

  std::vector<double> ibdXsecX_, ibdXsecY_;

  double evaluate() const override;
  int getAnalyticalIntegral(RooArgSet &allVars, RooArgSet &analVars,
                            const char *rangeName = 0) const override;
  double analyticalIntegral(int code, const char *rangeName = 0) const override;

  void loadFromTGraph(const TGraph *grp, std::vector<double> &xx, std::vector<double> &yy);
  double interpolate(const double x, const std::vector<double> &xx,
                     const std::vector<double> &yy) const;

private:
  ClassDef(NuOscIBDPdf, 1) // Your description goes here...
};

#endif
