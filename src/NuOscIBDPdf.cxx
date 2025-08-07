#include "Riostream.h"

#include "NuOscIBDPdf.h"
#include "RooAbsCategory.h"
#include "RooAbsReal.h"
// #include <math.h>
#include "Math/SpecFunc.h"
#include "TMath.h"

#include <algorithm>

ClassImp(NuOscIBDPdf);

NuOscIBDPdf::NuOscIBDPdf(const char *name, const char *title, RooAbsReal &x, RooAbsReal &l,
                         RooAbsReal &sin13, RooAbsReal &dm31, RooAbsReal &sin14, RooAbsReal &dm41,
                         const RooArgList &elemFracs, const std::vector<const TGraph *> elemSpects,
                         const TGraph *grpXsec)
    : RooAbsPdf(name, title), x_("x", "x", this, x), l_("l", "l", this, l),
      sin13_("sin13", "sin13", this, sin13), dm31_("dm31", "dm31", this, dm31),
      sin14_("sin14", "sin14", this, sin14), dm41_("dm41", "dm41", this, dm41) {
  assert(elemFracs.getSize() == elemSpects.size());
  elemFracs_.add(elemFracs);

  // Load neutrino energy spectra for each fuel component.
  for (int i = 0; i < elemFracs_.getSize(); ++i) {
    elemSpectsX_.push_back({});
    elemSpectsY_.push_back({});
    loadFromTGraph(elemSpects[i], elemSpectsX_[i], elemSpectsY_[i]);
  }

  // Load the IBD cross-section curve.
  loadFromTGraph(grpXsec, ibdXsecX_, ibdXsecY_);

  // Build a unified set of energy bin edges from all spectra and the cross section.
  std::vector<double> allEdges;
  for (const auto &xs : elemSpectsX_) {
    allEdges.insert(allEdges.end(), xs.begin(), xs.end());
  }
  allEdges.insert(allEdges.end(), ibdXsecX_.begin(), ibdXsecX_.end());
  std::sort(allEdges.begin(), allEdges.end());
  if (!allEdges.empty())
    xEdges_.push_back(allEdges.front());
  const double tol = 1e-7;
  for (auto x : allEdges) {
    if (std::abs(x - xEdges_.back()) > tol)
      xEdges_.push_back(x);
  }
}

NuOscIBDPdf::NuOscIBDPdf(const NuOscIBDPdf &other, const char *name)
    : RooAbsPdf(other, name), x_("x", this, other.x_), l_("l", this, other.l_),
      sin13_("sin13", this, other.sin13_), dm31_("dm31", this, other.dm31_),
      sin14_("sin14", this, other.sin14_), dm41_("dm41", this, other.dm41_),
      elemSpectsX_(other.elemSpectsX_), elemSpectsY_(other.elemSpectsY_),
      ibdXsecX_(other.ibdXsecX_), ibdXsecY_(other.ibdXsecY_), xEdges_(other.xEdges_) {
  elemFracs_.add(other.elemFracs_);
}

void NuOscIBDPdf::loadFromTGraph(const TGraph *grp, std::vector<double> &xx,
                                 std::vector<double> &yy) {
  if (!grp) {
    xx = {{0, 1}};
    yy = {{0, 0}};
    return;
  }

  xx.clear();
  yy.clear();

  for (int i = 0, n = grp->GetN(); i < n; ++i) {
    xx.push_back(grp->GetX()[i]);
    yy.push_back(grp->GetY()[i]);
  }
}

double NuOscIBDPdf::interpolate(const double x, const std::vector<double> &xx,
                                const std::vector<double> &yy) const {
  auto itr = std::lower_bound(xx.begin(), xx.end(), x);
  if (itr == xx.begin() || itr == xx.end()) {
    if (x == xx.front())
      return yy.front(); // at boundary
    else
      return 0; // outside the defined range
  }

  const size_t idx = itr - xx.begin();
  const double x2 = *itr, x1 = *(itr - 1);
  const double y2 = yy[idx], y1 = yy[idx - 1];

  return y1 + (y2 - y1) / (x2 - x1) * (x - x1);
}

int NuOscIBDPdf::getAnalyticalIntegral(RooArgSet &allVars, RooArgSet &analVars,
                                       const char * /*rangeName*/) const {
  if (matchArgs(allVars, analVars, x_))
    return 1;
  return 0;
}

double NuOscIBDPdf::analyticalIntegral(int code, const char *rangeName) const {
  R__ASSERT(code == 1);

  const double l = l_->getVal();

  const double sin13 = sin13_->getVal();
  const double dm31 = dm31_->getVal();
  const double sin14 = sin14_->getVal();
  const double dm41 = dm41_->getVal();

  const double k31 = 1.27 * dm31 * l;
  const double k41 = 1.27 * dm41 * l;

  std::vector<double> elemFracs;
  for (int i = 0; i < elemFracs_.getSize(); ++i) {
    const RooAbsReal &elemFrac = static_cast<const RooAbsReal &>(elemFracs_[i]);
    elemFracs.push_back(elemFrac.getVal());
  }

  double sumW = 0;
  for (int i = 0, nn = xEdges_.size() - 1; i < nn; ++i) {
    const double e0 = xEdges_[i], e1 = xEdges_[i + 1];

    double f0 = 0, f1 = 0;
    for (int j = 0; j < elemFracs.size(); ++j) {
      f0 += elemFracs[j] * interpolate(e0, elemSpectsX_[j], elemSpectsY_[j]);
      f1 += elemFracs[j] * interpolate(e1, elemSpectsX_[j], elemSpectsY_[j]);
    }
    const double s0 = interpolate(e0, ibdXsecX_, ibdXsecY_);
    const double s1 = interpolate(e1, ibdXsecX_, ibdXsecY_);

    sumW += subIntegral(e0, e1, f0, f1, s0, s1, sin13, k31, sin14, k41);
  }

  return sumW;
}

double NuOscIBDPdf::subIntegral(const double e0, const double e1, const double f0, const double f1,
                                const double s0, const double s1, const double s13,
                                const double k31, const double s14, const double k41) const {
  const double dE = e1 - e0;
  if (e0 <= 0 || dE <= 0)
    return 0;

  const double dFdE = (f1 - f0) / dE;
  const double dSdE = (s1 - s0) / dE;

  const double a3 = (dFdE * dSdE) / 3;
  const double b2 = (f0 * dSdE + dFdE * s0 - 2 * dFdE * dSdE * e0) / 2;
  const double c1 = (f0 - dFdE * e0) * (s0 - dSdE * e0);

  const double kk31e0 = 2 * k31 / e0;
  const double si31e0 = ROOT::Math::sinint(kk31e0);
  const double ci31e0 = kk31e0 <= 0 ? 0 : ROOT::Math::cosint(kk31e0);
  const double sin31e0 = std::sin(kk31e0), cos31e0 = std::cos(kk31e0);

  const double kk31e1 = 2 * k31 / e1;
  const double si31e1 = ROOT::Math::sinint(kk31e1);
  const double ci31e1 = kk31e1 <= 0 ? 0 : ROOT::Math::cosint(kk31e1);
  const double sin31e1 = std::sin(kk31e1), cos31e1 = std::cos(kk31e1);

  const double kk41e0 = 2 * k41 / e0;
  const double si41e0 = ROOT::Math::sinint(kk41e0);
  const double ci41e0 = kk41e0 <= 0 ? 0 : ROOT::Math::cosint(kk41e0);
  const double sin41e0 = std::sin(kk41e0), cos41e0 = std::cos(kk41e0);

  const double kk41e1 = 2 * k41 / e1;
  const double si41e1 = ROOT::Math::sinint(kk41e1);
  const double ci41e1 = kk41e1 <= 0 ? 0 : ROOT::Math::cosint(kk41e1);
  const double sin41e1 = std::sin(kk41e1), cos41e1 = std::cos(kk41e1);

  // Envelope term
  const double int0 = e0 * (a3 * e0 * e0 + b2 * e0 + c1);
  const double int1 = e1 * (a3 * e1 * e1 + b2 * e1 + c1);

  // Disappearance term for $\nu_3$
  const double int31A0 = 2 * k31 * k31 * k31 * si31e0 +
                         (k31 * k31 * e0 - e0 * e0 * e0 / 2) * cos31e0 + e0 * e0 * e0 / 2 +
                         k31 / 2 * e0 * e0 * sin31e0;
  const double int31B0 =
      -2 * k31 * k31 * ci31e0 + k31 * e0 * sin31e0 - e0 * e0 / 2 * cos31e0 + e0 * e0 / 2;
  const double int31C0 = -k31 * si31e0 + e0 / 2 * (1 - cos31e0);

  const double int31A1 = 2 * k31 * k31 * k31 * si31e1 +
                         (k31 * k31 * e1 - e1 * e1 * e1 / 2) * cos31e1 + e1 * e1 * e1 / 2 +
                         k31 / 2 * e1 * e1 * sin31e1;
  const double int31B1 =
      -2 * k31 * k31 * ci31e1 + k31 * e1 * sin31e1 - e1 * e1 / 2 * cos31e1 + e1 * e1 / 2;
  const double int31C1 = -k31 * si31e1 + e1 / 2 * (1 - cos31e1);

  // Disappearance term for $\nu_4$
  const double int41A0 = 2 * k41 * k41 * k41 * si41e0 +
                         (k41 * k41 * e0 - e0 * e0 * e0 / 2) * cos41e0 + e0 * e0 * e0 / 2 +
                         k41 / 2 * e0 * e0 * sin41e0;
  const double int41B0 =
      -2 * k41 * k41 * ci41e0 + k41 * e0 * sin41e0 - e0 * e0 / 2 * cos41e0 + e0 * e0 / 2;
  const double int41C0 = -k41 * si41e0 + e0 / 2 * (1 - cos41e0);

  const double int41A1 = 2 * k41 * k41 * k41 * si41e1 +
                         (k41 * k41 * e1 - e1 * e1 * e1 / 2) * cos41e1 + e1 * e1 * e1 / 2 +
                         k41 / 2 * e1 * e1 * sin41e1;
  const double int41B1 =
      -2 * k41 * k41 * ci41e1 + k41 * e1 * sin41e1 - e1 * e1 / 2 * cos41e1 + e1 * e1 / 2;
  const double int41C1 = -k41 * si41e1 + e1 / 2 * (1 - cos41e1);

  // Sum up values
  const double intEnv = int1 - int0;
  const double int31 =
      a3 * (int31A1 - int31A0) + b2 * (int31B1 - int31B0) + c1 * (int31C1 - int31C0);
  const double int41 =
      a3 * (int41A1 - int41A0) + b2 * (int41B1 - int41B0) + c1 * (int41C1 - int41C0);

  const double cq14 = TMath::Sq(1+std::sqrt(1-s14))/4; // cos^4(theta) = ((1 +- sqrt(1-sin^2(2theta)))/2)^2
  return (intEnv - s13 * cq14 * int31 - s14 * int41) / dE;
}

double NuOscIBDPdf::evaluate() const {
  const double x = x_->getVal();
  const double l = l_->getVal();
  if (x <= 0)
    return 0; // safeguard for unphysical energy

  const double sin14 = sin14_->getVal();
  const double dm41 = dm41_->getVal();
  const double sinD14 = std::sin(dm41 * 1.27 * l / x); // 1.27 = 1/(4\hbar c)
  const double prob14 = 1 - sin14 * sinD14 * sinD14;

  const double sin13 = sin13_->getVal();
  const double dm31 = dm31_->getVal();
  const double sinD13 = std::sin(dm31 * 1.27 * l / x); // 1.27 = 1/(4\hbar c)
  const double cq14 = TMath::Sq(1+std::sqrt(1-sin14))/4; // cos^4(theta) = ((1 +- sqrt(1-sin^2(2theta)))/2)^2
  const double prob13 = prob13 - sin13 * cq14 * sinD13 * sinD13;

  double spec = 0;
  for (int i = 0; i < elemFracs_.getSize(); ++i) {
    const RooAbsReal &elemFrac = static_cast<const RooAbsReal &>(elemFracs_[i]);
    spec += elemFrac.getVal() * interpolate(x, elemSpectsX_[i], elemSpectsY_[i]);
  }
  const double xsec = interpolate(x, ibdXsecX_, ibdXsecY_);

  return spec * xsec * prob14;
}
