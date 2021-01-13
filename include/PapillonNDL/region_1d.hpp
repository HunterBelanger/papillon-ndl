/*
 * Copyright 2021, Hunter Belanger
 *
 * hunter.belanger@gmail.com
 *
 * Ce logiciel est régi par la licence CeCILL soumise au droit français et
 * respectant les principes de diffusion des logiciels libres. Vous pouvez
 * utiliser, modifier et/ou redistribuer ce programme sous les conditions
 * de la licence CeCILL telle que diffusée par le CEA, le CNRS et l'INRIA
 * sur le site "http://www.cecill.info".
 *
 * En contrepartie de l'accessibilité au code source et des droits de copie,
 * de modification et de redistribution accordés par cette licence, il n'est
 * offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
 * seule une responsabilité restreinte pèse sur l'auteur du programme,  le
 * titulaire des droits patrimoniaux et les concédants successifs.
 *
 * A cet égard  l'attention de l'utilisateur est attirée sur les risques
 * associés au chargement,  à l'utilisation,  à la modification et/ou au
 * développement et à la reproduction du logiciel par l'utilisateur étant
 * donné sa spécificité de logiciel libre, qui peut le rendre complexe à
 * manipuler et qui le réserve donc à des développeurs et des professionnels
 * avertis possédant  des  connaissances  informatiques approfondies.  Les
 * utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
 * logiciel à leurs besoins dans des conditions permettant d'assurer la
 * sécurité de leurs systèmes et ou de leurs données et, plus généralement,
 * à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.
 *
 * Le fait que vous puissiez accéder à cet en-tête signifie que vous avez
 * pris connaissance de la licence CeCILL, et que vous en avez accepté les
 * termes.
 *
 * */
#ifndef PAPILLON_NDL_REGION_1D_H
#define PAPILLON_NDL_REGION_1D_H

#include <PapillonNDL/interpolation.hpp>
#include <PapillonNDL/tabulated_1d.hpp>
#include <memory>

namespace pndl {

class Region1D : public Tabulated1D {
 public:
  Region1D(const std::vector<double>& i_x, const std::vector<double>& i_y);
  virtual ~Region1D() = default;

  // Methods required by Function1D
  virtual double operator()(double x) const override = 0;
  virtual double integrate(double x_low, double x_hi) const override = 0;

  // Methods required by Tabulated1D
  std::vector<uint32_t> breakpoints() const override final;
  std::vector<Interpolation> interpolation() const override = 0;
  std::vector<double> x() const override final { return x_; }
  std::vector<double> y() const override final { return y_; }

  // Extra Methods
  size_t size() const { return x_.size(); }
  double min_x() const { return x_.front(); }
  double max_x() const { return x_.back(); }

 protected:
  std::vector<double> x_;
  std::vector<double> y_;
};

template<class I>
class Interp1D : public Region1D {
 public:
  Interp1D(const std::vector<double>& i_x, const std::vector<double>& i_y): Region1D(i_x, i_y), interpolator() {
    // Ensure x and y grids are valid
    interpolator.verify_x_grid(x_.cbegin(), x_.cend());
    interpolator.verify_y_grid(y_.cbegin(), y_.cend());
  }
  ~Interp1D() = default;

  double operator()(double x) const override final {
    if (x <= min_x())
      return y_.front();
    else if (x >= max_x())
      return y_.back();

    // Get bounding x1 < x < x2
    auto low_it = std::lower_bound(x_.begin(), x_.end(), x);
    low_it--;

    auto hi_it = low_it;
    hi_it++;

    size_t i = low_it - x_.begin();

    double x1 = *low_it;
    double x2 = *hi_it;
    double y1 = y_[i];
    double y2 = y_[i + 1];

    return interpolator.interpolate(x, x1, y1, x2, y2);
  }

  double integrate(double x_low, double x_hi) const override final {
    // Integration may only be carried out over the function's valid domain
    if (x_low <= min_x())
      x_low = min_x();
    else if (x_low >= max_x())
      x_low = max_x();

    if (x_hi >= max_x())
      x_hi = max_x();
    else if (x_hi <= min_x())
      x_hi = min_x();

    // Get iterator for lower bound of first interval
    auto low_it = std::lower_bound(x_.begin(), x_.end(), x_low);
    if (*low_it > x_low) low_it--;

    double integral = 0.;
    double x_low_lim = x_low;
    double x_upp_lim = x_hi;
    bool integrating = true;
    while (integrating) {
      auto hi_it = low_it;
      hi_it++;

      size_t i = low_it - x_.begin();

      double x1 = *low_it;
      double x2 = *hi_it;
      double y1 = y_[i];
      double y2 = y_[i + 1];

      if (x_low_lim < x1) x_low_lim = x1;
      if (x_upp_lim > x2) x_upp_lim = x2;

      integral += interpolator.integrate(x_low_lim, x_upp_lim, x1, y1, x2, y2);

      if (x_upp_lim == x_hi)
        integrating = false;
      else {
        x_low_lim = x_upp_lim;
        x_upp_lim = x_hi;
        low_it++;
      }
    }

    return integral;
  }

  std::vector<Interpolation> interpolation() const override final {
    return {interpolator.interpolation};
  }

 private:
  I interpolator;
};

std::shared_ptr<Region1D> build_Region1D(const std::vector<double>& x, const std::vector<double>& y, Interpolation interp);

// This operator overload is provided only to accomodate the
// std::lower_bound algorithm, in the MultiRegion1D class
bool operator<(const Region1D& R, const double& X);
bool operator<(const std::shared_ptr<Region1D>& R, const double& X);

}  // namespace pndl

#endif
