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

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/interpolation.hpp>
#include <PapillonNDL/tabulated_1d.hpp>
#include <memory>
#include <variant>

namespace pndl {

/**
 * @brief Implementation of a Tabulated1D which has only one interpolation
 *        retion.
 */
class Region1D : public Tabulated1D {
 public:
  /**
   * @param i_x Vector of all x points.
   * @param i_y Vector of all y points.
   * @param interp Interpolation method used for all points.
   */
  Region1D(const std::vector<double>& i_x, const std::vector<double>& i_y,
           Interpolation interp);

  // Methods required by Function1D
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

    std::size_t i = low_it - x_.begin();

    double x1 = *low_it;
    double x2 = *hi_it;
    double y1 = y_[i];
    double y2 = y_[i + 1];

    auto doInterp = [&x, &x1, &x2, &y1, &y2](auto& interp) {
      return interp.interpolate(x, x1, y1, x2, y2);
    };
    return std::visit(doInterp, interpolator);
  }

  double integrate(double x_low, double x_hi) const override final {
    bool inverted = x_low > x_hi;
    if (inverted) {
      double x_low_tmp = x_low;
      x_low = x_hi;
      x_hi = x_low_tmp;
    }

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

      std::size_t i = low_it - x_.begin();

      double x1 = *low_it;
      double x2 = *hi_it;
      double y1 = y_[i];
      double y2 = y_[i + 1];

      if (x_low_lim < x1) x_low_lim = x1;
      if (x_upp_lim > x2) x_upp_lim = x2;

      auto doIntegrl = [&x_low_lim, &x_upp_lim, &x1, &x2, &y1,
                        &y2](auto& interp) {
        return interp.integrate(x_low_lim, x_upp_lim, x1, y1, x2, y2);
      };
      integral += std::visit(doIntegrl, interpolator);

      // integral += interpolator.integrate(x_low_lim, x_upp_lim, x1, y1, x2,
      // y2);

      if (x_upp_lim == x_hi)
        integrating = false;
      else {
        x_low_lim = x_upp_lim;
        x_upp_lim = x_hi;
        low_it++;
      }
    }

    if (inverted) integral *= -1.;

    return integral;
  }

  // Methods required by Tabulated1D
  std::vector<uint32_t> breakpoints() const override final {
    return {static_cast<uint32_t>(x_.size())};
  }
  std::vector<Interpolation> interpolation() const override final {
    return {interpolation_};
  }
  std::vector<double> x() const override final { return x_; }
  std::vector<double> y() const override final { return y_; }

  /**
   * @brief Returns the number of (x,y) pairs.
   */
  std::size_t size() const { return x_.size(); }

  /**
   * @brief Returns the lowest x value.
   */
  double min_x() const { return x_.front(); }

  /**
   * @brief Returns the highest x value.
   */
  double max_x() const { return x_.back(); }

 private:
  std::vector<double> x_;
  std::vector<double> y_;
  Interpolation interpolation_;
  std::variant<Histogram, LinLin, LinLog, LogLin, LogLog> interpolator;
};

// This operator overload is provided only to accomodate the
// std::lower_bound algorithm, in the MultiRegion1D class
bool operator<(const Region1D& R, const double& X);

}  // namespace pndl

#endif
