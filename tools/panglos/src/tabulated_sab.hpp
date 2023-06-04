/*
 * Papillon Nuclear Data Library
 * Copyright 2021-2023, Hunter Belanger
 *
 * hunter.belanger@gmail.com
 *
 * This file is part of the Papillon Nuclear Data Library (PapillonNDL).
 *
 * PapillonNDL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PapillonNDL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PapillonNDL. If not, see <https://www.gnu.org/licenses/>.
 *
 * */
#ifndef PANGLOS_TABULATED_SAB_H
#define PANGLOS_TABULATED_SAB_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <ENDFtk/section/7.hpp>

#include "interpolator.hpp"
#include "sab.hpp"
#include "short_collision_time_sab.hpp"
using namespace njoy::ENDFtk;

#include <vector>

/**
 * @brief This class respresents a tabulated \f$S(\alpha,\beta)\f$ function.
 **/
class TabulatedSab : public Sab {
 public:
  /**
   * @param T Temperature of the scattering law to be read, in K.
   * @param Teff Effective temperature for the Short Collision Time
   *        approximation, in K.
   * @param A Atomic weight ratio of the primary scattering nuclide.
   **/
  TabulatedSab(double T, double Teff, double A) : Sab(T, A), sct_(T, Teff, A) {}

  /**
   * @brief Returns a reference to the ShortCollisionTimeSab which is used for
   *        \f$\alpha\f$ and \f$\beta\f$ values outside the grid.
   */
  const ShortCollisionTimeSab& sct() const { return sct_; }

  /**
   * @brief If ture, then only positive values are stored in the \f$\beta\f$
   *        grid, and negative values of \f$\beta\f$ can be evaluated by using
   *        the absolute value.
   **/
  bool symmetric() const { return symmetric_; }

  /**
   * @brief Returns a reference to the \f$\beta\f$ grid.
   */
  const std::vector<double>& beta() const { return beta_; }

  /**
   * @brief Returns a reference to the \f$\beta\f$ interpolation boundaries.
   */
  const std::vector<long>& beta_boundaries() const { return beta_bounds_; }

  /**
   * @brief Returns a reference to the \f$\beta\f$ Interpolator instances.
   */
  const std::vector<Interpolator>& beta_interpolators() const {
    return beta_interps_;
  }

  /**
   * @brief Returns a reference to the \f$\alpha\f$ grid.
   */
  const std::vector<double>& alpha() const { return alpha_; }

 protected:
  std::vector<double> beta_;
  std::vector<long> beta_bounds_;
  std::vector<Interpolator> beta_interps_;
  std::vector<double> alpha_;
  ShortCollisionTimeSab sct_;
  bool symmetric_;
};

#endif
