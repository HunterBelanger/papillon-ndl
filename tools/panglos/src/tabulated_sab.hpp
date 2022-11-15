/*
 * Papillon Nuclear Data Library
 * Copyright 2021-2022, Hunter Belanger
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
#include <ndarray.hpp>

#include "interpolator.hpp"
#include "sab.hpp"
#include "short_collision_time_sab.hpp"
using namespace njoy::ENDFtk;

#include <vector>

/**
 * @brief This class respresents a tabulated \f$S(\alpha,\beta)\f$ function, as
 *        is typically provided in File 7 Section 4 of ENDF evaulations.
 **/
class TabulatedSab : public Sab {
 public:
  /**
   * @param TSL The TabulatedFunctions from ENDFtk, containing the tabulated
   *        \f$S(\alpha,\beta)\f$.
   * @param indx_T The index corresponding to temperature T, for the data in
   *        the TSL.
   * @param T Temperature of the scattering law to be read, in K.
   * @param Teff Effective temperature for the Short Collision Time
   *        approximation, in K.
   * @param LAT Flag indicating that the \f$\alpha\f$ and \f$\beta\f$ grids
   *            are calculated at room temperature (when LAT = 1).
   * @param LASYM Flag indicating that the scattering law is asymmetric in
   *        \f$\beta\f$, and negative values are explicitly tabulated (when
   *        LASYM = 1).
   * @param LLN Flag indicating that \f$\ln(S)\f$ is stored in the tabulation,
   *        instead of \f$S\f$ directly (when LLN = 1).
   **/
  TabulatedSab(section::Type<7, 4>::TabulatedFunctions& TSL, std::size_t indx_T,
               double T, double Teff, double A, int LAT, int LASYM, int LLN);

  double operator()(double a, double b) const override final;

  double integrate_alpha(double a_low, double a_hi,
                         double b) const override final;

  double integrate_exp_beta(double E, double b_low,
                            double b_hi) const override final;

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

 private:
  std::vector<double> beta_;
  std::vector<long> beta_bounds_;
  std::vector<Interpolator> beta_interps_;
  std::vector<double> alpha_;
  std::vector<long> alpha_bounds_;
  std::vector<Interpolator> alpha_interps_;
  NDArray<double> data_;
  ShortCollisionTimeSab sct_;
  bool symmetric_;
};

#endif
