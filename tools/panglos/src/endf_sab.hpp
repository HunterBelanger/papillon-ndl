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
#ifndef PANGLOS_ENDF_SAB_H
#define PANGLOS_ENDF_SAB_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <ENDFtk/section/7.hpp>

#include "interpolator.hpp"
#include "tabulated_sab.hpp"
using namespace njoy::ENDFtk;

#include <ndarray.hpp>
#include <vector>

/**
 * @brief This class respresents a \f$S(\alpha,\beta)\f$ function which is
 *        tabulated as provided in File 7 Section 4 of ENDF evaulations.
 **/
class ENDFSab : public TabulatedSab {
 public:
  /**
   * @param TSL The TabulatedFunctions from ENDFtk, containing the tabulated
   *        \f$S(\alpha,\beta)\f$.
   * @param indx_T The index corresponding to temperature T, for the data in
   *        the TSL.
   * @param T Temperature of the scattering law to be read, in K.
   * @param Teff Effective temperature for the Short Collision Time
   *        approximation, in K.
   * @param A Atomic weight ratio of the primary scattering nuclide.
   * @param LAT Flag indicating that the \f$\alpha\f$ and \f$\beta\f$ grids
   *            are calculated at room temperature (when LAT = 1).
   * @param LASYM Flag indicating that the scattering law is asymmetric in
   *        \f$\beta\f$, and negative values are explicitly tabulated (when
   *        LASYM = 1).
   * @param LLN Flag indicating that \f$\ln(S)\f$ is stored in the tabulation,
   *        instead of \f$S\f$ directly (when LLN = 1).
   **/
  ENDFSab(section::Type<7, 4>::TabulatedFunctions& TSL, std::size_t indx_T,
          double T, double Teff, double A, int LAT, int LASYM, int LLN);

  double operator()(double a, double b) const override final;

  double integrate_alpha(double a_low, double a_hi,
                         double b) const override final;

  double integrate_exp_beta(double E, double b_low,
                            double b_hi) const override final;

  /**
   * @brief Returns a reference to the \f$\alpha\f$ interpolation boundaries.
   */
  const std::vector<long>& alpha_boundaries() const { return alpha_bounds_; }

  /**
   * @brief Returns a reference to the \f$\alpha\f$ Interpolator instances.
   */
  const std::vector<Interpolator>& alpha_interpolators() const {
    return alpha_interps_;
  }

  /**
   * @brief Returns a reference to the tabulated S values.
   */
  const NDArray<double>& data() const { return data_; }

 private:
  std::vector<long> alpha_bounds_;
  std::vector<Interpolator> alpha_interps_;
  NDArray<double> data_;
};

#endif
