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
#ifndef PAPILLON_NDL_LEGENDRE_H
#define PAPILLON_NDL_LEGENDRE_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/angle_law.hpp>
#include <vector>

namespace pndl {

/**
 * @brief Angular distribution represented by a series of Legendre
 *        polynomials.
 */
class Legendre : public AngleLaw {
 public:
  /**
   * @brief Constructs an isotropic distribution with a 0th order
   *        Legendre series.
   */
  Legendre();

  /**
   * @brief Constructs a Legendre series from coefficients, starting
   *        with the 1st order moment.
   * @param a Vector containing the Legendre moments, from the 1st order
   *          up. These coefficients will be multiplied by (2l + 1)/2.
   */
  Legendre(const std::vector<double>& a);

  double sample_mu(const std::function<double()>& rng) const override final;

  double pdf(double mu) const override final;

  /**
   * @brief Sets the coefficient for the Legendre moment l.
   * @param l Legendre order to set.
   * @param a Coefficient to be used for Legendre order l.
   *          Will be multiplied by (2l + 1)/2.
   */
  void set_moment(std::size_t l, double a);

  /**
   * @brief Returns a reference to the vector of all Legendre coefficients.
   */
  const std::vector<double>& coefficients() const { return a_; }

 private:
  std::vector<double> a_;
  double pdf_max_;

  // Computes the constants in Ci_, H_, and cdf_
  void find_max();

  // Checks that the distribution is positive over the intervale [-1,1].
  bool positive() const;
};

}  // namespace pndl

#endif
