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
#ifndef PAPILLON_NDL_ANGLE_LAW_H
#define PAPILLON_NDL_ANGLE_LAW_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <functional>
#include <memory>

namespace pndl {

/**
 * @brief Interface to represent an angular distribution for a single energy.
 */
class AngleLaw : public std::enable_shared_from_this<AngleLaw> {
 public:
  virtual ~AngleLaw() = default;

  /**
   * @brief Samples a scattering cosine from the distribution.
   * @param xi Random variable from the interval [0,1).
   * @param rng Random number generator function.
   */
  virtual double sample_mu(const std::function<double()>& rng) const = 0;

  /**
   * @brief Returns the PDF for the desired scattering cosine.
   * @param mu Scatter cosnine at which to evaluate the PDF.
   */
  virtual double pdf(double mu) const = 0;
};

}  // namespace pndl

#endif
