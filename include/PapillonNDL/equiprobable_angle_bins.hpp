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
#ifndef PAPILLON_NDL_EQUIPROBABLE_ANGLE_BINS_H
#define PAPILLON_NDL_EQUIPROBABLE_ANGLE_BINS_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/angle_law.hpp>

namespace pndl {

/**
 * @brief Angular distribution represented as equiprobable cosine bins.
 */
class EquiprobableAngleBins : public AngleLaw {
 public:
  /**
   * @param ace ACE file to take data from.
   * @param i Starting index of distribution in the XSS array.
   */
  EquiprobableAngleBins(const ACE& ace, std::size_t i);

  /**
   * @param bounds Vector of 33 bin bounds.
   */
  EquiprobableAngleBins(const std::vector<double>& bounds);
  ~EquiprobableAngleBins() = default;

  double sample_mu(const std::function<double()>& rng) const override final;

  double pdf(double mu) const override final;

  /**
   * @brief Returns the number of bin boundaries (number of bins + 1);
   */
  std::size_t size() const { return bounds_.size(); }

  /**
   * @brief Returns the vector with the bin boundaries.
   */
  const std::vector<double>& bin_bounds() const { return bounds_; }

 private:
  std::vector<double> bounds_;

  static constexpr size_t NBOUNDS = 33;
  static constexpr size_t NBINS = 32;
  static constexpr double P_BIN = 1. / 32.;
};

}  // namespace pndl

#endif
