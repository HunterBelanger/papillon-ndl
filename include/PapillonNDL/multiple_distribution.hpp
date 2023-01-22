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
#ifndef PAPILLON_NDL_MULTIPLE_DISTRIBUTION_H
#define PAPILLON_NDL_MULTIPLE_DISTRIBUTION_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/angle_energy.hpp>
#include <PapillonNDL/tabulated_1d.hpp>
#include <vector>

namespace pndl {

/**
 * @brief A dsitribution which is composed of mutliple possible
 *        distributions, each with a tabulated probability.
 */
class MultipleDistribution : public AngleEnergy {
 public:
  MultipleDistribution(
      const std::vector<std::shared_ptr<AngleEnergy>>& distributions,
      const std::vector<std::shared_ptr<Tabulated1D>>& probabilities);

  AngleEnergyPacket sample_angle_energy(
      double E_in, const std::function<double()>& rng) const override final;

  std::optional<double> angle_pdf(double E_in, double mu) const override final;

  std::optional<double> pdf(double E_in, double mu,
                            double E_out) const override final;

  /**
   * @brief Returns the number of distributions for the reaction.
   */
  std::size_t size() const { return distributions_.size(); }

  /**
   * @brief Returns the ith distribution for the reaction.
   * @param i Index of distribution to fetch.
   */
  const AngleEnergy& distribution(std::size_t i) const {
    return *distributions_[i];
  }

  /**
   * @brief Returns the ith distribution's probability function.
   * @param i Index of distribution to fetch.
   */
  const Tabulated1D& probability(std::size_t i) const {
    return *probabilities_[i];
  }

 private:
  std::vector<std::shared_ptr<AngleEnergy>> distributions_;
  std::vector<std::shared_ptr<Tabulated1D>> probabilities_;
};

}  // namespace pndl

#endif
