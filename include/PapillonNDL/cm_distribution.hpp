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
#ifndef PAPILLON_NDL_CM_DISTRIBUTION_H
#define PAPILLON_NDL_CM_DISTRIBUTION_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/angle_energy.hpp>
#include <PapillonNDL/frame.hpp>
#include <memory>
#include <vector>

namespace pndl {

/**
 * @brief A dsitribution for which the data is provided in the center of
 *        mass frame.
 */
class CMDistribution : public AngleEnergy {
 public:
  /**
   * @param A Atomic weight ratio of the nuclide.
   * @param Q The Q-value of the reaction.
   * @param distribution Pointer to the distribution object in the center of
   * mass frame.
   *
   */
  CMDistribution(double A, double Q, std::shared_ptr<AngleEnergy> distribution)
      : awr_(A), q_(Q), distribution_(distribution) {}

  AngleEnergyPacket sample_angle_energy(
      double E_in, std::function<double()> rng) const override final;

  std::optional<double> angle_pdf(double E_in, double mu) const override final;

  std::optional<double> pdf(double E_in, double mu,
                            double E_out) const override final;
  /**
   * @brief Returns the distribution in the Center of Mass frame.
   */
  const AngleEnergy& distribution() const { return *distribution_; }

  /**
   * @brief Returns the nuclide Atomic Weight Ratio.
   */
  double awr() const { return awr_; }

  /**
   * @brief Returns the Q-value of the reaction.
   */
  double q() const { return q_; }

 private:
  double awr_, q_;
  std::shared_ptr<AngleEnergy> distribution_;
};

}  // namespace pndl

#endif
