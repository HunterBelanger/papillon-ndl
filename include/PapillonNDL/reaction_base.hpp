/*
 * Papillon Nuclear Data Library
 * Copyright 2021, Hunter Belanger
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
#ifndef PAPILLON_NDL_REACTION_BASE_H
#define PAPILLON_NDL_REACTION_BASE_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/angle_energy.hpp>
#include <PapillonNDL/frame.hpp>
#include <PapillonNDL/function_1d.hpp>
#include <memory>

namespace pndl {

/**
 * @brief Holds the non temperature-dependent info and product distributions
 *        for a single MT.
 */
class ReactionBase {
 public:
  ReactionBase(const ReactionBase&) = default;

  /**
   * @brief Returns the MT of the reaction.
   */
  uint32_t mt() const { return mt_; }

  /**
   * @brief Returns the Q-value of the reaction.
   */
  double q() const { return q_; }

  /**
   * @brief Returns the threshold energy for the reaction.
   */
  double threshold() const { return threshold_; }

  /**
   * @brief Returns the function for the reaction yield.
   */
  const Function1D& yield() const { return *yield_; }

  /**
   * @brief Samples and angle and energy from the neutron reaction
   *        product distribution.
   * @param E_in Incident energy in MeV.
   * @param rng Random number generation function.
   */
  AngleEnergyPacket sample_neutron_angle_energy(
      double E_in, std::function<double()> rng) const {
    if (E_in < threshold_) return {0., 0.};

    return neutron_distribution_->sample_angle_energy(E_in, rng);
  }

  /**
   * @brief Returns the distribution for neutron reaction products.
   */
  const AngleEnergy& neutron_distribution() const {
    return *neutron_distribution_;
  }

 protected:
  uint32_t mt_;
  double q_;
  double awr_;
  double threshold_;
  std::shared_ptr<Function1D> yield_;
  std::shared_ptr<AngleEnergy> neutron_distribution_;

  /**
   * @param ace ACE file to take reaction from.
   * @param indx Reaction index in the MT array.
   */
  ReactionBase(const ACE& ace, std::size_t indx);

  // Private helper methods
  void load_neutron_distributions(
      const ACE& ace, std::size_t indx,
      std::vector<std::shared_ptr<AngleEnergy>>& distributions,
      std::vector<std::shared_ptr<Tabulated1D>>& probabilities);
};

}  // namespace pndl

#endif
