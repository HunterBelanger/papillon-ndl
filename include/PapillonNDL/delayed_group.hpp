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
#ifndef PAPILLON_NDL_DELAYED_GROUP_H
#define PAPILLON_NDL_DELAYED_GROUP_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/energy_law.hpp>
#include <PapillonNDL/tabulated_1d.hpp>
#include <memory>

// The delayed group numbers start at g = 1, and go up !
// g = 0 would be the prompt neutorns.

namespace pndl {

/**
 * @brief Contains data for a delayed neutron group.
 */
class DelayedGroup {
 public:
  /**
   * @param ace ACE file to take delayed neutron data from.
   * @param i Index to the beinning of the delayed group data
   *          in the XSS block.
   * @param g Delayed group index.
   */
  DelayedGroup(const ACE& ace, std::size_t i, std::size_t g);
  ~DelayedGroup() = default;

  /**
   * @brief Returns the decay constant for the group in inverse seconds.
   */
  double decay_constant() const { return decay_constant_; }

  /**
   * @brief Returns the Tabulated1D function for the probability
   *        of selecting the delayed group for a given energy.
   */
  const Tabulated1D& probability() const { return *probability_; }

  /**
   * @brief Samples and energy from the delayed group distribution.
   * @param E Incident energy.
   * @param rng Random number generation function.
   */
  double sample_energy(double E, std::function<double()> rng) const {
    return energy_->sample_energy(E, rng);
  }

  /**
   * @brief Returns the EnergyLaw for the group.
   */
  const EnergyLaw& energy() const { return *energy_; }

 private:
  double decay_constant_;
  std::shared_ptr<Tabulated1D> probability_;
  std::shared_ptr<EnergyLaw> energy_;
};

}  // namespace pndl

#endif
