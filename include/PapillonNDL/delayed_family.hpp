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
#ifndef PAPILLON_NDL_DELAYED_FAMILY_H
#define PAPILLON_NDL_DELAYED_FAMILY_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/energy_law.hpp>
#include <PapillonNDL/tabulated_1d.hpp>
#include <memory>

// The delayed family numbers start at g = 1, and go up !
// g = 0 would be the prompt neutorns.

namespace pndl {

/**
 * @brief Contains data for a delayed neutron family.
 */
class DelayedFamily {
 public:
  /**
   * @param ace ACE file to take delayed neutron data from.
   * @param i Index to the beinning of the delayed family data
   *          in the XSS block.
   * @param g Delayed family index.
   */
  DelayedFamily(const ACE& ace, std::size_t i, std::size_t g);
  ~DelayedFamily() = default;

  /**
   * @brief Returns the decay constant for the family in inverse seconds.
   */
  double decay_constant() const { return decay_constant_; }

  /**
   * @brief Returns the Tabulated1D function for the probability
   *        of selecting the delayed family for a given energy.
   */
  const Tabulated1D& probability() const { return *probability_; }

  /**
   * @brief Samples and energy from the delayed family distribution.
   * @param E Incident energy.
   * @param rng Random number generation function.
   */
  double sample_energy(double E, std::function<double()> rng) const {
    return energy_->sample_energy(E, rng);
  }

  /**
   * @brief Returns the EnergyLaw for the family.
   */
  const EnergyLaw& energy() const { return *energy_; }

 private:
  double decay_constant_;
  std::shared_ptr<Tabulated1D> probability_;
  std::shared_ptr<EnergyLaw> energy_;
};

}  // namespace pndl

#endif
