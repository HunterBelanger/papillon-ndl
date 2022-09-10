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
#ifndef PAPILLON_NDL_DISCRETE_COSINES_ENERGIES_H
#define PAPILLON_NDL_DISCRETE_COSINES_ENERGIES_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/angle_energy.hpp>

namespace pndl {

/**
 * @brief Class which represents equiprobably and skewed discrete energy and
 *        discrete cosine distributions for incoherent inelastic scattering.
 */
class DiscreteCosinesEnergies : public AngleEnergy {
 public:
  /**
   * @param ace ACE file which contains thermal scattering law.
   */
  DiscreteCosinesEnergies(const ACE& ace);

  /**
   * @brief Struct to contain a discrete outgoing energy, with its
   *        associated discrete cosines.
   */
  struct DiscreteEnergy {
    double energy;               /**< Discrete outgoing energy */
    std::vector<double> cosines; /**< Discrete cosines */
  };

  AngleEnergyPacket sample_angle_energy(
      double E_in, std::function<double()> rng) const override final;

  std::optional<double> angle_pdf(double E_in, double mu) const override final;

  std::optional<double> pdf(double E_in, double mu,
                            double E_out) const override final;

  /**
   * @brief Returns true if the outgoing energies are skewed.
   */
  bool skewed() const { return skewed_; }

  /**
   * @brief Returns vector to the incoming energy grid.
   */
  const std::vector<double>& incoming_energy() const {
    return incoming_energy_;
  }

  /**
   * @brief Returns the vector of outgoing energies for all incoming energies.
   */
  const std::vector<std::vector<DiscreteEnergy>>& outgoing_energies() const {
    return outgoing_energies_;
  }

 private:
  std::vector<double> incoming_energy_;
  std::vector<std::vector<DiscreteEnergy>> outgoing_energies_;
  uint32_t Noe, Nmu;
  bool skewed_;
};

}  // namespace pndl

#endif
