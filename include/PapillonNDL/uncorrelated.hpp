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
#ifndef PAPILLON_NDL_UNCORRELATED_H
#define PAPILLON_NDL_UNCORRELATED_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/angle_distribution.hpp>
#include <PapillonNDL/angle_energy.hpp>
#include <PapillonNDL/energy_law.hpp>
#include <memory>

namespace pndl {

/**
 * @brief AngleDistribution - EnergyLaw pair for secondary distributions
 *        where the angle and energy are not correlated.
 */
class Uncorrelated : public AngleEnergy {
 public:
  /**
   * @param angle Angle distribution for all incoming energies.
   * @param energy Shared pointer to the energy distribution.
   */
  Uncorrelated(const AngleDistribution& angle,
               std::shared_ptr<EnergyLaw> energy);
  ~Uncorrelated() = default;

  AngleEnergyPacket sample_angle_energy(
      double E_in, const std::function<double()>& rng) const override final;

  std::optional<double> angle_pdf(double E_in, double mu) const override final;

  std::optional<double> pdf(double E_in, double mu,
                            double E_out) const override final;

  /**
   * @brief Returns a pointer to the angular distribution.
   */
  const AngleDistribution& angle() const { return angle_; }

  /**
   * @brief Returns a pointer to the energy distribution.
   */
  const EnergyLaw& energy() const { return *energy_; }

 private:
  AngleDistribution angle_;
  std::shared_ptr<EnergyLaw> energy_;
};

}  // namespace pndl

#endif
