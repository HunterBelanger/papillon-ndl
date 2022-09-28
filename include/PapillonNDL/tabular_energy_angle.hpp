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
#ifndef PAPILLON_NDL_TABULAR_ENERGY_ANGLE_H
#define PAPILLON_NDL_TABULAR_ENERGY_ANGLE_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/angle_energy.hpp>
#include <PapillonNDL/energy_angle_table.hpp>

namespace pndl {

/**
 * @brief A product Angle-Energy distribution where the angle and energy
 *        PDFs and CDFs are tabulated for different incoming energies.
 */
class TabularEnergyAngle : public AngleEnergy {
 public:
  /**
   * @param ace ACE file to take data from.
   * @param i Starting index of distribution in the XSS array.
   * @param JED Relative index for finding energy and angle distributions.
   */
  TabularEnergyAngle(const ACE& ace, std::size_t i, std::size_t JED);

  /**
   * @param incoming_energy Incoming energy grid.
   * @param tables vector of EnergyAngleTable for each point in the
   *               incoming energy grid.
   */
  TabularEnergyAngle(const std::vector<double>& incoming_energy,
                     const std::vector<EnergyAngleTable>& tables);
  ~TabularEnergyAngle() = default;

  AngleEnergyPacket sample_angle_energy(
      double E_in, const std::function<double()>& rng) const override final;

  std::optional<double> angle_pdf(double E_in, double mu) const override final;

  std::optional<double> pdf(double E_in, double mu,
                            double E_out) const override final;

  /**
   * @brief Returns a vector to the grid of incoming energies.
   */
  const std::vector<double>& incoming_energy() const {
    return incoming_energy_;
  }

  /**
   * @brief Returns the ith incoming energy in MeV.
   * @param i Index to the incoming energy grid.
   */
  double incoming_energy(std::size_t i) const { return incoming_energy_[i]; }

  /**
   * @brief Returns an EnergyAngleTable which contains the distributions
   *        for the ith incoming energy.
   * @param i Index to the incoming energy.
   */
  const EnergyAngleTable& table(std::size_t i) const { return tables_[i]; }

  /**
   * @brief Returns the number of incoming energy points.
   */
  std::size_t size() const { return incoming_energy_.size(); }

 private:
  std::vector<double> incoming_energy_;
  std::vector<EnergyAngleTable> tables_;
};

}  // namespace pndl

#endif
